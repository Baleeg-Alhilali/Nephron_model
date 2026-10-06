/*
 * Copyright (c) 2010-2017:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "navier_stokes_fv.h"

#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/fvho_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesFV<TDomain>::NavierStokesFV(const char* functions,
                                          const char* subsets)
: IncompressibleNavierStokesBase<TDomain>(functions, subsets)
{
	init();
};

template<typename TDomain>
NavierStokesFV<TDomain>::NavierStokesFV(const std::vector<std::string>& vFct,
                                          const std::vector<std::string>& vSubset)
: IncompressibleNavierStokesBase<TDomain>(vFct, vSubset)
{
	init();
};


template<typename TDomain>
void NavierStokesFV<TDomain>::init()
{
//	check number of functions
	if(this->num_fct() != dim+1)
		UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");
//	register imports
	this->register_import(m_imSource);
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensitySCVF);
	this->register_import(m_imDensitySCVFp);
	this->register_import(m_imDensitySCV);

	m_bQuadOrderUserDef = false;

	m_imSource.set_rhs_part();
	m_imDensitySCV.set_mass_part();

	//	default value for density
	base_type::set_density(1.0);

	//	ensure that we do not use the virtual assembling functions at all
	this->clear_add_fct();
}

template<typename TDomain>
void NavierStokesFV<TDomain>::set_quad_order(size_t order)
{
	m_quadOrder = order;
	m_bQuadOrderUserDef = true;
}

template<typename TDomain>
void NavierStokesFV<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	if(bNonRegularGrid)
		UG_THROW("NavierStokes: only regular grid implemented.");

//	check number
	if(vLfeID.size() != dim+1)
		UG_THROW("NavierStokes: Needa exactly "<<dim+1<<" functions");

	for(int d = 1; d < dim; ++d)
		if(vLfeID[0] != vLfeID[d])
			UG_THROW("NavierStokes: trial spaces for velocity expected to be"
					" identical for all velocity components.");

//	remember lfeID;
	m_vLFEID = vLfeID[0];
	m_pLFEID = vLfeID[dim];
	if(!m_bQuadOrderUserDef) m_quadOrder = std::max(m_vLFEID.order(), m_pLFEID.order()) + 1;

	//	update assemble functions
	register_all_funcs(m_vLFEID, m_pLFEID);
}

template<typename TDomain>
void NavierStokesFV<TDomain>::
set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > data)
{
	m_imKinViscosity.set_data(data);
}

template<typename TDomain>
void NavierStokesFV<TDomain>::
set_density(SmartPtr<CplUserData<number, dim> > data)
{
	m_imDensitySCVF.set_data(data);
	m_imDensitySCVFp.set_data(data);
	m_imDensitySCV.set_data(data);
}

template<typename TDomain>
void NavierStokesFV<TDomain>::
set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
{
	m_imSource.set_data(data);
}

////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
//	check, that kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("NavierStokes: Kinematic Viscosity has not been set, but is required.");

//	check, that Density has been set
	if(!m_imDensitySCVF.data_given())
		UG_THROW("NavierStokes: Density has not been set, but is required.");

//	check, that Density has been set
	if(!m_imDensitySCV.data_given())
		UG_THROW("NavierStokes: Density has not been set, but is required.");

	VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
	PGeom& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
	try{
		vgeo.update_local(roid, m_vLFEID);
		pgeo.update_local(roid, m_pLFEID);
	}
	UG_CATCH_THROW("NavierStokes: Cannot update Finite Volume Geometry.");

	static const int refDim = TElem::dim;
	{
		const MathVector<dim>* vSCVFip = vgeo.scvf_global_ips();
		const size_t numSCVFip = vgeo.num_scvf_ips();
		const MathVector<dim>* vSCVip = vgeo.scv_global_ips();
		const size_t numSCVip = vgeo.num_scv_ips();

		m_imKinViscosity.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imDensitySCVF.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imDensitySCV.template set_local_ips<refDim>(vSCVip,numSCVip);
		m_imSource.template set_local_ips<refDim>(vSCVip,numSCVip);
	}

	{
		m_imDensitySCVFp.template set_local_ips<refDim>(pgeo.scvf_local_ips(),
		                                               pgeo.num_scvf_ips());

		const LocalShapeFunctionSet<dim>& rVTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid,m_vLFEID);
		const MathVector<dim>* PLocIP = pgeo.scvf_local_ips();

		m_vvVShape.resize(pgeo.num_scvf_ips());
		for(size_t ip = 0; ip < m_vvVShape.size(); ++ip){
			m_vvVShape[ip].resize(rVTrialSpace.num_sh());
			for(size_t sh = 0; sh < rVTrialSpace.num_sh(); ++sh){
				m_vvVShape[ip][sh] = rVTrialSpace.shape(sh, PLocIP[ip]);
			}
		}

		const LocalShapeFunctionSet<dim>& rPTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, m_pLFEID);
		const MathVector<dim>* VLocIP = vgeo.scvf_local_ips();

		m_vvPShape.resize(vgeo.num_scvf_ips());
		for(size_t ip = 0; ip < m_vvPShape.size(); ++ip){
			m_vvPShape[ip].resize(rPTrialSpace.num_sh());
			for(size_t sh = 0; sh < rPTrialSpace.num_sh(); ++sh){
				m_vvPShape[ip][sh] = rPTrialSpace.shape(sh, VLocIP[ip]);
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
fsh_elem_loop()
{}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
	VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
	PGeom& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
	try{
		vgeo.update(elem, vCornerCoords, &(this->subset_handler()));
	    pgeo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("NavierStokes: Cannot update Finite Volume Geometry.");

//	set global positions for imports
	const MathVector<dim>* vSCVFip = vgeo.scvf_global_ips();
	const size_t numSCVFip = vgeo.num_scvf_ips();
	const MathVector<dim>* vSCVip = vgeo.scv_global_ips();
	const size_t numSCVip = vgeo.num_scv_ips();
	m_imKinViscosity.	set_global_ips(vSCVFip, numSCVFip);
	m_imDensitySCVF.	set_global_ips(vSCVFip, numSCVFip);
	m_imDensitySCV.		set_global_ips(vSCVip, numSCVip);
	m_imSource.			set_global_ips(vSCVip, numSCVip);
	m_imDensitySCVFp.	set_global_ips(pgeo.scvf_global_ips(), pgeo.num_scvf_ips());
}

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
	const PGeom& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);

	////////////////////////////////////////////////////
	////////////////////////////////////////////////////
	// Momentum Equation (conservation of momentum)
	////////////////////////////////////////////////////
	////////////////////////////////////////////////////

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t s = 0, ip = 0; s < vgeo.num_scvf(); ++s){
		const typename VGeom::SCVF& scvf = vgeo.scvf(s);
		for(size_t i = 0; i < scvf.num_ip(); ++i, ++ip){

			////////////////////////////////////////////////////
			// Diffusive Term (Momentum Equation)
			////////////////////////////////////////////////////
			const number scale = -1.0 * m_imKinViscosity[ip]
								 * m_imDensitySCVF[ip] * scvf.weight(i);

		// 	loop shape functions
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			{

			// 	Compute flux derivative at IP
				const number flux_sh =  scale * VecDot(scvf.global_grad(i, sh), scvf.normal());

			// 	Add flux derivative  to local matrix
				for(int d1 = 0; d1 < dim; ++d1)
				{
					J(d1, scvf.from(), d1, sh) += flux_sh;
					J(d1, scvf.to()  , d1, sh) -= flux_sh;
				}

				if(!m_bLaplace)
				{
					for(int d1 = 0; d1 < dim; ++d1)
						for(int d2 = 0; d2 < dim; ++d2)
						{
							const number flux2_sh = scale * scvf.global_grad(i, sh)[d1]
													* scvf.normal()[d2];
							J(d1, scvf.from(), d2, sh) += flux2_sh;
							J(d1, scvf.to()  , d2, sh) -= flux2_sh;
						}
				}
			}

			////////////////////////////////////////////////////
			// Convective Term (Momentum Equation)
			////////////////////////////////////////////////////

			if (!m_bStokes) {
				// 	Interpolate Velocity at ip
				MathVector<dim> Vel(0.0);
				for(int d = 0; d < dim; ++d)
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						Vel[d] += u(d, sh) * scvf.shape(i, sh);

				//	compute product of stabilized vel and normal
				const number prod = VecProd(Vel, scvf.normal())
									* m_imDensitySCVF[ip] * scvf.weight(i);

				// linearization of u^T n u w.r.t second u (i.e. keeping first as old iterate)
				for (int d1 = 0; d1 < dim; ++d1){
					for (size_t sh = 0; sh < scvf.num_sh(); ++sh){
						J(d1, scvf.from(), d1, sh) += prod * scvf.shape(i, sh);
						J(d1, scvf.to(), d1, sh)   -= prod * scvf.shape(i, sh);
					}
				}

				if(m_bFullNewtonFactor){
					// linearization of u^T n u w.r.t first u (i.e. keeping second as old iterate)
					for (int d1 = 0; d1 < dim; ++d1){
						for (size_t sh = 0; sh < scvf.num_sh(); ++sh){

							const number prod =  m_bFullNewtonFactor * scvf.shape(i, sh)
												 * m_imDensitySCVF[ip] * scvf.weight(i);
							for (int d2 = 0; d2 < dim; ++d2){
								J(d1, scvf.from(), d2, sh) += Vel[d1] * scvf.normal()[d2] * prod;
								J(d1, scvf.to(), d2, sh)   -= Vel[d1] * scvf.normal()[d2] * prod;
							}
						}
					}
				}
			}

			////////////////////////////////////////////////////
			// Pressure Term (Momentum Equation)
			////////////////////////////////////////////////////

		//	Add flux derivative for local matrix
			for(size_t sh = 0; sh < m_vvPShape[ip].size(); ++sh)
				for(int d1 = 0; d1 < dim; ++d1)
				{
					const number flux_sh = m_vvPShape[ip][sh]
					                       * scvf.normal()[d1]
	                                       * scvf.weight(i);
					J(d1, scvf.from(), _P_, sh) += flux_sh;
					J(d1, scvf.to()  , _P_, sh) -= flux_sh;
				}
		} // end ip
	} // end scvf


	////////////////////////////////////////////////////
	////////////////////////////////////////////////////
	// Continuity Equation (conservation of mass)
	////////////////////////////////////////////////////
	////////////////////////////////////////////////////

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t s = 0, ip = 0; s < pgeo.num_scvf(); ++s){
		const typename PGeom::SCVF& scvf = pgeo.scvf(s);
		for(size_t i = 0; i < scvf.num_ip(); ++i, ++ip){
			for(size_t sh = 0; sh < m_vvVShape[ip].size(); ++sh)
			{
			//	compute flux at ip
				const number contFlux = m_vvVShape[ip][sh]
										* m_imDensitySCVFp[ip]
										* scvf.weight(i);

			//	Add contributions to local defect
				for(int d1 = 0; d1 < dim; ++d1){
					J(_P_, scvf.from(), d1, sh) += contFlux * scvf.normal()[d1];
					J(_P_, scvf.to()  , d1, sh) -= contFlux * scvf.normal()[d1];
				}
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
	const PGeom& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);


	////////////////////////////////////////////////////
	////////////////////////////////////////////////////
	// Momentum Equation (conservation of momentum)
	////////////////////////////////////////////////////
	////////////////////////////////////////////////////

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t s = 0, ip = 0; s < vgeo.num_scvf(); ++s){
		const typename VGeom::SCVF& scvf = vgeo.scvf(s);
		for(size_t i = 0; i < scvf.num_ip(); ++i, ++ip){

			////////////////////////////////////////////////////
			// Diffusive Term (Momentum Equation)
			////////////////////////////////////////////////////

		// 	1. Interpolate Functional Matrix of velocity at ip
			MathMatrix<dim, dim> gradVel;
			for(int d1 = 0; d1 < dim; ++d1){
				for(int d2 = 0; d2 <dim; ++d2){
					gradVel(d1, d2) = 0.0;
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						gradVel(d1, d2) += u(d1, sh) * scvf.global_grad(i, sh)[d2];
				}
			}

		//	2. Compute flux
			MathVector<dim> diffFlux;

		//	Add (\nabla u) \cdot \vec{n}
			MatVecMult(diffFlux, gradVel, scvf.normal());

		//	Add (\nabla u)^T \cdot \vec{n}
			if(!m_bLaplace)
				TransposedMatVecMultAdd(diffFlux, gradVel, scvf.normal());

		//	scale by viscosity
			VecScale(diffFlux, diffFlux, (-1.0) * m_imKinViscosity[ip]
			                             * m_imDensitySCVF[ip] * scvf.weight(i));

		//	3. Add flux to local defect
			for(int d1 = 0; d1 < dim; ++d1)
			{
				d(d1, scvf.from()) += diffFlux[d1];
				d(d1, scvf.to()  ) -= diffFlux[d1];
			}

			////////////////////////////////////////////////////
			// Convective Term (Momentum Equation)
			////////////////////////////////////////////////////

			if (! m_bStokes) {
			//	compute ip velocity
				MathVector<dim> Vel(0.0);
				for(int d1 = 0; d1 < dim; ++d1)
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						Vel[d1] += u(d1, sh) * scvf.shape(i, sh);

			//	compute product of standard velocity and normal
				const number prod = VecProd(Vel, scvf.normal())
									* m_imDensitySCVF[ip] * scvf.weight(i);

			//	Add contributions to local velocity components
				for(int d1 = 0; d1 < dim; ++d1)
				{
					d(d1, scvf.from()) += Vel[d1] * prod;
					d(d1, scvf.to()  ) -= Vel[d1] * prod;
				}
			}

			////////////////////////////////////////////////////
			// Pressure Term (Momentum Equation)
			////////////////////////////////////////////////////

		//	1. Interpolate pressure at ip
			number pressure = 0.0;
			for(size_t sh = 0; sh < m_vvPShape[ip].size(); ++sh)
				pressure += m_vvPShape[ip][sh] * u(_P_, sh);
			pressure *= scvf.weight(i);

		//	2. Add contributions to local defect
			for(int d1 = 0; d1 < dim; ++d1)
			{
				d(d1, scvf.from()) += pressure * scvf.normal()[d1];
				d(d1, scvf.to()  ) -= pressure * scvf.normal()[d1];
			}
		}
	}

//	interpolate velocity at ip with standard lagrange interpolation
	std::vector<MathVector<dim> > PStdVel(m_vvVShape.size());
	for(size_t ip = 0; ip < m_vvVShape.size(); ++ip){
		VecSet(PStdVel[ip], 0.0);
		for(size_t sh = 0; sh < m_vvVShape[ip].size(); ++sh)
			for(int d1 = 0; d1 < dim; ++d1)
				PStdVel[ip][d1] += u(d1, sh) * m_vvVShape[ip][sh];
	}

	////////////////////////////////////////////////////
	////////////////////////////////////////////////////
	// Continuity Equation (conservation of mass)
	////////////////////////////////////////////////////
	////////////////////////////////////////////////////

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t s = 0, ip = 0; s < pgeo.num_scvf(); ++s){
		const typename PGeom::SCVF& scvf = pgeo.scvf(s);
		for(size_t i = 0; i < scvf.num_ip(); ++i, ++ip){

			//	compute flux at ip
			const number contFlux = VecProd(PStdVel[ip], scvf.normal())
									* m_imDensitySCVFp[ip] * scvf.weight(i);

			//	Add contributions to local defect
			d(_P_, scvf.from()) += contFlux;
			d(_P_, scvf.to()  ) -= contFlux;
		}
	}
}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);

// 	loop Sub Control Volumes (SCV)
	for(size_t s = 0, ip = 0; s < vgeo.num_scv(); ++s)
	{
	// 	get current SCV
		const typename VGeom::SCV& scv = vgeo.scv(s);

	// 	get associated node
		const int co = scv.node_id();

	//	loop shapes
		for(size_t sh = 0; sh < scv.num_sh(); ++sh)
		{
		//	reset integral
			number integral = 0;

		//	loop integration points
			for(size_t i = 0; i < scv.num_ip(); ++i)
			{
				integral += scv.shape(i, sh) * scv.weight(i)
							* m_imDensitySCV[ip+i];
			}

		// 	loop velocity components
			for(int d1 = 0; d1 < dim; ++d1)
			{
			// 	Add to local matrix
				J(d1, co, d1, sh) += integral;
			}
		}

	//	increase ip offset
		ip += scv.num_ip();
	}
}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);

// 	loop Sub Control Volumes (SCV)
	for(size_t s = 0, ip = 0; s < vgeo.num_scv(); ++s)
	{
	// 	get current SCV
		const typename VGeom::SCV& scv = vgeo.scv(s);

	// 	get associated node
		const int co = scv.node_id();

	//	loop integration points
		for(size_t i = 0; i < scv.num_ip(); ++i, ++ip)
		{
			MathVector<dim> Vel; VecSet(Vel, 0.0);
			for(size_t sh = 0; sh < scv.num_sh(); ++sh)
				for(int d1 = 0; d1 < dim; ++d1)
					Vel[d1] += u(d1, sh) * scv.shape(i, sh);

		// 	Add to local defect
			for(int d1 = 0; d1 < dim; ++d1)
				d(d1, co) +=  Vel[d1] * m_imDensitySCV[ip]
                              * scv.weight(i);
		}
	}
}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	if zero data given, return
	if(!m_imSource.data_given()) return;
	const VGeom& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t s = 0, ip = 0; s < vgeo.num_scv(); ++s){
		const typename VGeom::SCV& scv = vgeo.scv(s);
		for(size_t i = 0; i < scv.num_ip(); ++i, ++ip){
			const int co = scv.node_id();

			for(int d1 = 0; d1 < dim; ++d1)
			{
				d(d1, co) += m_imSource[ip][d1] * scv.weight(i);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void NavierStokesFV<Domain1d>::
register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID)
{
	UG_THROW("Not implemented.");
}
#endif

#ifdef UG_DIM_2
template<>
void NavierStokesFV<Domain2d>::
register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID)
{
	typedef DimFVGeometry<dim> FVGeom;
	register_func<Triangle, FVGeom, FVGeom >();
	register_func<Quadrilateral, FVGeom, FVGeom >();
}
#endif

#ifdef UG_DIM_3
template<>
void NavierStokesFV<Domain3d>::
register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID)
{
	typedef DimFVGeometry<dim> FVGeom;
	register_func<Tetrahedron, FVGeom, FVGeom >();
	register_func<Prism, FVGeom, FVGeom >();
	register_func<Hexahedron, FVGeom, FVGeom >();
}
#endif

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFV<TDomain>::register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, VGeom, PGeom>);
	this->set_prep_elem_fct(	 	id, &T::template prep_elem<TElem, VGeom, PGeom>);
	this->set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop<TElem, VGeom, PGeom>);
	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, VGeom, PGeom>);
	this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, VGeom, PGeom>);
	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, VGeom, PGeom>);
	this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, VGeom, PGeom>);
	this->set_add_rhs_elem_fct(	id, &T::template add_rhs_elem<TElem, VGeom, PGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_2
template class NavierStokesFV<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesFV<Domain3d>;
#endif

} // namespace NavierStokes
} // namespace ug
