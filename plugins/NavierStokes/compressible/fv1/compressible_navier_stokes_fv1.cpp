/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
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

/*
 *      main parts are copied from the discretization of the incompressible Navier-Stokes Equations
 *      of Andreas Vogel and Christian Wehner
 */

#include "compressible_navier_stokes_fv1.h"

#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
CompressibleNavierStokesFV1<TDomain>::CompressibleNavierStokesFV1(const char* functions,
                                          const char* subsets)
: CompressibleNavierStokesBase<TDomain>(functions, subsets)
{
	init();
};

template<typename TDomain>
CompressibleNavierStokesFV1<TDomain>::CompressibleNavierStokesFV1(const std::vector<std::string>& vFct,
                                          const std::vector<std::string>& vSubset)
: CompressibleNavierStokesBase<TDomain>(vFct, vSubset)
{
	init();
};


template<typename TDomain>
void CompressibleNavierStokesFV1<TDomain>::init()
{
//	check number of functions
	if(this->num_fct() != dim+2)
		UG_THROW("Wrong number of functions: The ElemDisc 'CompressibleNavierStokes'"
					   " needs exactly "<<dim+2<<" symbolic function.");

	m_maxVel = 0.0; m_maxPressure = 0.0; m_maxDensity = 0.0; m_refLength = 0.0;
	m_refMachNrSq = 0.0; m_numeratorOfRefReynoldsNr = 0.0;

//	register imports
	this->register_import(m_imSourceSCV);
	this->register_import(m_imSourceSCVF);
	this->register_import(m_imKinViscosity);
	this->register_import(m_imAdiabaticIndex);

	m_imSourceSCV.set_rhs_part();
	m_imSourceSCVF.set_rhs_part();

	//	update assemble functions
	register_all_funcs(false);
}

template<typename TDomain>
void CompressibleNavierStokesFV1<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	if(bNonRegularGrid)
		UG_THROW("CompressibleNavierStokes: only regular grid implemented.");

//	check number
	if(vLfeID.size() != dim+2)
		UG_THROW("CompressibleNavierStokes: Need exactly "<<dim+2<<" functions");

	for(int d = 0; d <= dim+1; ++d)
		if(vLfeID[d].type() != LFEID::LAGRANGE || vLfeID[d].order() != 1)
			UG_THROW("CompressibleNavierStokes: 'fv1' expects Lagrange P1 trial space "
					"for velocity,pressure and density.");

	//	update assemble functions
	register_all_funcs(false);
}

template<typename TDomain>
void CompressibleNavierStokesFV1<TDomain>::
set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > data)
{
	m_imKinViscosity.set_data(data);
}

template<typename TDomain>
void CompressibleNavierStokesFV1<TDomain>::
set_adiabatic_index(SmartPtr<CplUserData<number, dim> > data)
{
	m_imAdiabaticIndex.set_data(data);
}

template<typename TDomain>
void CompressibleNavierStokesFV1<TDomain>::
set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
{
	m_imSourceSCV.set_data(data);
	m_imSourceSCVF.set_data(data);
}



////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void CompressibleNavierStokesFV1<TDomain>::
prep_timestep_elem(const number time, const LocalVector& u, GridObject* elem,
		const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
/*	TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("CompressibleNavierStokes::prep_timestep_elem:"
						" Cannot update Finite Volume Geometry.");*/

	const size_t numVertices = static_cast<TElem*>(elem)->num_vertices();

	//	compute norm of velocity-vector, norm of density
	//	and norm of pressure at every element corner
	for (size_t co = 0; co < numVertices; ++co)
	{
		MathVector<dim> vel;
		for(int d1 = 0; d1 < dim; ++d1)
			vel[d1] = u(d1,co);

		number velNorm = VecTwoNorm(vel);
		number pressure = u(_P_,co);
		number density = u(_Rho_,co);

		//	update the maximal values
		if (velNorm > m_maxVel) 		m_maxVel = velNorm;
		if (pressure > m_maxPressure) 	m_maxPressure = pressure;
		if (density > m_maxDensity) 	m_maxDensity = density;
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void CompressibleNavierStokesFV1<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
// 	Only first order implementation
	if(!(TFVGeom::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
						" Geometry set.");

//	check, that convective upwinding has been set
	if(m_spConvUpwind.invalid())
		UG_THROW("Upwinding for convective Term in Momentum eq. not set.");

//	init convection stabilization for element type
	m_spConvUpwind->template set_geometry_type<TFVGeom >();

//	check, that kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("CompressibleNavierStokes::prep_elem_loop:"
						" Kinematic Viscosity has not been set, but is required.");

//	check, that adiabatic index has been set
	if(!m_imAdiabaticIndex.data_given())
		UG_THROW("CompressibleNavierStokes::prep_elem_loop:"
						" Adiabatic Index has not been set, but is required.");

//	set local positions for imports
	if(!TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;
		TFVGeom& geo = GeomProvider<TFVGeom>::get();
		const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
		const size_t numSCVFip = geo.num_scvf_ips();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();

		m_imKinViscosity.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imAdiabaticIndex.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imSourceSCV.template set_local_ips<refDim>(vSCVip,numSCVip);
		m_imSourceSCVF.template set_local_ips<refDim>(vSCVFip,numSCVFip);
	}

//	check, that a reference length and a reference density have been set
	if (m_maxDensity <= 0.0)
		UG_THROW("CompressibleNavierStokes::prep_elem_loop:"
						" maximal density is " <<m_maxDensity<< " and needs to be positive!");
	if (m_refLength <= 0.0)
		UG_THROW("CompressibleNavierStokes::prep_elem_loop:"
						" reference length is " <<m_refLength<< " and needs to be positive!");

//	set reference Mach-number squared and reference reynolds number
	m_refMachNrSq = m_maxDensity * m_maxVel * m_maxVel / m_maxPressure;
	m_numeratorOfRefReynoldsNr = m_maxVel * m_refLength * m_maxDensity;
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void CompressibleNavierStokesFV1<TDomain>::
fsh_elem_loop()
{}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void CompressibleNavierStokesFV1<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
	TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("CompressibleNavierStokes::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	set local positions for imports
	if(TFVGeom::usesHangingNodes)
	{
	//	request ip series
		static const int refDim = TElem::dim;
		const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
		const size_t numSCVFip = geo.num_scvf_ips();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();
		m_imKinViscosity.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imAdiabaticIndex.template set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imSourceSCV.template set_local_ips<refDim>(vSCVip,numSCVip);
		m_imSourceSCVF.template set_local_ips<refDim>(vSCVFip,numSCVFip);
	}

//	set global positions for imports
	const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
	const size_t numSCVFip = geo.num_scvf_ips();
	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	const size_t numSCVip = geo.num_scv_ips();
	m_imKinViscosity.set_global_ips(vSCVFip, numSCVFip);
	m_imAdiabaticIndex.set_global_ips(vSCVFip, numSCVFip);
	m_imSourceSCV.set_global_ips(vSCVip, numSCVip);
	m_imSourceSCVF.set_global_ips(vSCVFip, numSCVFip);
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void CompressibleNavierStokesFV1<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

	// 	TODO: NEEDS TO BE IMPLEMENTED!

/*
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	check for source term to pass to the stabilization
	const DataImport<MathVector<dim>, dim>* pSource = NULL;
	if(m_imSourceSCVF.data_given())	pSource = &m_imSourceSCVF;

//	check for solutions to pass to stabilization in time-dependent case
	const LocalVector *pSol = &u, *pOldSol = NULL;
	number dt = 0.0;
	if(this->is_time_dependent())
	{
	//	get and check current and old solution
		const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
		if(vLocSol->size() != 2)
			UG_THROW("NavierStokes::add_jac_A_elem: "
							" Stabilization needs exactly two time points.");

	//	remember local solutions
		pSol = &vLocSol->solution(0);
		pOldSol = &vLocSol->solution(1);
		dt = vLocSol->time(0) - vLocSol->time(1);
	}

//	interpolate velocity at ip with standard lagrange interpolation
	static const size_t numSCVF = TFVGeom::numSCVF;
	MathVector<dim> StdVel[numSCVF];
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		VecSet(StdVel[ip], 0.0);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			for(int d1 = 0; d1 < dim; ++d1)
				StdVel[ip][d1] += u(d1, sh) * scvf.shape(sh);
	}

//	compute stabilized velocities and shapes for continuity equation
	m_spStab->update(&geo, *pSol, StdVel, m_bStokes, m_imKinViscosity, pSource, pOldSol, dt);

	if (! m_bStokes) // no convective terms in the Stokes eq. => no upwinding
	{
	//	compute upwind shapes
		if(m_spConvUpwind.valid())
			m_spConvUpwind->update(&geo, StdVel);
	}

//	get a const (!!) reference to the stabilization
	const INavierStokesFV1Stabilization<dim>& stab = *m_spStab;
	const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// 	get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	// 	loop shape functions
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
		{
			////////////////////////////////////////////////////
			////////////////////////////////////////////////////
			// Momentum Equation (conservation of momentum)
			////////////////////////////////////////////////////
			////////////////////////////////////////////////////

			////////////////////////////////////////////////////
			// Diffusive Term (Momentum Equation)
			////////////////////////////////////////////////////

		// 	Compute flux derivative at IP
			const number flux_sh =  -1.0 * m_imKinViscosity[ip]
									* VecDot(scvf.global_grad(sh), scvf.normal());

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
						const number flux2_sh = -1.0 * m_imKinViscosity[ip] * m_imDensitySCVF[ip]
												* scvf.global_grad(sh)[d1] * scvf.normal()[d2];
						J(d1, scvf.from(), d2, sh) += flux2_sh;
						J(d1, scvf.to()  , d2, sh) -= flux2_sh;
					}
			}

			////////////////////////////////////////////////////
			// Pressure Term (Momentum Equation)
			////////////////////////////////////////////////////

		//	Add flux derivative for local matrix
			for(int d1 = 0; d1 < dim; ++d1)
			{
				const number flux_sh = scvf.shape(sh) * scvf.normal()[d1];
				J(d1, scvf.from(), _P_, sh) += flux_sh;
				J(d1, scvf.to()  , _P_, sh) -= flux_sh;
			}

			////////////////////////////////////////////////////
			// Convective Term (Momentum Equation)
			////////////////////////////////////////////////////

			// note: StdVel will be the advecting velocity (not upwinded)
			//       UpwindVel/StabVel will be the transported velocity

			if (! m_bStokes) // no convective terms in the Stokes equation
			{
			//	compute upwind velocity
				if(m_spConvUpwind.invalid())
					UG_THROW("Cannot find upwind for convective term.");
				MathVector<dim> UpwindVel;
				UpwindVel = upwind.upwind_vel(ip, u, StdVel);

			//	Mach-number blend
				number w = 1.0;
				if(m_bMachNrBlend)
					w = mach_number_blending(UpwindVel, geo, ip, StdVel[ip], m_imKinViscosity[ip]);

			//	compute product of stabilized vel and normal
				const number prod = VecProd(StdVel[ip], scvf.normal()) * m_imDensitySCVF[ip];

			///////////////////////////////////
			//	Add fixpoint linearization
			///////////////////////////////////

			//	Stabilization used as upwind

			//	Upwind used as upwind
				if(m_spConvUpwind.valid())
				{
					number convFlux_vel = upwind.upwind_shape_sh(ip, sh);

				//	in some cases (e.g. PositiveUpwind, RegularUpwind) the upwind
				//	velocity in an ip depends also on the upwind velocity in
				//	other ips. This is reflected by the fact, that the ip
				//	shapes are non-zero. In that case, we can interpolate an
				//	approximate upwind only from the corner velocities by using
				//	u_up = \sum shape_co U_co + \sum shape_ip \tilde{u}_ip
				//	     = \sum shape_co U_co + \sum \sum shape_ip norm_shape_co|_ip * U_co
					if(m_spConvUpwind->non_zero_shape_ip())
					{
						for(size_t ip2 = 0; ip2 < geo.num_scvf(); ++ip2)
						{
							const typename TFVGeom::SCVF& scvf2 = geo.scvf(ip2);
							convFlux_vel += scvf2.shape(sh) * upwind.upwind_shape_ip(ip, ip2);
						}
					}

					convFlux_vel *= prod * w;

					for(int d1 = 0; d1 < dim; ++d1)
					{
						J(d1, scvf.from(), d1, sh) += convFlux_vel;
						J(d1, scvf.to()  , d1, sh) -= convFlux_vel;
					}
				}

			//	derivative due to Mach-number blending
				if(m_bMachNrBlend)
				{
					const number convFluxPe = prod * (1.0-w) * scvf.shape(sh);
					for(int d1 = 0; d1 < dim; ++d1)
					{
						J(d1, scvf.from(), d1, sh) += convFluxPe;
						J(d1, scvf.to()  , d1, sh) -= convFluxPe;
					}
				}

			/////////////////////////////////////////
			//	Add full jacobian (remaining part)
			/////////////////////////////////////////

			//	Add remaining term for exact jacobian
				if(m_bFullNewtonFactor)
				{
				//	Stabilization used as upwind

				//	Upwind used as upwind
					if(m_spConvUpwind.valid())
					{
					//	loop defect components
						for(int d1 = 0; d1 < dim; ++d1)
							for(int d2 = 0; d2 < dim; ++d2)
							{
							//	derivatives w.r.t. velocity
								number prod_vel = w * upwind.upwind_shape_sh(ip,sh)
													* scvf.normal()[d2] * m_imDensitySCVF[ip];

								J(d1, scvf.from(), d2, sh) += prod_vel * UpwindVel[d1];
								J(d1, scvf.to()  , d2, sh) -= prod_vel * UpwindVel[d1];
							}
					}

				//	derivative due to Mach-number blending
					if(m_bMachNrBlend)
					{
						for(int d1 = 0; d1 < dim; ++d1)
							for(int d2 = 0; d2 < dim; ++d2)
							{
								const number convFluxPe = UpwindVel[d1] * (1.0-w)
														  * scvf.shape(sh)
														  * scvf.normal()[d2]
														  * m_imDensitySCVF[ip];
								J(d1, scvf.from(), d2, sh) += convFluxPe;
								J(d1, scvf.to()  , d2, sh) -= convFluxPe;
							}
					}
				} // end exact jacobian part

			} // end of if (! m_bStokes) for the convective terms

			////////////////////////////////////////////////////
			////////////////////////////////////////////////////
			// Continuity Equation (conservation of mass)
			////////////////////////////////////////////////////
			////////////////////////////////////////////////////

		//	Add derivative of stabilized flux w.r.t velocity comp to local matrix
			if(stab.vel_comp_connected())
			{
				for(int d1 = 0; d1 < dim; ++d1)
				{
					number contFlux_vel = 0.0;
					for(int d2 = 0; d2 < dim; ++d2)
						contFlux_vel += stab.stab_shape_vel(ip, d2, d1, sh)
										* scvf.normal()[d2] * m_imDensitySCVF[ip];

					J(_P_, scvf.from(), d1, sh) += contFlux_vel;
					J(_P_, scvf.to()  , d1, sh) -= contFlux_vel;
				}
			}
			else
			{
				for(int d1 = 0; d1 < dim; ++d1)
				{
					const number contFlux_vel = stab.stab_shape_vel(ip, d1, d1, sh)
													* scvf.normal()[d1] * m_imDensitySCVF[ip];

					J(_P_, scvf.from(), d1, sh) += contFlux_vel;
					J(_P_, scvf.to()  , d1, sh) -= contFlux_vel;
				}
			}

		//	Add derivative of stabilized flux w.r.t pressure to local matrix
			number contFlux_p = 0.0;
			for(int d1 = 0; d1 < dim; ++d1)
				contFlux_p += stab.stab_shape_p(ip, d1, sh) * scvf.normal()[d1] * m_imDensitySCVF[ip];

			J(_P_, scvf.from(), _P_, sh) += contFlux_p;
			J(_P_, scvf.to()  , _P_, sh) -= contFlux_p;
		}
	}*/
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void CompressibleNavierStokesFV1<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implemented
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	check for source term to pass to the stabilization
	//const DataImport<MathVector<dim>, dim>* pSource = NULL;
	//if(m_imSourceSCVF.data_given())	pSource = &m_imSourceSCVF;

//	interpolate velocity at ip with standard lagrange interpolation
	static const size_t numSCVF = TFVGeom::numSCVF;
	MathVector<dim> StdVel[numSCVF];
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		VecSet(StdVel[ip], 0.0);

		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			for(int d1 = 0; d1 < dim; ++d1)
				StdVel[ip][d1] += u(d1, sh) * scvf.shape(sh);
	}

//	compute upwind shapes
	if(m_spConvUpwind.valid())
		m_spConvUpwind->update(&geo, StdVel);

//	get a const (!!) reference to the upwind scheme
	const INavierStokesUpwind<dim>& upwind = *m_spConvUpwind;

// 	loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// 	get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// Momentum Equation (conservation of momentum)
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////

	//	\todo: we could also add all fluxes at once in order to save several
	//			accesses to the local defect, implement and loose clarity?

		////////////////////////////////////////////////////
		// Diffusive Term (Momentum Equation)
		////////////////////////////////////////////////////

	// 	1. a) Interpolate Functional Matrix of velocity at ip
		MathMatrix<dim, dim> gradVel;
		for(int d1 = 0; d1 < dim; ++d1)
			for(int d2 = 0; d2 <dim; ++d2)
			{
			//	sum up contributions of each shape
				gradVel(d1, d2) = 0.0;
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					gradVel(d1, d2) += scvf.global_grad(sh)[d2]
					                    * u(d1, sh);
			}

	//	b) Interpolate divergence of velocity at ip
		number divVel = 0.0;
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			for(int d1 = 0; d1 < dim; ++d1)
				divVel += scvf.global_grad(sh)[d1] * u(d1, sh);


	//	2. Compute flux = (1.0/Re) \cdot \tau \cdot \vec{n}
		MathVector<dim> diffFlux;

	//	Add (\nabla u) \cdot \vec{n}
		MatVecMult(diffFlux, gradVel, scvf.normal());

	//	Add (\nabla u)^T \cdot \vec{n}
		// TODO: m_bLaplace applicable here?!
		if(true) //!m_bLaplace)
			TransposedMatVecMultAdd(diffFlux, gradVel, scvf.normal());

	//	scale by viscosity and subtract volumetric part of the shear stress tensor
		// diffFlux = 	- kinVisco * rho * ((\nabla u) + (\nabla u)^T) \cdot \vec{n}
		//				+ 2.0/3.0 \delta_{ij} * kinVisco * rho * div(u)

		// TODO: use the right viscosity here! the density has its own equation now!!!
		//	vorher: VecScale(diffFlux, diffFlux, (-1.0) * m_imKinViscosity[ip] * m_imDensitySCVF[ip]);
		VecScaleAdd(diffFlux, (-1.0) * m_imKinViscosity[ip], diffFlux,
				(2.0/3.0) * m_imKinViscosity[ip] * divVel, scvf.normal());

	//	diffusive coefficient: 1.0 / Re
	//	(Re: Reynolds number = vel_{ref} * L_{ref} * density_{ref} / viscosity)
		number diffCoeff = m_imKinViscosity[ip] / m_numeratorOfRefReynoldsNr;

	//	3. Add flux to local defect
		for(int d1 = 0; d1 < dim; ++d1)
		{
			d(d1, scvf.from()) += diffCoeff * diffFlux[d1];
			d(d1, scvf.to()  ) -= diffCoeff * diffFlux[d1];
		}

		////////////////////////////////////////////////////
		// Convective Term (Momentum Equation)
		////////////////////////////////////////////////////

		//	interpolate density at ip
		number density = 0.0;
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			density += scvf.shape(sh) * u(_Rho_, sh);

		// note: StdVel will be the advecting velocity (not upwinded)
		//       UpwindVel will be the transported velocity

		//	find the upwind velocity at ip
		if(m_spConvUpwind.invalid())
			UG_THROW("Cannot find upwind for convective term.");

		MathVector<dim> UpwindVel;
		UpwindVel = upwind.upwind_vel(ip, u, StdVel);

	//	Mach-number Blending
		if(m_bMachNrBlend)
			mach_number_blending(UpwindVel, geo, ip, StdVel[ip], m_imKinViscosity[ip]);

	//	compute product of standard velocity and normal
		const number prod = VecProd(StdVel[ip], scvf.normal()) * density;

	//	Add contributions to local velocity components
		for(int d1 = 0; d1 < dim; ++d1)
		{
			d(d1, scvf.from()) += UpwindVel[d1] * prod;
			d(d1, scvf.to()  ) -= UpwindVel[d1] * prod;
		}

		////////////////////////////////////////////////////
		// Pressure Term (Momentum Equation)
		////////////////////////////////////////////////////

	//	1. Interpolate pressure at ip
		number pressure = 0.0;
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			pressure += scvf.shape(sh) * u(_P_, sh);
		const number pressureFac = 1.0 / (m_imAdiabaticIndex[ip] * m_refMachNrSq);

	//	2. Add contributions to local defect
		for(int d1 = 0; d1 < dim; ++d1)
		{
			d(d1, scvf.from()) += pressureFac * pressure * scvf.normal()[d1];
			d(d1, scvf.to()  ) -= pressureFac * pressure * scvf.normal()[d1];
		}

		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// Continuity Equation (conservation of mass)
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////

	//	compute flux at ip
		const number contFlux = VecProd(StdVel[ip], scvf.normal()) * density;

	//	Add contributions to local defect
		d(_Rho_, scvf.from()) += contFlux;
		d(_Rho_, scvf.to()  ) -= contFlux;

		////////////////////////////////////////////////////
		////////////////////////////////////////////////////
		// Energy Equation
		// (conservation of energy for an ideal polytropic gas)
		////////////////////////////////////////////////////
		////////////////////////////////////////////////////

	//	compute pressure flux at ip
		//	TODO: use correct velocity here!
		const number pressFlux = VecProd(StdVel[ip], scvf.normal()) * pressure;

	//	compute 'density-velocity' flux at ip
		const number densityVelFlux = 0.5 * VecProd(StdVel[ip], StdVel[ip])
				* VecProd(StdVel[ip], scvf.normal()) * density;
		const number densityVelFac = (1.0 - m_imAdiabaticIndex[ip]) * m_refMachNrSq;

	//	compute shear-stress flux at ip
		MathMatrix<dim, dim> shearStressTensor;
		for(int d1 = 0; d1 < dim; ++d1){
			for(int d2 = 0; d2 < dim; ++d2){
				shearStressTensor(d1, d2) = m_imKinViscosity[ip] * (gradVel(d1, d2) + gradVel(d2, d1));
			}
			shearStressTensor(d1, d1) -= (2.0/3.0) * m_imKinViscosity[ip] * divVel;
		}

		MathVector<dim> shearStressVelProd;
		MatVecMult(shearStressVelProd, shearStressTensor, StdVel[ip]);

		const number shearStressFlux = VecProd(shearStressVelProd, scvf.normal());
		const number shearStressFac = densityVelFac * diffCoeff;

	//	Add contributions to local defect
		d(_P_, scvf.from()) += pressFlux + densityVelFac * densityVelFlux
				+ shearStressFac * shearStressFlux;
		d(_P_, scvf.to()  ) -= pressFlux + densityVelFac * densityVelFlux
				+ shearStressFac * shearStressFlux;

	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void CompressibleNavierStokesFV1<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

	// 	TODO: NEEDS TO BE IMPLEMENTED!

/*
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int sh = scv.node_id();

	// 	loop velocity components
		for(int d1 = 0; d1 < dim; ++d1)
		{
		// 	Add to local matrix
			J(d1, sh, d1, sh) += scv.volume() * m_imDensitySCV[ip];
		}
	}*/
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void CompressibleNavierStokesFV1<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

	number massMatCoeff = 1.0 / sqrt(m_refMachNrSq);

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int sh = scv.node_id();

	//	contribution to the mass-term of the continuity equation
		d(_Rho_, sh) += massMatCoeff * u(_Rho_, sh) * scv.volume();

	//	contribution to the mass-term of the momentum equation
		MathVector<dim> vel;
		// 	loop velocity components
		for(int d1 = 0; d1 < dim; ++d1)
		{
		// 	Add to local matrix
			d(d1, sh) += massMatCoeff * u(_Rho_, sh) * u(d1, sh) * scv.volume();
			vel[d1] = u(d1, sh);
		}

	//	contribution to the mass-term of the energy-conserving equation
		d(_P_, sh) += massMatCoeff * scv.volume() * ( u(_P_, sh) / m_imAdiabaticIndex[ip]
		              + m_refMachNrSq * (m_imAdiabaticIndex[ip] - 1.0)
		              * 0.5 * u(_P_, sh) * VecProd(vel, vel) );
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void CompressibleNavierStokesFV1<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// 	TODO: NEEDS TO BE IMPLEMENTED!


// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

//	if zero data given, return
	if(!m_imSourceSCV.data_given()) return;

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int sh = scv.node_id();

	// 	Add to local rhs
		for(int d1 = 0; d1 < dim; ++d1)
			d(d1, sh) += m_imSourceSCV[ip][d1] * scv.volume();
	}
}

template<typename TDomain>
template<typename TFVGeom>
inline
number
CompressibleNavierStokesFV1<TDomain>::
mach_number_blending(MathVector<dim>& UpwindVel, const TFVGeom& geo, size_t ip,
             const MathVector<dim>& StdVel, number kinVisco)
{
	// TODO: distinguish between the treatment of the different equations (s. Gordner p. 80)!!!

	const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
//	compute peclet number
	number Pe = VecProd(StdVel, scvf.normal())/VecTwoNormSq(scvf.normal())
	 * VecDistance(geo.corners() [scvf.to()], geo.corners() [scvf.from()]) / kinVisco;

//	compute weight
	const number Pe2 = Pe * Pe;
	const number w = Pe2 / (5.0 + Pe2);

//	compute upwind vel
	VecScaleAdd(UpwindVel, w, UpwindVel, (1.0-w), StdVel);

	return w;
}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void CompressibleNavierStokesFV1<Domain1d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
	}
	else
	{
		UG_THROW("CompressibleNavierStokesFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_2
template<>
void CompressibleNavierStokesFV1<Domain2d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
	}
	else
	{
		UG_THROW("CompressibleNavierStokesFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_3
template<>
void CompressibleNavierStokesFV1<Domain3d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
		register_func<Prism, FV1Geometry<Prism, dim> >();
		register_func<Pyramid, FV1Geometry<Pyramid, dim> >();
		register_func<Hexahedron, FV1Geometry<Hexahedron, dim> >();
	}
	else
	{
		UG_THROW("CompressibleNavierStokesFV1: Hanging Nodes not implemented.")
	}
}
#endif

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void
CompressibleNavierStokesFV1<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->clear_add_fct(id);
	this->set_prep_timestep_elem_fct(id, &T::template prep_timestep_elem<TElem, TFVGeom>);
	this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 	id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(	id, &T::template add_rhs_elem<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_2
template class CompressibleNavierStokesFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class CompressibleNavierStokesFV1<Domain3d>;
#endif

} // namespace NavierStokes
} // namespace ug
