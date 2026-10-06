/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
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

#include "navier_stokes_fe.h"

#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
NavierStokesFE<TDomain>::NavierStokesFE(const char* functions,
                                          const char* subsets)
: IncompressibleNavierStokesBase<TDomain>(functions, subsets)
{
	init();
};

template<typename TDomain>
NavierStokesFE<TDomain>::NavierStokesFE(const std::vector<std::string>& vFct,
                                          const std::vector<std::string>& vSubset)
: IncompressibleNavierStokesBase<TDomain>(vFct, vSubset)
{
	init();
};


template<typename TDomain>
void NavierStokesFE<TDomain>::init()
{
//	check number of functions
	if(this->num_fct() != dim+1)
		UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");

	m_stabParam = 0.0;
	m_bQuadOrderUserDef = false;

//	register imports
	this->register_import(m_imSource);
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensity);

	m_imSource.set_rhs_part();
	m_imDensity.set_mass_part();

	//	default value for density
	base_type::set_density(1.0);

	// use only non-virtual assembling functions
	this->clear_add_fct();
}

template<typename TDomain>
void NavierStokesFE<TDomain>::set_quad_order(size_t order)
{
	m_quadOrder = order;
	m_bQuadOrderUserDef = true;
}


template<typename TDomain>
void NavierStokesFE<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	if(bNonRegularGrid)
		UG_THROW("NavierStokes: only implemented for regular grids.");

//	check number
	if(vLfeID.size() != dim+1)
		UG_THROW("NavierStokes: Needs exactly "<<dim+1<<" functions.");

	for(int d = 1; d < dim; ++d)
		if(vLfeID[0] != vLfeID[d])
			UG_THROW("NavierStokes: trial spaces for velocity expected to be"
					" identical for all velocity components.");

//	remember lfeID;
	m_vLFEID = vLfeID[0];
	m_pLFEID = vLfeID[dim];

	if(!m_bQuadOrderUserDef) m_quadOrder = 2*m_vLFEID.order()+1;

	//	update assemble functions
	register_all_funcs(m_vLFEID, m_pLFEID, m_quadOrder);
}

template<typename TDomain>
void NavierStokesFE<TDomain>::
set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > data)
{
	m_imKinViscosity.set_data(data);
}

template<typename TDomain>
void NavierStokesFE<TDomain>::
set_density(SmartPtr<CplUserData<number, dim> > data)
{
	m_imDensity.set_data(data);
}

template<typename TDomain>
void NavierStokesFE<TDomain>::
set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
{
	m_imSource.set_data(data);
}

////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
//	check, that kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("NavierStokes: Kinematic Viscosity has not been set, but is required.");

//	check, that Density has been set
	if(!m_imDensity.data_given())
		UG_THROW("NavierStokes: Density has not been set, but is required.");

	DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
	DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
	try{
		vgeo.update_local(roid, m_vLFEID, m_quadOrder);
		pgeo.update_local(roid, m_pLFEID, m_quadOrder);
	}
	UG_CATCH_THROW("NavierStokes: Cannot update Finite Element Geometry.");

//	set local positions for imports
	m_imKinViscosity.set_local_ips(vgeo.local_ips(), vgeo.num_ip());
	m_imDensity. set_local_ips(vgeo.local_ips(), vgeo.num_ip());
	m_imSource. set_local_ips(vgeo.local_ips(), vgeo.num_ip());
}

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
fsh_elem_loop()
{}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
	m_pElem = elem;

// 	Update Geometry for this element
	DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
	DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
	try{
		vgeo.update(elem, vCornerCoords, m_vLFEID, m_quadOrder);
	    pgeo.update(elem, vCornerCoords, m_pLFEID, m_quadOrder);
	}
	UG_CATCH_THROW("NavierStokes: Cannot update Finite Element Geometry.");

//	set global positions for imports
	m_imKinViscosity.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
	m_imDensity.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
	m_imSource.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
}

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
	const DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);

	for (size_t ip = 0; ip < vgeo.num_ip(); ++ip){

		////////////////////////////////////////////////////
		// Diffusive Term (Momentum Equation)
		////////////////////////////////////////////////////

		const number scale = m_imKinViscosity[ip] * m_imDensity[ip] * vgeo.weight(ip);
		for (int vdim = 0; vdim < dim; ++vdim){
			for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){

				for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
					for (int udim = 0; udim < dim; ++udim) {

						J(vdim, vsh, vdim, ush) +=  scale
											       * vgeo.global_grad(ip, ush)[udim]
											       * vgeo.global_grad(ip, vsh)[udim];

						if(!m_bLaplace) {
							J(vdim, vsh, udim, ush) +=  scale
												       * vgeo.global_grad(ip, ush)[udim]
												       * vgeo.global_grad(ip, vsh)[udim];
						}
					}
				}
			}
		}

		////////////////////////////////////////////////////
		// Convective Term (Momentum Equation)
		////////////////////////////////////////////////////

		if (!m_bStokes) {
			// 	Interpolate Velocity at ip
			MathVector<dim> Vel;
			for(int d = 0; d < dim; ++d){
				Vel[d] = 0.0;
				for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
					Vel[d] += u(d, sh) * vgeo.shape(ip, sh);
			}

			// linearization of u \nabla u w.r.t second u (i.e. keeping first as old iterate)
			for (int vdim = 0; vdim < dim; ++vdim){
				for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){

					for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
						for (int udim = 0; udim < dim; ++udim) {

							J(vdim, vsh, vdim, ush) += m_imDensity[ip]
							                           * Vel[udim]
													   * vgeo.global_grad(ip, ush)[udim]
							                   		   * vgeo.shape(ip, vsh)
							                   		   * vgeo.weight(ip);
						}
					}
				}
			}

			if(m_bFullNewtonFactor){

				// 	Interpolate Functional Matrix of velocity at ip
				MathMatrix<dim, dim> gradVel;
				for(int d1 = 0; d1 < dim; ++d1){
					for(int d2 = 0; d2 <dim; ++d2){
						gradVel(d1, d2) = 0.0;
						for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
							gradVel(d1, d2) += u(d1, sh) * vgeo.global_grad(ip, sh)[d2];
					}
				}

				// linearization of u \nabla u w.r.t first u (i.e. keeping second as old iterate)
				for (int vdim = 0; vdim < dim; ++vdim){
					for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){

						for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
							for (int udim = 0; udim < dim; ++udim) {

								J(vdim, vsh, udim, ush) += m_bFullNewtonFactor * m_imDensity[ip]
								                           * gradVel(vdim, udim)
								                   		   * vgeo.shape(ip, ush)
														   * vgeo.shape(ip, vsh)
														   * vgeo.weight(ip);
							}
						}
					}
				}
			}
		}

		////////////////////////////////////////////////////
		// Pressure Term (Momentum Equation)
		////////////////////////////////////////////////////

		for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
			for (int vdim = 0; vdim < dim; ++vdim){
				for (size_t psh = 0; psh < pgeo.num_sh(); ++psh)
					J(vdim, vsh, _P_, psh) -= pgeo.shape(ip, psh)
								   * vgeo.global_grad(ip, vsh)[vdim]
								   * vgeo.weight(ip);
			}
		}

		////////////////////////////////////////////////////
		// Continuity Equation (conservation of mass)
		////////////////////////////////////////////////////

		for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
			for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
				for (int udim = 0; udim < dim; ++udim) {
						J(_P_, psh, udim, ush) +=
										m_imDensity[ip]
								       * vgeo.global_grad(ip, ush)[udim]
						               * pgeo.shape(ip, psh)
						               * vgeo.weight(ip);
				}
			}
		}
	}

	// stabilization
	if(m_stabParam != 0)
	{
		const number scale = m_stabParam * ElementDiameterSq<GridObject, TDomain>(*m_pElem, *this->domain());
		for (size_t ip = 0; ip < pgeo.num_ip(); ++ip){
			for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
				for (size_t psh2 = 0; psh2 < pgeo.num_sh(); ++psh2){
					J(_P_, psh,_P_, psh2) += scale
								   * VecDot(pgeo.global_grad(ip, psh), pgeo.global_grad(ip, psh2))
								   * pgeo.weight(ip);
				}
			}
		}
	}

}

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
	const DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);

	// loop integration points, note: pgeo and vgeo have same ip
	for (size_t ip = 0; ip < vgeo.num_ip(); ++ip){

		// 	Interpolate Functional Matrix of velocity at ip
		MathMatrix<dim, dim> gradVel;
		for(int d1 = 0; d1 < dim; ++d1){
			for(int d2 = 0; d2 <dim; ++d2){
				gradVel(d1, d2) = 0.0;
				for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
					gradVel(d1, d2) += u(d1, sh) * vgeo.global_grad(ip, sh)[d2];
			}
		}

		number divu = 0.0;
		for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
			for (int udim = 0; udim < dim; ++udim) {
				divu += u(udim, ush) * vgeo.global_grad(ip, ush)[udim];
			}
		}

		////////////////////////////////////////////////////
		// Diffusive Term (Momentum Equation)
		////////////////////////////////////////////////////

		const number scale = m_imKinViscosity[ip] * m_imDensity[ip] * vgeo.weight(ip);
		for (int vdim = 0; vdim < dim; ++vdim){
			for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
				for (int udim = 0; udim < dim; ++udim) {

					d(vdim, vsh) +=  scale * gradVel(vdim, udim)
								   * vgeo.global_grad(ip, vsh)[udim];

					if(!m_bLaplace) {
						d(vdim, vsh) +=  scale * gradVel(udim, vdim)
									   * vgeo.global_grad(ip, vsh)[udim];
					}
				}
			}
		}

		////////////////////////////////////////////////////
		// Convective Term (Momentum Equation)
		////////////////////////////////////////////////////

		if (!m_bStokes) {
			// 	Interpolate Velocity at ip
			MathVector<dim> Vel;
			for(int d1 = 0; d1 < dim; ++d1){
				Vel[d1] = 0.0;
				for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
					Vel[d1] += u(d1, sh) * vgeo.shape(ip, sh);
			}

			MathVector<dim> convFlux;
			MatVecMult(convFlux, gradVel, Vel);

			for (int vdim = 0; vdim < dim; ++vdim){
				for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
					d(vdim, vsh) +=  m_imDensity[ip] * convFlux[vdim]
								   * vgeo.shape(ip, vsh) * vgeo.weight(ip);
				}
			}

			// \todo: if density == constant, this will always be zero
			//		  therefore it is excluded
			if(false){
				for (int vdim = 0; vdim < dim; ++vdim){
					for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
						d(vdim, vsh) +=  m_imDensity[ip] * Vel[vdim] * divu
									   * vgeo.shape(ip, vsh) * vgeo.weight(ip);
					}
				}
			}
		}

		////////////////////////////////////////////////////
		// Pressure Term (Momentum Equation)
		////////////////////////////////////////////////////

		number pressure = 0.0;
		for (size_t psh = 0; psh < pgeo.num_sh(); ++psh)
			pressure += u(_P_, psh) * pgeo.shape(ip, psh);

		for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
			for (int vdim = 0; vdim < dim; ++vdim){
				d(vdim, vsh) -= pressure
							   * vgeo.global_grad(ip, vsh)[vdim]
							   * vgeo.weight(ip);
			}
		}

		////////////////////////////////////////////////////
		// Continuity Equation (conservation of mass)
		////////////////////////////////////////////////////

		for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
			d(_P_, psh) += divu * m_imDensity[ip] * pgeo.shape(ip, psh)* vgeo.weight(ip);
		}
	}

	// stabilization
	if(m_stabParam != 0)
	{
		const number scale = m_stabParam * ElementDiameterSq<GridObject, TDomain>(*m_pElem, *this->domain());
		for (size_t ip = 0; ip < pgeo.num_ip(); ++ip){

			MathVector<dim> pressGrad; VecSet(pressGrad, 0.0);
			for (size_t psh = 0; psh < pgeo.num_sh(); ++psh)
				VecScaleAppend(pressGrad, u(_P_, psh), pgeo.global_grad(ip, psh));

			for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
				d(_P_, psh) += scale
							   * VecDot(pgeo.global_grad(ip, psh), pressGrad)
							   * pgeo.weight(ip);
			}
		}

		// \todo: only 2nd order accurate. For higher order the laplace of the
		//			velocity test space is needed.
	}
}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	request geometry
	const DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);

	//	loop integration points
	for(size_t ip = 0; ip < vgeo.num_ip(); ++ip){
		for(size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
			for(size_t ush= 0; ush < vgeo.num_sh(); ++ush){

				number val = m_imDensity[ip]
				             * vgeo.shape(ip, vsh) *vgeo.shape(ip, ush)
				             * vgeo.weight(ip);

				for (int vdim = 0; vdim < dim; ++vdim)
					J(vdim, vsh, vdim, ush) += val;
			}
		}
	}
}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	request geometry
	const DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);

	//	loop integration points
	for(size_t ip = 0; ip < vgeo.num_ip(); ++ip){
		for (int vdim = 0; vdim < dim; ++vdim){

			//	compute value of current solution at ip
			number u_dim_ip = 0.0;
			for(size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh)
				u_dim_ip += u(vdim,vsh) * vgeo.shape(ip, vsh);

			//	add density
			u_dim_ip *= m_imDensity[ip];

			//	loop test spaces
			for(size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
				d(vdim, vsh) +=  u_dim_ip * vgeo.shape(ip, vsh) * vgeo.weight(ip);
			}
		}
	}
}


template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	//	if zero data given, return
	if(!m_imSource.data_given()) return;

	const DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);

	for(size_t ip = 0; ip < vgeo.num_ip(); ++ip){
		for(size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
			for (int vdim = 0; vdim < dim; ++vdim){
				d(vdim, vsh) += m_imSource[ip][vdim] * vgeo.shape(ip, vsh) * vgeo.weight(ip) * m_imDensity[ip];
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void NavierStokesFE<Domain1d>::
register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder)
{
	UG_THROW("Not implemented.");
}
#endif

#ifdef UG_DIM_2
template<>
void NavierStokesFE<Domain2d>::
register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder)
{
	typedef DimFEGeometry<dim> FVGeom;
	register_func<Triangle, FVGeom, FVGeom >();
	register_func<Quadrilateral, FVGeom, FVGeom >();
}
#endif

#ifdef UG_DIM_3
template<>
void NavierStokesFE<Domain3d>::
register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder)
{
	typedef DimFEGeometry<dim> FVGeom;
	register_func<Tetrahedron, FVGeom, FVGeom >();
	register_func<Prism, FVGeom, FVGeom >();
	register_func<Hexahedron, FVGeom, FVGeom >();
}
#endif

template<typename TDomain>
template<typename TElem, typename VGeom, typename PGeom>
void NavierStokesFE<TDomain>::register_func()
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
template class NavierStokesFE<Domain2d>;
#endif
#ifdef UG_DIM_3
template class NavierStokesFE<Domain3d>;
#endif

} // namespace NavierStokes
} // namespace ug
