/*
 * Copyright (c) 2012-2014:  G-CSC, Goethe University Frankfurt
 * Author: Christian Wehner
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

#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

#include "symmetric_boundary_fvcr.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace NavierStokes{

template<typename TDomain>
void CRNavierStokesSymBC<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	if(bNonRegularGrid)
		UG_THROW("CRNavierStokesSymBC: only regular grid implemented.");

//	check number
	if(vLfeID.size() != dim+1)
		UG_THROW("CRNavierStokesSymBC: Needs exactly "<<dim+1<<" functions.");
}

/**
 * converts the subset names where the BC is imposed to the corresponding subset
 * indices (i.e. m_vScheduledBndSubSets -> m_vBndSubSetIndex):
 */
template<typename TDomain>
void
CRNavierStokesSymBC<TDomain>::
extract_scheduled_data()
{
//	clear all extracted data
	m_vBndSubSetIndex.clear();

//	loop all scheduled subsets
	for(size_t i = 0; i < m_vScheduledBndSubSets.size(); ++i)
	{
	//	create Subset Group
		SubsetGroup subsetGroup;

	//	convert strings
		try{
			subsetGroup = this->approx_space()->subset_grp_by_name(m_vScheduledBndSubSets[i].c_str());
		}UG_CATCH_THROW("'CRNavierStokesSymBC:extract_scheduled_data':"
						" Subsets '" <<m_vScheduledBndSubSets[i].c_str() <<"' not"
						" all contained in ApproximationSpace.");

	//	get subsethandler
		const ISubsetHandler& rSH = *this->function_pattern()->subset_handler();

	// 	loop subsets
		for(size_t si = 0; si < subsetGroup.size(); ++si)
		{
		//	get subset index
			const int subsetIndex = subsetGroup[si];

		//	check that subsetIndex is valid
			if(subsetIndex < 0 || subsetIndex >= rSH.num_subsets())
			{
				UG_LOG("ERROR in 'CRNavierStokesSymBC:extract_scheduled_data':"
						" Invalid subset Index " << subsetIndex <<
						". (Valid is 0, .. , " << rSH.num_subsets() <<").\n");
				return;
			}

		// save the index
			m_vBndSubSetIndex.push_back(subsetIndex);
		}
	}
}

/**
 * The add method for the boundary subsets:
 */
template<typename TDomain>
void
CRNavierStokesSymBC<TDomain>::
add
(
	const char* subsets // string with the ','-separated names of the subsets
)
{
	m_vScheduledBndSubSets.push_back(subsets);
}

/**
 * Prepares the element loop for a given element type: computes the FV-geo, ...
 * Note that there are separate loops for every type of the grid elements.
 */
template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void CRNavierStokesSymBC<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
//	register subsetIndex at Geometry
	static TFVGeom<TElem, dim>& geo = GeomProvider<TFVGeom<TElem,dim> >::get();

// 	Only first order implementation
	if(!(TFVGeom<TElem, dim>::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
						" Geometry set.");

//	check if kinematic Viscosity has been set
	if(!m_imKinViscosity.data_given())
		UG_THROW("CRNavierStokesSymBC::prep_elem_loop:"
						" Kinematic Viscosity has not been set, but is required.\n");

//	check if Density has been set
	if(!m_imDensity.data_given())
		UG_THROW("CRNavierStokesSymBC::prep_elem_loop:"
						" Density has not been set, but is required.\n");

//	extract indices of boundary
	extract_scheduled_data();

//	request the subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
		geo.add_boundary_subset(*subsetIter);
}

/**
 * Finalizes the element loop for a given element type.
 */
template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void CRNavierStokesSymBC<TDomain>::
fsh_elem_loop()
{
	static TFVGeom<TElem, dim>& geo = GeomProvider<TFVGeom<TElem,dim> >::get();

//	remove the bnd subsets
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
		geo.remove_boundary_subset(*subsetIter);
}


/**
 * General initializations of a given grid element for the assembling.
 */
template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void CRNavierStokesSymBC<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
	static TFVGeom<TElem, dim>& geo = GeomProvider<TFVGeom<TElem,dim> >::get();
	try{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("CRNavierStokesSymBC::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	find and set the local and the global positions of the IPs for imports
	typedef typename TFVGeom<TElem, dim>::BF BF;
	typename std::vector<int>::const_iterator subsetIter;

	m_vLocIP.clear(); m_vGloIP.clear();
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
		const int bndSubset = *subsetIter;
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);
		for(size_t i = 0; i < vBF.size(); ++i)
		{
			m_vLocIP.push_back(vBF[i].local_ip());
			m_vGloIP.push_back(vBF[i].global_ip());
		}
	}
	// REMARK: The loop above determines the ordering of the integration points:
	// The "outer ordering" corresponds to the ordering of the subsets in
	// m_vBndSubSetIndex, and "inside" of this ordering, the ip's are ordered
	// according to the order of the boundary faces in the FV geometry structure.

	m_imKinViscosity.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imKinViscosity.set_global_ips(&m_vGloIP[0], m_vGloIP.size());

	m_imDensity.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
	m_imDensity.set_global_ips(&m_vGloIP[0], m_vGloIP.size());
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void CRNavierStokesSymBC<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = GeomProvider<TFVGeom<TElem,dim> >::get();
	typedef typename TFVGeom<TElem, dim>::BF BF;

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	size_t ip = 0;
	for(subsetIter = m_vBndSubSetIndex.begin();
		subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;

	//	get the list of the ip's:
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);

		for (size_t i=0;i < vBF.size();i++)
		{
			BF bf=vBF[i];
			// 	loop shape functions
			for(size_t sh = 0; sh < bf.num_sh(); ++sh)
			{
				////////////////////////////////////////////////////
				// Add normal flux
				////////////////////////////////////////////////////

			// 	Compute flux derivative at IP
				const number flux_sh =  -1.0 * m_imKinViscosity[ip] * m_imDensity[ip]
										* VecDot(bf.global_grad(sh), bf.normal());

			// 	Add flux derivative  to local matrix
				for(int d1 = 0; d1 < dim; ++d1)
				{
					J(d1, bf.node_id(), d1, sh) += flux_sh*bf.normal()[d1];
				}
				if(!m_spMaster->laplace())
				{
					for(int d1 = 0; d1 < dim; ++d1)
						for(int d2 = 0; d2 < dim; ++d2)
						{
							const number flux2_sh = -1.0 * m_imKinViscosity[ip] * m_imDensity[ip]
													* bf.global_grad(sh)[d1] * bf.normal()[d2];
							J(d1, bf.node_id(), d2, sh) += flux2_sh*bf.normal()[d1];
						}
				}

			}// end of loop shape functions

			/// pressure term
			for(int d1 = 0; d1 < dim; ++d1)
			{
				// pressure is constant over element
				J(d1, bf.node_id(), _P_, 0) += bf.normal()[d1];
			}
			// enforce un = 0 in bip by subtracting previously assembled flux
			for(int d1 = 0; d1 < dim; ++d1)
			{
				J(_P_, 0 , d1, bf.node_id()) -= bf.normal()[d1];
			}
		}
	}
}

template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void CRNavierStokesSymBC<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implemented
	UG_ASSERT((TFVGeom<TElem, dim>::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom<TElem, dim>& geo = GeomProvider<TFVGeom<TElem,dim> >::get();
	typedef typename TFVGeom<TElem, dim>::BF BF;

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	size_t ip = 0;
	for(subsetIter = m_vBndSubSetIndex.begin();
		subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;

	//	get the list of the ip's:
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);

		for (size_t i=0;i < vBF.size();i++)
		{
			BF bf = vBF[i];
		// 	1. Interpolate Functional Matrix of velocity at ip
			MathMatrix<dim, dim> gradVel;
			for(int d1 = 0; d1 < dim; ++d1)
				for(int d2 = 0; d2 <dim; ++d2)
				{
				//	sum up contributions of each shape
					gradVel(d1, d2) = 0.0;
					for(size_t sh = 0; sh < bf.num_sh(); ++sh)
						gradVel(d1, d2) += bf.global_grad(sh)[d2]
			                 * u(d1, sh);
				}

		//	2. Compute flux
			MathVector<dim> diffFlux;

		//	Add (\nabla u) \cdot \vec{n}
			MatVecMult(diffFlux, gradVel, bf.normal());

		//	Add (\nabla u)^T \cdot \vec{n}
			if(!m_spMaster->laplace())
				TransposedMatVecMultAdd(diffFlux, gradVel, bf.normal());

		//	scale by viscosity
			VecScale(diffFlux, diffFlux, (-1.0) * m_imKinViscosity[ip]* m_imDensity[ip]);

			for(int d1 = 0; d1 < dim; ++d1)
			{
				d(d1, bf.node_id()) += diffFlux[d1]*bf.normal()[d1];
			}
			// pressure term
			for(int d1 = 0; d1 < dim; ++d1)
			{
				d(d1, bf.node_id()) += u(_P_, 0) * bf.normal()[d1];
			}
			// un = 0 in bip
			for(int d1 = 0; d1 < dim; ++d1)
			{
				d(_P_, 0 ) -= bf.normal()[d1] * u(d1,bf.node_id());
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
CRNavierStokesSymBC<TDomain>::
CRNavierStokesSymBC(SmartPtr< IncompressibleNavierStokesBase<TDomain> > spMaster)
: IElemDisc<TDomain>(spMaster->symb_fcts(), spMaster->symb_subsets()), m_spMaster (spMaster)
{
//	check number of functions
	if(this->num_fct() != dim+1)
		UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");

//	yet no boundary subsets
	m_vBndSubSetIndex.clear ();

//	register imports
	this->register_import(m_imKinViscosity);
	this->register_import(m_imDensity);

//	initialize the imports from the master discretization
	m_imKinViscosity.set_data(spMaster->kinematic_viscosity ());
	m_imDensity.set_data(spMaster->density ());

//	register assemble functions
	this->register_all_funcs(false);
}


////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void CRNavierStokesSymBC<Domain1d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		UG_THROW("Crouxeiz-Raviart only senseful in dimension >= 2");
	}
	else
	{
		UG_THROW("CRNavierStokesSymBC: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_2
template<>
void CRNavierStokesSymBC<Domain2d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Triangle, CRFVGeometry >();
		register_func<Quadrilateral, CRFVGeometry >();
	}
	else
	{
		UG_THROW("CRNavierStokesSymBC: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_3
template<>
void CRNavierStokesSymBC<Domain3d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Tetrahedron, CRFVGeometry >();
		register_func<Prism, CRFVGeometry>();
		register_func<Pyramid, CRFVGeometry >();
		register_func<Hexahedron, CRFVGeometry >();
	}
	else
	{
		UG_THROW("CRNavierStokesSymBC: Hanging Nodes not implemented.")
	}
}
#endif


template<typename TDomain>
template<typename TElem, template <class Elem, int WorldDim> class TFVGeom>
void
CRNavierStokesSymBC<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	   	id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct(	id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(	id, &T::template add_rhs_elem<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class CRNavierStokesSymBC<Domain1d>;
#endif
#ifdef UG_DIM_2
template class CRNavierStokesSymBC<Domain2d>;
#endif
#ifdef UG_DIM_3
template class CRNavierStokesSymBC<Domain3d>;
#endif

} // namespace NavierStokes
} // end namespace ug
