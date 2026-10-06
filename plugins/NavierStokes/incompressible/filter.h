/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__NAVIER_STOKES__INCOMPRESSIBLE__FILTER__
#define __H__UG__NAVIER_STOKES__INCOMPRESSIBLE__FILTER__

#include <vector>
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_grid/tools/periodic_boundary_manager.h"
#include "lib_grid/lg_base.h"
#include "common/profiler/profiler.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "fvcr/disc_constraint_fvcr.h"

namespace ug{
namespace NavierStokes{

template <typename TGridFunction>
class WallObject{
public:
	inline int si(){return m_subset_index;}
	inline number coord(){return m_coord;}
	inline size_t direction(){return m_direction;}
	number dist(MathVector<TGridFunction::dim> co){
		return std::abs(co[m_direction]-m_coord);
	}
	WallObject(SmartPtr<TGridFunction> u,size_t direction,number coord,const char* subset){
		m_subset_index = u->subset_id_by_name(subset);
		m_direction=direction;
		m_coord=coord;
	}
private:
	int m_subset_index;
	size_t m_direction;
	number m_coord;
};

template <int dim,typename elem_type,typename TGridFunction>
void copyAttachmentToGridFunction(SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<elem_type,Attachment<MathVector<dim> > >& aaU){
	/// element iterator
	typedef typename TGridFunction::template traits<elem_type>::const_iterator ElemIterator;
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	//	get domain
	domain_type& domain = *u->domain().get();
	//	create Multiindex
	std::vector<DoFIndex> multInd;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		ElemIterator iter = u->template begin<elem_type>(si);
		ElemIterator iterEnd = u->template end<elem_type>(si);
		for(  ;iter !=iterEnd; ++iter)
		{
			elem_type* elem = *iter;
			for (int d=0;d<dim;d++){
				u->dof_indices(elem, d, multInd);
				DoFRef(*u,multInd[0])=aaU[elem][d];
			}
		}
	}
}

template <int dim,typename elem_type,typename TGridFunction>
void copyGridFunctionToAttachment(PeriodicAttachmentAccessor<elem_type,Attachment<MathVector<dim> > >& aaU,SmartPtr<TGridFunction> u){
	/// element iterator
	typedef typename TGridFunction::template traits<elem_type>::const_iterator ElemIterator;
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	//	get domain
	domain_type& domain = *u->domain().get();
	//	create Multiindex
	std::vector<DoFIndex> multInd;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		ElemIterator iter = u->template begin<elem_type>(si);
		ElemIterator iterEnd = u->template end<elem_type>(si);
		for(  ;iter !=iterEnd; ++iter)
		{
			elem_type* elem = *iter;
			for (int d=0;d<dim;d++){
				u->dof_indices(elem, d, multInd);
				aaU[elem][d]=DoFRef(*u,multInd[0]);
			}
		}
	}
}
	
template <typename elem_type,typename attachment_type,typename grid_type>
void initAttachment(PeriodicAttachmentAccessor<elem_type,attachment_type >& accessor,
					attachment_type& attachment,
					grid_type& grid
					){
	grid.template attach_to<elem_type>(attachment);
	accessor.access(grid,attachment);
};
		
template <typename TGridFunction>
class FilterBaseClass{
public:
	/// dimension
	static const size_t dim = TGridFunction::dim;

    /// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	/// side type
	typedef typename elem_type::side side_type;

	// vertex type
	typedef Vertex vertex_type;

	// data types used in attachments
	typedef MathVector<dim> type0;
	typedef number type1;
	typedef MathSymmetricMatrix<dim> type2;

	// filter for side attachment
	virtual void apply(PeriodicAttachmentAccessor<side_type,Attachment<type0> >& aaUHat,
			PeriodicAttachmentAccessor<side_type,Attachment<type0> >& aaU){}

	virtual void apply(PeriodicAttachmentAccessor<side_type,Attachment<type1> >& aaUHat,
			PeriodicAttachmentAccessor<side_type,Attachment<type1> >& aaU){}
	
	virtual void apply(PeriodicAttachmentAccessor<side_type,Attachment<type2> >& aaUHat,
			PeriodicAttachmentAccessor<side_type,Attachment<type2> >& aaU){}

	// filter for crouzeix-raviart type grid function
	virtual void apply(PeriodicAttachmentAccessor<side_type,Attachment<type0> >& aaUHat,
			SmartPtr<TGridFunction> u){}

	// filter for vertex attachment
	virtual void apply(PeriodicAttachmentAccessor<vertex_type,Attachment<type0> >& aaUHat,
					   PeriodicAttachmentAccessor<vertex_type,Attachment<type0> >& aaU){}

	virtual void apply(PeriodicAttachmentAccessor<vertex_type,Attachment<type1> >& aaUHat,
					   PeriodicAttachmentAccessor<vertex_type,Attachment<type1> >& aaU){}
	
	virtual void apply(PeriodicAttachmentAccessor<vertex_type,Attachment<type2> >& aaUHat,
					   PeriodicAttachmentAccessor<vertex_type,Attachment<type2> >& aaU){}

	// filter for general grid function, result stored in nodes
	virtual void apply(PeriodicAttachmentAccessor<vertex_type,Attachment<type0> >& aaUHat,
					   SmartPtr<TGridFunction> u){}

	// filter grid function
	virtual void apply(SmartPtr<TGridFunction> u){}
	
	virtual number width(vertex_type* vrt){return -1;}
	
	virtual number width(side_type* side){return -1;}
	
	virtual void compute_filterwidth(){}
	
	FilterBaseClass(){}
	virtual ~FilterBaseClass(){}
};

// wrapper class for filtering
template <typename TImpl,typename TGridFunction>
class FilterImplBaseClass
: public FilterBaseClass<TGridFunction>
{
	public:
	/// dimension
	static const size_t dim = TGridFunction::dim;
	
    /// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;
	
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	/// side type
	typedef typename elem_type::side side_type;
	
	// vertex type
	typedef Vertex vertex_type;
	
	// data types
	typedef typename FilterBaseClass<TGridFunction>::type0 type0;
	typedef typename FilterBaseClass<TGridFunction>::type1 type1;
	typedef typename FilterBaseClass<TGridFunction>::type2 type2;
	
	// filter for side attachment
	void apply(PeriodicAttachmentAccessor<side_type,Attachment<type0> >& aaUHat,
			   PeriodicAttachmentAccessor<side_type,Attachment<type0> >& aaU){
		getImpl().apply_filter(aaUHat,SPNULL,aaU);
	}
	
	void apply(PeriodicAttachmentAccessor<side_type,Attachment<type1> >& aaUHat,
			   PeriodicAttachmentAccessor<side_type,Attachment<type1> >& aaU){
		getImpl().apply_filter(aaUHat,SPNULL,aaU);
	}
	
	void apply(PeriodicAttachmentAccessor<side_type,Attachment<type2> >& aaUHat,
			   PeriodicAttachmentAccessor<side_type,Attachment<type2> >& aaU){
		getImpl().apply_filter(aaUHat,SPNULL,aaU);
	}
	
	// filter for crouzeix-raviart type grid function
	void apply(PeriodicAttachmentAccessor<side_type,Attachment<type0> >& aaUHat,
			SmartPtr<TGridFunction> u){
			UG_LOG("son class\n");
			getImpl().apply_filter(aaUHat,u,aaUHat);
	}
	
	// filter for vertex attachment
	void apply(PeriodicAttachmentAccessor<vertex_type,Attachment<type0> >& aaUHat,
			   PeriodicAttachmentAccessor<vertex_type,Attachment<type0> >& aaU){
		getImpl().apply_filter(aaUHat,SPNULL,aaU);
	}
	
	void apply(PeriodicAttachmentAccessor<vertex_type,Attachment<type1> >& aaUHat,
			   PeriodicAttachmentAccessor<vertex_type,Attachment<type1> >& aaU){
		getImpl().apply_filter(aaUHat,SPNULL,aaU);
	}
	
	void apply(PeriodicAttachmentAccessor<vertex_type,Attachment<type2> >& aaUHat,
			   PeriodicAttachmentAccessor<vertex_type,Attachment<type2> >& aaU){
		getImpl().apply_filter(aaUHat,SPNULL,aaU);
	}
	
	// filter for general grid function, result stored in nodes
	void apply(PeriodicAttachmentAccessor<vertex_type,Attachment<type0> >& aaUHat,
			SmartPtr<TGridFunction> u){
			if (u->local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
				getImpl().apply_filter(aaUHat,u);
			} else {
				getImpl().apply_filter(aaUHat,u,aaUHat);
			}
	}

	template <typename TElem>
	void apply_(SmartPtr<TGridFunction> u){
		// define attachment types
		Attachment<MathVector<dim> > aUHat;
		PeriodicAttachmentAccessor<TElem,Attachment<MathVector<dim> > > acUHat;

		domain_type& domain = *u->domain().get();
		grid_type& grid = *domain.grid();

		// attach
		grid.template attach_to<TElem>(aUHat);
		acUHat.access(grid,aUHat);
		getImpl().apply_filter(acUHat,u,acUHat);
		copyAttachmentToGridFunction<dim,TElem,TGridFunction>(u,acUHat);
		grid.template detach_from<TElem>(aUHat);
	}

	// filter grid function
	void apply(SmartPtr<TGridFunction> u){
		if (u->local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
			apply_<side_type>(u);
		} else {
			apply_<vertex_type>(u);
		}
	}
	
	template<typename VType,typename TElem>
	void copyWallData(PeriodicAttachmentAccessor<TElem,Attachment<VType> >& acUHat,SmartPtr<TGridFunction> u,std::vector<WallObject<TGridFunction> > walls){
		typedef typename TGridFunction::domain_type TDomain;
		typedef typename TDomain::grid_type TGrid;
		TDomain domain = *u->domain().get();
		typedef typename TGridFunction::template traits<TElem>::const_iterator ElemIterator;
		PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
		std::vector<DoFIndex> dofInd;
		for (int j=0;j<(int)walls.size();j++){
			int si = walls[j].si();
			ElemIterator elemIter = u->template begin<TElem>(si);
			ElemIterator elemIterEnd = u->template end<TElem>(si);
			for(  ;elemIter !=elemIterEnd; elemIter++)
			{
				TElem* elem = *elemIter;
				if (pbm && pbm->is_slave(elem)) continue;
				VType localValue;
				for	(int d=0;d<TGridFunction::dim;d++){
					u->dof_indices(elem, d, dofInd);
					assignVal(localValue,d,DoFRef(*u,dofInd[0]));
				}
				acUHat[elem]=localValue;
			}
		}
	}

	template<typename VType,typename TElem>
	void copyWallData(PeriodicAttachmentAccessor<TElem,Attachment<VType> >& acUHat,PeriodicAttachmentAccessor<TElem,Attachment<VType> >&  acU,
				  SmartPtr<TGridFunction> u,std::vector<WallObject<TGridFunction> > walls){
		typedef typename TGridFunction::domain_type TDomain;
		typedef typename TDomain::grid_type TGrid;
		TDomain domain = *u->domain().get();
		typedef typename TGridFunction::template traits<TElem>::const_iterator ElemIterator;
		PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
		std::vector<DoFIndex> dofInd;
		for (int j=0;j<(int)walls.size();j++){
			int si = walls[j].si();
			ElemIterator elemIter = u->template begin<TElem>(si);
			ElemIterator elemIterEnd = u->template end<TElem>(si);
			for(  ;elemIter !=elemIterEnd; elemIter++)
			{
				TElem* elem = *elemIter;
				if (pbm && pbm->is_slave(elem)) continue;
				acUHat[elem]=acU[elem];
			}
		}
	}
	
	void assignVal(number& v,size_t ind,number value){
		v=value;	
	}
	
	void assignVal(MathVector<2>& v,size_t ind,number value){
		v[ind]=value;
	}
	
	void assignVal(MathVector<3>& v,size_t ind,number value){
		v[ind]=value;
	}
	
	void assignVal(MathSymmetricMatrix<2>& v,size_t ind,number value){
		v[ind]=value;
	}
	
	void assignVal(MathSymmetricMatrix<3>& v,size_t ind,number value){
		v[ind]=value;
	}
	
	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};
	
template <int dim>
struct Region{
	size_t firstIndex;
	MathVector<dim> periodicOffset;
};
	
// Box filter with constant filter width
template <typename TGridFunction>
class ConstantBoxFilter
	: public FilterImplBaseClass<ConstantBoxFilter<TGridFunction>,TGridFunction>
{
public:
	// base type
	typedef FilterImplBaseClass<ConstantBoxFilter<TGridFunction>,TGridFunction> base_type;
	
	// dimension
	static const size_t dim = TGridFunction::dim;
	
    // element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;
	
	//	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
	// side type
	typedef typename elem_type::side side_type;
	
	// vertex type
	typedef Vertex vertex_type;
	
	//	grid type
	typedef typename domain_type::grid_type grid_type;	
	
	//	position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// side iterator
	typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;

	/// vertex iterator
	typedef typename TGridFunction::template traits<vertex_type>::const_iterator VertexIterator;

	// accessor types
	typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	typedef PeriodicAttachmentAccessor<vertex_type,ANumber > aVertexNumber;

	// filter for crouzeix-raviart type data
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaUHat,
			   SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaU);
			
	// filter for fv1 type data
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaUHat,
					  SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaU){
							UG_THROW("not implemented");
						};

	// filter for crouzeix-raviart grid function with output in node attachment
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaUHat,
			   SmartPtr<TGridFunction> u){
							UG_THROW("not implemented");
						};
	
	// compute average mesh size
	number compute_average_element_size(SmartPtr<TGridFunction> u);
	
	inline number width(vertex_type* vrt){ return m_width; }
	
	inline number width(side_type* side){ return m_width; }
	
	void add_wall(WallObject<TGridFunction> w){
		m_walls.push_back(w);
	}
	
	// constructors	
	ConstantBoxFilter(SmartPtr<TGridFunction> u,number width){
		m_width = width;
		m_uInfo = u;
		domain_type& domain = *u->domain().get();
		if (u->local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
			domain.grid()->template attach_to<side_type>(m_aSideVolume);
			m_acSideVolume.access(*domain.grid(),m_aSideVolume);
		} else {
			domain.grid()->template attach_to<vertex_type>(m_aVertexVolume);
			m_acSideVolume.access(*domain.grid(),m_aVertexVolume);
		}
	}
	
	ConstantBoxFilter(SmartPtr<TGridFunction> u){
		m_uInfo = u;
		domain_type& domain = *u->domain().get();
		if (u->local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
			domain.grid()->template attach_to<side_type>(m_aSideVolume);
			m_acSideVolume.access(*domain.grid(),m_aSideVolume);
		} else {
			domain.grid()->template attach_to<vertex_type>(m_aVertexVolume);
			m_acSideVolume.access(*domain.grid(),m_aVertexVolume);
		}
		m_width = pow(compute_average_element_size(u),(number)1.0/dim);
	}
	
	~ConstantBoxFilter(){
		domain_type& domain = *m_uInfo->domain().get();
		domain.grid()->template detach_from<side_type>(m_aSideVolume);
		domain.grid()->template detach_from<vertex_type>(m_aVertexVolume);
	}
	
private:
	template <typename VType>
	void collectSides(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acUHat,
					  PeriodicAttachmentAccessor<side_type,Attachment<number> >& acVolume,
					  std::vector< MathVector<dim> >& coord,
					  VType values[DimFV1Geometry<dim>::maxNumSCV],
					  number volumes[DimFV1Geometry<dim>::maxNumSCV],
					  std::vector<side_type*>& sides
					  );
	
	//  volume attachment
	ANumber m_aSideVolume;
	aSideNumber m_acSideVolume;
	ANumber m_aVertexVolume;
	aVertexNumber m_acVertexVolume;

	number m_width;
	
	position_accessor_type m_posAcc;

	// grid function
	SmartPtr<TGridFunction> m_uInfo;
	
	std::vector<WallObject<TGridFunction> > m_walls;
	
	using base_type::assignVal;
	using base_type::copyWallData;
};
	
// Box filter with constant filter width
template <typename TGridFunction>
class VariableBoxFilter
: public FilterImplBaseClass<VariableBoxFilter<TGridFunction>,TGridFunction>
{
public:
	// base type
	typedef FilterImplBaseClass<VariableBoxFilter<TGridFunction>,TGridFunction> base_type;
	
	// dimension
	static const size_t dim = TGridFunction::dim;
		
	// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;
		
	//	domain type
	typedef typename TGridFunction::domain_type domain_type;
		
	// side type
	typedef typename elem_type::side side_type;
		
	// vertex type
	typedef Vertex vertex_type;
		
	//	grid type
	typedef typename domain_type::grid_type grid_type;	
		
	//	position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;
		
	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
		
	/// side iterator
	typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;
		
	/// vertex iterator
	typedef typename TGridFunction::template traits<vertex_type>::const_iterator VertexIterator;
		
	// accessor types
	typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	typedef PeriodicAttachmentAccessor<vertex_type,ANumber > aVertexNumber;
		
	// filter for crouzeix-raviart type data
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaUHat,
					  SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaU);
		
	// filter for fv1 type data
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaUHat,
					  SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaU){
			UG_THROW("not implemented");
	};
		
	// filter for crouzeix-raviart grid function with output in node attachment
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaUHat,
					  SmartPtr<TGridFunction> u){
			UG_THROW("not implemented");
	};
		
	// compute average mesh size
	number compute_average_element_size(SmartPtr<TGridFunction> u);
		
	inline number width(vertex_type* vrt){ 
		std::vector<DoFIndex> multInd; 
		m_width->inner_dof_indices(vrt,0, multInd);
		return DoFRef(*m_width,multInd[0]); 
	}
	
	inline number width(side_type* side){
		std::vector<DoFIndex> multInd; 
		m_width->inner_dof_indices(side,0, multInd);
		return DoFRef(*m_width,multInd[0]);
	}
	
	template <typename TElem>
	void init_width(number width){
		domain_type& domain = *m_width->domain().get();
		std::vector<DoFIndex> multInd; 
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
		{
			//	get iterators
			ElemIterator iter = m_width->template begin<elem_type>(si);
			ElemIterator iterEnd = m_width->template end<elem_type>(si);
			
			//	loop elements of dimension
			for(  ;iter !=iterEnd; ++iter)
			{
				//	get Elem
				elem_type* elem = *iter;
				m_width->inner_dof_indices(elem,0, multInd);
				DoFRef(*m_width,multInd[0])=width;
			}
		}
	}
	
	template <typename TElem>
	number max_width(number width){
		domain_type& domain = *m_width->domain().get();
		std::vector<DoFIndex> multInd; 
		number maxwidth = 0;
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
		{
			//	get iterators
			ElemIterator iter = m_width->template begin<elem_type>(si);
			ElemIterator iterEnd = m_width->template end<elem_type>(si);
			
			//	loop elements of dimension
			for(  ;iter !=iterEnd; ++iter)
			{
				//	get Elem
				elem_type* elem = *iter;
				m_width->inner_dof_indices(elem,0, multInd);
				maxwidth=std::max(DoFRef(*m_width,multInd[0]),maxwidth);
			}
		}
		m_maxwidth = maxwidth;
		return maxwidth;
	}
	
	void add_wall(WallObject<TGridFunction> w){
		m_walls.push_back(w);
	}
	
	template <typename VType,typename TElem>
	void check_volume_sizes(PeriodicAttachmentAccessor<TElem,Attachment<VType> >& acUHat,SmartPtr<TGridFunction> u){
		domain_type& domain = *m_width->domain().get();
		
	}

	// constructors	
	VariableBoxFilter(SmartPtr<TGridFunction> u,SmartPtr<TGridFunction> fwidth,bool initWidth=false){
		m_width = fwidth;
		m_uInfo = u;
		domain_type& domain = *u->domain().get();
		if (u->local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
			domain.grid()->template attach_to<side_type>(m_aSideVolume);
			m_acSideVolume.access(*domain.grid(),m_aSideVolume);
		} else {
			domain.grid()->template attach_to<vertex_type>(m_aVertexVolume);
			m_acSideVolume.access(*domain.grid(),m_aVertexVolume);
		}
		if (initWidth){
			number width = pow(compute_average_element_size(u),(number)1.0/dim);
			if (m_width->local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
				init_width<side_type>(width);
			}
			if (m_width->local_finite_element_id(0) == LFEID(LFEID::LAGRANGE, dim, 1)){
				init_width<vertex_type>(width);
			}
		}
	}
		
	~VariableBoxFilter(){
		domain_type& domain = *m_uInfo->domain().get();
		domain.grid()->template detach_from<side_type>(m_aSideVolume);
		domain.grid()->template detach_from<vertex_type>(m_aVertexVolume);
	}
		
private:
	template <typename VType>
	void collectSides(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acUHat,
					  PeriodicAttachmentAccessor<side_type,Attachment<number> >& acVolume,
					  std::vector< MathVector<dim> >& coord,
					  VType values[DimFV1Geometry<dim>::maxNumSCV],
					  number volumes[DimFV1Geometry<dim>::maxNumSCV],
					  std::vector<side_type*>& sides
					);
		
	//  volume attachment
	ANumber m_aSideVolume;
	aSideNumber m_acSideVolume;
	ANumber m_aVertexVolume;
	aVertexNumber m_acVertexVolume;
		
	position_accessor_type m_posAcc;
		
	// grid function
	SmartPtr<TGridFunction> m_uInfo;
		
	// filter width
	SmartPtr<TGridFunction> m_width;
	
	number m_maxwidth;
	
	std::vector<WallObject<TGridFunction> > m_walls;
	
	using base_type::assignVal;
	using base_type::copyWallData;
};
	


// Box filter with local filter volume given by vertex-centered FV (fv1 geometry)
template <typename TGridFunction>
class FV1BoxFilter
	: public FilterImplBaseClass<FV1BoxFilter<TGridFunction>,TGridFunction>
{
public:
	// base type
	typedef FilterImplBaseClass<FV1BoxFilter<TGridFunction>,TGridFunction> base_type;

	// dimension
	static const size_t dim = TGridFunction::dim;

	// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	//	domain type
	typedef typename TGridFunction::domain_type domain_type;

	// side type
	typedef typename elem_type::side side_type;

	// vertex type
	typedef Vertex vertex_type;

	//	grid type
	typedef typename domain_type::grid_type grid_type;

	//	position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	// side iterator
	typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;

	// vertex iterator
	typedef typename TGridFunction::template traits<vertex_type>::const_iterator VertexIterator;

	typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	typedef PeriodicAttachmentAccessor<vertex_type,ANumber > aVertexNumber;

	// filter for crouzeix-raviart type data
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaUHat,
					   SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaU);

	// filter for fv1 type data
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaUHat,
					   SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaU);

	// filter for crouzeix-raviart grid function with output in node attachment
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaUHat,
			   SmartPtr<TGridFunction> u)
						{
							UG_THROW("not implemented");
						};

	// return filter width
	inline number width(vertex_type* vrt){ return pow(m_acVertexVolume[vrt],(number)1.0/dim); }

	inline number width(side_type* side){ return pow(m_acSideVolume[side],(number)1.0/dim); }
	
	void compute_filterwidth(){
		if (m_uInfo->local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
			compute_filterwidth_fvcr();
		} else {
			compute_filterwidth_fv1();
		}
	}
	
	void compute_filterwidth_fv1();
	
	void compute_filterwidth_fvcr();
	
	void add_wall(WallObject<TGridFunction> w){
		m_walls.push_back(w);
	}

	// constructor
	FV1BoxFilter(SmartPtr<TGridFunction> u){
		m_uInfo = u;
		domain_type& domain = *u->domain().get();
		if (u->local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
			domain.grid()->template attach_to<side_type>(m_aSideVolume);
			m_acSideVolume.access(*domain.grid(),m_aSideVolume);
		} else {
			domain.grid()->template attach_to<vertex_type>(m_aVertexVolume);
			m_acSideVolume.access(*domain.grid(),m_aVertexVolume);
		}
	}

	// destructor
	~FV1BoxFilter(){
		domain_type& domain = *m_uInfo->domain().get();
		domain.grid()->template detach_from<side_type>(m_aSideVolume);
		domain.grid()->template detach_from<vertex_type>(m_aVertexVolume);
	}

private:
	//  volume attachment
	ANumber m_aSideVolume;
	aSideNumber m_acSideVolume;
	ANumber m_aVertexVolume;
	aVertexNumber m_acVertexVolume;

	// grid function
	SmartPtr<TGridFunction> m_uInfo;
	
	std::vector<WallObject<TGridFunction> > m_walls;
	
	using base_type::assignVal;
	using base_type::copyWallData;
};


// Box filter with local filter volume given by vertex-centered FV (fv1 geometry)
template <typename TGridFunction>
class FVCRBoxFilter
	: public FilterImplBaseClass<FVCRBoxFilter<TGridFunction>,TGridFunction>
{
public:
	// base type
	typedef FilterImplBaseClass<FVCRBoxFilter<TGridFunction>,TGridFunction> base_type;
	
	// dimension
	static const size_t dim = TGridFunction::dim;

	// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	//	domain type
	typedef typename TGridFunction::domain_type domain_type;

	// side type
	typedef typename elem_type::side side_type;

	// vertex type
	typedef Vertex vertex_type;

	//	grid type
	typedef typename domain_type::grid_type grid_type;

	//	position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	// side iterator
	typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;

	// vertex iterator
	typedef typename TGridFunction::template traits<vertex_type>::const_iterator VertexIterator;

	typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	typedef PeriodicAttachmentAccessor<vertex_type,ANumber > aVertexNumber;

	// filter for crouzeix-raviart type data
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaUHat,
					   SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaU);

	// filter for fv1 type data
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaUHat,
					   SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaU)
							{
								UG_THROW("not implemented");
							};

	// filter for crouzeix-raviart grid function with output in node attachment
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaUHat,
						SmartPtr<TGridFunction> u)
						{
							UG_THROW("not implemented");
						};
						
	void add_wall(WallObject<TGridFunction> w){
		m_walls.push_back(w);
	}

	// return filter width
	inline number width(vertex_type* vrt){ UG_THROW("not implemented"); return 0; }

	inline number width(side_type* side){ return pow(m_acSideVolume[side],(number)1.0/dim); }
	
	void compute_filterwidth(){
		if (m_uInfo->local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
			compute_filterwidth_fvcr();
		} else {
			compute_filterwidth_fv1();
		}
	}
	
	void compute_filterwidth_fv1(){
		UG_LOG("not implemented\n");
	}
	
	void compute_filterwidth_fvcr();

	// constructors
	FVCRBoxFilter(SmartPtr<TGridFunction> u){
		m_uInfo = u;
		domain_type& domain = *u->domain().get();
		if (u->local_finite_element_id(0) != LFEID(LFEID::CROUZEIX_RAVIART, dim, 1))
			UG_THROW("Only implemented for Crouzeix-Raviart type functions.");
		domain.grid()->template attach_to<side_type>(m_aSideVolume);
		m_acSideVolume.access(*domain.grid(),m_aSideVolume);
	}

	~FVCRBoxFilter(){
		domain_type& domain = *m_uInfo->domain().get();
		domain.grid()->template detach_from<side_type>(m_aSideVolume);
	}

	private:
	//  volume attachment
	ANumber m_aSideVolume;
	aSideNumber m_acSideVolume;
	//ANumber m_aVertexVolume;
	//aVertexNumber m_acVertexVolume;

	// grid function
	SmartPtr<TGridFunction> m_uInfo;
	
	std::vector<WallObject<TGridFunction> > m_walls;
	
	using base_type::assignVal;
	using base_type::copyWallData;
};
	
// Box filter with local filter volume given by vertex-centered FV (fv1 geometry)
template <typename TGridFunction>
class ElementBoxFilter
: public FilterImplBaseClass<ElementBoxFilter<TGridFunction>,TGridFunction>
{
public:
	// base type
	typedef FilterImplBaseClass<ElementBoxFilter<TGridFunction>,TGridFunction> base_type;
	
	// dimension
	static const size_t dim = TGridFunction::dim;
		
	// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;
		
	//	domain type
	typedef typename TGridFunction::domain_type domain_type;
		
	// side type
	typedef typename elem_type::side side_type;
		
	// vertex type
	typedef Vertex vertex_type;
		
	//	grid type
	typedef typename domain_type::grid_type grid_type;
		
	//	position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;
		
	// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;
		
	// side iterator
	typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;
		
	// vertex iterator
	typedef typename TGridFunction::template traits<vertex_type>::const_iterator VertexIterator;
		
	typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	typedef PeriodicAttachmentAccessor<vertex_type,ANumber > aVertexNumber;
		
	// filter for crouzeix-raviart type data
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaUHat,
						SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaU);
		
	// filter for fv1 type data
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaUHat,
					  SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaU);
		
	// filter for crouzeix-raviart grid function with output in node attachment
	template <typename VType>
	void apply_filter(PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& aaUHat,
						SmartPtr<TGridFunction> u)
	{
		UG_THROW("not implemented");
	};
		
	// return filter width
	inline number width(vertex_type* vrt){ return pow(m_acVertexVolume[vrt],(number)1.0/dim); }
		
	inline number width(side_type* side){ return pow(m_acSideVolume[side],(number)1.0/dim); }
		
	void compute_filterwidth(){
		if (m_uInfo->local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
			compute_filterwidth_fvcr();
		} else {
			compute_filterwidth_fv1();
		}
	}
	
	void compute_filterwidth_fv1();
	
	void compute_filterwidth_fvcr();
	
	void add_wall(WallObject<TGridFunction> w){
		m_walls.push_back(w);
	}
	
	// constructors
	ElementBoxFilter(SmartPtr<TGridFunction> u){
		m_uInfo = u;
		domain_type& domain = *u->domain().get();
		if (u->local_finite_element_id(0) == LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
			domain.grid()->template attach_to<side_type>(m_aSideVolume);
			m_acSideVolume.access(*domain.grid(),m_aSideVolume);
		} else {
			domain.grid()->template attach_to<side_type>(m_aVertexVolume);
			m_acSideVolume.access(*domain.grid(),m_aVertexVolume);
		}
	}

	~ElementBoxFilter(){
		domain_type& domain = *m_uInfo->domain().get();
		domain.grid()->template detach_from<side_type>(m_aSideVolume);
		domain.grid()->template detach_from<vertex_type>(m_aVertexVolume);
	}
		
private:
	//  volume attachment
	ANumber m_aSideVolume;
	aSideNumber m_acSideVolume;
	ANumber m_aVertexVolume;
	aVertexNumber m_acVertexVolume;
		
	// grid function
	SmartPtr<TGridFunction> m_uInfo;
	
	std::vector<WallObject<TGridFunction> > m_walls;
	
	using base_type::assignVal;
	using base_type::copyWallData;
};
	
} // end namespace NavierStokes
} // end namespace ug

#include "filter_impl.h"

#endif /* __H__UG__LIB_DISC__NAVIER_STOKES__INCOMPRESSIBLE__FILTER__ */
