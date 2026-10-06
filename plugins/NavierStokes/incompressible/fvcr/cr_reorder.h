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

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FVCR__CR_REORDER__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FVCR__CR_REORDER__

#include <vector>
#include "lib_disc/function_spaces/approximation_space.h"
#include "common/profiler/profiler.h"
#include <boost/graph/adjacency_list.hpp>
#include "lib_disc/dof_manager/dof_distribution.h"

namespace ug{

/**
 * This class is used to provide an ordering for indices. The ordering relation
 * is based on the connectivity-degree, i.e. on the number of connections the
 * index has. The indices with less connections are ordered first.
 */
struct CompareDeg {
///	constructor, passing field with connections for each index
	CompareDeg(const std::vector<size_t>& vInfo) : m_degree(vInfo) {}

///	comparison operator
	bool operator() (size_t i,size_t j)
	{
		return (m_degree[i] < m_degree[j]);
	}

private:
///	storage field of connections of each index
	const std::vector<size_t>& m_degree;
};

void computeDegree(std::vector<size_t>& degree,std::vector<std::vector<size_t> >& vvConnections,size_t minpind);

// compute adjacency graph on Crouzeix-Raviart-elements
template <typename TElem, typename TGridFunction>
void cr_get_connections(std::vector<std::vector<size_t> >& vvConnection,size_t& minpind, DoFDistribution& dd,TGridFunction& u){
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	world dimension
	static const int dim = domain_type::dim;

	typename DoFDistribution::traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	LocalIndices ind;

	size_t num_dd_ind = dd.num_indices();

	//	clear neighbors and resize
	vvConnection.clear();
	vvConnection.resize(num_dd_ind);

	minpind=num_dd_ind;

	for (int si=0;si<u.num_subsets();si++){
		iterBegin = dd.begin<TElem>(si);
		iterEnd = dd.end<TElem>(si);
		// 	check if at least one element exist
		if(iterBegin == iterEnd) continue;
		// 	Loop over all elements
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get Element
			TElem* elem = *iter;
		// 	get global indices
			dd.indices(elem, ind, false);
			size_t num_sides = elem->num_sides();
			size_t num_ind = num_sides * dim + 1;
			std::vector<size_t> globind(num_ind);
		 // velocity component indices
			for (size_t side=0;side<num_sides;side++)
				for (int d=0;d<dim;d++)
					globind[dim*side+d]=ind.index(d,side);
		 // pressure component index
			size_t pind = ind.index(dim,0);
			globind[num_ind-1]=pind;
			if (pind<minpind) minpind=pind;
			//	add couplings to adjacency graph
			std::vector<size_t>::iterator it;
			for (size_t i=0;i<num_ind;i++){
				size_t iIndex = globind[i];
				for (size_t j=0;j<num_ind;j++){
					size_t jIndex = globind[j];
					//	add connection (iIndex->jIndex)
					it = find(vvConnection[iIndex].begin(), vvConnection[iIndex].end(), jIndex);
					if(it == vvConnection[iIndex].end()) vvConnection[iIndex].push_back(jIndex);
					//	add the opposite direction (adjInd -> index)
					it = find(vvConnection[jIndex].begin(), vvConnection[jIndex].end(), iIndex);
					if(it == vvConnection[jIndex].end()) vvConnection[jIndex].push_back(iIndex);
				}
			}
		}
		std::vector<size_t> degree;
		computeDegree(degree,vvConnection,minpind);
		//	Sort neighbours by degreee
		CompareDeg myCompDegree(degree);
		for (size_t i=0;i<vvConnection.size();i++){
			std::sort(vvConnection[i].begin(),vvConnection[i].end(),myCompDegree);
		};
	}
}
	
// extract p graph from connections, two p nodes are defined adjacent if they share a common velocity neighbour
template <class TGraph>
void extractPGraph(TGraph& G,std::vector<std::vector<size_t> >& vvConnection,size_t pmin,size_t dim,bool directed){
	for (size_t i=0;i<pmin;i=i+dim){
		size_t s=vvConnection[i].size();
		int ind1 = vvConnection[i][s-2]-pmin;
		if (ind1<0) continue;
		int ind2 = vvConnection[i][s-1]-pmin;
		boost::add_edge(ind1,ind2,G);
		if (directed==true) boost::add_edge(ind2,ind1,G);
	}	
}
	
// extract v graph from connections
template <class TGraph>
void extractVGraph(TGraph& G,std::vector<std::vector<size_t> >& vvConnection,size_t pmin,size_t dim,bool directed){
	for (size_t i=0;i<pmin;i=i+dim){
		size_t ind1 = (size_t)round((number)i/dim);
		size_t s = vvConnection[i].size();
		for (size_t j=0;j<s-1;j=j+2){
			size_t ind2 = (size_t)round((number)vvConnection[i][j]/dim);
			if (dim*ind2>=pmin) break;
			if (ind1==ind2) continue;
			boost::add_edge(ind1,ind2,G);
			if (directed==true) boost::add_edge(ind2,ind1,G);
		}
	}
}

void CRCuthillMcKee(std::vector<size_t>& newIndex,std::vector<std::vector<size_t> >& vvConnection,size_t dim,
						 size_t minpind,bool bReverse,bool bseparate,bool orderp,bool orderv);
	
/// orders the all DofDistributions of the ApproximationSpace using Cuthill-McKee
template <typename TDomain,typename TGridFunction>
void CROrderCuthillMcKee(ApproximationSpace<TDomain>& approxSpace,TGridFunction& u,
                       bool bReverse,bool bseparate,bool orderp,bool orderv)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
	///	world dimension
	static const int dim = domain_type::dim;
	
	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;
	
	std::vector<std::vector<size_t> > vvConnection;
	size_t minpind;
	
	for (int si=0;si<u.num_subsets();++si){
		if (u.num_fct(si)<2){
			if (u.num_fct(si)!=dim+1){
				UG_THROW("Wrong number of components in approximation space. Needed " << dim+1 << ", given " << u.num_fct(si) << ".\n");
			}
		}
		for (int d=0;d<dim;d++){
			if (u.local_finite_element_id(d) != LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Crouzeix-Raviart type.");
			}
		}
		if (u.local_finite_element_id(dim) != LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0)){
			UG_THROW("Component dim in approximation space must be of piecewise constant type.");
		}
	}
	
	std::vector<SmartPtr<DoFDistribution> > vDD = approxSpace.dof_distributions();

	//	order levels
	for(size_t i = 0; i < vDD.size(); ++i){
			std::vector<size_t> newIndex;
			cr_get_connections<elem_type,TGridFunction>(vvConnection,minpind,*vDD[i],u);
			CRCuthillMcKee(newIndex,vvConnection,dim,minpind,bReverse,bseparate,orderp,orderv);
			vDD[i]->permute_indices(newIndex);
	}
}

void CRSloan(std::vector<size_t>& newIndex,std::vector<std::vector<size_t> >& vvConnection,size_t dim,
						 size_t minpind,bool bseparate,bool orderp,bool orderv);
	
/// orders the all DofDistributions of the ApproximationSpace using Sloan	
template <typename TDomain,typename TGridFunction>
void CROrderSloan(ApproximationSpace<TDomain>& approxSpace,TGridFunction& u,
                       bool bseparate,bool orderp,bool orderv)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
	///	world dimension
	static const int dim = domain_type::dim;
	
	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;
	
	std::vector<std::vector<size_t> > vvConnection;
	size_t minpind;
	
	for (int si=0;si<u.num_subsets();++si){
		if (u.num_fct(si)<2){
			if (u.num_fct(si)!=dim+1){
				UG_THROW("Wrong number of components in approximation space. Needed " << dim+1 << ", given " << u.num_fct(si) << ".\n");
			}
		}
		for (int d=0;d<dim;d++){
			if (u.local_finite_element_id(d) != LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Crouzeix-Raviart type.");
			}
		}
		if (u.local_finite_element_id(dim) != LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0)){
			UG_THROW("Component dim in approximation space must be of piecewise constant type.");
		}
	}
	
	std::vector<SmartPtr<DoFDistribution> > vDD = approxSpace.dof_distributions();
	for(size_t i = 0; i < vDD.size(); ++i){
			std::vector<size_t> newIndex;
			cr_get_connections<elem_type,TGridFunction>(vvConnection,minpind,*vDD[i],u);
			CRSloan(newIndex,vvConnection,dim,minpind,bseparate,orderp,orderv);
			vDD[i]->permute_indices(newIndex);
	}
}

void CRKing(std::vector<size_t>& newIndex,std::vector<std::vector<size_t> >& vvConnection,size_t dim,
						 size_t minpind,bool bReverse,bool bseparate,bool orderp,bool orderv);

/// orders the all DofDistributions of the ApproximationSpace using Cuthill-McKee
template <typename TDomain,typename TGridFunction>
void CROrderKing(ApproximationSpace<TDomain>& approxSpace,TGridFunction& u,
                       bool bReverse,bool bseparate,bool orderp,bool orderv)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
	///	world dimension
	static const int dim = domain_type::dim;
	
	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;
	
	std::vector<std::vector<size_t> > vvConnection;
	size_t minpind;
	
	for (int si=0;si<u.num_subsets();++si){
		if (u.num_fct(si)<2){
			if (u.num_fct(si)!=dim+1){
				UG_THROW("Wrong number of components in approximation space. Needed " << dim+1 << ", given " << u.num_fct(si) << ".\n");
			}
		}
		for (int d=0;d<dim;d++){
			if (u.local_finite_element_id(d) != LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Crouzeix-Raviart type.");
			}
		}
		if (u.local_finite_element_id(dim) != LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0)){
			UG_THROW("Component dim in approximation space must be of piecewise constant type.");
		}
	}

	std::vector<SmartPtr<DoFDistribution> > vDD = approxSpace.dof_distributions();
	for(size_t i = 0; i < vDD.size(); ++i){
			std::vector<size_t> newIndex;
			cr_get_connections<elem_type,TGridFunction>(vvConnection,minpind,*vDD[i],u);
			CRKing(newIndex,vvConnection,dim,minpind,bReverse,bseparate,orderp,orderv);
			vDD[i]->permute_indices(newIndex);
	}
}

void CRMinimumDegree(std::vector<size_t>& newIndex,std::vector<std::vector<size_t> >& vvConnection,size_t dim,
						 size_t minpind,bool bseparate,bool orderp,bool orderv);
	
/// orders the all DofDistributions of the ApproximationSpace using minimum degree algorithm	
template <typename TDomain,typename TGridFunction>
void CROrderMinimumDegree(ApproximationSpace<TDomain>& approxSpace,TGridFunction& u,
                       bool bseparate,bool orderp,bool orderv)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	
	///	world dimension
	static const int dim = domain_type::dim;
	
	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;
	
	std::vector<std::vector<size_t> > vvConnection;
	size_t minpind;
	
	for (int si=0;si<u.num_subsets();++si){
		if (u.num_fct(si)<2){
			if (u.num_fct(si)!=dim+1){
				UG_THROW("Wrong number of components in approximation space. Needed " << dim+1 << ", given " << u.num_fct(si) << ".\n");
			}
		}
		for (int d=0;d<dim;d++){
			if (u.local_finite_element_id(d) != LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Crouzeix-Raviart type.");
			}
		}
		if (u.local_finite_element_id(dim) != LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0)){
			UG_THROW("Component dim in approximation space must be of piecewise constant type.");
		}
	}

	std::vector<SmartPtr<DoFDistribution> > vDD = approxSpace.dof_distributions();
	for(size_t i = 0; i < vDD.size(); ++i){
			std::vector<size_t> newIndex;
			cr_get_connections<elem_type,TGridFunction>(vvConnection,minpind,*vDD[i],u);
			CRMinimumDegree(newIndex,vvConnection,dim,minpind,bseparate,orderp,orderv);
			vDD[i]->permute_indices(newIndex);
	}
}

// original CR cuthill McKee without boost (faster)
void ComputeCRCuthillMcKeeOrder(std::vector<size_t>& vNewIndex,
                              std::vector<std::vector<size_t> >& vvConnections,
                              size_t minpind,
                              bool bReverse);

/// orders the all DofDistributions of the ApproximationSpace using Cuthill-McKee
template <typename TDomain,typename TGridFunction>
void OrderCRCuthillMcKee(ApproximationSpace<TDomain>& approxSpace,TGridFunction& u,
                       bool bReverse)
{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	world dimension
	static const int dim = domain_type::dim;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	for (int si=0;si<u.num_subsets();++si){
		if (u.num_fct(si)<2){
			if (u.num_fct(si)!=dim+1){
				UG_THROW("Wrong number of components in approximation space. Needed " << dim+1 << ", given " << u.num_fct(si) << ".\n");
			}
		}
		for (int d=0;d<dim;d++){
			if (u.local_finite_element_id(d) != LFEID(LFEID::CROUZEIX_RAVIART, dim, 1)){
				UG_THROW("Component " << d << " in approximation space must be of Crouzeix-Raviart type.");
			}
		}
		if (u.local_finite_element_id(dim) != LFEID(LFEID::PIECEWISE_CONSTANT, dim, 0)){
			UG_THROW("Component dim in approximation space must be of piecewise constant type.");
		}
	}

	std::vector<SmartPtr<DoFDistribution> > vDD = approxSpace.dof_distributions();
	for(size_t i = 0; i < vDD.size(); ++i){
			//	get adjacency graph
			std::vector<std::vector<size_t> > vvConnection;
			size_t minpind;
			cr_get_connections<elem_type,TGridFunction>(vvConnection,minpind,*vDD[i],u);
			//	get mapping for cuthill-mckee order
			std::vector<size_t> vNewIndex;
			ComputeCRCuthillMcKeeOrder(vNewIndex,vvConnection,minpind,bReverse);
			//	reorder indices
			vDD[i]->permute_indices(vNewIndex);
		}
}

} // end namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FVCR__CR_REORDER__ */
