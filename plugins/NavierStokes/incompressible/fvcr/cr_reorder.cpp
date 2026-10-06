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

#include <vector>
#include "lib_disc/function_spaces/approximation_space.h"
#include "common/profiler/profiler.h"
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "boost/graph/graph_utility.hpp"
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/sloan_ordering.hpp>
#include <boost/graph/king_ordering.hpp>
#include <boost/graph/minimum_degree_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/math/special_functions/round.hpp>
#include "lib_disc/dof_manager/dof_distribution.h"
#include "cr_reorder.h"

namespace ug{
void computeNewIndicesFromPIndexset(std::vector<std::vector<size_t> >& vvConnection,std::vector<size_t>& newIndex,size_t pmin,
		const std::vector<int>& inv_perm,bool bseparate){
	size_t n = vvConnection.size();
	newIndex.resize(n);
	size_t pn = n - pmin;
	size_t undefined = n+1;
	for (size_t i=0;i<n;i++){
		newIndex[i] = undefined;
	}
	// group p indices behind v indices
	if (bseparate==true){
		size_t vindex=0;
		for (size_t i=0;i<pn;i++){
			size_t pind = inv_perm[i]+pmin;
			newIndex[pind]=pmin+i;
			size_t isize = vvConnection[pind].size()-1;
			for (size_t jj=0;jj<isize;jj++){
				size_t j = vvConnection[pind][jj];
				if (newIndex[j]==undefined){
					newIndex[j]=vindex;
					vindex++;
				}
			}
		}
	// mixed ordering, p indices inserted directly behind associated v nodes
	} else {
		size_t index=0;
		for (size_t i=0;i<pn;i++){
			size_t pind = inv_perm[i]+pmin;
			size_t isize = vvConnection[pind].size()-1;
			for (size_t jj=0;jj<isize;jj++){
				size_t j = vvConnection[pind][jj];
				if (newIndex[j]==undefined){
					newIndex[j]=index;
					index++;
				}
			}
			newIndex[pind]=index;
			index++;
		}
	}
}

	
void computeNewIndicesFromVIndexset(std::vector<std::vector<size_t> >& vvConnection,std::vector<size_t>& newIndex,size_t pmin,size_t dim,
		const std::vector<int>& inv_perm,bool bseparate,size_t strategy){
	size_t n = vvConnection.size();
	newIndex.resize(n);
	size_t vn = boost::math::iround((number)1.0/dim*pmin);
	size_t undefined = n+1;
	for (size_t i=0;i<n;i++){
		newIndex[i] = undefined;
	}
	if (bseparate==true){
		size_t pindex=0;
		// group p behind using first given v node
		if (strategy==1){
			for (size_t i=0;i<vn;i++){
				size_t vind = dim*inv_perm[i];
				for (size_t d=0;d<dim;d++) newIndex[vind+d]=dim*i+d;
				size_t isize = vvConnection[vind].size();
				for (size_t jj=isize-2;jj<isize;jj++){
					size_t j = vvConnection[vind][jj];
					if (j<pmin) continue;
					if (newIndex[j]==undefined){
						newIndex[j]=pindex+pmin;
						pindex++;
					}
				}
			}
		}
		// group p behind when all associated v nodes are handled
		if (strategy==2){
			for (size_t i=0;i<vn;i++){
				size_t vind = dim*inv_perm[i];
				for (size_t d=0;d<dim;d++) newIndex[vind+d]=dim*i+d;
				size_t isize = vvConnection[vind].size();
				for (size_t jj=isize-2;jj<isize;jj++){
					size_t j = vvConnection[vind][jj];
					if (j<pmin) continue;
					if (newIndex[j]==undefined){
						bool giveindex=true;
						size_t jsize = vvConnection[j].size()-1;
						for (size_t kk=0;kk<jsize;kk=kk+dim){
							size_t k = vvConnection[j][kk];
							if (newIndex[k]==n+1){
								giveindex=false;
								break;
							}
						}
						if (giveindex==true){
							newIndex[j]=pindex+pmin;
							pindex++;
						}
					}
				}
			}
		}
	// separate==false mix v and p indices, p is inserted when all associated v nodes are handled
	} else {
		size_t index=0;
		for (size_t i=0;i<vn;i++){
			size_t vind = dim*inv_perm[i];
			for (size_t d=0;d<dim;d++) newIndex[vind+d]=index+d;
			index=index+dim;
			size_t isize = vvConnection[vind].size();
			for (size_t jj=isize-2;jj<isize;jj++){
				size_t j = vvConnection[vind][jj];
				if (j<pmin) continue;
				if (newIndex[j]==undefined){
					bool giveindex=true;
					size_t jsize = vvConnection[j].size()-1;
					for (size_t kk=0;kk<jsize;kk=kk+dim){
						size_t k = vvConnection[j][kk];
						if (newIndex[k]==undefined){
							giveindex=false;
							break;
						}
					}
					if (giveindex==true){
						newIndex[j]=index;
						index++;
					}
				}
			}
		}
	}
}

void computeDegree(std::vector<size_t>& degree,std::vector<std::vector<size_t> >& vvConnections,size_t minpind){
	size_t n=vvConnections.size();
	degree.resize(n);
	for (size_t i=0;i<minpind;i++){
		degree[i]=vvConnections[i].size();
	}
	for (size_t i=minpind;i<n;i++){
		degree[i]=0;
		for (size_t j=0;j<vvConnections[i].size();j++){
			degree[i]+=degree[j];
		}
	}
}

/// orders the dof distribution using Cuthill-McKee
void CRCuthillMcKee(std::vector<size_t>& newIndex,std::vector<std::vector<size_t> >& vvConnection,size_t dim,
						 size_t minpind,bool bReverse,bool bseparate,bool orderp,bool orderv)
{
	PROFILE_FUNC();
		
	size_t n = vvConnection.size();
	size_t pnr = n-minpind;
	size_t vnr = boost::math::iround((number)minpind/2.0);
		
	typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
	boost::property<boost::vertex_color_t, boost::default_color_type,
	boost::property<boost::vertex_degree_t,int> > >
	Graph;
	std::vector<int> inv_permP(0);
	std::vector<int> inv_permV(0);
	newIndex.clear();
	newIndex.resize(n);

	if (orderp==true){
		Graph pGraph(pnr);
		inv_permP.resize(pnr);
		extractPGraph(pGraph,vvConnection,minpind,dim,false);
		if (bReverse==true)
			cuthill_mckee_ordering(pGraph, inv_permP.rbegin(), boost::get(boost::vertex_color, pGraph),
							   boost::make_degree_map(pGraph));
		else
			cuthill_mckee_ordering(pGraph, inv_permP.begin(), boost::get(boost::vertex_color, pGraph),
								   boost::make_degree_map(pGraph));
		if (orderv==false)
			computeNewIndicesFromPIndexset(vvConnection,newIndex,minpind,inv_permP,bseparate);
	}
	
	if (orderv==true){
		Graph vGraph(vnr);
		inv_permV.resize(vnr);
		extractVGraph(vGraph,vvConnection,minpind,dim,false);
		if (bReverse==true)
			cuthill_mckee_ordering(vGraph, inv_permV.rbegin(), boost::get(boost::vertex_color, vGraph),
							   boost::make_degree_map(vGraph));
		else 
			cuthill_mckee_ordering(vGraph, inv_permV.begin(), boost::get(boost::vertex_color, vGraph),
								   boost::make_degree_map(vGraph));
		size_t strategy = 1;
		if (orderp==false)
			computeNewIndicesFromVIndexset(vvConnection,newIndex,minpind,dim,inv_permV,bseparate,strategy);
	}

	if ((orderp==true)&&(orderv==true)){
		// order from p and v index set
		for (size_t i=0;i<(number)1.0/dim*minpind;i++){
			for (size_t d=0;d<dim;d++) newIndex[dim*inv_permV[i]+d]=dim*i+d;
		}
		for (size_t i=0;i<pnr;i++){
			newIndex[inv_permP[i]+minpind]=i+minpind;
		}
	}
}

/// orders the dof distribution using Sloan
void CRSloan(std::vector<size_t>& newIndex,std::vector<std::vector<size_t> >& vvConnection,size_t dim,
						 size_t minpind,bool bseparate,bool orderp,bool orderv)
{
	PROFILE_FUNC();

	size_t n = vvConnection.size();
	size_t pnr = n-minpind;
	size_t vnr = boost::math::iround((number)minpind/2.0);

	newIndex.clear();
	newIndex.resize(n);

	typedef boost::adjacency_list<
		boost::setS,
		boost::vecS,
		boost::undirectedS,
		boost::property<
		boost::vertex_color_t,
		boost::default_color_type,
		boost::property<
		boost::vertex_degree_t,
		int,
		boost::property<
		boost::vertex_priority_t,
		double > > > > Graph;

	typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef boost::property_map<Graph,boost::vertex_degree_t>::type deg_type;
	typedef boost::property_map<Graph, boost::vertex_index_t>::type map_type;

	std::vector<int> inv_permP(0);
	std::vector<int> inv_permV(0);

	if (orderp==true){
		Graph G(pnr);
		inv_permP.resize(pnr);
		extractPGraph(G,vvConnection,minpind,dim,false);

		//Creating two iterators over the vertices
		boost::graph_traits<Graph>::vertex_iterator ui, ui_end;

		//Creating a property_map with the degrees of the degrees of each vertex
		deg_type deg = get(boost::vertex_degree, G);
		for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
			deg[*ui] = degree(*ui, G);

		//Creating a property_map for the indices of a vertex
		//map_type index_map = get(boost::vertex_index, G);

		//Setting the start node
		Vertex s = vertex(0, G);
		int ecc;   //defining a variable for the pseudoperipheral radius

		//Calculating the pseudoeperipheral node and radius
		Vertex e = pseudo_peripheral_pair(G, s, ecc, boost::get(boost::vertex_color, G), boost::get(boost::vertex_degree, G) );

		//Sloan ordering
		sloan_ordering(G, s, e, inv_permP.begin(), get(boost::vertex_color, G),
				   get(boost::vertex_degree, G), get(boost::vertex_priority, G));

		if (orderv==false)
			computeNewIndicesFromPIndexset(vvConnection,newIndex,minpind,inv_permP,bseparate);
	}

	if (orderv==true){
		Graph G(vnr);
		inv_permV.resize(vnr);
		extractVGraph(G,vvConnection,minpind,dim,false);

		//Creating two iterators over the vertices
		boost::graph_traits<Graph>::vertex_iterator ui, ui_end;

		//Creating a property_map with the degrees of the degrees of each vertex
		deg_type deg = get(boost::vertex_degree, G);
		for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
			deg[*ui] = degree(*ui, G);

		//Creating a property_map for the indices of a vertex
		//map_type index_map = get(boost::vertex_index, G);

		//Setting the start node
		Vertex s = vertex(0, G);
		int ecc;   //defining a variable for the pseudoperipheral radius

		//Calculating the pseudoeperipheral node and radius
		Vertex e = pseudo_peripheral_pair(G, s, ecc, boost::get(boost::vertex_color, G), boost::get(boost::vertex_degree, G) );

		//Sloan ordering
		sloan_ordering(G, s, e, inv_permV.begin(), get(boost::vertex_color, G),
				   get(boost::vertex_degree, G), get(boost::vertex_priority, G));

		size_t strategy = 1;
		if (orderp==false)
			computeNewIndicesFromVIndexset(vvConnection,newIndex,minpind,dim,inv_permV,bseparate,strategy);
	}

	if ((orderp==true)&&(orderv==true)){
		// order from p and v index set
		for (size_t i=0;i<(number)1.0/dim*minpind;i++){
			for (size_t d=0;d<dim;d++) newIndex[dim*inv_permV[i]+d]=dim*i+d;
		}
		for (size_t i=0;i<pnr;i++){
			newIndex[inv_permP[i]+minpind]=i+minpind;
		}
	}
}

/// orders the dof distribution using King
void CRKing(std::vector<size_t>& newIndex,std::vector<std::vector<size_t> >& vvConnection,size_t dim,
						 size_t minpind,bool bReverse,bool bseparate,bool orderp,bool orderv)
{
	PROFILE_FUNC();

	size_t n = vvConnection.size();
	size_t pnr = n-minpind;
	size_t vnr = boost::math::iround((number)minpind/2.0);

	newIndex.clear();
	newIndex.resize(n);

	typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
		boost::property<boost::vertex_color_t, boost::default_color_type,
		boost::property<boost::vertex_degree_t,int> > >
	Graph;

	std::vector<int> inv_permP(0);
	std::vector<int> inv_permV(0);

	if (orderp==true){
		Graph pGraph(pnr);
		inv_permP.resize(pnr);
		extractPGraph(pGraph,vvConnection,minpind,dim,false);
		if (bReverse==true)
			king_ordering(pGraph, inv_permP.rbegin());
		else
			king_ordering(pGraph, inv_permP.begin());

		if (orderv==false)
			computeNewIndicesFromPIndexset(vvConnection,newIndex,minpind,inv_permP,bseparate);
	}

	if (orderv==true){
		Graph vGraph(vnr);
		inv_permV.resize(vnr);
		extractVGraph(vGraph,vvConnection,minpind,dim,false);
		if (bReverse==true)
			king_ordering(vGraph, inv_permV.rbegin());
		else
			king_ordering(vGraph, inv_permV.begin());

		size_t strategy = 1;
		if (orderp==false)
			computeNewIndicesFromVIndexset(vvConnection,newIndex,minpind,dim,inv_permV,bseparate,strategy);
	}

	if ((orderp==true)&&(orderv==true)){
		// order from p and v index set
		for (size_t i=0;i<(number)1.0/dim*minpind;i++){
			for (size_t d=0;d<dim;d++) newIndex[dim*inv_permV[i]+d]=dim*i+d;
		}
		for (size_t i=0;i<pnr;i++){
			newIndex[inv_permP[i]+minpind]=i+minpind;
		}
	}

}

/// orders the dof distribution using minimum degree algorithm
void CRMinimumDegree(std::vector<size_t>& newIndex,std::vector<std::vector<size_t> >& vvConnection,size_t dim,
						 size_t minpind,bool bseparate,bool orderp,bool orderv)
{
	PROFILE_FUNC();

	size_t n = vvConnection.size();
	size_t pnr = n-minpind;
	size_t vnr = boost::math::iround((number)minpind/2.0);

	newIndex.clear();
	newIndex.resize(n);

	typedef boost::adjacency_list<boost::setS, boost::vecS, boost::directedS>  Graph;

	typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef boost::property_map<Graph, boost::vertex_index_t>::type map_type;

	std::vector<int> inv_permP(0);
	std::vector<int> inv_permV(0);

	int delta=0;

	if (orderp==true){
		Graph G(pnr);

		inv_permP.resize(pnr,0);
		std::vector<int> perm(pnr,0);
		std::vector<int> degree(pnr,0);
		std::vector<int> supernode_sizes(pnr,1);

		extractPGraph(G,vvConnection,minpind,dim,true);

		// Creating a property_map for the indices of a vertex
		map_type index_map = get(boost::vertex_index, G);

		// boost minimum degree algorithm does not work for 1 edge graph use king instead
		if (pnr<=2)
			king_ordering(G, inv_permP.rbegin());
		else
			minimum_degree_ordering
			(G,
			boost::make_iterator_property_map(&degree[0],index_map, degree[0]),
			perm.begin(),
			inv_permP.begin(),
			boost::make_iterator_property_map(&supernode_sizes[0],index_map, supernode_sizes[0]),
			delta, index_map);

		if (orderv==false)
			computeNewIndicesFromPIndexset(vvConnection,newIndex,minpind,inv_permP,bseparate);
	}

	if (orderv==true){
		Graph G(vnr);
		inv_permV.resize(vnr,0);
		std::vector<int> perm(vnr,0);
		std::vector<int> degree(vnr,0);
		std::vector<int> supernode_sizes(vnr,1);

		extractVGraph(G,vvConnection,minpind,dim,true);

		// Creating a property_map for the indices of a vertex
		map_type index_map = get(boost::vertex_index, G);

		minimum_degree_ordering
		(G,
		boost::make_iterator_property_map(&degree[0],index_map, degree[0]),
		perm.begin(),
		inv_permV.begin(),
		boost::make_iterator_property_map(&supernode_sizes[0],index_map, supernode_sizes[0]),
		delta, index_map);

		size_t strategy = 1;
		if (orderp==false)
			computeNewIndicesFromVIndexset(vvConnection,newIndex,minpind,dim,inv_permV,bseparate,strategy);
	}

	if ((orderp==true)&&(orderv==true)){
		// order from p and v index set
		for (size_t i=0;i<(number)1.0/dim*minpind;i++){
			for (size_t d=0;d<dim;d++) newIndex[dim*inv_permV[i]+d]=dim*i+d;
		}
		for (size_t i=0;i<pnr;i++){
			newIndex[inv_permP[i]+minpind]=i+minpind;
		}
	}
}

// Cuthill-McKee order on Crouzeix-Raviart elements
void ComputeCRCuthillMcKeeOrder(std::vector<size_t>& vNewIndex,
                              std::vector<std::vector<size_t> >& vvConnections,
                              size_t minpind,
                              bool bReverse){
	size_t numInd = vvConnections.size();
	vNewIndex.resize(numInd);
	size_t numP = numInd-minpind;
	std::vector<std::vector<size_t> > pConnections(numP);
	std::vector<size_t> degree(numInd);
	std::vector<size_t>::iterator it;
	// go over all velocity nodes and fill degree vector (number of connections)
	// set up p adjacency graph, two pressure nodes are seen as associated
	// when they have a common associated velocity node
	for (size_t i=0;i<minpind;i=i+2){
		size_t localNumInd=vvConnections[i].size();
		degree[i]=localNumInd;
		degree[i+1]=localNumInd;
		vNewIndex[i]=numInd;
		vNewIndex[i+1]=numInd;
		size_t localNumP=0;
		size_t associatedP[2];
		for (size_t j=0;j<localNumInd;j++){
			if (vvConnections[i][j]>=minpind){
				associatedP[localNumP]=vvConnections[i][j]-minpind;
				localNumP++;
			}
		}
		if (localNumP>1){
			size_t pi = associatedP[0];
			size_t pj = associatedP[1];
			it = find(pConnections[pi].begin(), pConnections[pi].end(), pj);
			if(it == pConnections[pi].end()) pConnections[pi].push_back(pj+minpind);
			//	add the opposite direction (adjInd -> index)
			it = find(pConnections[pj].begin(), pConnections[pj].end(), pi);
			if(it == pConnections[pj].end()) pConnections[pj].push_back(pi+minpind);
		}
	}
	// go over all pressure nodes, fill degree vector (sum of number of connections of all associated velocity nodes)
	// compute minimum degree pressure node
	size_t minpdeg = numInd;
	size_t minpdegind = minpind;
	for (size_t i=minpind;i<numInd;i++){
		degree[i]=0;
		vNewIndex[i]=numInd;
		size_t localNumInd=vvConnections[i].size();
		for (size_t j=0;j<localNumInd;j++){
			if (vvConnections[i][j]!=i) degree[i]+=degree[vvConnections[i][j]];
		}
		if (degree[i]<minpdeg){
			minpdeg=degree[i];
			minpdegind=i;
		}
	}
	// set up list of pressure nodes, first entry is minimum degree index
	std::vector<size_t> plist(1);
	plist[0]=minpdegind;
	size_t count=0;
	//	Sort neighbours by degreee
	CompareDeg myCompDegree(degree);
	// go over pressure nodes and sort
	for (size_t j=0;j<numP;j++){
		size_t i=plist[j];
		// sort adjacent indices by degree (remark: by definition the degree of velocity indices is always lower than
		// the associated pressure index degree)
		std::sort(	vvConnections[i].begin(),
				          	vvConnections[i].end(), myCompDegree);
		// assign new index
		for (size_t k=0;k<vvConnections[i].size();k++){
			if (vNewIndex[vvConnections[i][k]]<numInd) continue;
			vNewIndex[vvConnections[i][k]]=count;
			count++;
		}
		size_t pind=i-minpind;
		// sort associated p nodes (using p adjacency list)
		std::sort(	pConnections[pind].begin(),
	          	pConnections[pind].end(), myCompDegree);
		// check if already in p list, if not insert
		for (size_t k=0;k<pConnections[pind].size();k++){
			if (vNewIndex[pConnections[pind][k]]<numInd) continue;
			bool found=false;
			for (size_t iiiiii=j;iiiiii<plist.size();iiiiii++){
				if (plist[iiiiii]==pConnections[pind][k]){
					found=true;
					break;
				}
			}
			if (found==false){
				plist.push_back(pConnections[pind][k]);
			}
		}
	}
}

} // end namespace ug
