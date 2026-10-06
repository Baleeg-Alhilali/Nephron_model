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

#ifndef __H__UG__NAVIER_STOKES__INCOMPRESSIBLE__FILTER__IMPL__
#define __H__UG__NAVIER_STOKES__INCOMPRESSIBLE__FILTER__IMPL__

#include "common/util/provider.h"

namespace ug{
namespace NavierStokes{

template<typename VType,typename TElem,typename TDomain,typename TGridFunction>
void averageByVolume(PeriodicAttachmentAccessor<TElem,Attachment<VType> >& acUHat,
					 PeriodicAttachmentAccessor<TElem,Attachment<number> >& acVolume,TDomain& domain,SmartPtr<TGridFunction> uInfo,std::vector<WallObject<TGridFunction> > walls){
	typedef typename TDomain::grid_type TGrid;
	typedef typename TGridFunction::template traits<TElem>::const_iterator ElemIterator;
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		for (int j=0;j<(int)walls.size();j++) if ((int)si == (int)walls[j].si()) continue;
		ElemIterator elemIter = uInfo->template begin<TElem>(si);
		ElemIterator elemIterEnd = uInfo->template end<TElem>(si);
		for(  ;elemIter !=elemIterEnd; elemIter++)
		{
			TElem* elem = *elemIter;
			if (pbm && pbm->is_slave(elem)) continue;
			acUHat[elem]/=(number)acVolume[elem];
		}
	}
}

template <typename vector_t>
inline
typename vector_t::value_type
VecDistanceMaxNorm(const vector_t& v1, const vector_t& v2)
{
	typename vector_t::value_type m = std::abs(v1[0]-v2[0]);
	typedef typename vector_t::size_type size_type;

	for(size_type i = 1; i < v1.size(); ++i)
	{
		m = std::max(m, std::abs(v1[i]-v2[i]) );
	}
	return m;
}
	
// compute average element size of grid
template<typename TGridFunction>
number ConstantBoxFilter<TGridFunction>::compute_average_element_size(SmartPtr<TGridFunction> u){
	typedef typename TGridFunction::const_element_iterator const_iterator;
	typedef typename TGridFunction::element_type element_type;
	typename TGridFunction::domain_type::position_accessor_type& aaPos
					= u->domain()->position_accessor();
	number averageElemSize = 0;
	std::vector<MathVector<dim> > vCorner;
	size_t nrOfElem = 0;
	const_iterator iter = m_uInfo->template begin<elem_type>();
	const_iterator iterEnd = m_uInfo->template end<elem_type>();
	for(; iter != iterEnd; ++iter)
	{
		elem_type* elem = *iter;
		ReferenceObjectID roid = elem->reference_object_id();
		CollectCornerCoordinates(vCorner, *elem, aaPos);
		averageElemSize += ElementSize<dim>(roid, &vCorner[0]);
		nrOfElem++;
	}
	return averageElemSize/(number)nrOfElem;
}

template <typename TGridFunction>
template <typename VType>
void ConstantBoxFilter<TGridFunction>::collectSides(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acUHat,
														PeriodicAttachmentAccessor<side_type,Attachment<number> >& acVolume,
														std::vector< MathVector<dim> >& coord,
														VType values[DimFV1Geometry<dim>::maxNumSCV],
														number volumes[DimFV1Geometry<dim>::maxNumSCV],
														std::vector<side_type*>& sides
														){
	domain_type& domain = *m_uInfo->domain().get();
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();
	side_secure_container assoSides;
	size_t startIndex = sides.size()-1;
	// handle start side
	side_type* s = sides[startIndex];
	MathVector<dim> bary = posAcc[s->vertex(0)];
	for (size_t k=1;k<s->num_vertices();k++){
		bary += posAcc[s->vertex(k)];
	}
	bary /= s->num_vertices();
	number width=0.5*m_width;
	for (size_t w=0;w<m_walls.size();w++){
		width = std::min(width,m_walls[w].dist(bary));
	}
	for (size_t k=0;k<coord.size();k++){
		number dist = VecDistanceMaxNorm(coord[k],bary);
		if (dist<width){ 
			if ((!pbm)||(pbm && !pbm->is_slave(s))){
				acVolume[s]+=volumes[k];
				acUHat[s]+=values[k];
			}
		}
	}
	for (size_t i=startIndex;i<sides.size();i++){
		for (size_t v=0;v<sides[i]->num_vertices();v++){
			domain.grid()->associated_elements(assoSides,sides[i]->vertex(v));
			for (size_t j=0;j<assoSides.size();j++){
				s = assoSides[j];
				// check if new side
				bool isNew = true;
				for (int k=sides.size()-1;k>=0;k--){
					if (sides[k]==s){
						isNew=false;
						break;
					}
				}
				if (isNew==false) continue;
				// compute dof coord
				bary = posAcc[s->vertex(0)];
				for (size_t k=1;k<s->num_vertices();k++){
					bary += posAcc[s->vertex(k)];
				}
				bary /= s->num_vertices();
				number width=0.5*m_width;
				for (size_t w=0;w<m_walls.size();w++){
					width = std::min(width,m_walls[w].dist(bary));
				}
				// check if inside filter region
				bool inside=false;
				for (size_t k=0;k<coord.size();k++){
					number dist = VecDistanceMaxNorm(coord[k],bary);
					if (dist<0.5*m_width) inside=true;
					if (dist<width){
						if ((!pbm)||(pbm && !pbm->is_slave(s))){
							acVolume[s]+=volumes[k];
							acUHat[s]+=values[k];
						}
					} 
					
				}
				// if inside add to list
				if (inside==true) sides.push_back(s);
			}
		}
	}
}
	
// constant box filter for side data
template <typename TGridFunction>
template <typename VType>
void ConstantBoxFilter<TGridFunction>::apply_filter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acUHat,
					  SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acU){
		bool useGridFunction;
		if (u==SPNULL) useGridFunction = false; else useGridFunction = true;
		// set attachments to zero
		domain_type& domain = *m_uInfo->domain().get();
		SetAttachmentValues(m_acSideVolume, domain.grid()->template begin<side_type>(), domain.grid()->template end<side_type>(), 0);
		SetAttachmentValues(acUHat, domain.grid()->template begin<side_type>(), domain.grid()->template end<side_type>(), 0);

		std::vector<DoFIndex> multInd;

		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];

		//	get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

		typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
		side_secure_container elemsides;

		DimFV1Geometry<dim> geo;
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
		{
			//	get iterators
			ElemIterator iter = m_uInfo->template begin<elem_type>(si);
			ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

			//	loop elements of dimension
			for(  ;iter !=iterEnd; ++iter)
			{
				//	get Elem
				elem_type* elem = *iter;

				//	get corners of element
				size_t noc=elem->num_vertices();
				for(size_t i = 0; i < noc; ++i)
					coCoord[i] = posAcc[elem->vertex(i)];

				// get sides
				domain.grid()->associated_elements_sorted(elemsides, static_cast<elem_type*>(elem) );

				//	evaluate finite volume geometry
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

				//	reference object id
				ReferenceObjectID roid = elem->reference_object_id();

				//	get trial space
				const LocalShapeFunctionSet<dim>& rTrialSpace =
				LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART,dim, 1));

				VType values[DimFV1Geometry<dim>::maxNumSCV];
				number volumes[DimFV1Geometry<dim>::maxNumSCV];
				std::vector< MathVector<dim> > baryCoords(geo.num_scv());
				std::vector<side_type*> sides;
				sides.reserve(20);
				std::vector< Region<dim> > regions(1);
				regions[0].firstIndex=0;
				regions[0].periodicOffset=0;


				for (size_t i=0;i<geo.num_scv();i++){
					typename DimFV1Geometry<dim>::SCV scv = geo.scv(i);
					baryCoords.resize(geo.num_scv());
					MathVector<dim> scvLocalBary = scv.local_corner(0);
					baryCoords[i] = scv.global_corner(0);
					for (size_t co=1;co<scv.num_corners();co++){
						scvLocalBary+=scv.local_corner(co);
						baryCoords[i]+=scv.global_corner(co);
					}
					scvLocalBary/=scv.num_corners();
					baryCoords[i]/=scv.num_corners();

					std::vector<number> vShape;
					rTrialSpace.shapes(vShape, scvLocalBary);
				   //MathVector<dim> localValue = 0;

					size_t nos=elemsides.size();
					// sum up shapes
					values[i] = 0;
					for (size_t sh=0;sh<nos;sh++){
						VType localValue;
						if (useGridFunction)
						for (size_t d=0;d<dim;d++){
							u->dof_indices(elemsides[sh], d, multInd);
							assignVal(localValue,d,DoFRef(*u,multInd[0]));
						}
						else
							localValue = acU[elemsides[sh]];
						localValue *= vShape[sh];
						values[i] += localValue;
					}
					values[i] *= scv.volume();
					volumes[i] = scv.volume();
				}
				// check neighborhood for sides in the filter region
				sides.push_back(elemsides[0]);
				collectSides(acUHat,m_acSideVolume,baryCoords,values,volumes,sides);
				PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
				// check for periodic sides and add corresponding side
				if (pbm){
					for (size_t i=0;i<sides.size();i++){
						side_type* s = sides[i];
						if ( (pbm->is_slave(s)) || (pbm->is_master(s)) ) {
							side_type* sNew;
							if (pbm->is_slave(s))
								sNew = pbm->master(s);
							else
								sNew = (*pbm->slaves(s))[0];
							bool isNew = true;
							for (size_t j=0;j<sides.size();j++){
								if (sides[j]==sNew){
									isNew=false;
									break;
								}
							}
							if (isNew==true){
								sides.push_back(sNew);
								// identify base region
								size_t regInd;
								if (regions.size()==0)
									regInd=0;
								else{
									for (regInd=regions.size()-1; ;regInd--){
										if (i>regions[regInd].firstIndex) break;
										if (regInd==0) break;
									}
								}
								// compute periodic offset
								MathVector<dim> sReg0Co = posAcc[s->vertex(0)];
								MathVector<dim> sReg1Co = posAcc[sNew->vertex(0)];
								for (size_t k=1;k<sNew->num_vertices();k++){
									sReg0Co += posAcc[s->vertex(k)];
									sReg1Co += posAcc[sNew->vertex(k)];
								}
								sReg0Co/=s->num_vertices();
								sReg1Co/=sNew->num_vertices();
								MathVector<dim> pOffset;
								for (size_t k=0;k<dim;k++){
									pOffset[k] = regions[regInd].periodicOffset[k]+sReg0Co[k]-sReg1Co[k];
								}
								// add region
								size_t rsize = regions.size();
								regions.resize(rsize+1);
								regions[rsize].firstIndex=sides.size()-1;
								regions[rsize].periodicOffset=pOffset;
								std::vector< MathVector<dim> > baryCoordsNew(geo.num_scv());
								for (size_t k=0;k<baryCoordsNew.size();k++){
									for (size_t d0=0;d0<dim;d0++){
										baryCoordsNew[k][d0] = baryCoords[k][d0]-pOffset[d0];
									}
								}
								// find associated sides
								collectSides(acUHat,m_acSideVolume,baryCoordsNew,values,volumes,sides);
							}
						}
					}
				}
			}
		}
		averageByVolume(acUHat,m_acSideVolume,domain,m_uInfo,m_walls);
		if (useGridFunction) copyWallData(acUHat,u,m_walls);
		else copyWallData(acUHat,acU,m_uInfo,m_walls);
		if (dim==2) constrainingSideAveraging<TGridFunction,side_type,ConstrainingEdge,VType>(acUHat,m_uInfo);
		else {
			constrainingSideAveraging<TGridFunction,side_type,ConstrainingTriangle,VType>(acUHat,m_uInfo);
			constrainingSideAveraging<TGridFunction,side_type,ConstrainingQuadrilateral,VType>(acUHat,m_uInfo);
		}
};
	
// compute average element size of grid
template<typename TGridFunction>
number VariableBoxFilter<TGridFunction>::compute_average_element_size(SmartPtr<TGridFunction> u){
	typedef typename TGridFunction::const_element_iterator const_iterator;
	typedef typename TGridFunction::element_type element_type;
	typename TGridFunction::domain_type::position_accessor_type& aaPos
	= u->domain()->position_accessor();
	number averageElemSize = 0;
	std::vector<MathVector<dim> > vCorner;
	size_t nrOfElem = 0;
	const_iterator iter = m_uInfo->template begin<elem_type>();
	const_iterator iterEnd = m_uInfo->template end<elem_type>();
	for(; iter != iterEnd; ++iter)
	{
		elem_type* elem = *iter;
		ReferenceObjectID roid = elem->reference_object_id();
		CollectCornerCoordinates(vCorner, *elem, aaPos);
		averageElemSize += ElementSize<dim>(roid, &vCorner[0]);
		nrOfElem++;
	}
	return averageElemSize/(number)nrOfElem;
}
	
template <typename TGridFunction>
template <typename VType>
void VariableBoxFilter<TGridFunction>::collectSides(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acUHat,
														PeriodicAttachmentAccessor<side_type,Attachment<number> >& acVolume,
														std::vector< MathVector<dim> >& coord,
														VType values[DimFV1Geometry<dim>::maxNumSCV],
														number volumes[DimFV1Geometry<dim>::maxNumSCV],
														std::vector<side_type*>& sides
														){
	domain_type& domain = *m_uInfo->domain().get();
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
	//	get position accessor
	std::vector<DoFIndex> multInd; 
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();
	side_secure_container assoSides;
	size_t startIndex = sides.size()-1;
	// handle start side
	side_type* s = sides[startIndex];
	MathVector<dim> bary = posAcc[s->vertex(0)];
	for (size_t k=1;k<s->num_vertices();k++){
		bary += posAcc[s->vertex(k)];
	}
	bary /= s->num_vertices();
	m_width->inner_dof_indices(s,0, multInd);
	number localwidth=0.5*DoFRef(*m_width,multInd[0]);
	for (size_t w=0;w<m_walls.size();w++){
		localwidth = std::min(localwidth,m_walls[w].dist(bary));
	}
	for (size_t k=0;k<coord.size();k++){
		number dist = VecDistanceMaxNorm(coord[k],bary);
		if (dist<0.5*localwidth){ 
			if ((!pbm)||(pbm && !pbm->is_slave(s))){
				acVolume[s]+=volumes[k];
				acUHat[s]+=values[k];
			}
		}
	}
	for (size_t i=startIndex;i<sides.size();i++){
		for (size_t v=0;v<sides[i]->num_vertices();v++){
			domain.grid()->associated_elements(assoSides,sides[i]->vertex(v));
			for (size_t j=0;j<assoSides.size();j++){
				s = assoSides[j];
				// check if new side
				bool isNew = true;
				for (int k=sides.size()-1;k>=0;k--){
					if (sides[k]==s){
						isNew=false;
						break;
					}
				}
				if (isNew==false) continue;
				// compute dof coord
				bary = posAcc[s->vertex(0)];
				for (size_t k=1;k<s->num_vertices();k++){
					bary += posAcc[s->vertex(k)];
				}
				bary /= s->num_vertices();
				// get local filterwidth
				m_width->inner_dof_indices(s,0, multInd);
				localwidth=0.5*DoFRef(*m_width,multInd[0]);
				for (size_t w=0;w<m_walls.size();w++){
					localwidth = std::min(localwidth,m_walls[w].dist(bary));
				}
				// check if inside filter region
				bool inside=false;
				for (size_t k=0;k<coord.size();k++){
					number dist = VecDistanceMaxNorm(coord[k],bary);
					if (dist<0.5*localwidth){ 
						if ((!pbm)||(pbm && !pbm->is_slave(s))){
							acVolume[s]+=volumes[k];
							acUHat[s]+=values[k];
						}
					} 
					if (dist<0.5*m_maxwidth) inside=true;
				}
				// if inside add to list
				if (inside==true) sides.push_back(s);
			}
		}
	}
}
	
// variable box filter for side data
template <typename TGridFunction>
template <typename VType>
void VariableBoxFilter<TGridFunction>::apply_filter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acUHat,
													SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acU){
	bool useGridFunction;
	if (u==SPNULL) useGridFunction = false; else useGridFunction = true;
	// set attachments to zero
	domain_type& domain = *m_uInfo->domain().get();
	SetAttachmentValues(m_acSideVolume, domain.grid()->template begin<side_type>(), domain.grid()->template end<side_type>(), 0);
	SetAttachmentValues(acUHat, domain.grid()->template begin<side_type>(), domain.grid()->template end<side_type>(), 0);
		
	std::vector<DoFIndex> multInd;
		
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		
	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();
		
	typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
	side_secure_container elemsides;
		
	DimFV1Geometry<dim> geo;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);
			
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
				
			//	get corners of element
			size_t noc=elem->num_vertices();
			for(size_t i = 0; i < noc; ++i)
				coCoord[i] = posAcc[elem->vertex(i)];
				
			// get sides
			domain.grid()->associated_elements_sorted(elemsides, static_cast<elem_type*>(elem) );
				
			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				
			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();
				
			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART,dim, 1));
				
			VType values[DimFV1Geometry<dim>::maxNumSCV];
			number volumes[DimFV1Geometry<dim>::maxNumSCV];
			std::vector< MathVector<dim> > baryCoords(geo.num_scv());
			std::vector<side_type*> sides;
			sides.reserve(20);
			std::vector< Region<dim> > regions(1);
			regions[0].firstIndex=0;
			regions[0].periodicOffset=0;
				
				
			for (size_t i=0;i<geo.num_scv();i++){
				typename DimFV1Geometry<dim>::SCV scv = geo.scv(i);
				baryCoords.resize(geo.num_scv());
				MathVector<dim> scvLocalBary = scv.local_corner(0);
				baryCoords[i] = scv.global_corner(0);
				for (size_t co=1;co<scv.num_corners();co++){
					scvLocalBary+=scv.local_corner(co);
					baryCoords[i]+=scv.global_corner(co);
				}
				scvLocalBary/=scv.num_corners();
				baryCoords[i]/=scv.num_corners();
					
				std::vector<number> vShape;
				rTrialSpace.shapes(vShape, scvLocalBary);
				//MathVector<dim> localValue = 0;
					
				size_t nos=elemsides.size();
				// sum up shapes
				values[i] = 0;
				for (size_t sh=0;sh<nos;sh++){
					VType localValue;
					if (useGridFunction)
						for (size_t d=0;d<dim;d++){
							u->dof_indices(elemsides[sh], d, multInd);
							assignVal(localValue,d,DoFRef(*u,multInd[0]));
						}
					else
						localValue = acU[elemsides[sh]];
					localValue *= vShape[sh];
					values[i] += localValue;
				}
				values[i] *= scv.volume();
				volumes[i] = scv.volume();
			}
			// check neighborhood for sides in the filter region
			sides.push_back(elemsides[0]);
			collectSides(acUHat,m_acSideVolume,baryCoords,values,volumes,sides);
			PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
			// check for periodic sides and add corresponding side
			if (pbm){
				for (size_t i=0;i<sides.size();i++){
					side_type* s = sides[i];
					if ( (pbm->is_slave(s)) || (pbm->is_master(s)) ) {
						side_type* sNew;
						if (pbm->is_slave(s))
							sNew = pbm->master(s);
						else
							sNew = (*pbm->slaves(s))[0];
						bool isNew = true;
						for (size_t j=0;j<sides.size();j++){
							if (sides[j]==sNew){
								isNew=false;
								break;
							}
						}
						if (isNew==true){
							sides.push_back(sNew);
							// identify base region
							size_t regInd;
							if (regions.size()==0)
								regInd=0;
							else{
								for (regInd=regions.size()-1; ;regInd--){
									if (i>regions[regInd].firstIndex) break;
									if (regInd==0) break;
								}
							}
							// compute periodic offset
							MathVector<dim> sReg0Co = posAcc[s->vertex(0)];
							MathVector<dim> sReg1Co = posAcc[sNew->vertex(0)];
							for (size_t k=1;k<sNew->num_vertices();k++){
								sReg0Co += posAcc[s->vertex(k)];
								sReg1Co += posAcc[sNew->vertex(k)];
							}
							sReg0Co/=s->num_vertices();
							sReg1Co/=sNew->num_vertices();
							MathVector<dim> pOffset;
							for (size_t k=0;k<dim;k++){
								pOffset[k] = regions[regInd].periodicOffset[k]+sReg0Co[k]-sReg1Co[k];
							}
							// add region
							size_t rsize = regions.size();
							regions.resize(rsize+1);
							regions[rsize].firstIndex=sides.size()-1;
							regions[rsize].periodicOffset=pOffset;
							std::vector< MathVector<dim> > baryCoordsNew(geo.num_scv());
							for (size_t k=0;k<baryCoordsNew.size();k++){
								for (size_t d0=0;d0<dim;d0++){
									baryCoordsNew[k][d0] = baryCoords[k][d0]-pOffset[d0];
								}
							}
							// find associated sides
							collectSides(acUHat,m_acSideVolume,baryCoordsNew,values,volumes,sides);
						}
					}
				}
			}
		}
	}
	averageByVolume(acUHat,m_acSideVolume,domain,m_uInfo,m_walls);
	if (useGridFunction) copyWallData(acUHat,u,m_walls);
	else copyWallData(acUHat,acU,m_uInfo,m_walls);
	if (dim==2) constrainingSideAveraging<TGridFunction,side_type,ConstrainingEdge,VType>(acUHat,m_uInfo);
	else {
		constrainingSideAveraging<TGridFunction,side_type,ConstrainingTriangle,VType>(acUHat,m_uInfo);
		constrainingSideAveraging<TGridFunction,side_type,ConstrainingQuadrilateral,VType>(acUHat,m_uInfo);
	}
};

// fv1 Filter for vertex data
template <typename TGridFunction>
template <typename VType>
void FV1BoxFilter<TGridFunction>::apply_filter(PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& acUHat,
			   SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& acU){
	bool useGridFunction = true;
	if (u==SPNULL) useGridFunction = false;
	
	// set attachments to zero
	domain_type& domain = *m_uInfo->domain().get();
	SetAttachmentValues(m_acVertexVolume, domain.grid()->template begin<vertex_type>(), domain.grid()->template end<vertex_type>(), 0);
	SetAttachmentValues(acUHat, domain.grid()->template begin<vertex_type>(), domain.grid()->template end<vertex_type>(), 0);
	
	std::vector<DoFIndex> multInd;

	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	DimFV1Geometry<dim> geo;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;

			//	get corners of element
			size_t noc=elem->num_vertices();
			for(size_t i = 0; i < noc; ++i)
				coCoord[i] = posAcc[elem->vertex(i)];

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE,dim, 1));

			for (size_t i=0;i<geo.num_scv();i++){
				typename DimFV1Geometry<dim>::SCV scv = geo.scv(i);
				MathVector<dim> scvLocalBary = scv.local_corner(0);
				for (size_t co=1;co<scv.num_corners();co++){
					scvLocalBary+=scv.local_corner(co);
				}
				scvLocalBary/=scv.num_corners();
				size_t node = scv.node_id();

				std::vector<number> vShape;
				rTrialSpace.shapes(vShape, scvLocalBary);
				// sum up shapes
				VType value;
				value = 0;
				for (size_t sh=0;sh<noc;sh++){
					VType localValue;
					if (useGridFunction)
						for (size_t d=0;d<dim;d++){
							u->dof_indices(elem->vertex(sh), d, multInd);
							assignVal(localValue,d,DoFRef(*u,multInd[0]));
						}
					else
						localValue = acU[elem->vertex(sh)];
					localValue *= vShape[sh];
					value += localValue;
				}
				value *= scv.volume();
				m_acVertexVolume[elem->vertex(node)]  += scv.volume();
				acUHat[elem->vertex(node)] += value;
			}
		}
	}
	averageByVolume(acUHat,m_acVertexVolume,domain,m_uInfo,m_walls);
	if (useGridFunction) copyWallData(acUHat,u,m_walls);
	else copyWallData(acUHat,acU,m_uInfo,m_walls);
}
	
// fv1 Filter for vertex data
template <typename TGridFunction>
void FV1BoxFilter<TGridFunction>::compute_filterwidth_fv1()
{
	// set attachments to zero
	domain_type& domain = *m_uInfo->domain().get();
	SetAttachmentValues(m_acVertexVolume, domain.grid()->template begin<vertex_type>(), domain.grid()->template end<vertex_type>(), 0);
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		
	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();
		
	DimFV1Geometry<dim> geo;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);
			
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
				
			//	get corners of element
			size_t noc=elem->num_vertices();
			for(size_t i = 0; i < noc; ++i)
				coCoord[i] = posAcc[elem->vertex(i)];
				
			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				
			for (size_t i=0;i<geo.num_scv();i++){
				typename DimFV1Geometry<dim>::SCV scv = geo.scv(i);
				size_t node = scv.node_id();
				m_acVertexVolume[elem->vertex(node)]  += scv.volume();
			}
		}
	}
}
	
// fv1 Filter for side data
template <typename TGridFunction>
void FV1BoxFilter<TGridFunction>::compute_filterwidth_fvcr()
{
	domain_type& domain = *m_uInfo->domain().get();
	SetAttachmentValues(m_acSideVolume, domain.grid()->template begin<side_type>(), domain.grid()->template end<side_type>(), 0);
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		
	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();
		
	typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
	side_secure_container sides;
		
	DimFV1Geometry<dim> geo;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);
			
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
				
			//	get corners of element
			size_t noc=elem->num_vertices();
			for(size_t i = 0; i < noc; ++i)
				coCoord[i] = posAcc[elem->vertex(i)];
				
			// get sides
			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );
				
			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
			
			for (size_t i=0;i<geo.num_scv();i++){
				typename DimFV1Geometry<dim>::SCV scv = geo.scv(i);
				size_t node = scv.node_id();
				// find sides associated with scv
				for (size_t j=0;j<sides.size();j++){
					side_type* sidej = sides[j];
					for (size_t sco=0;sco<sidej->num_vertices();sco++){
						if (VecDistance(coCoord[node],posAcc[sidej->vertex(sco)])<1e-12){
							m_acSideVolume[sidej]  += scv.volume();
						}
					}
				}
			}
		}
	}
}
	
	
// fv1 Filter for side data
template <typename TGridFunction>
template <typename VType>
void FV1BoxFilter<TGridFunction>::apply_filter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acUHat,
			   SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acU){
	bool useGridFunction = true;
	if (u==SPNULL) useGridFunction = false;
	
	// set attachments to zero
	domain_type& domain = *m_uInfo->domain().get();
	SetAttachmentValues(m_acSideVolume, domain.grid()->template begin<side_type>(), domain.grid()->template end<side_type>(), 0);
	SetAttachmentValues(acUHat, domain.grid()->template begin<side_type>(), domain.grid()->template end<side_type>(), 0);
	
	std::vector<DoFIndex> multInd;

	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
	side_secure_container sides;

	DimFV1Geometry<dim> geo;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;

			//	get corners of element
			size_t noc=elem->num_vertices();
			for(size_t i = 0; i < noc; ++i)
				coCoord[i] = posAcc[elem->vertex(i)];

			// get sides
			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART,dim, 1));

			for (size_t i=0;i<geo.num_scv();i++){
				typename DimFV1Geometry<dim>::SCV scv = geo.scv(i);
				MathVector<dim> scvLocalBary = scv.local_corner(0);
				for (size_t co=1;co<scv.num_corners();co++){
					scvLocalBary+=scv.local_corner(co);
				}
				scvLocalBary/=scv.num_corners();
				size_t node = scv.node_id();

				std::vector<number> vShape;
				rTrialSpace.shapes(vShape, scvLocalBary);
//				MathVector<dim> localValue = 0;

				size_t nos=sides.size();
				// sum up shapes
				VType value;
				value = 0;
				for (size_t sh=0;sh<nos;sh++){
					VType localValue;
					if (useGridFunction)
					for (size_t d=0;d<dim;d++){
						u->dof_indices(sides[sh], d, multInd);
						assignVal(localValue,d,DoFRef(*u,multInd[0]));
					}
					else
						localValue = acU[sides[sh]];
					localValue *= vShape[sh];
					value += localValue;
				}
				value *= scv.volume();
				// find sides associated with scv
				for (size_t j=0;j<nos;j++){
					side_type* sidej = sides[j];
					for (size_t sco=0;sco<sidej->num_vertices();sco++){
						if (VecDistance(coCoord[node],posAcc[sidej->vertex(sco)])<1e-12){
							m_acSideVolume[sidej]  += scv.volume();
							acUHat[sidej] += value;
						}
					}
				}
			}
		}
	}
	averageByVolume(acUHat,m_acSideVolume,domain,m_uInfo,m_walls);
	if (useGridFunction) copyWallData(acUHat,u,m_walls);
	else copyWallData(acUHat,acU,m_uInfo,m_walls);
	if (dim==2) constrainingSideAveraging<TGridFunction,side_type,ConstrainingEdge,VType>(acUHat,m_uInfo);
	else {
		constrainingSideAveraging<TGridFunction,side_type,ConstrainingTriangle,VType>(acUHat,m_uInfo);
		constrainingSideAveraging<TGridFunction,side_type,ConstrainingQuadrilateral,VType>(acUHat,m_uInfo);
	}
}

// fvcr Filter for side data
template <typename TGridFunction>
template <typename VType>
void FVCRBoxFilter<TGridFunction>::apply_filter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acUHat,
			   SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acU){
	bool useGridFunction = true;
	if (u==SPNULL) useGridFunction = false;
	
	// set attachments to zero
	domain_type& domain = *m_uInfo->domain().get();
	SetAttachmentValues(m_acSideVolume, domain.grid()->template begin<side_type>(), domain.grid()->template end<side_type>(), 0);
	SetAttachmentValues(acUHat, domain.grid()->template begin<side_type>(), domain.grid()->template end<side_type>(), 0);
	
	std::vector<DoFIndex> multInd;

	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
	side_secure_container sides;

	DimCRFVGeometry<dim> geo;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;

			//	get corners of element
			size_t noc=elem->num_vertices();
			for(size_t i = 0; i < noc; ++i)
				coCoord[i] = posAcc[elem->vertex(i)];

			// get sides
			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());

			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART,dim, 1));

			for (size_t i=0;i<geo.num_scv();i++){
				typename DimCRFVGeometry<dim>::SCV scv = geo.scv(i);
				MathVector<dim> scvLocalBary = scv.local_corner(0);
				for (size_t co=1;co<scv.num_corners();co++){
					scvLocalBary+=scv.local_corner(co);
				}
				scvLocalBary/=scv.num_corners();
				size_t sidei = scv.node_id();

				std::vector<number> vShape;
				rTrialSpace.shapes(vShape, scvLocalBary);

				size_t nos=sides.size();
				// sum up shapes
				VType value;
				value = 0;
				for (size_t sh=0;sh<nos;sh++){
					VType localValue;
					if (useGridFunction)
					for (size_t d=0;d<dim;d++){
						u->dof_indices(sides[sh], d, multInd);
						assignVal(localValue,d,DoFRef(*u,multInd[0]));
					}
					else
						localValue = acU[sides[sh]];
					localValue *= vShape[sh];
					value += localValue;
				}
				value *= scv.volume();
				m_acSideVolume[sides[sidei]]  += scv.volume();
				acUHat[sides[sidei]] += value;
			}
		}
	}
	averageByVolume(acUHat,m_acSideVolume,domain,m_uInfo,m_walls);
	if (useGridFunction) copyWallData(acUHat,u,m_walls);
	else copyWallData(acUHat,acU,m_uInfo,m_walls);
	if (dim==2) constrainingSideAveraging<TGridFunction,side_type,ConstrainingEdge,VType>(acUHat,m_uInfo);
	else {
		constrainingSideAveraging<TGridFunction,side_type,ConstrainingTriangle,VType>(acUHat,m_uInfo);
		constrainingSideAveraging<TGridFunction,side_type,ConstrainingQuadrilateral,VType>(acUHat,m_uInfo);
	}
}
	
template <typename TGridFunction>
void FVCRBoxFilter<TGridFunction>::compute_filterwidth_fvcr()
{
	// set attachments to zero
	domain_type& domain = *m_uInfo->domain().get();
	SetAttachmentValues(m_acSideVolume, domain.grid()->template begin<side_type>(), domain.grid()->template end<side_type>(), 0);
				
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		
	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();
		
	typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
	side_secure_container sides;
		
	DimCRFVGeometry<dim> geo;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);
			
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
				
			//	get corners of element
			size_t noc=elem->num_vertices();
			for(size_t i = 0; i < noc; ++i)
				coCoord[i] = posAcc[elem->vertex(i)];
				
			// get sides
			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );
				
			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				
			for (size_t i=0;i<geo.num_scv();i++){
				typename DimCRFVGeometry<dim>::SCV scv = geo.scv(i);
				size_t sidei = scv.node_id();
				m_acSideVolume[sides[sidei]]  += scv.volume();
			}
		}
	}
}
	
template <typename TGridFunction>
void ElementBoxFilter<TGridFunction>::compute_filterwidth_fv1(){
	// set attachments to zero
	domain_type& domain = *m_uInfo->domain().get();
	SetAttachmentValues(m_acVertexVolume, domain.grid()->template begin<vertex_type>(), domain.grid()->template end<vertex_type>(), 0);
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		
	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();
		
	DimFV1Geometry<dim> geo;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);
			
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
				
			//	get corners of element
			size_t noc=elem->num_vertices();
			for(size_t i = 0; i < noc; ++i)
				coCoord[i] = posAcc[elem->vertex(i)];
				
			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				
			for (size_t i=0;i<geo.num_scv();i++){
				typename DimFV1Geometry<dim>::SCV scv = geo.scv(i);
				for (size_t node=0;node<elem->num_vertices();node++){
					m_acVertexVolume[elem->vertex(node)]  += scv.volume();
				}
			}
		}
	}
}
	
// Element filter for vertex data
template <typename TGridFunction>
template <typename VType>
void ElementBoxFilter<TGridFunction>::apply_filter(PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& acUHat,
												   SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<vertex_type,Attachment<VType> >& acU){
	bool useGridFunction = true;
	if (u==SPNULL) useGridFunction = false;
	
	// set attachments to zero
	domain_type& domain = *m_uInfo->domain().get();
	SetAttachmentValues(m_acVertexVolume, domain.grid()->template begin<vertex_type>(), domain.grid()->template end<vertex_type>(), 0);
	SetAttachmentValues(acUHat, domain.grid()->template begin<vertex_type>(), domain.grid()->template end<vertex_type>(), 0);
	
	std::vector<DoFIndex> multInd;
		
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		
	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();
		
	DimFV1Geometry<dim> geo;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);
			
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
				
			//	get corners of element
			size_t noc=elem->num_vertices();
			for(size_t i = 0; i < noc; ++i)
				coCoord[i] = posAcc[elem->vertex(i)];
				
			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				
			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();
				
			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE,dim, 1));
				
			for (size_t i=0;i<geo.num_scv();i++){
				typename DimFV1Geometry<dim>::SCV scv = geo.scv(i);
				MathVector<dim> scvLocalBary = scv.local_corner(0);
				for (size_t co=1;co<scv.num_corners();co++){
					scvLocalBary+=scv.local_corner(co);
				}
				scvLocalBary/=scv.num_corners();
					
				std::vector<number> vShape;
				rTrialSpace.shapes(vShape, scvLocalBary);
				// sum up shapes
				VType value;
				value = 0;
				for (size_t sh=0;sh<noc;sh++){
					VType localValue;
					if (useGridFunction)
						for (size_t d=0;d<dim;d++){
							u->dof_indices(elem->vertex(sh), d, multInd);
							assignVal(localValue,d,DoFRef(*u,multInd[0]));
						}
					else
						localValue = acU[elem->vertex(sh)];
					localValue *= vShape[sh];
					value += localValue;
				}
				value *= scv.volume();
				for (size_t node=0;node<elem->num_vertices();node++){
					m_acVertexVolume[elem->vertex(node)]  += scv.volume();
					acUHat[elem->vertex(node)] += value;
				}
			}
		}
	}
	averageByVolume(acUHat,m_acVertexVolume,domain,m_uInfo,m_walls);
	if (useGridFunction) copyWallData(acUHat,u,m_walls);
	else copyWallData(acUHat,acU,m_uInfo,m_walls);
}
	
// fvcr Filter for side data
template <typename TGridFunction>
void ElementBoxFilter<TGridFunction>::compute_filterwidth_fvcr(){
	// set attachments to zero
	domain_type& domain = *m_uInfo->domain().get();
	SetAttachmentValues(m_acSideVolume, domain.grid()->template begin<side_type>(), domain.grid()->template end<side_type>(), 0);
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		
	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();
		
	typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
	side_secure_container sides;
		
	DimCRFVGeometry<dim> geo;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);
			
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
				
			//	get corners of element
			size_t noc=elem->num_vertices();
			for(size_t i = 0; i < noc; ++i)
				coCoord[i] = posAcc[elem->vertex(i)];
				
			// get sides
			domain.grid()->associated_elements(sides, static_cast<elem_type*>(elem) );
				
			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				
			for (size_t i=0;i<geo.num_scv();i++){
				typename DimCRFVGeometry<dim>::SCV scv = geo.scv(i);
				// sum up shapes
				for (size_t j=0;j<sides.size();j++){
					m_acSideVolume[sides[j]]  += scv.volume();
				}
			}
		}
	}
}
	

// fvcr Filter for side data
template <typename TGridFunction>
template <typename VType>
void ElementBoxFilter<TGridFunction>::apply_filter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acUHat,
													SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<side_type,Attachment<VType> >& acU){
	bool useGridFunction = true;
	if (u==SPNULL) useGridFunction = false;
	
	// set attachments to zero
	domain_type& domain = *m_uInfo->domain().get();
	SetAttachmentValues(m_acSideVolume, domain.grid()->template begin<side_type>(), domain.grid()->template end<side_type>(), 0);
	SetAttachmentValues(acUHat, domain.grid()->template begin<side_type>(), domain.grid()->template end<side_type>(), 0);
	
	std::vector<DoFIndex> multInd;
		
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	
	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();
		
	typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
	side_secure_container sides;
		
	DimCRFVGeometry<dim> geo;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);
			
		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
				
			//	get corners of element
			size_t noc=elem->num_vertices();
			for(size_t i = 0; i < noc; ++i)
				coCoord[i] = posAcc[elem->vertex(i)];
				
			// get sides
			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );
				
			//	evaluate finite volume geometry
			geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				
			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();
				
			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART,dim, 1));
				
			for (size_t i=0;i<geo.num_scv();i++){
				typename DimCRFVGeometry<dim>::SCV scv = geo.scv(i);
				MathVector<dim> scvLocalBary = scv.local_corner(0);
				for (size_t co=1;co<scv.num_corners();co++){
					scvLocalBary+=scv.local_corner(co);
				}
				scvLocalBary/=scv.num_corners();
					
				std::vector<number> vShape;
				rTrialSpace.shapes(vShape, scvLocalBary);
					
				size_t nos=sides.size();
				// sum up shapes
				VType value;
				value = 0;
				for (size_t sh=0;sh<nos;sh++){
					VType localValue;
					if (useGridFunction)
						for (size_t d=0;d<dim;d++){
							u->dof_indices(sides[sh], d, multInd);
							assignVal(localValue,d,DoFRef(*u,multInd[0]));
						}
					else
						localValue = acU[sides[sh]];
					localValue *= vShape[sh];
					value += localValue;
				}
				value *= scv.volume();
				for (size_t j=0;j<sides.size();j++){
					m_acSideVolume[sides[j]]  += scv.volume();
					acUHat[sides[j]] += value;
				}
			}
		}
	}
	averageByVolume(acUHat,m_acSideVolume,domain,m_uInfo,m_walls);
	if (useGridFunction) copyWallData(acUHat,u,m_walls);
	else copyWallData(acUHat,acU,m_uInfo,m_walls);
	if (dim==2) constrainingSideAveraging<TGridFunction,side_type,ConstrainingEdge,VType>(acUHat,m_uInfo);
	else {
		constrainingSideAveraging<TGridFunction,side_type,ConstrainingTriangle,VType>(acUHat,m_uInfo);
		constrainingSideAveraging<TGridFunction,side_type,ConstrainingQuadrilateral,VType>(acUHat,m_uInfo);
	}
}
	
	
	
} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__INCOMPRESSIBLE__FILTER__IMPL__ */
