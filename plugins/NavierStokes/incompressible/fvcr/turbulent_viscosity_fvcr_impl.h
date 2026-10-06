/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_
#define __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_

namespace ug{
namespace NavierStokes{

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::transferToLowerLevels(aSideNumber& aaData,ApproximationSpace<domain_type>& approximationSpace){
	for(int si = 0; si < approximationSpace.domain()->subset_handler()->num_subsets(); ++si){
		// transfer to lower levels, averaging over child edges (2d) / child faces (3d)
		for (size_t lev=approximationSpace.num_levels()-2;(int)lev>=0;lev--){
			const DoFDistribution& lDD = *approximationSpace.dof_distribution(GridLevel(lev, GridLevel::LEVEL));
			const MultiGrid& grid = *lDD.multi_grid();
			typedef typename DoFDistribution::traits<side_type>::const_iterator coarseLevelSideIter;
			coarseLevelSideIter clsIter, clsIterEnd;
			clsIter = lDD.template begin<side_type>(si);
			clsIterEnd = lDD.template end<side_type>(si);
			for (;clsIter != clsIterEnd;clsIter++){
				side_type* elem = *clsIter;
				size_t numChildren = grid.num_children<side_type>(elem);
				number avgValue=0;
				for (size_t i=0;i<numChildren;i++){
					avgValue += aaData[grid.get_child<side_type>(elem, i)];
				}
				avgValue/=numChildren;
				aaData[elem] = avgValue;
			}
			if (lev==0) break;
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::fillAttachment(aSideDimVector& aaU,SmartPtr<TGridFunction> u){
	//	get domain
	domain_type& domain = *u->domain().get();
	//	create Multiindex
	std::vector<DoFIndex> multInd;
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		ElemIterator iter = u->template begin<side_type>(si);
		ElemIterator iterEnd = u->template end<side_type>(si);
		for(  ;iter !=iterEnd; ++iter)
		{
			side_type* side = iter;
			for (int d=0;d<dim;d++){
				u->dof_indices(side, d, multInd);
				aaU[side][d]=DoFRef(*u,multInd[0]);
			}
		}
	}
}

void assignValue(number& v,size_t ind,number value){
	v=value;	
}
	
void assignValue(MathVector<2>& v,size_t ind,number value){
	v[ind]=value;
}
	
void assignValue(MathVector<3>& v,size_t ind,number value){
	v[ind]=value;
}
	
void assignValue(MathSymmetricMatrix<2>& v,size_t ind,number value){
	v[ind]=value;
}
	
void assignValue(MathSymmetricMatrix<3>& v,size_t ind,number value){
	v[ind]=value;
}

// go over all elements, interpolate data to barycenter, average by multiplying with corresponding element volume and deviding by complete adjacent element volume
// parameters: filtered values, filter volume (computed in function), original values
template <typename TData, int dim, typename TImpl,typename TGridFunction>
template <typename VType>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::elementFilter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaUHat,aSideNumber& aaVol,
																			 SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<side_type,Attachment<VType> >* aaU){
	bool useGridFunction = true;
	if (u==SPNULL) useGridFunction = false;
	
	std::vector<DoFIndex> multInd;
	
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	DimCRFVGeometry<dim> geo;

	// set attachment values to zero
	SetAttachmentValues(aaUHat , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);
	SetAttachmentValues(aaVol , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	// assemble deformation tensor fluxes
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
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			MathVector<dim> bary,localBary;
			bary = 0;

			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				bary += coCoord[i];
			};
			bary /= numVertices;

			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();

			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);

			map.global_to_local(localBary,bary);

			typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
			side_secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );
			if (m_bAdaptive){
				std::vector<DoFIndex> ind;
				m_uInfo->dof_indices(elem,0,ind,true,true);
				get_constrained_sides_cr<side_type,side_secure_container,TGridFunction>(sides,*m_uInfo,ind);
			}

			//	memory for shapes
			std::vector<number> vShape;

			//	evaluate finite volume geometry
			if (!m_bAdaptive)
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
			else
				geo.update_hanging(elem, &(coCoord[0]), domain.subset_handler().get());

			rTrialSpace.shapes(vShape, localBary);

			size_t nofsides = geo.num_scv();

			VType value;
			value = 0;
			number elementVolume = 0;
			for (size_t i=0;i<nofsides;i++){
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(i);
				size_t s = scv.node_id();
				VType localValue;
				if (useGridFunction)
					for (int d=0;d<dim;d++){
						u->dof_indices(sides[s], d, multInd);
						assignValue(localValue,d,DoFRef(*u,multInd[0]));
					}
				else
					localValue = (*aaU)[sides[s]];
				localValue *= vShape[s];
				value += localValue;
				elementVolume += scv.volume();
			}
			value *= elementVolume;
			for (size_t s=0;s<nofsides;s++){
				aaUHat[sides[s]] += value;
				aaVol[sides[s]] += elementVolume;
			}
		}
	}
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		// no filtering in boundary nodes
		if ((m_turbZeroSg.size()!=0) && (m_turbZeroSg.contains(si)==true)){
			for(  ;sideIter !=sideIterEnd; sideIter++)
			{
				side_type* side = *sideIter;
				if (pbm && pbm->is_slave(side)) continue;
				if (useGridFunction){
					for (int d=0;d<dim;d++){
						u->dof_indices(side, d, multInd);
						assignValue(aaUHat[side],d,DoFRef(*u,multInd[0]));
					}
				}
				else
					aaUHat[side] = (*aaU)[side];
			}
		} else {
			for(  ;sideIter !=sideIterEnd; sideIter++)
			{
				side_type* side = *sideIter;
				if (pbm && pbm->is_slave(side)) continue;
				// constrained edges
				if (aaVol[side]==0) continue;
				aaUHat[side]/=(number)aaVol[side];
			}
		}
	}	
	if (m_bAdaptive){ 
		if (dim==2) constrainingSideAveraging<TGridFunction,side_type,ConstrainingEdge,VType>(aaUHat,m_uInfo);
		else {
			constrainingSideAveraging<TGridFunction,side_type,ConstrainingTriangle,VType>(aaUHat,m_uInfo);
			constrainingSideAveraging<TGridFunction,side_type,ConstrainingQuadrilateral,VType>(aaUHat,m_uInfo);
		}
	}
}

// go over all elements, interpolate data to scv barycenter, average by multiplying with corresponding scv volume and deviding by volume of complete control volume
template <typename TData, int dim, typename TImpl,typename TGridFunction>
template <typename VType>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::scvFilter(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaUHat,aSideNumber& aaVol,
																		 SmartPtr<TGridFunction> u,PeriodicAttachmentAccessor<side_type,Attachment<VType> >* aaU){
	bool useGridFunction = true;
	if (u==NULL) useGridFunction = false;

	std::vector<DoFIndex> multInd;
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	DimCRFVGeometry<dim> geo;
	
	// set attachment values to zero
	SetAttachmentValues(aaUHat , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);
	SetAttachmentValues(aaVol , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);
	
	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
	
	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();
	
	// assemble deformation tensor fluxes
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
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			MathVector<dim> bary,localBary;
			bary = 0;
			
			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
				bary += coCoord[i];
			};
			bary /= numVertices;
			
			//	reference object id
			ReferenceObjectID roid = elem->reference_object_id();
			
			//	get trial space
			const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));
			
			//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, coCoord);
			
			map.global_to_local(localBary,bary);
			
			typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
			side_secure_container sides;
			
			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");
			
			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );
			
			//	evaluate finite volume geometry
			if (!m_bAdaptive){
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
			} else {
				geo.update_hanging(elem, &(coCoord[0]), domain.subset_handler().get());
				if (geo.num_constrained_dofs()>0){
					m_uInfo->dof_indices(elem,0,multInd,true,true);
					get_constrained_sides_cr<side_type,side_secure_container,TGridFunction>(sides,*m_uInfo,multInd);
				}
			}
			
			size_t nNaturalSides = geo.num_sh();
			size_t nSides = geo.num_scv();
			
			MathVector<dim> scvBary,scvLocalBary;
			for (size_t i=0;i<nSides;i++){
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(i);
				size_t s = geo.node_id();
				scvBary = 0;
				// compute barycenter of scv (average of side corner nodes + element bary)
				for (size_t j=0;j<scv.num_corners();j++){
					scvBary += scv.global_corner(j);
				}
				scvBary/=(number)scv.num_corners();
				map.global_to_local(scvLocalBary,scvBary);
				//	memory for shapes
				std::vector<number> vShape;
				rTrialSpace.shapes(vShape, scvLocalBary);
				VType localValue = 0;
				for (size_t j=0;j<nNaturalSides;j++){
					if (useGridFunction){
						VType uValue;
						for (int d=0;d<dim;d++){
							u->dof_indices(sides[j], d, multInd);
							assignValue(uValue,d,DoFRef(*u,multInd[0]));
						}
						localValue += vShape[j]*uValue;
					} else
						localValue += vShape[j]*aaU[sides[j]];
				}
				localValue *= scv.volume();
				aaVol[sides[s]]  += scv.volume();
				aaUHat[sides[s]] += localValue;
			}
		}
	}
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		// no filtering in boundary nodes
		if ((m_turbZeroSg.size()!=0) && (m_turbZeroSg.contains(si)==true)){
			for(  ;sideIter !=sideIterEnd; sideIter++)
			{
				side_type* side = *sideIter;
				if (pbm && pbm->is_slave(side)) continue;
				aaUHat[side] = aaU[side];
			}
		} else {
			for(  ;sideIter !=sideIterEnd; sideIter++)
			{
				side_type* side = *sideIter;
				if (pbm && pbm->is_slave(side)) continue;
				aaUHat[side]/=(number)aaVol[side];
			}
		}
	}
	if (m_bAdaptive){
		if (dim==2) constrainingSideAveraging<TGridFunction,side_type,ConstrainingEdge,VType>(aaUHat,m_uInfo);
		else {
			constrainingSideAveraging<TGridFunction,side_type,ConstrainingTriangle,VType>(aaUHat,m_uInfo);
			constrainingSideAveraging<TGridFunction,side_type,ConstrainingQuadrilateral,VType>(aaUHat,m_uInfo);
		}
	}
}
	
template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::assembleDeformationTensor(aSideTensor& aaDefTensor,aSideNumber& aaVol,SmartPtr<TGridFunction> u,aSideDimVector* aaU){
	bool useGridFunction = true;
	if (u==SPNULL) useGridFunction = false;

	//	get domain
	domain_type& domain = *m_uInfo->domain().get();
	// get periodic boundary manager
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();

	//	create Multiindex
	std::vector<DoFIndex> multInd;

	DimCRFVGeometry<dim> geo;

	// add boundary subsets to enforce boundary subset computations in geo.update()
	for(size_t i = 0; i < m_turbZeroSg.size(); ++i){
		geo.add_boundary_subset(m_turbZeroSg[i]);
	}

	// set attachment values to zero
	SetAttachmentValues(aaDefTensor , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);
	SetAttachmentValues(aaVol , m_uInfo->template begin<side_type>(), m_uInfo->template end<side_type>(), 0);

	//	coord and vertex array
	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	//	get position accessor
	typedef typename domain_type::position_accessor_type position_accessor_type;
	const position_accessor_type& posAcc = domain.position_accessor();

	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		//	get iterators
		ElemIterator iter = m_uInfo->template begin<elem_type>(si);
		ElemIterator iterEnd = m_uInfo->template end<elem_type>(si);

		//	loop elements of dimension
		for(  ;iter !=iterEnd; ++iter)
		{
			//	get Elem
			elem_type* elem = *iter;
			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			for(size_t i = 0; i < numVertices; ++i){
				vVrt[i] = elem->vertex(i);
				coCoord[i] = posAcc[vVrt[i]];
			};

			static const size_t MaxNumSidesOfElem = 10;

			typedef MathVector<dim> MVD;
			std::vector<MVD> uValue(MaxNumSidesOfElem);
			MVD ipVelocity;

			typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
			side_secure_container sides;

			UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

			domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

			//	evaluate finite volume geometry
			if (!m_bAdaptive){
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
			} else {
				geo.update_hanging(elem, &(coCoord[0]), domain.subset_handler().get());
				if (geo.num_constrained_dofs()>0){
					m_uInfo->dof_indices(elem,0,multInd,true,true);
					get_constrained_sides_cr<side_type,side_secure_container,TGridFunction>(sides,*m_uInfo,multInd);
				}
			}

			size_t nNaturalSides = geo.num_sh();
			size_t nSides = geo.num_scv();
			size_t nip = geo.num_scvf();

			if (useGridFunction){
				for (size_t s=0;s < nNaturalSides;s++){
					for (int d=0;d<dim;d++){
						m_uInfo->dof_indices(sides[s], d, multInd);
						uValue[s][d]=DoFRef(*u,multInd[0]);
					}
				}
			} else
			{
				for (size_t s=0;s < nNaturalSides;s++)
					uValue[s]=(*aaU)[sides[s]];
			}

			for (size_t i=0;i < nSides;i++)
			{
				const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(i);
				size_t s = scv.node_id();
				aaVol[sides[s]] += scv.volume();
			}

			for (size_t ip=0;ip<nip;ip++){
				// 	get current SCVF
				ipVelocity = 0;
				const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
				for (size_t s=0;s < nNaturalSides;s++){
					for (int d=0;d<dim;d++){
					    ipVelocity[d] += scvf.shape(s)*uValue[s][d];
					};
				};
				dimMat ipDefTensorFlux;
				ipDefTensorFlux = 0;
				size_t ind = 0;
				for (int i=0;i<dim;i++){
					for (int j=i;j<dim;j++){
						ipDefTensorFlux[ind]= 0.5 * (ipVelocity[i] * scvf.normal()[j] + ipVelocity[j] * scvf.normal()[i]);
						ind++;
					}
				}
				aaDefTensor[sides[scvf.from()]]+=ipDefTensorFlux;
				aaDefTensor[sides[scvf.to()]]-=ipDefTensorFlux;
			}
			for(size_t sgi = 0; sgi < m_turbZeroSg.size(); ++sgi){
				const size_t sgsi=m_turbZeroSg[sgi];
				if (geo.num_bf(sgsi) == 0) continue;
				for(size_t bfi = 0; bfi < geo.num_bf(sgsi); ++bfi){
					const typename DimCRFVGeometry<dim>::BF& bf = geo.bf(sgsi, bfi);
					const size_t sideID = bf.node_id();
					size_t ind = 0;
					for (int i=0;i<dim;i++){
						for (int j=i;j<dim;j++){
							// bf ip and u position are identical for CR-FV-Geometry
							//for debug UG_LOG("[" << i << "," << j << "]" <<  0.5 * (uValue[sideID][i] * bf.normal()[j] + uValue[sideID][j] * bf.normal()[i]) << "\n");
							aaDefTensor[sides[sideID]][ind] -= 0.5 * (uValue[sideID][i] * bf.normal()[j] + uValue[sideID][j] * bf.normal()[i]);
							ind++;
						}
					}
				}
			}
		}
	}
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		//for debug UG_LOG("|||||||||||||||||||||||||||||||||||\n");
		//for debug UG_LOG("si=" << si << "\n");
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			aaDefTensor[side]/=(number)aaVol[side];
			//for debug UG_LOG("$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
		//	MathVector<dim> posCo;
		//				for (int d=0;d<dim;d++)
		//								posCo[d] = 0.5*posAcc[side->vertex(0)][d] + 0.5*posAcc[side->vertex(1)][d];
						//for debug UG_LOG(" c=" << posCo << "\n");
		//for debug				for (int d1=0;d1<dim;d1++)
		//for debug					for (int d2=0;d2<dim;d2++)
									//for debug UG_LOG(" tensor(" << d1 << "," << d2 << ")=" << aaDefTensor[side][d1][d2] << "\n");
			//for debug UG_LOG(" norm=" << FNorm(aaDefTensor[side]) << "\n");
		}
	}
	// handle constrained sides
	if (m_bAdaptive){
		if (dim==2) constrainingSideAveraging<TGridFunction,side_type,ConstrainingEdge,dimMat>(aaDefTensor,m_uInfo);
		else {
			constrainingSideAveraging<TGridFunction,side_type,ConstrainingTriangle,dimMat>(aaDefTensor,m_uInfo);
			constrainingSideAveraging<TGridFunction,side_type,ConstrainingQuadrilateral,dimMat>(aaDefTensor,m_uInfo);
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::scaleTensorByNorm(aSideTensor& aaTensor){
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	// get periodic boundary manager
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	// average
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			//UG_LOG(aaTensor[side] << "\n");
			//UG_LOG(FNorm(aaTensor[side]) << "\n");
			//UG_LOG("--------------------\n");
			//for debug UG_LOG("&&&&&&&& norm = " << FNorm(aaTensor[side]) << "\n");
			aaTensor[side].scale_by_fnorm();
		}
	}
}

template <typename TData, int dim, typename TImpl,typename TGridFunction>
void StdTurbulentViscosityData<TData,dim,TImpl,TGridFunction>::addUiUjTerm(aSideTensor& aaResult,const number factor,
																		   SmartPtr<TGridFunction> u,aSideDimVector* aaU){
	bool useGridFunction = true;
	if (u==SPNULL) useGridFunction = false;
	
	//	get domain of grid function
	domain_type& domain = *m_uInfo->domain().get();
	// get periodic boundary manager
	PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
	//	create Multiindex
	std::vector<DoFIndex> multInd;
		
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
		SideIterator sideIter = m_uInfo->template begin<side_type>(si);
		SideIterator sideIterEnd = m_uInfo->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			side_type* side = *sideIter;
			if (pbm && pbm->is_slave(side)) continue;
			dimMat Tij;
			MathVector<dim> uValue;
			if (useGridFunction)
				for (int d=0;d<dim;d++){
					u->dof_indices(side, d, multInd);
					uValue[d]=DoFRef(*u,multInd[0]);
				}
			else 
				for (int d=0;d<dim;d++) uValue[d] = (*aaU)[side][d];
			size_t ind=0;
			for (int d1=0;d1 < dim;d1++)
				for (int d2=d1;d2 < dim;d2++){
					Tij[ind] = uValue[d1]*uValue[d2];
					ind++;
				}
			Tij*=factor;
			aaResult[side]+=Tij;
		}
	}
}

template<typename TGridFunction>
void CRSmagorinskyTurbViscData<TGridFunction>::update(){
	//	get domain of grid function
	domain_type& domain = *m_u->domain().get();
	SetAttachmentValues(m_acTurbulentViscosity, m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);

	//	coord and vertex array
	//	MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
	//	Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

	// assemble deformation tensor fluxes
	this->assembleDeformationTensor(m_acDeformation,m_acVolume,m_u,NULL);
	// compute turbulent viscosity , loop over sides
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		if ((m_turbZeroSg.size()!=0) && (m_turbZeroSg.contains(si)==true)) continue;
		SideIterator sideIter = m_u->template begin<side_type>(si);
		SideIterator sideIterEnd = m_u->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++)
		{
			//	get Elem
			side_type* side = *sideIter;
			if (m_pbm && m_pbm->is_slave(side)) continue;
			number delta = m_acVolume[side];
			// for possible other choices of delta see Fr�hlich p 160
			delta = pow(delta,(number)1.0/(number)dim);
			m_acTurbulentViscosity[side] = m_c*m_c * delta*delta * m_acDeformation[side].fnorm();
		}
	}
	if (m_bAdaptive){
		if (dim==2) constrainingSideAveraging<TGridFunction,side_type,ConstrainingEdge,number>(m_acTurbulentViscosity,m_u);
		else {
			constrainingSideAveraging<TGridFunction,side_type,ConstrainingTriangle,number>(m_acTurbulentViscosity,m_u);
			constrainingSideAveraging<TGridFunction,side_type,ConstrainingQuadrilateral,number>(m_acTurbulentViscosity,m_u);
		}
	}
	// transfer attachment data to lower levels
	this->transferToLowerLevels(m_acTurbulentViscosity,*m_spApproxSpace);
}

template<typename TGridFunction>
void CRDynamicTurbViscData<TGridFunction>::update(){
	//	get domain of grid function
	domain_type& domain = *m_u->domain().get();

	//	get position accessor
	// for debug 
	typedef typename domain_type::position_accessor_type position_accessor_type;
	// for debug 
	const position_accessor_type& posAcc = domain.position_accessor();

	// initialize attachment values with 0
//	SetAttachmentValues(m_acDeformation , m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acTurbulentViscosity, m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
//	SetAttachmentValues(m_acVolume,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acTurbulentCNew,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
//	SetAttachmentValues(m_acVolumeHat,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
//	SetAttachmentValues(m_acUHat,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
//	SetAttachmentValues(m_acDeformationHat,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
//	SetAttachmentValues(m_acLij,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);
	SetAttachmentValues(m_acMij,m_grid->template begin<side_type>(), m_grid->template end<side_type>(), 0);

	// compute Lij term \hat{u_i u_j} - \hat{u_i} \hat{u_j}
	// \hat{u}
	this->elementFilter(m_acUHat,m_acVolumeHat,m_u,&m_acUHat);
	// use Mij attachment to store first Lij part
	// u_i u_j
	this->addUiUjTerm(m_acMij,1.0,m_u,NULL);
	// \hat{u_i u_j}
	this->elementFilter(m_acLij,m_acVolumeHat,SPNULL,&m_acMij);
	// \hat{u_i u_j} - \hat{u_i} \hat{u_j}
	this->addUiUjTerm(m_acLij,-1.0,SPNULL,&m_acUHat);

	// Mij term
	// first term |\hat{S}| \hat{S}
	// assemble \hat{S} using \hat{u}
	this->assembleDeformationTensor(m_acDeformationHat,m_acVolume,SPNULL,&m_acUHat);
	// normalize \hat{S}
	this->scaleTensorByNorm(m_acDeformationHat);
	// Mij second term \hat{|S|S}
	// compute S
	this->assembleDeformationTensor(m_acDeformation,m_acVolume,m_u,NULL);
	// compute |S| S
	this->scaleTensorByNorm(m_acDeformation);
	// filter |S| S
	//for debug UG_LOG("------------------------------------------------------\n");
	this->elementFilter(m_acMij,m_acVolumeHat,SPNULL,&m_acDeformation);

	// complete Mij term computation by scaling and adding the two terms,
	// solve the local least squares problem and compute local c and local turbulent viscosity
	for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
	{
		if ((m_turbZeroSg.size()!=0) && (m_turbZeroSg.contains(si)==true)) continue;
		SideIterator sideIter = m_u->template begin<side_type>(si);
		SideIterator sideIterEnd = m_u->template end<side_type>(si);
		for(  ;sideIter !=sideIterEnd; sideIter++){
			side_type* side = *sideIter;
			if (m_pbm && m_pbm->is_slave(side)) continue;
			// use fixed ratio \hat{delta} / delta or compute it from filter width
			number kappa;
			if (m_fixedRatio) kappa = m_kappa;
			else {
				number delta = m_acVolume[side];
				delta = pow(delta,(number)1.0/(number)dim);
				number deltaHat = m_acVolumeHat[side];
				deltaHat = pow(deltaHat,(number)1.0/(number)dim);
				kappa = deltaHat /(number) delta;
			}
			m_acDeformationHat[side] *= -kappa*kappa;
			m_acMij[side] += m_acDeformationHat[side];
			// compute local c
			// solve least squares problem
			number c = 0;
			for (int d1=0;d1<dim;d1++)
				for (int d2=0;d2<dim;d2++)
					c += m_acLij[side](d1,d2)*m_acMij[side](d1,d2);
			number denom=0;
			//for debug UG_LOG("c=" << c << "\n");
			for (int d1=0;d1<dim;d1++)
				for (int d2=0;d2<dim;d2++)
					denom += m_acMij[side](d1,d2)*m_acMij[side](d1,d2);
			if (denom>1e-15)
				c/=(number)(2.0*denom);
			else c=0;

//			UG_LOG("c=" << c << "\n");

			if (m_spaceFilter==false){
				if (m_timeFilter==false){
					m_acTurbulentViscosity[side] = c * m_acDeformation[side].fnorm();
					UG_LOG("nu_t = " << m_acTurbulentViscosity[side]  << " c = " << c << " kappa = " << kappa << " co=[" << 0.5*(posAcc[side->vertex(0)][0] + posAcc[side->vertex(1)][0]) << "," << 0.5*(posAcc[side->vertex(0)][1] + posAcc[side->vertex(1)][1]) << "]\n");
				} else {
					m_acTurbulentC[side]= (m_timeFilterEps * c + (1-m_timeFilterEps)*m_acTurbulentC[side]);
					m_acTurbulentViscosity[side] = m_acTurbulentC[side] * m_acDeformation[side].fnorm();
				}
				if (m_acTurbulentViscosity[side]+m_viscosityNumber<m_small) m_acTurbulentViscosity[side] = m_viscosityNumber + m_small;			}
			else{
				// store c in viscosity attachment
				m_acTurbulentViscosity[side] = c;
			}
		}
	}
	if (m_spaceFilter==true){
		// filter c
		if (m_timeFilter==false)
			// c has been stored in viscosity attachment
			this->elementFilter(m_acTurbulentC,m_acVolumeHat,SPNULL,&m_acTurbulentViscosity);
		else
			this->elementFilter(m_acTurbulentCNew,m_acVolumeHat,SPNULL,&m_acTurbulentViscosity);
		// compute turbulent viscosity
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
		{
			//for debug UG_LOG("si = " << si << "\n");
			if ((m_turbZeroSg.size()!=0) && (m_turbZeroSg.contains(si)==true)) continue;
			SideIterator sideIter = m_u->template begin<side_type>(si);
			SideIterator sideIterEnd = m_u->template end<side_type>(si);
			for(  ;sideIter !=sideIterEnd; sideIter++){
				side_type* side = *sideIter;
				if (m_pbm && m_pbm->is_slave(side)) continue;
				if (m_timeFilter==true)
					// time averaging, note that c has been stored in m_acVolumeHat
					m_acTurbulentC[side]= (m_timeFilterEps * m_acTurbulentCNew[side] + (1-m_timeFilterEps)*m_acTurbulentC[side]);
				m_acTurbulentViscosity[side] = m_acTurbulentC[side] * m_acDeformation[side].fnorm();
				UG_LOG("nu_t = " << m_acTurbulentViscosity[side]  << " c = " << m_acTurbulentC[side] << " co=[" << 0.5*(posAcc[side->vertex(0)][0] + posAcc[side->vertex(1)][0]) << "," << 0.5*(posAcc[side->vertex(0)][1] + posAcc[side->vertex(1)][1]) << "]\n");
				if (m_acTurbulentViscosity[side]+m_viscosityNumber<m_small) m_acTurbulentViscosity[side] = m_viscosityNumber+m_small;
				// for debug 
				// UG_LOG("nu_t = " << m_acTurbulentViscosity[side]  << " c = " << m_acTurbulentC[side] << " co=[" << 0.5*(posAcc[side->vertex(0)][0] + posAcc[side->vertex(1)][0]) << "," << 0.5*(posAcc[side->vertex(0)][1] + posAcc[side->vertex(1)][1]) << "]\n");
			}
		}
	}
	if (m_bAdaptive){
		if (dim==2) constrainingSideAveraging<TGridFunction,side_type,ConstrainingEdge,number>(m_acTurbulentViscosity,m_u);
		else {
			constrainingSideAveraging<TGridFunction,side_type,ConstrainingTriangle,number>(m_acTurbulentViscosity,m_u);
			constrainingSideAveraging<TGridFunction,side_type,ConstrainingQuadrilateral,number>(m_acTurbulentViscosity,m_u);
		}
	}
	// transfer attachment data to lower levels
	this->transferToLowerLevels(m_acTurbulentViscosity,*m_spApproxSpace);
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES_TURBULENT_VISCOSITY_DATA_IMPL_H_ */
