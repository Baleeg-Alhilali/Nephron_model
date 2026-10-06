/*
 * Copyright (c) 2010-2013:  G-CSC, Goethe University Frankfurt
 * Author: Josef Dubsky, Andreas Vogel, Raphael Prohl
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


#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__DIFFUSION_LENGTH__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__DIFFUSION_LENGTH__

// other ug4 modules
#include "common/common.h"

namespace ug{
namespace NavierStokes{

/**
 *
 * \param[out]		DiffLengthSqInv		Inverse of squared Diffusion Length for each integration point
 * \param[in]		geo					Finite Volume Geometry
 */
template <typename TFVGeometry>
void NSDiffLengthFivePoint(number DiffLengthSqInv[], const TFVGeometry& geo)
{
    //	dimension of element
	static const size_t dim = TFVGeometry::dim;

    // Check the SCVF - number of IPs (order of elements)
    UG_ASSERT(geo.scvf(0).num_ip() == 1, "Only implemented for first order.");

	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	//	get SubControlVolumeFace
		const typename TFVGeometry::SCVF& scvf = geo.scvf(i);

	//	get associated SubControlVolumes
		const typename TFVGeometry::SCV& scvFrom = geo.scv(scvf.from());
		const typename TFVGeometry::SCV& scvTo = geo.scv(scvf.to());

	//	Norm of Normal to SCVF squared
		const number normNormalSq = VecTwoNormSq(scvf.normal());

	//	average squared size of associated SCV
		number areaSCVSq = 0.5 * (scvFrom.volume() + scvTo.volume());
		areaSCVSq *= areaSCVSq;

		if(dim == 2)
		{
			DiffLengthSqInv[i] = 2.0 * normNormalSq/ areaSCVSq + 8.0 / normNormalSq;
		}
		else if (dim == 3)
		{
		//	Distance between edge midpoint and center of element, that are part
		//	of the SCVF
			const number distSq = VecDistanceSq(scvf.global_corner(0),
												scvf.global_corner(2));

			DiffLengthSqInv[i] = 2.0 * normNormalSq/ areaSCVSq + 8.0 * distSq / normNormalSq;
		}
		else
			UG_THROW("NSDiffLengthFivePoint not implemented for dimension "<<dim);
	}
}

template <typename TFVGeometry>
void NSDiffLengthRaw(number DiffLengthSqInv[], const TFVGeometry& geo)
{
    //	dimension of element
	static const size_t dim = TFVGeometry::dim;

    // Check the SCVF - number of IPs (order of elements)
    UG_ASSERT(geo.scvf(0).num_ip() == 1, "Only implemented for first order.");


	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	//	get SubControlVolumeFace
		const typename TFVGeometry::SCVF& scvf = geo.scvf(i);

	//	get associated SubControlVolumes
		const typename TFVGeometry::SCV& scvFrom = geo.scv(scvf.from());
		const typename TFVGeometry::SCV& scvTo = geo.scv(scvf.to());

	//	Norm of Normal to SCVF squared
		const number normNormalSq = VecTwoNormSq(scvf.normal());

	//	average squared size of associated SCV
		number areaSCVSq = 0.5 * (scvFrom.volume() + scvTo.volume());
		areaSCVSq *= areaSCVSq;

		if(dim == 2)
		{
			DiffLengthSqInv[i] = 1. / (0.5*areaSCVSq/normNormalSq + 3.0*normNormalSq/8.);
		}
		else if (dim == 3)
		{
		//	Distance between edge midpoint and center of element, that are part
		//	of the SCVF
			const number distSq = VecDistanceSq(scvf.global_corner(0),
												scvf.global_corner(2));

			DiffLengthSqInv[i] = 1. / (0.5*areaSCVSq/normNormalSq + 3.0*distSq/8.);
		}
		else
			UG_THROW("NSDiffLengthRaw not implemented for dimension "<<dim);
    }
}

template <typename TFVGeometry>
void NSDiffLengthCor(number DiffLengthSqInv[], const TFVGeometry& geo)
{
    //	dimension of element
	static const size_t dim = TFVGeometry::dim;

    // Check the SCVF - number of IPs (order of elements)
    UG_ASSERT(geo.scvf(0).num_ip() == 1, "Only implemented for first order.");

    // compute min, max and avg of ||n||^2
    number minNormSq = std::numeric_limits<number>::max();
    number minDistSq = std::numeric_limits<number>::max();
    number avgNormSq = 0.0;
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	//	get SubControlVolumeFace
		const typename TFVGeometry::SCVF& scvf = geo.scvf(i);

	//	Norm of Normal to SCVF squared
		const number normNormalSq = VecTwoNormSq(scvf.normal());

	//	get min
		if(normNormalSq < minNormSq) minNormSq = normNormalSq;

	//	sum for average
		avgNormSq += normNormalSq;

	//	Distance between edge midpoint and center of element, that are part
	//	of the SCVF
		if(dim == 3)
		{
			const number distSq = VecDistanceSq(scvf.global_corner(0),
											scvf.global_corner(2));

		//	get min
			if(distSq < minDistSq) minDistSq = distSq;
		}
	}

	//	devide by num scvf
	avgNormSq /= geo.num_scvf();

	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
	//	get associated SubControlVolumes
		const typename TFVGeometry::SCVF& scvf = geo.scvf(i);
		const typename TFVGeometry::SCV& scvFrom = geo.scv(scvf.from());
		const typename TFVGeometry::SCV& scvTo = geo.scv(scvf.to());

	//	average squared size of associated SCV
		number areaSCVSq = 0.5 * (scvFrom.volume() + scvTo.volume());
		areaSCVSq *= areaSCVSq;

		if(dim == 2)
		{
			DiffLengthSqInv[i] = 2.0 * minNormSq/ areaSCVSq + 8.0 / (3.0*avgNormSq);
		}
		else if (dim == 3)
		{
			DiffLengthSqInv[i] = 2.0 * minNormSq/ areaSCVSq + 8.0 * minDistSq / (3.0*avgNormSq);
		}
		else
			UG_THROW("NSDiffLengthCor not implemented for dimension "<<dim);
	}
}

} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__DIFFUSION_LENGTH__ */
