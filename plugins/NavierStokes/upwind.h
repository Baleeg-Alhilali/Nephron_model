/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__NAVIER_STOKES__UPWIND__
#define __H__UG__NAVIER_STOKES__UPWIND__

//#define UG_NSUPWIND_ASSERT(cond, exp)
// include define below to assert arrays used in stabilization
#define UG_NSUPWIND_ASSERT(cond, exp) UG_ASSERT((cond), exp)

#include <vector>
#include <string>

#include "upwind_interface.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"
#include "lib_disc/spatial_disc/disc_util/hfvcr_geom.h"

namespace ug{
namespace NavierStokes{

/////////////////////////////////////////////////////////////////////////////
// No Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesNoUpwind
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<FV1Geometry, dim, NavierStokesNoUpwind<dim> >,
  public NavierStokesUpwindRegister<CRFVGeometry, dim, NavierStokesNoUpwind<dim> >,
  public NavierStokesUpwindRegister<HCRFVGeometry, dim, NavierStokesNoUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesNoUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		static void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);

	///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const CRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);

		///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const HCRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);
};

/////////////////////////////////////////////////////////////////////////////
// Full Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesFullUpwind
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<FV1Geometry, dim, NavierStokesFullUpwind<dim> >,
  public NavierStokesUpwindRegister<CRFVGeometry, dim, NavierStokesFullUpwind<dim> >,
  public NavierStokesUpwindRegister<HCRFVGeometry, dim, NavierStokesFullUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesFullUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		static void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);

	///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const CRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);

	///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const HCRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);
};


////////////////////////////////////////////////////////////////////////////////////
// Weighted Upwind
// upwinding between full and no upwind
// shapes computed as m_weight*no_upwind_shape + (1-m_weight)*full_upwind_shape
////////////////////////////////////////////////////////////////////////////////////
/*
template <int dim>
class NavierStokesWeightedUpwind
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<CRFVGeometry, dim, NavierStokesWeightedUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesWeightedUpwind(number weight = 0.5) : m_weight(weight)
		{
			this->set_shape_ip_flag(false);
		}

		void set_weight(number weight) {m_weight = weight;}

	///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const CRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);

	protected:
		number m_weight;
};
*/
/////////////////////////////////////////////////////////////////////////////
// Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesSkewedUpwind
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<FV1Geometry, dim, NavierStokesSkewedUpwind<dim> >,
  public NavierStokesUpwindRegister<CRFVGeometry, dim, NavierStokesSkewedUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesSkewedUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		static void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);

	///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const CRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);
};

/////////////////////////////////////////////////////////////////////////////
// Linear Profile Skewed Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesLinearProfileSkewedUpwind
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<FV1Geometry, dim, NavierStokesLinearProfileSkewedUpwind<dim> >,
  public NavierStokesUpwindRegister<CRFVGeometry, dim, NavierStokesLinearProfileSkewedUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesLinearProfileSkewedUpwind(){this->set_shape_ip_flag(false);}

	///	update of values for FV1Geometry
		template <typename TElem>
		static void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);

	///	update of values for CRFVGeometry
		template <typename TElem>
		static void compute(const CRFVGeometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);
};


/////////////////////////////////////////////////////////////////////////////
// Positive Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesPositiveUpwind
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<FV1Geometry,dim, NavierStokesPositiveUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesPositiveUpwind(){this->set_shape_ip_flag(true);}

	///	update of values for FV1Geometry
		template <typename TElem>
		static void compute(const FV1Geometry<TElem, dim>* geo,
		             const MathVector<dim> vIPVel[maxNumSCVF],
		             number vUpShapeSh[maxNumSCVF][maxNumSH],
		             number vUpShapeIp[maxNumSCVF][maxNumSCVF],
		             number vConvLength[maxNumSCVF]);
};

/////////////////////////////////////////////////////////////////////////////
// Regular Upwind
/////////////////////////////////////////////////////////////////////////////

template <int dim>
class NavierStokesRegularUpwind
: public INavierStokesUpwind<dim>,
  public NavierStokesUpwindRegister<FV1Geometry, dim, NavierStokesRegularUpwind<dim> >
{
	public:
		typedef INavierStokesUpwind<dim> base_type;
		static const size_t maxNumSCVF = base_type::maxNumSCVF;
		static const size_t maxNumSH = base_type::maxNumSH;

	public:
	///	constructor
		NavierStokesRegularUpwind(){this->set_shape_ip_flag(true);}

	///	update of values for FV1Geometry
		template <typename TElem>
		static void compute(const FV1Geometry<TElem, dim>* geo,
					 const MathVector<dim> vIPVel[maxNumSCVF],
					 number vUpShapeSh[maxNumSCVF][maxNumSH],
					 number vUpShapeIp[maxNumSCVF][maxNumSCVF],
					 number vConvLength[maxNumSCVF]);
};


} // namespace NavierStokes
} // end namespace ug

#endif /* __H__UG__NAVIER_STOKES__UPWIND__ */
