/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
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


#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

#include "../register_navier_stokes.h"

#include "compressible_navier_stokes_base.h"

#include "fv1/register_fv1.h"


using namespace std;
using namespace ug::bridge;

namespace ug{
namespace NavierStokes{

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct FunctionalityComp
{

/**
 * Function called for the registration of Domain and Algebra dependent parts
 * of the plugin. All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{}

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{}

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	const int dim = TDomain::dim;
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

//	Compressible Navier-Stokes Base
	{
		typedef CompressibleNavierStokesBase<TDomain> T;
		typedef NavierStokesBase<TDomain> TBase;
		string name = string("CompressibleNavierStokesBase").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_adiabatic_index", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adiabatic_index), "", "AdiabaticIndex")
			.add_method("set_adiabatic_index", static_cast<void (T::*)(number)>(&T::set_adiabatic_index), "", "AdiabaticIndex")
#ifdef UG_FOR_LUA
			.add_method("set_adiabatic_index", static_cast<void (T::*)(const char*)>(&T::set_adiabatic_index), "", "AdiabaticIndex")
#endif
			.add_method("set_mach_number_blend", &T::set_mach_number_blend, "", "Set mach number blending");
		reg.add_class_to_group(name, "CompressibleNavierStokesBase", tag);
	}
}

/**
 * Function called for the registration of Dimension dependent parts
 * of the plugin. All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{}

}; // end Functionality
} // end namespace NavierStokes

/**
 * This function is called when the plugin is loaded.
 */
void Init___CompressibleNavierStokes(Registry* reg, string grp)
{
	grp.append("SpatialDisc/CompressibleNavierStokes/");
	typedef NavierStokes::FunctionalityComp Functionality;

	try{
		//RegisterDimension2d3dDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dDependent<Functionality>(*reg,grp);
//		RegisterAlgebraDependent<Functionality>(*reg,grp);
		//RegisterDomain2d3dAlgebraDependent<Functionality>(*reg,grp);

		Init___CompressibleNavierStokes___FV1(reg, grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}


}// namespace ug
