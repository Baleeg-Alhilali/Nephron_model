/*
 * Copyright (c) 2014:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_algebra_dependent.h"
#include "bridge/util_domain_algebra_dependent.h"

// preconditioner
#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/interface/constrained_linear_iterator.h"
#include "bridge/util_overloaded.h"
#include "super_lu.h"

using namespace std;

namespace ug{

namespace bridge{
namespace SuperLUBridge{


struct Functionality
{


/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra, typename TRegistry=ug::bridge::Registry>
static void Algebra(TRegistry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

//	typedefs for this algebra
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

	// 	SuperLU Solver
	{
		typedef SuperLUSolver<TAlgebra> T;
		typedef IExternalSolver<TAlgebra> TBase1;
		typedef ILinearOperatorInverse<vector_type> TBase2;
		string name = string("SuperLU").append(suffix);
		reg.template add_class_<T,TBase1,TBase2>(name, grp, "SuperLU")
			.add_constructor()
			.add_method("print_stat", &T::print_stat)
			.add_method("equil", &T::equil)
			.add_method("col_perm_natural", &T::col_perm_natural)
			.add_method("col_perm_mdo_ATA", &T::col_perm_mdo_ATA)
			.add_method("col_perm_mdo_AT_plus_A", &T::col_perm_mdo_AT_plus_A)
			.add_method("col_perm_approximate", &T::col_perm_approximate)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SuperLU", tag);
	}
}

template <typename TDomain, typename TAlgebra, typename TRegistry=ug::bridge::Registry>
static void DomainAlgebra(TRegistry& reg, string grp)
{
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	// constrained SuperLU
	{
		typedef ConstrainedLinearIterator<TDomain, TAlgebra, SuperLUSolver<TAlgebra> > T;
		typedef SuperLUSolver<TAlgebra> TBase;
		string name = string("SuperLU_c").append(suffix);
		reg.template add_class_<T,TBase>(name, grp, "SuperLU solver respecting constraints")
		   .template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("domain discretization")
			.add_method("set_time", &T::set_time, "", "time")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SuperLU_c", tag);
	}
}


}; // end Functionality

// end group precond_bridge
/// \}

}// end Preconditioner

/// \addtogroup precond_bridge
template <typename TRegistry=ug::bridge::Registry>
void RegisterBridge_SuperLU(TRegistry& reg, string grp)
{
	grp.append("/Algebra/ExternalSolvers/");
	typedef SuperLUBridge::Functionality Functionality;

	try{
#ifndef UG_USE_PYBIND11
		RegisterAlgebraDependent<Functionality>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(reg,grp);
#else
		RegisterAlgebraDependent<Functionality, TRegistry>(reg,grp);
		RegisterDomainAlgebraDependent<Functionality, TRegistry>(reg,grp);
#endif
	}
	UG_REGISTRY_CATCH_THROW(grp);
}


// The following functions are the exported
extern "C" void
InitUGPlugin_SuperLU6(Registry* reg, string grp)
{
	RegisterBridge_SuperLU(*reg, grp);
}


#ifdef UG_USE_PYBIND11 // Expose for pybind11.
namespace SuperLUBridge{
	void InitUGPlugin(ug::pybind::Registry* reg, string grp)
	{
		RegisterBridge_SuperLU<ug::pybind::Registry>(*reg, grp);
	}
}
#endif


} // namespace bridge
} // namespace ug
