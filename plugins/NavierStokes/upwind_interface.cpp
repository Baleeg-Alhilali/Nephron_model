/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

// for minimum
#include <limits>
#include <algorithm>
#include <locale>

#include "upwind.h"

namespace ug {
namespace NavierStokes{

template <int dim>
SmartPtr<INavierStokesUpwind<dim> > CreateNavierStokesUpwind(const std::string& name)
{
	std::string n = TrimString(name);
	std::transform(n.begin(), n.end(), n.begin(), ::tolower);

	if(n == "no") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesNoUpwind<dim>());
	if(n == "full") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesFullUpwind<dim>());
//	if(n == "weighted") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesWeightedUpwind<dim>());
	if(n == "skewed") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesSkewedUpwind<dim>());
	if(n == "linearprofileskewed" ||
	   n == "lps") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesLinearProfileSkewedUpwind<dim>());
	if(n == "positive" ||
	   n == "pos") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesPositiveUpwind<dim>());
	if(n == "regular" ||
	   n == "reg") return SmartPtr<INavierStokesUpwind<dim> >(new NavierStokesRegularUpwind<dim>());

	UG_THROW("NavierStokes: upwind type '"<<name<<"' not found. Options are: "
	         "no, full, weighted, skewed, LinearProfileSkewed (lps), positive (pos), regular (reg)");
}

////////////////////////////////////////////////////////////////////////////////
//	explicit instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_2
template SmartPtr<INavierStokesUpwind<2> >CreateNavierStokesUpwind<2>(const std::string& name);
#endif
#ifdef UG_DIM_3
template SmartPtr<INavierStokesUpwind<3> >CreateNavierStokesUpwind<3>(const std::string& name);
#endif

} // namespace NavierStokes
} // end namespace ug
