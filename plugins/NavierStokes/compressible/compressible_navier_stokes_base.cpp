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

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

#include "compressible_navier_stokes_base.h"
#include <vector>
#include <string>

namespace ug{
namespace NavierStokes{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
CompressibleNavierStokesBase<TDomain>::CompressibleNavierStokesBase(const char* functions,
                                    const char* subsets)
: NavierStokesBase<TDomain>(functions, subsets),
  m_bMachNrBlend(false)
{
};

template<typename TDomain>
CompressibleNavierStokesBase<TDomain>::CompressibleNavierStokesBase(const std::vector<std::string>& vFct,
                                    const std::vector<std::string>& vSubset)
: NavierStokesBase<TDomain>(vFct, vSubset),
  m_bMachNrBlend(false)
{
};

/////////// adiabatic Index

template<typename TDomain>
void CompressibleNavierStokesBase<TDomain>::
set_adiabatic_index(number val)
{
	set_adiabatic_index(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void CompressibleNavierStokesBase<TDomain>::
set_adiabatic_index(const char* fctName)
{
	set_adiabatic_index(LuaUserDataFactory<number, dim>::create(fctName));
}
#endif


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class CompressibleNavierStokesBase<Domain1d>;
#endif
#ifdef UG_DIM_2
template class CompressibleNavierStokesBase<Domain2d>;
#endif
#ifdef UG_DIM_3
template class CompressibleNavierStokesBase<Domain3d>;
#endif

} // namespace NavierStokes
} // end namespace ug




