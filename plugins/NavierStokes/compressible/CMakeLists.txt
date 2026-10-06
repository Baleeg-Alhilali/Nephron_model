# Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
# Author: Raphael Prohl
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

################################################################################
# COMPRESSIBLE_NAVIER_STOKES
#
# The compressible-navier-stokes plugin.
################################################################################

cmake_minimum_required(VERSION 2.6)

project(COMPRESSIBLE_NAVIER_STOKES)

# include the definitions and dependencies for ug-plugins.
include("../../../../cmake/ug_plugin_includes.cmake")


# sources
set(SOURCES	../navier_stokes_base.cpp
			../upwind_interface.cpp
			../upwind.cpp			
			
			
			compressible_navier_stokes_base.cpp
			
			fv1/compressible_navier_stokes_fv1.cpp
			fv1/register_fv1.cpp
			
			
			../register_navier_stokes.cpp
			compressible_navier_stokes_plugin.cpp)



if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/navier_stokes/compressible ${SOURCES})
else(buildEmbeddedPlugins)
	add_library(CompressibleNavierStokes SHARED ${SOURCES})
	target_link_libraries (CompressibleNavierStokes ug4)
endif(buildEmbeddedPlugins)
