# Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
# Author: Martin Rupp
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
#      SuperLU plugin
#----------------------------
# information about SuperLU is here http://crd-legacy.lbl.gov/~xiaoye/SuperLU
# if you are using SuperLU, remember to CITE it!
# bibtex information is here
#     http://crd-legacy.lbl.gov/~xiaoye/SuperLU/referencing.html
#
# ug4 plugin by Martin Rupp
# how to use:
#     cmake -DSuperLU6=ON ..
# this will install SuperLU with uginstall to ~/local/SuperLU
# this is tested for Macs and Cekon
# if this does not work, you can supply your own SuperLU path with
#     cmake -DSuperLU6=ON -DSUPERLU_PATH=/path/to/my/superlu .. 
################################################################################



message(STATUS "# SuperLU6 - BEGIN ##############################################################################")
set(pluginName	SuperLU)
set(SOURCES 	super_lu.h super_lu.cpp super_lu_bridge.h super_lu_bridge.cpp  )

cmake_minimum_required(VERSION 2.8.12...4.4)
project(UG_PLUGIN_${pluginName})

# Include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

# Include SuperLU.
option (USE_INTERNAL_SUPERLU "Use internal superlu shipped with this plugin" ON) 
message(STATUS "SuperLU6-Info: using internal version ${USE_INTERNAL_SUPERLU}")


# Clear internal variables
SET(SLU_LIB SLU_LIB-NOTFOUND)

if (USE_INTERNAL_SUPERLU)

	# Option 0: Install SuperLU shipped with this version.
	# SET (BUILD_SHARED_LIBS OFF)
	SET(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
	SET(enable_double ON CACHE BOOL "Enable double")
	SET(enable_single OFF CACHE BOOL "Disable single")
	SET(enable_complex OFF CACHE BOOL "Disable complex")
	SET(enable_complex16 OFF CACHE BOOL "Disable complex16")
	
	SET(enable_tests OFF CACHE BOOL "Disable tests" FORCE)
	SET(enable_examples OFF CACHE BOOL "Disable examples" FORCE)
	SET(enable_doc OFF CACHE BOOL "Disable doc" FORCE)
	
	add_subdirectory(external/superlu)
	
	# The static library has to be build with posititon independent code.
	set_property(TARGET superlu PROPERTY POSITION_INDEPENDENT_CODE ON)
	
	# Internal variables
	SET (SLU_LIB superlu)
	SET(INTERNAL_SUPERLU_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/superlu)
	SET(SUPERLU_INCLUDE_DIRS "${INTERNAL_SUPERLU_PATH}/SRC/")

else (USE_INTERNAL_SUPERLU)


	# Option 1: Try to find SuperLU using PkgConfig.
	find_package(PkgConfig)
	if (PkgConfig_FOUND)
			pkg_check_modules(SUPERLU superlu)
			message(STATUS "SuperLU-Info: ${SUPERLU_PATH} ${SUPERLU_INCLUDE_DIRS} ${SUPERLU_VERSION}")
	endif (PkgConfig_FOUND)


	# Option 2: Detect path for SuperLU-Library
	if(SUPERLU_PATH)
	    # the user can supply a SUPER_LU path for his own compiled SuperLU   
	    set(INTERNAL_SUPERLU_PATH "${SUPERLU_PATH}")
	    message(STATUS "SuperLU-Info: SuperLU is ON, using user provided path SUPERLU_PATH=${INTERNAL_SUPERLU_PATH}")
	endif(SUPERLU_PATH)


	

	# we only support version 4.3
	#find_library(SLU_LIB NAMES superlu_4.3 superlu_4.2 superlu_4.1 superlu_4.0 superlu_5.1.1 superlu_5.1 
	find_library(SLU_LIB NAMES superlu_4.3 superlu
			 PATHS ${INTERNAL_SUPERLU_PATH} PATH_SUFFIXES lib
			 NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)



endif (USE_INTERNAL_SUPERLU)



# Check, if internal variables are set
if ("${SLU_LIB}" STREQUAL "SLU_LIB-NOTFOUND")
	message(FATAL_ERROR "SuperLU library was not found in ${INTERNAL_SUPERLU_PATH}")
endif()



# Check variables.
message(STATUS "SuperLU6-Info: version  ${SUPERLU_VERSION}")  
message(STATUS "SuperLU6-Info: includes ${SUPERLU_INCLUDE_DIRS}") 
message(STATUS "SuperLU6-Info: library  ${SLU_LIB}") 


################################################################################
# Classic binding (static or dynamic plugin).
################################################################################
if(NOT USE_PYBIND11)
# Add to build process 
if(buildEmbeddedPlugins)
	ExportSources(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
    ExportDependencies("${SLU_LIB}")
    #if(EXISTS "${INTERNAL_SUPERLU_PATH}/lib/blas.a")
    # 	ExportDependencies("${INTERNAL_SUPERLU_PATH}/lib/blas.a")
    #endif(EXISTS "${INTERNAL_SUPERLU_PATH}/lib/blas.a")
    ExportIncludes(${SUPERLU_INCLUDE_DIRS})
else(buildEmbeddedPlugins)
    include_directories(${SUPERLU_INCLUDE_DIRS})
	add_library(SuperLU6 SHARED ${SOURCES})	
	target_link_libraries(SuperLU6 ug4)
	target_link_libraries(SuperLU6 ${SLU_LIB})
	
	# CPack specific
	set_target_properties(SuperLU6 PROPERTIES INSTALL_RPATH "$ORIGIN/../../lib/")
	install(TARGETS SuperLU6 LIBRARY DESTINATION bin/plugins COMPONENT plugins)
endif(buildEmbeddedPlugins)
endif(NOT USE_PYBIND11)

################################################################################
# Python binding (static plugin, dynamic python interface).
################################################################################
if(USE_PYBIND11)

	include_directories(${SUPERLU_INCLUDE_DIRS})
	
	###############################################
	# 0) Let's try to locate ug_4 static library....
	###############################################
	
	set(UG4_STATIC_LIBRARY ug4_s)
	
	if(buildDynamicLib)
		set(UG4_LIBRARY ug4)
	else(buildDynamicLib)
		set(UG4_LIBRARY ug4_s)
	endif(buildDynamicLib)
	
	if (BUILD_WITH_PIP)
		get_libug4_from_pip(PYUG4_STATIC_LIBRARY ${UG4_LIBRARY})
		message(STATUS "Library found: ${PYUG4_STATIC_LIBRARY}")
		set(UG4_LIBRARY ${PYUG4_STATIC_LIBRARY})
	endif()
	message(STATUS "Linking with: ${UG4_LIBRARY}")
	

	###############################################
	# 1) Create (static) library for plugin.
	###############################################
	add_library(SuperLU6 STATIC ${SOURCES})


	# IMPORTANT: ensures install RPATH is used during build
	set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

	# avoid adding build directory paths into binaries
	set(CMAKE_SKIP_BUILD_RPATH FALSE)
	set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
	set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

	
	###############################################
 	# 2) Create (dynamic) python library module.
 	# First argumentin 'python_add_library' must match module name in PYBIND11_MODULE call!
	###############################################
	SET(myPluginSources superlu_pybind.cpp)
	SET(myLibraries ug4 SuperLU6)
	python_add_library(pysuperlu MODULE ${myPluginSources} ${SOURCES} WITH_SOABI)
	target_link_libraries (pysuperlu PRIVATE SuperLU6 ${UG4_LIBRARY})
	
	# Link ug4 dynamically.
	if(APPLE)
    	set_target_properties(pysuperlu PROPERTIES 
    		INSTALL_RPATH "@loader_path/../lib")
	elseif(LINUX AND NOT APPLE)
		set_target_properties(pysuperlu PROPERTIES 
			INSTALL_RPATH "$ORIGIN:$ORIGIN/../lib"
			)
	endif()

	install(TARGETS pysuperlu LIBRARY DESTINATION ug4py COMPONENT pymodules)
endif(USE_PYBIND11)

message(STATUS "# SuperLU6 - END ##############################################################################
")
