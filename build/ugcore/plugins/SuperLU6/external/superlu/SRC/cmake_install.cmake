# Install script for directory: /Users/alhilaba/UG4_promesh_ogrid/plugins/SuperLU6/external/superlu/SRC

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/alhilaba/UG4_promesh_ogrid/lib/libsuperlu.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsuperlu.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsuperlu.a")
    execute_process(COMMAND "/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsuperlu.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/Users/alhilaba/UG4_promesh_ogrid/plugins/SuperLU6/external/superlu/SRC/supermatrix.h"
    "/Users/alhilaba/UG4_promesh_ogrid/plugins/SuperLU6/external/superlu/SRC/slu_Cnames.h"
    "/Users/alhilaba/UG4_promesh_ogrid/plugins/SuperLU6/external/superlu/SRC/slu_dcomplex.h"
    "/Users/alhilaba/UG4_promesh_ogrid/plugins/SuperLU6/external/superlu/SRC/slu_scomplex.h"
    "/Users/alhilaba/UG4_promesh_ogrid/plugins/SuperLU6/external/superlu/SRC/slu_util.h"
    "/Users/alhilaba/UG4_promesh_ogrid/plugins/SuperLU6/external/superlu/SRC/superlu_enum_consts.h"
    "/Users/alhilaba/UG4_promesh_ogrid/build/ugcore/plugins/SuperLU6/external/superlu/SRC/superlu_config.h"
    "/Users/alhilaba/UG4_promesh_ogrid/plugins/SuperLU6/external/superlu/SRC/slu_ddefs.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/superlu/superluTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/superlu/superluTargets.cmake"
         "/Users/alhilaba/UG4_promesh_ogrid/build/ugcore/plugins/SuperLU6/external/superlu/SRC/CMakeFiles/Export/29844f99febfd7af969da705a88a5078/superluTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/superlu/superluTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/superlu/superluTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/superlu" TYPE FILE FILES "/Users/alhilaba/UG4_promesh_ogrid/build/ugcore/plugins/SuperLU6/external/superlu/SRC/CMakeFiles/Export/29844f99febfd7af969da705a88a5078/superluTargets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/superlu" TYPE FILE FILES "/Users/alhilaba/UG4_promesh_ogrid/build/ugcore/plugins/SuperLU6/external/superlu/SRC/CMakeFiles/Export/29844f99febfd7af969da705a88a5078/superluTargets-release.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/superlu" TYPE FILE FILES
    "/Users/alhilaba/UG4_promesh_ogrid/build/ugcore/plugins/SuperLU6/external/superlu/SRC/superluConfig.cmake"
    "/Users/alhilaba/UG4_promesh_ogrid/build/ugcore/plugins/SuperLU6/external/superlu/SRC/superluConfigVersion.cmake"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/alhilaba/UG4_promesh_ogrid/build/ugcore/plugins/SuperLU6/external/superlu/SRC/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
