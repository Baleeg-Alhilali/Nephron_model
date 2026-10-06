# Install script for directory: /Users/alhilaba/UG4_promesh_ogrid/plugins/ConvectionDiffusion

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

if(CMAKE_INSTALL_COMPONENT STREQUAL "plugins" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin/plugins" TYPE SHARED_LIBRARY FILES "/Users/alhilaba/UG4_promesh_ogrid/bin/plugins/libConvectionDiffusion.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/libConvectionDiffusion.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/libConvectionDiffusion.dylib")
    execute_process(COMMAND "/usr/bin/install_name_tool"
      -id "libConvectionDiffusion.dylib"
      -change "/Users/alhilaba/UG4_promesh_ogrid/lib/libCompileInfo.dylib" "libCompileInfo.dylib"
      -change "/Users/alhilaba/UG4_promesh_ogrid/lib/libug4.dylib" "libug4.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/libConvectionDiffusion.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/alhilaba/UG4_promesh_ogrid/ugcore/cmake/../../lib"
      -add_rpath "$ORIGIN/../../lib/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/libConvectionDiffusion.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/plugins/libConvectionDiffusion.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "plugins" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/alhilaba/UG4_promesh_ogrid/build/ugcore/plugins/ConvectionDiffusion/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
