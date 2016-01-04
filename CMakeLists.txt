# Copyright Télécom ParisTech (2015).
# Copyright CEA (2014).
# Copyright Université Paris XI (2014).
#
# Contributor: Yann Leprince <yann.leprince@ylep.fr>.
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved. This file is offered as-is,
# without any warranty.

cmake_minimum_required(VERSION 2.6)
find_package(brainvisa-cmake REQUIRED)
BRAINVISA_PROJECT(CXX)

option(BUILD_SHARED_LIBS "Build highres-cortex as a shared library" ON)

BRAINVISA_FIND_PACKAGE(aims-free 4.5 REQUIRED)
include("${AIMS-FREE_USE_FILE}")
BRAINVISA_DEPENDENCY(DEV DEPENDS aims-free DEV ">= 4.5")
BRAINVISA_DEPENDENCY(RUN DEPENDS aims-free RUN ">= 4.5")


# Do not use BoostConfig.cmake from boost-cmake, because its behaviour may be
# different from regular FindBoost.cmake.
set(Boost_NO_BOOST_CMAKE ON)
# All suitable Boost versions should be included here, if they are not already
# contained in the FindBoost.cmake module provided with CMake. Hence, new
# versions should be added here as they are released...
set(Boost_ADDITIONAL_VERSIONS
  "1.60.0" "1.60" "1.59.0" "1.59" "1.58.0" "1.58" "1.57.0" "1.57"
  "1.56.0" "1.56" "1.55.0" "1.55" "1.54.0" "1.54" "1.53.0" "1.53" "1.52.0"
  "1.52" "1.51.0" "1.51" "1.50.0" "1.50" "1.49.0" "1.49")
find_package(Boost 1.49 REQUIRED)
# The version of Boost is not checked by find_package if it is in the cache.
if(${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION} VERSION_LESS 1.49)
    message(FATAL_ERROR
      "Unsuitable Boost version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}."
      "${Boost_SUBMINOR_VERSION} is already in the CMake cache. "
      "Please set the BOOST_ROOT variable to a directory containing "
      "Boost 1.49 or newer (under BOOST_ROOT/include).")
endif()
BRAINVISA_DEPENDENCY(DEV DEPENDS libboost DEV ">= 1.49")
# No run-time dependency as we use header-only Boost libraries


find_package(GSL REQUIRED)
BRAINVISA_DEPENDENCY(DEV DEPENDS libgsl DEV)
BRAINVISA_DEPENDENCY(RUN DEPENDS libgsl RUN)

# OpenMP is optional
find_package(OpenMP)

find_package(Doxygen)

find_package(python REQUIRED)
BRAINVISA_DEPENDENCY(DEV DEPENDS python DEV ">= 2.6; << 3.0")
BRAINVISA_DEPENDENCY(RUN RECOMMENDS python RUN ">= 2.6; << 3.0")
# No need for find_package(NumPy) as that is for C headers
BRAINVISA_DEPENDENCY(RUN RECOMMENDS python-numpy RUN)

# For fix_cortex_topology
BRAINVISA_DEPENDENCY(RUN SUGGESTS morphologist-private RUN ">= 4.5")


BRAINVISA_COPY_PYTHON_DIRECTORY("${CMAKE_CURRENT_SOURCE_DIR}/python"
                                 ${PROJECT_NAME})

BRAINVISA_COPY_DIRECTORY("${CMAKE_CURRENT_SOURCE_DIR}/bin"
                         bin
                         ${PROJECT_NAME})

add_subdirectory(src)
