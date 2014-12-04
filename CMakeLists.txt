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
BRAINVISA_PROJECT()

BRAINVISA_FIND_PACKAGE(aims-free REQUIRED)
include("${AIMS-FREE_USE_FILE}")

find_package(Boost 1.49 REQUIRED)
if(${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION} VERSION_LESS 1.49)
    message(FATAL_ERROR "Boost version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION} found - need at least 1.49")
endif()

BRAINVISA_DEPENDENCY(DEV DEPENDS libboost DEV)

find_package(GSL REQUIRED)

find_package(python REQUIRED)

BRAINVISA_COPY_PYTHON_DIRECTORY("${CMAKE_CURRENT_SOURCE_DIR}/python"
                                 ${PROJECT_NAME})

add_subdirectory(src)
