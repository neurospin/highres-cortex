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

set(_headers
  cortex.hh
  field.hh
  front.hh
  advection.hh
  advection.tcc
  cortex_advection.hh
  propagate_along_field.hh
  propagate_along_field.tcc
  label_volume.hh
  label_volume.tcc
  iterative_region_merger.hh
  iterative_region_merger.tcc
  cortex_column_region_quality.hh
  cortex_column_region_quality.tcc
  laplace_solver.hh
  laplace_solver.tcc
  upwinding.hh
  volume_util.hh
  volume_util.tcc
)
set(_sources
  field.cc
  front.cc
  advection.cc
  cortex_advection.cc
  propagate_along_field.cc
  label_volume.cc
  iterative_region_merger.cc
  cortex_column_region_quality.cc
  laplace_solver.cc
  upwinding.cc
  volume_util.cc
)

BRAINVISA_COPY_AND_INSTALL_HEADERS(_headers ${PROJECT_NAME} headersTarget)

include_directories(${GSL_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

add_library(${PROJECT_NAME} ${_sources})
add_dependencies(${PROJECT_NAME} ${headersTarget})

if(OPENMP_FOUND)
  set_source_files_properties(${_sources} PROPERTIES
    COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
  # This is needed for GCC, untested with any other compiler.
  set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS ${OpenMP_CXX_FLAGS})
endif()

set_property(TARGET ${PROJECT_NAME} PROPERTY VERSION ${${PROJECT_NAME}_VERSION})
target_link_libraries(${PROJECT_NAME} ${AIMS-FREE_LIBRARIES} ${GSL_LIBRARIES})
BRAINVISA_INSTALL(TARGETS ${PROJECT_NAME}
                  DESTINATION lib
                  COMPONENT ${PROJECT_NAME})

if(DOXYGEN_FOUND)
  set( DOXYFILE_PREDEFINED "${AIMS_DEFINITIONS}")
  # Enable cross-linking to AIMS documentation (only works for in-tree build)
  set(AIMS-FREE_DOC_DIR "${AIMS-FREE_BINARY_DIR}/share/doc")
  set(AIMS-FREE_VERSION_MAJMIN
    "${aims-free_VERSION_MAJOR}.${aims-free_VERSION_MINOR}")
  set(SOMA-IO_DOC_DIR "${SOMA-IO_BINARY_DIR}/share/doc")
  set(SOMA-IO_VERSION_MAJMIN "${soma-io_VERSION_MAJOR}.${soma-io_VERSION_MINOR}")

  # TODO: change this when new website gets in production
  set(BV_DOC_URL "http://brainvisa.fr/" CACHE STRING
    "Doxygen documentation links to the BrainVISA documentation at this URL")

  set(DOXYFILE_TAGFILES "${AIMS-FREE_DOC_DIR}/aimsdata-${AIMS-FREE_VERSION_MAJMIN}/doxygen/aimsdata.tag=${BV_DOC_URL}aimsdata-${AIMS-FREE_VERSION_MAJMIN}/doxygen ${AIMS-FREE_DOC_DIR}/aimsalgo-${AIMS-FREE_VERSION_MAJMIN}/doxygen/aimsalgo.tag=${BV_DOC_URL}aimsalgo-${AIMS-FREE_VERSION_MAJMIN}/doxygen ${AIMS-FREE_DOC_DIR}/cartodata-${AIMS-FREE_VERSION_MAJMIN}/doxygen/cartodata.tag=${BV_DOC_URL}cartodata-${AIMS-FREE_VERSION_MAJMIN}/doxygen ${SOMA-IO_DOC_DIR}/cartobase-${SOMA-IO_VERSION_MAJMIN}/doxygen/cartobase.tag=${BV_DOC_URL}cartobase-${SOMA-IO_VERSION_MAJMIN}/doxygen")

  BRAINVISA_GENERATE_DOXYGEN_DOC(_headers
                                 INPUT_PREFIX "${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}")

  # Fix warning for out-of-tree build
  cmake_policy(SET CMP0046 OLD)
  add_dependencies(${PROJECT_NAME}-doxygen aimsdata-doxygen aimsalgo-doxygen cartobase-doxygen cartodata-doxygen)
endif(DOXYGEN_FOUND)
