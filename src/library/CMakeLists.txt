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
  field.hh
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
)
set(_sources
  field.cc
  advection.cc
  cortex_advection.cc
  propagate_along_field.cc
  label_volume.cc
  iterative_region_merger.cc
  cortex_column_region_quality.cc
)

BRAINVISA_COPY_AND_INSTALL_HEADERS(_headers ${PROJECT_NAME} headersTarget)

include_directories(${GSL_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

add_library(${PROJECT_NAME} SHARED ${_sources})
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

find_package(Doxygen)
if(DOXYGEN_FOUND)
  # Enable cross-linking to AIMS documentation
  set(aims-free_version "${aims-free_VERSION_MAJOR}.${aims-free_VERSION_MINOR}")
  set(DOXYFILE_TAGFILES "${CMAKE_BINARY_DIR}/share/doc/aimsdata-${aims-free_version}/doxygen/aimsdata.tag=../../aimsdata-${aims-free_version}/doxygen  ${CMAKE_BINARY_DIR}/share/doc/graph-${aims-free_version}/doxygen/graph.tag=../../graph-${aims-free_version}/doxygen  ${CMAKE_BINARY_DIR}/share/doc/cartobase-${aims-free_version}/doxygen/cartobase.tag=../../cartobase-${aims-free_version}/doxygen  ${CMAKE_BINARY_DIR}/share/doc/cartodata-${aims-free_version}/doxygen/cartodata.tag=../../cartodata-${aims-free_version}/doxygen  ${CMAKE_BINARY_DIR}/share/doc/pyaims-${aims-free_version}/doxygen/pyaims.tag=../../pyaims-${aims-free_version}/doxygen")

  BRAINVISA_GENERATE_DOXYGEN_DOC(_headers
                                 INPUT_PREFIX "${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}")
  add_dependencies(${PROJECT_NAME}-doxygen aimsdata-doxygen graph-doxygen cartobase-doxygen cartodata-doxygen pyaims-doxygen)
endif(DOXYGEN_FOUND)
