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

set(_commands
  ylAdvectEuclidean
  ylAdvectTubes
  ylPropagateAlongField
  ylLabelEachVoxel
  ylLaplacian
  ylMergeCortexColumnRegions
  ylMakeTraversePseudoAreaMap
)

include_directories("${CMAKE_BINARY_DIR}/include")

foreach(_command ${_commands})
  add_executable(${_command} ${_command}.cc)
  target_link_libraries(${_command} ${PROJECT_NAME} ${AIMS-FREE_LIBRARIES})
  BRAINVISA_INSTALL(TARGETS ${_command}
                    DESTINATION bin
                    COMPONENT ${PROJECT_NAME})
endforeach()

set(OpenMP_commands
  ylLaplacian
)
if(OPENMP_FOUND)
  set_source_files_properties(${OpenMP_commands} PROPERTIES
    COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
  # This is needed for GCC, untested with any other compiler.
  set_target_properties(${OpenMP_commands} PROPERTIES
    LINK_FLAGS ${OpenMP_CXX_FLAGS})
endif()
