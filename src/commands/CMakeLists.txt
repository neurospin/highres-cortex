# Copyright Forschungszentrum Jülich GmbH (2016).
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

set(OpenMP_commands
  ylLaplacian
  ylIsoCurvature
)
set(_commands
  ${OpenMP_commands}
  ylAdvectEuclidean
  ylAdvectPath
  ylAdvectTubes
  ylAdvectValues
  ylPropagateAlongField
  ylLabelEachVoxel
  ylMergeCortexColumnRegions
  ylMakeTraversePseudoAreaMap
  ylUpwindDistance
)

include_directories("${CMAKE_BINARY_DIR}/include")

foreach(_command ${_commands})
  BRAINVISA_ADD_EXECUTABLE(${_command} ${_command}.cc)
  target_link_libraries(${_command} ${PROJECT_NAME} ${AIMS_AIMSDATA_LIBRARIES})
  BRAINVISA_INSTALL(TARGETS ${_command}
                    DESTINATION bin
                    COMPONENT ${PROJECT_NAME})
endforeach()

if(OPENMP_FOUND)
  set_source_files_properties(${OpenMP_commands} PROPERTIES
    COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
  if(OpenMP_CXX_LIBRARIES)
    foreach(_command ${OpenMP_commands})
      target_link_libraries(${_command} ${OpenMP_CXX_LIBRARIES})
    endforeach()
  else()
    # Old versions of FindOpenMP.cmake (before CMake 3.9) do not define
    # _LIBRARIES. This works with GCC and clang.
    set_target_properties(${OpenMP_commands} PROPERTIES
      LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endif()
