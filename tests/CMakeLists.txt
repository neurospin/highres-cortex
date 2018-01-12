# Copyright Forschungszentrum Jülich GmbH (2018).
# Copyright Forschungszentrum Jülich GmbH (2017).
#
# Contributor: Yann Leprince <yann.leprince@ylep.fr>.
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved. This file is offered as-is,
# without any warranty.

# Determine if we are compiling in the BrainVISA build tree. If not, set up
# TEST_ENV to point to a wrapper script which sets the environment.
file(RELATIVE_PATH
  PATH_TO_BRAINVISA_TREE
  "${CMAKE_BINARY_DIR}"
  "${brainvisa-cmake_DIR}/../../.."
  )
if(PATH_TO_BRAINVISA_TREE STREQUAL "")
  set(TEST_ENV "")
else()
  configure_file(test_env.sh.in test_env.sh @ONLY)
  set(TEST_ENV "sh" "-e" "${CMAKE_CURRENT_BINARY_DIR}/test_env.sh")
endif()

# These tests need /bin/bash, exclude them from Windows cross-builds
if(NOT WIN32)
  brainvisa_add_test(
    NAME highres_cortex_all_scripts
    COMMAND ${TEST_ENV} "${CMAKE_CURRENT_SOURCE_DIR}/test_all_scripts.sh"
    )
endif()

if(NOT DEFINED PYTHON_HOST_EXECUTABLE)
  # Compatibility with BrainVISA <= 4.5
  set(PYTHON_HOST_EXECUTABLE "${PYTHON_EXECUTABLE}")
endif()

brainvisa_add_test(
  NAME highres_cortex_capsul
  COMMAND ${TEST_ENV} "${PYTHON_HOST_EXECUTABLE}" -m highres_cortex.test.test_capsul
  )
