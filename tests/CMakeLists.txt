# Copyright Forschungszentrum Jülich GmbH (2017, 2018).
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
  "${AIMS-FREE_BINARY_DIR}"
  )

set(TEST_ENV "")
if(NOT PATH_TO_BRAINVISA_TREE STREQUAL "")
  message(STATUS "Detected a build outside of the BrainVISA tree (PATH_TO_BRAINVISA_TREE=${PATH_TO_BRAINVISA_TREE}). Activating the test_env.sh wrapper.")
  if(WIN32)
    message(WARNING "Tests for ${PROJECT_NAME} will fail on Windows when building outside of the BrainVISA tree. You are welcome to contribute a test_env.bat wrapper.")
  else()
    configure_file(test_env.sh.in test_env.sh @ONLY)
    set(TEST_ENV "sh" "-e" "${CMAKE_CURRENT_BINARY_DIR}/test_env.sh")
  endif()
endif()

# These tests need /bin/bash, exclude them from Windows cross-builds
if(NOT WIN32)
  # We use an extra "sh -e" prefix here, because the test is called through
  # bv_env_test, which calls the *basename* of the first argument as the
  # command. Without this trick, the basename is "test_all_scripts.sh", which
  # cannot be found on the PATH.
  brainvisa_add_test(
    NAME highres_cortex_example_scripts
    COMMAND ${TEST_ENV} "sh" "-e" "${CMAKE_CURRENT_SOURCE_DIR}/test_all_scripts.sh"
    )
endif()

brainvisa_add_test(
  NAME highres_cortex_capsul
  COMMAND ${TEST_ENV} "${PYTHON_EXECUTABLE_NAME}" -m highres_cortex.test.test_capsul
  )
