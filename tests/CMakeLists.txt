# Copyright Forschungszentrum Jülich GmbH (2017).
#
# Contributor: Yann Leprince <yann.leprince@ylep.fr>.
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved. This file is offered as-is,
# without any warranty.

configure_file(test_env.sh.in test_env.sh @ONLY)
set(TEST_ENV "${CMAKE_CURRENT_BINARY_DIR}/test_env.sh")

brainvisa_add_test(
  NAME highres_cortex_all_scripts
  COMMAND /bin/sh -e "${TEST_ENV}" "${CMAKE_CURRENT_SOURCE_DIR}/test_all_scripts.sh"
  )

if(NOT DEFINED PYTHON_HOST_EXECUTABLE)
  # Compatibility with BrainVISA <= 4.5
  set(PYTHON_HOST_EXECUTABLE "${PYTHON_EXECUTABLE}")
endif()
brainvisa_add_test(
  NAME highres_cortex_capsul
  COMMAND "${PYTHON_HOST_EXECUTABLE}" -m highres_cortex.test.capsul_tests
)
