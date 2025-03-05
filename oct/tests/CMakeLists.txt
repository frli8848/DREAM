cmake_minimum_required(VERSION 3.10)

message (STATUS "Current CTest bin dir: ${CMAKE_CURRENT_BINARY_DIR}")
message (STATUS "Current CTest PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

#
# Platform specific settings - typically include directories and link libraries.
#
#include("${CMAKE_HOST_SYSTEM_NAME}.cmake")

#
# User specific settings - e.g. to override platform file locations or play with compiler options.
#
#include("user.cmake" OPTIONAL)

#
# Octave tests
#

add_test(test_octave_dreamline_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_dreamline)
add_test(test_octave_dreamcirc_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_dreamcirc)
add_test(test_octave_dreamcirc_f_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_dreamcirc_f)
add_test(test_octave_dreamrect_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_dreamrect)
add_test(test_octave_dreamrect_f_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_dreamrect_f)
add_test(test_octave_dreamsphere_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_dreamsphere)

add_test(test_octave_dream_arr_cylind_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_dream_arr_cylind)
add_test(test_octave_dream_arr_circ_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_dream_arr_circ)
add_test(test_octave_dream_arr_rect_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_dream_arr_rect)
add_test(test_octave_dreamcylind_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_dreamcylind)
add_test(test_octave_dream_arr_annu_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_dream_arr_annu)

add_test(test_octave_conv_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_conv)

add_test(test_octave_das_saft_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_das_saft)
add_test(test_octave_das_tfm_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_das_tfm)
add_test(test_octave_das_tfm_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_das_rca)

if (USE_OPENCL)
  add_test(test_octave_opencl_run octave --path "${PROJECT_SOURCE_DIR}/../tests/m_files/" --eval test_opencl)
endif (USE_OPENCL)
