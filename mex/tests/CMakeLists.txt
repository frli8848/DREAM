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

add_test(test_matlab_dreamline_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_dreamline)
add_test(test_matlab_dreamcirc_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_dreamcirc)
add_test(test_matlab_dreamcirc_f_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_dreamcirc_f)
add_test(test_matlab_dreamrect_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_dreamrect)
add_test(test_matlab_dreamrect_f_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_dreamrect_f)
add_test(test_matlab_dreamsphere_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_dreamsphere)

add_test(test_matlab_dream_arr_cylind_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_dream_arr_cylind)
add_test(test_matlab_dream_arr_circ_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_dream_arr_circ)
add_test(test_matlab_dream_arr_rect_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_dream_arr_rect)
add_test(test_matlab_dreamcylind_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_dreamcylind)
add_test(test_matlab_dream_arr_annu_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_dream_arr_annu)

# These fails on MATLAB due to FFTW issues unless one do something like:
#
# $ export LD_PRELOAD=/usr/lib64/libfftw3.so

add_test(test_matlab_conv_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_conv)

add_test(test_matlab_das_saft_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_das_saft)
add_test(test_matlab_das_tfm_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_das_tfm)
add_test(test_matlab_das_tfm_run matlab -sd "${PROJECT_SOURCE_DIR}/../tests/m_files/" -batch test_das_rca)
