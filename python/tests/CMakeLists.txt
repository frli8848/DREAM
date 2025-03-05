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

add_test(test_python_dreamline_run python "${PROJECT_SOURCE_DIR}/tests/test_dreamline.py")
add_test(test_python_dreamcirc_run python "${PROJECT_SOURCE_DIR}/tests/test_dreamcirc.py")
add_test(test_python_dreamcirc_f_run python "${PROJECT_SOURCE_DIR}/tests/test_dreamcirc_f.py")
add_test(test_python_dreamrect_run python "${PROJECT_SOURCE_DIR}/tests/test_dreamrect.py")
add_test(test_python_dreamrect_f_run python "${PROJECT_SOURCE_DIR}/tests/test_dreamrect_f.py")
add_test(test_python_dreamsphere_run python "${PROJECT_SOURCE_DIR}/tests/test_dreamsphere.py")

add_test(test_python_dream_arr_cylind_run python "${PROJECT_SOURCE_DIR}/tests/test_dream_arr_cylind.py")
add_test(test_python_dream_arr_circ_run python "${PROJECT_SOURCE_DIR}/tests/test_dream_arr_circ.py")
add_test(test_python_dream_arr_rect_run python "${PROJECT_SOURCE_DIR}/tests/test_dream_arr_rect.py")
add_test(test_python_dreamcylind_run python "${PROJECT_SOURCE_DIR}/tests/test_dreamcylind.py")
add_test(test_python_dream_arr_annu_run python "${PROJECT_SOURCE_DIR}/tests/test_dream_arr_annu.py")

add_test(test_python_conv_run python "${PROJECT_SOURCE_DIR}/tests/test_conv.py")

add_test(test_python_das_saft_run python "${PROJECT_SOURCE_DIR}/tests/test_das_saft.py")
add_test(test_python_das_tfm_run python "${PROJECT_SOURCE_DIR}/tests/test_das_tfm.py")
add_test(test_python_das_rca_run python "${PROJECT_SOURCE_DIR}/tests/test_das_rca.py")
