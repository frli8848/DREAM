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

add_test(test_julia_dreamline_run julia "${PROJECT_SOURCE_DIR}/tests/test_dreamline.jl")
add_test(test_julia_dreamcirc_run julia "${PROJECT_SOURCE_DIR}/tests/test_dreamcirc.jl")
add_test(test_julia_dreamcirc_f_run julia "${PROJECT_SOURCE_DIR}/tests/test_dreamcirc_f.jl")
add_test(test_julia_dreamrect_run julia "${PROJECT_SOURCE_DIR}/tests/test_dreamrect.jl")
add_test(test_julia_dreamrect_f_run julia "${PROJECT_SOURCE_DIR}/tests/test_dreamrect_f.jl")
add_test(test_julia_dreamsphere_run julia "${PROJECT_SOURCE_DIR}/tests/test_dreamsphere.jl")

add_test(test_julia_dream_arr_cylind_run julia "${PROJECT_SOURCE_DIR}/tests/test_dream_arr_cylind.jl")
add_test(test_julia_dream_arr_circ_run julia "${PROJECT_SOURCE_DIR}/tests/test_dream_arr_circ.jl")
add_test(test_julia_dream_arr_rect_run julia "${PROJECT_SOURCE_DIR}/tests/test_dream_arr_rect.jl")
add_test(test_julia_dreamcylind_run julia "${PROJECT_SOURCE_DIR}/tests/test_dreamcylind.jl")
add_test(test_julia_dream_arr_annu_run julia "${PROJECT_SOURCE_DIR}/tests/test_dream_arr_annu.jl")

add_test(test_julia_conv_run julia "${PROJECT_SOURCE_DIR}/tests/test_conv.jl")

add_test(test_julia_das_saft_run julia "${PROJECT_SOURCE_DIR}/tests/test_das_saft.jl")
add_test(test_julia_das_tfm_run julia "${PROJECT_SOURCE_DIR}/tests/test_das_tfm.jl")
add_test(test_julia_das_rca_run julia "${PROJECT_SOURCE_DIR}/tests/test_das_rca.jl")
