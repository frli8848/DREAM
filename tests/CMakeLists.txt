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

# Check that we have FFTW installed
find_package (FFTW)  # Will define FFTW_INCLUDE_DIR, FFTW_LIBRARIES, FFTW_STATIC_LIBRARIES, and FFTW_FOUND

set (UNIT_TEST_LIST
  test_fftconv_p
)

set (DREAM_SRC
  "${PROJECT_SOURCE_DIR}/src/affinity.cc"
  "${PROJECT_SOURCE_DIR}/extra_src/fftconv.cc"
)



# Use the list to make executables and what to link with (the dependency)
foreach(NAME IN LISTS UNIT_TEST_LIST)

  add_executable(${NAME} src/${NAME}.cc ${DREAM_SRC})

  set_target_properties (${NAME} PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/include;"
    )

  target_link_libraries(${NAME} PUBLIC "${FFTW_LIBRARIES}")

endforeach()

#
# C++ tests
#

# fftconv_p
add_test(test_fftconv_p_run ${CMAKE_CURRENT_BINARY_DIR}/test_fftconv_p)

#
# Octave tests
#

#"C:\Program Files\GNU Octave\Octave-8.4.0\octave-launch.exe" --gui

if (BUILD_OCT)

  add_test(test_octave_dreamline_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_dreamline)
  add_test(test_octave_dreamcirc_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_dreamcirc)
  add_test(test_octave_dreamcirc_f_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_dreamcirc_f)
  add_test(test_octave_dreamrect_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_dreamrect)
  add_test(test_octave_dreamrect_f_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_dreamrect_f)
  add_test(test_octave_dreamsphere_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_dreamsphere)

  add_test(test_octave_dream_arr_cylind_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_dream_arr_cylind)
  add_test(test_octave_dream_arr_circ_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_dream_arr_circ)
  add_test(test_octave_dream_arr_rect_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_dream_arr_rect)
  add_test(test_octave_dreamcylind_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_dreamcylind)
  add_test(test_octave_dream_arr_annu_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_dream_arr_annu)

  add_test(test_octave_conv_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_conv)

  add_test(test_octave_das_saft_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_das_saft)
  add_test(test_octave_das_tfm_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_das_tfm)
  add_test(test_octave_das_tfm_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_das_rca)

  if (USE_OPENCL)
    add_test(test_octave_opencl_run octave --path "${PROJECT_SOURCE_DIR}/tests/m_files/" --eval test_opencl)
  endif (USE_OPENCL)

endif (BUILD_OCT)

#
# MATLAB tests
#

# NB. We need a relatively new MATLAB version that has the -batch flag.

if (BUILD_MEX)

  add_test(test_matlab_dreamline_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_dreamline)
  add_test(test_matlab_dreamcirc_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_dreamcirc)
  add_test(test_matlab_dreamcirc_f_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_dreamcirc_f)
  add_test(test_matlab_dreamrect_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_dreamrect)
  add_test(test_matlab_dreamrect_f_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_dreamrect_f)
  add_test(test_matlab_dreamsphere_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_dreamsphere)

  add_test(test_matlab_dream_arr_cylind_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_dream_arr_cylind)
  add_test(test_matlab_dream_arr_circ_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_dream_arr_circ)
  add_test(test_matlab_dream_arr_rect_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_dream_arr_rect)
  add_test(test_matlab_dreamcylind_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_dreamcylind)
  add_test(test_matlab_dream_arr_annu_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_dream_arr_annu)

  # These fails on MATLAB due to FFTW issues unless one do something like:
  #
  # $ export LD_PRELOAD=/usr/lib64/libfftw3.so

  add_test(test_matlab_conv_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_conv)

  add_test(test_matlab_das_saft_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_das_saft)
  add_test(test_matlab_das_tfm_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_das_tfm)
  add_test(test_matlab_das_tfm_run matlab -sd "${PROJECT_SOURCE_DIR}/tests/m_files/" -batch test_das_rca)

endif (BUILD_MEX)

#
# Python tests
#

if (BUILD_PYTHON)

  add_test(test_python_dreamline_run python "${PROJECT_SOURCE_DIR}/python/tests/test_dreamline.py")
  add_test(test_python_dreamcirc_run python "${PROJECT_SOURCE_DIR}/python/tests/test_dreamcirc.py")
  add_test(test_python_dreamcirc_f_run python "${PROJECT_SOURCE_DIR}/python/tests/test_dreamcirc_f.py")
  add_test(test_python_dreamrect_run python "${PROJECT_SOURCE_DIR}/python/tests/test_dreamrect.py")
  add_test(test_python_dreamrect_f_run python "${PROJECT_SOURCE_DIR}/python/tests/test_dreamrect_f.py")
  add_test(test_python_dreamsphere_run python "${PROJECT_SOURCE_DIR}/python/tests/test_dreamsphere.py")

  add_test(test_python_dream_arr_cylind_run python "${PROJECT_SOURCE_DIR}/python/tests/test_dream_arr_cylind.py")
  add_test(test_python_dream_arr_circ_run python "${PROJECT_SOURCE_DIR}/python/tests/test_dream_arr_circ.py")
  add_test(test_python_dream_arr_rect_run python "${PROJECT_SOURCE_DIR}/python/tests/test_dream_arr_rect.py")
  add_test(test_python_dreamcylind_run python "${PROJECT_SOURCE_DIR}/python/tests/test_dreamcylind.py")
  add_test(test_python_dream_arr_annu_run python "${PROJECT_SOURCE_DIR}/python/tests/test_dream_arr_annu.py")

  add_test(test_python_conv_run python "${PROJECT_SOURCE_DIR}/python/tests/test_conv.py")

  add_test(test_python_das_saft_run python "${PROJECT_SOURCE_DIR}/python/tests/test_das_saft.py")
  add_test(test_python_das_tfm_run python "${PROJECT_SOURCE_DIR}/python/tests/test_das_tfm.py")
  add_test(test_python_das_rca_run python "${PROJECT_SOURCE_DIR}/python/tests/test_das_rca.py")

endif (BUILD_PYTHON)

#
# Julia tests
#

if (BUILD_JULIA)

  add_test(test_julia_dreamline_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_dreamline.jl")
  add_test(test_julia_dreamcirc_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_dreamcirc.jl")
  add_test(test_julia_dreamcirc_f_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_dreamcirc_f.jl")
  add_test(test_julia_dreamrect_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_dreamrect.jl")
  add_test(test_julia_dreamrect_f_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_dreamrect_f.jl")
  add_test(test_julia_dreamsphere_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_dreamsphere.jl")

  add_test(test_julia_dream_arr_cylind_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_dream_arr_cylind.jl")
  add_test(test_julia_dream_arr_circ_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_dream_arr_circ.jl")
  add_test(test_julia_dream_arr_rect_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_dream_arr_rect.jl")
  add_test(test_julia_dreamcylind_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_dreamcylind.jl")
  add_test(test_julia_dream_arr_annu_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_dream_arr_annu.jl")

  add_test(test_julia_conv_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_conv.jl")

  add_test(test_julia_das_saft_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_das_saft.jl")
  add_test(test_julia_das_tfm_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_das_tfm.jl")
  add_test(test_julia_das_rca_run julia "${PROJECT_SOURCE_DIR}/julia/tests/test_das_rca.jl")

endif (BUILD_JULIA)
