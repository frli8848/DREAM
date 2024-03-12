#
# Copyright (C) 2023,2024 Fredrik Lingvall
#

project(py_dream)

# Hack to make cmake use the active python version on Gentoo
find_package(Python 3.11 COMPONENTS Interpreter Development REQUIRED EXACT)
#find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# Strip lead and trailing whitepasce
#string(STRIP "${OCTAVE_INCLUDE_DIRS}" OCTAVE_INCLUDE_DIRS)

set (DREAM_PYTHON_INCLUDE_DIRS
  ${PROJECT_BINARY_DIR};
  "${PROJECT_BINARY_DIR}/../opencl;"
  "${PROJECT_SOURCE_DIR}/../include;"
  "${PROJECT_SOURCE_DIR}/include;"
  ${FFTW_INCLUDE_DIR};
)
message (STATUS "DREAM_PYTHON_INCLUDE_DIRS: ${DREAM_PYTHON_INCLUDE_DIRS}")

if (pybind11_FOUND)

  set (DREAM_PYTHON_FLAGS "-DDREAM_PYTHON")

  message (STATUS "Python_INCLUDE_DIRS ${pybind11_INCLUDE_DIR}")
  message (STATUS "PYTHON_LD_FLAGS ${PYTHON_LD_FLAGS}")

  if (FFTW_FOUND)
    add_definitions( -DHAVE_FFTW )
    set (PYTHON_LD_FLAGS "${PYTHON_LD_FLAGS} -lfftw3")
  endif (FFTW_FOUND)

endif (pybind11_FOUND)

# Single element transducers

#
# dreamline
#

set (py_dreamline_SOURCE_FILES
  src/py_dreamline.cc
  ../src/dreamline.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(py_dreamline MODULE ${py_dreamline_SOURCE_FILES})
target_link_libraries(py_dreamline pybind11::module
  ${FFTW_LIBRARIES}
)

set_target_properties(py_dreamline PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamline"
)

#
# dreamrect
#

set (py_dreamrect_SOURCE_FILES
  src/py_dreamrect.cc
  ../src/dreamrect.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(py_dreamrect MODULE ${py_dreamrect_SOURCE_FILES})
target_link_libraries(py_dreamrect pybind11::module
  ${FFTW_LIBRARIES}
)

set_target_properties(py_dreamrect PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamrect"
)

#
# dreamrect_f
#

set (py_dreamrect_f_SOURCE_FILES
  src/py_dreamrect_f.cc
  ../src/dreamrect_f.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(py_dreamrect_f MODULE ${py_dreamrect_f_SOURCE_FILES})
target_link_libraries(py_dreamrect_f pybind11::module
  ${FFTW_LIBRARIES}
)

set_target_properties(py_dreamrect_f PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamrect_f"
)

#
# dreamcirc
#

set (py_dreamcirc_SOURCE_FILES
  src/py_dreamcirc.cc
  ../src/dreamcirc.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(py_dreamcirc MODULE ${py_dreamcirc_SOURCE_FILES})
target_link_libraries(py_dreamcirc pybind11::module
  ${FFTW_LIBRARIES}
)

set_target_properties(py_dreamcirc PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamcirc"
)

#
# dreamcirc_f
#

set (py_dreamcirc_f_SOURCE_FILES
  src/py_dreamcirc_f.cc
  ../src/dreamcirc_f.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(py_dreamcirc_f MODULE ${py_dreamcirc_f_SOURCE_FILES})
target_link_libraries(py_dreamcirc_f pybind11::module
  ${FFTW_LIBRARIES}
)

set_target_properties(py_dreamcirc_f PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamcirc_f"
)

#
# dreamsphere
#

set (py_dreamsphere_SOURCE_FILES
  src/py_dreamsphere.cc
  ../src/dreamsphere.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(py_dreamsphere MODULE ${py_dreamsphere_SOURCE_FILES})
target_link_libraries(py_dreamsphere pybind11::module
  ${FFTW_LIBRARIES}
)

set_target_properties(py_dreamsphere PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamsphere"
)

#
# dreamcylind
#

set (py_dreamcylind_SOURCE_FILES
  src/py_dreamcylind.cc
  ../src/dreamcylind.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(py_dreamcylind MODULE ${py_dreamcylind_SOURCE_FILES})
target_link_libraries(py_dreamcylind pybind11::module
  ${FFTW_LIBRARIES}
)

set_target_properties(py_dreamcylind PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamcylind"
)

# Arrays

#
# dream_arr_rect
#

set (py_dream_arr_rect_SOURCE_FILES
  src/py_dream_arr_rect.cc
  ../src/dream_arr_rect.cc
  ../src/dreamrect.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

#Python_add_library(py_dream_arr_rect ${py_dream_arr_rect_SOURCE_FILES})
add_library(py_dream_arr_rect MODULE ${py_dream_arr_rect_SOURCE_FILES})
target_link_libraries(py_dream_arr_rect pybind11::module
  ${FFTW_LIBRARIES}
)

set_target_properties(py_dream_arr_rect PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dream_arr_rect"
)

#
# dream_arr_circ
#

set (py_dream_arr_circ_SOURCE_FILES
  src/py_dream_arr_circ.cc
  ../src/dream_arr_circ.cc
  ../src/dreamcirc.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(py_dream_arr_circ MODULE ${py_dream_arr_circ_SOURCE_FILES})
target_link_libraries(py_dream_arr_circ pybind11::module
  ${FFTW_LIBRARIES}
)

set_target_properties(py_dream_arr_circ PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dream_arr_circ"
)

#
# dream_arr_cylind
#

set (py_dream_arr_cylind_SOURCE_FILES
  src/py_dream_arr_cylind.cc
  ../src/dream_arr_cylind.cc
  ../src/dreamcylind.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(py_dream_arr_cylind MODULE ${py_dream_arr_cylind_SOURCE_FILES})
target_link_libraries(py_dream_arr_cylind pybind11::module
  ${FFTW_LIBRARIES}
)

set_target_properties(py_dream_arr_cylind PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dream_arr_cylind"
)

#
# dream_arr_annu
#

set (py_dream_arr_annu_SOURCE_FILES
  src/py_dream_arr_annu.cc
  ../src/dream_arr_annu.cc
  ../src/dreamcirc.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(py_dream_arr_annu MODULE ${py_dream_arr_annu_SOURCE_FILES})
target_link_libraries(py_dream_arr_annu pybind11::module
  ${FFTW_LIBRARIES}
)

set_target_properties(py_dream_arr_annu PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dream_arr_annu"
)


#
# fftconv_p
#

set (py_fftconv_p_SOURCE_FILES
  src/py_fftconv_p.cc
  ../extra_src/fftconv_p.cc
  ../extra_src/fftconv.cc
  ../src/affinity.cc
  ../src/dream_error.cc
)

add_library(py_fftconv_p MODULE ${py_fftconv_p_SOURCE_FILES})
target_link_libraries(py_fftconv_p pybind11::module)

set_target_properties(py_fftconv_p PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "fftconv_p"
)

#
# das
#

# Single

set (py_das_f_SOURCE_FILES
  src/py_das.cc
  ../extra_src/das.cc
  ../src/affinity.cc
  ../src/dream_error.cc
)

if (OpenCL_FOUND)
  set (py_das_f_SOURCE_FILES ${py_das_f_SOURCE_FILES} "../opencl/cl_das.cc")
endif (OpenCL_FOUND)

add_library(py_das_f MODULE ${py_das_f_SOURCE_FILES})
target_link_libraries(py_das_f pybind11::module
  ${OpenCL_LIBRARIES}
)

set_target_properties(py_das_f PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS} -DUSE_FLOAT"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "das_f"
)

# Double

set (py_das_SOURCE_FILES
  src/py_das.cc
  ../extra_src/das.cc
  ../src/affinity.cc
  ../src/dream_error.cc
)

if (OpenCL_FOUND)
  set (py_das_SOURCE_FILES ${py_das_SOURCE_FILES} "../opencl/cl_das.cc")
endif (OpenCL_FOUND)

add_library(py_das MODULE ${py_das_SOURCE_FILES})
target_link_libraries(py_das pybind11::module
  ${OpenCL_LIBRARIES}
)

set_target_properties(py_das PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTHON_INCLUDE_DIRS}"
  INTERPROCEDURAL_OPTIMIZATION ON
  #CXX_VISIBILITY_PRESET ON
  #VISIBILITY_INLINES_HIDDEN ON
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "das"
)
