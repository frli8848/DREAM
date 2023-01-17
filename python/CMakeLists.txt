#
# Copyright (C) 2023 Fredrik Lingvall
#

project(py_dream)

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

find_package(FFTW)

# Strip lead and trailing whitepasce
#string(STRIP "${OCTAVE_INCLUDE_DIRS}" OCTAVE_INCLUDE_DIRS)

set (DREAM_PYTHON_INCLUDE_DIRS
  "${PROJECT_BINARY_DIR};"
  "${PROJECT_SOURCE_DIR}/../include;"
  "${FFTW_INCLUDE_DIR};"
  #"${OCTAVE_INCLUDE_DIRS}"
)
message (STATUS "DREAM_PYTHON_INCLUDE_DIRS: ${DREAM_PYTHON_INCLUDE_DIRS}")


if (pybind11_FOUND)

  set (DREAM_PYTHON_FLAGS "-DDREAM_PYTHON")

  message (STATUS "Python_INCLUDE_DIRS ${pybind11_INCLUDE_DIR}")
  message (STATUS "PYTHON_LD_FLAGS ${PYTHON_LD_FLAGS}")

  if (FFTW_FOUND)
    set (PYTHON_LD_FLAGS "${PYTHON_LD_FLAGS} -lfftw3")
  endif (FFTW_FOUND)

endif (pybind11_FOUND)

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
