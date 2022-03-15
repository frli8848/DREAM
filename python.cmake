#
# Copyright (C) 2022 Fredrik Lingvall
#

find_package (pybind11 CONFIG)

if (pybind11_FOUND)

  set (DREAM_PYTHON_FLAGS "-DDREAM_PYTHON")

  message (STATUS "Python_INCLUDE_DIRS ${pybind11_INCLUDE_DIR}")
  message (STATUS "Python_LIBRARIES  ${pybind11_LIBRARIES}")
  message (STATUS "PYTHON_LD_FLAGS ${PYTHON_LD_FLAGS}")
  if (FFTW_FOUND)
    set (PYTHON_LD_FLAGS "${PYTHON_LD_FLAGS} -lfftw3")
  endif (FFTW_FOUND)
endif (pybind11_FOUND)


#
# dream_arr_rect
#

set (py_dream_arr_rect_SOURCE_FILES
  python/py_dream_arr_rect.cc
  src/attenuation.cc
  src/affinity.cc
  src/arr_functions.cc
  src/dream_arr_rect.cc
  src/dreamrect.cc
  src/dream_error.cc
  )

add_library (py_dream_arr_rect MODULE
  ${py_dream_arr_rect_SOURCE_FILES}
  )

target_link_libraries (py_dream_arr_rect PUBLIC pybind11::module
  ${FFTW_LIBRARIES}
  )

set_target_properties (py_dream_arr_rect PROPERTIES
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_PYTHON_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_PYTNON_INCLUDE_DIRS}"
  LINK_FLAGS ${PYTHON_LD_FLAGS}
  OUTPUT_NAME "dream_arr_rect")
