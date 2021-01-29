#
# Copyright (C) 2016, 2019, 2021 Fredrik Lingvall
#

project (dream-toolbox)

cmake_minimum_required(VERSION 3.10)

# Add path to the various cmake files to detect Matlab, Octave, FTTW, etc.
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Version number.
set (DREAM_VERSION_MAJOR 3)
set (DREAM_VERSION_MINOR 0)
set (DREAM_VERSION_PATCH 0)

option(BUILD_MEX "Enable building of the mex-files." OFF)
option(BUILD_OCT "Enable building of the oct-files." ON)
option(BUILD_USERMAN "Enable building of the user manual." OFF)
option(USE_OPENCL "Enable OpenCL support." OFF)

if (USE_OPENCL)

  find_package (OpenCL 1.2)

  if (OpenCL_FOUND)
    add_definitions(-DUSE_OPENCL)
  endif (OpenCL_FOUND)

endif (USE_OPENCL)

if (BUILD_OCT OR BUILD_MEX)

  find_package (FFTW REQUIRED)  # Will define FFTW_INCLUDE_DIR and FFTW_LIBRARIES
  message (STATUS "${FFTW_LIBRARIES}")

endif (BUILD_OCT OR BUILD_MEX)

if (BUILD_OCT)

  # From http://www.coolprop.org/coolprop/wrappers/Octave/index.html
  find_package (Octave)

  # Octave oct flags.
  set (OCT_FLAGS "-DOCTAVE") # So that octave_idx_type is used for matrix/vector indexing.
  if (UNIX)
    set (OCT_LD_FLAGS "${CMAKE_LD_FLAGS} ${OCTAVE_LINK_FLAGS} -pthread") # We need phread to use pthread_setaffinity_np.
  else (UNIX)
    set (OCT_LD_FLAGS "${CMAKE_LD_FLAGS} ${OCTAVE_LINK_FLAGS}")
  endif (UNIX)
endif (BUILD_OCT)

if (BUILD_MEX)
  find_package (Matlab)
  message (STATUS "Matlab_MEX_EXTENSION = ${Matlab_MEX_EXTENSION}")
  #message (STATUS "MATLAB_INCLUDE_DIR = ${MATLAB_INCLUDE_DIR}")
  #message (STATUS "MATLAB_LIBARARIES = ${MATLAB_LIBRARIES}")

  if (Matlab_FOUND)
    # mex libs
    set (MEX_LD_FLAGS "${CMAKE_LD_FLAGS} ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY} -lstdc++ -lfftw3")
  endif (Matlab_FOUND)

endif (BUILD_MEX)

#
# Common flags for the oct and mex interfaces.
#
#
# NB We use C99 bool end rint etc., and we use C++11 threads.

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

#
# Include paths
#

include_directories (
  ${PROJECT_BINARY_DIR}
  ${PROJECT_SOURCE_DIR}/include
  ${FFTW_INCLUDE_DIR}
  if (OCTAVE_FOUND)
    ${OCTAVE_INCLUDE_DIRS}
  endif (OCTAVE_FOUND)
  if (Matlab_FOUND)
    ${Matlab_INCLUDE_DIRS}
  endif (Matlab_FOUND)
  )

add_definitions( -DUSE_FFTW )	# Build with FFTW support
add_definitions( -DPARALLEL )	# att.cc needs this.

#
# dreamline
#

# oct

if (OCTAVE_FOUND)

  set (oct_dreamline_SOURCE_FILES
    oct/oct_dreamline.cc
    src/att.cc
    src/affinity.cc
    src/dreamline.c
    src/fft_octave.cc
    src/dream_error.c
    )

  add_library (oct_dreamline MODULE
    ${oct_dreamline_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamline
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dreamline PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamline")

endif (OCTAVE_FOUND)

# mex

if (Matlab_FOUND)

  set (mex_dreamline_SOURCE_FILES
    mex/mex_dreamline.cc
    src/att.cc
    src/affinity.cc
    src/dreamline.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dreamline MODULE
    ${mex_dreamline_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamline
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamline PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamline")

endif (Matlab_FOUND)

#
# dreamrect
#

if (OCTAVE_FOUND)

  set (oct_dreamrect_SOURCE_FILES
    oct/oct_dreamrect.cc
    src/dreamrect.c
    src/att.cc
    src/affinity.cc
    src/fft_octave.cc
    src/dream_error.c
    )

  if (OpenCL_FOUND)
    set (oct_dreamrect_SOURCE_FILES ${oct_dreamrect_SOURCE_FILES} "opencl/cl_dreamrect.cc")
  endif (OpenCL_FOUND)

  add_library (oct_dreamrect MODULE
    ${oct_dreamrect_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamrect
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${OpenCL_LIBRARIES}
    )

  set_target_properties (oct_dreamrect PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamrect")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dreamrect_SOURCE_FILES
    mex/mex_dreamrect.cc
    src/att.cc
    src/affinity.cc
    src/dreamrect.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dreamrect MODULE
    ${mex_dreamrect_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamrect
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamrect PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamrect")

endif (Matlab_FOUND)

#
# dreamrect_f
#

if (OCTAVE_FOUND)

  set (oct_dreamrect_f_SOURCE_FILES
    oct/oct_dreamrect_f.cc
    src/att.cc
    src/affinity.cc
    src/dreamrect_f.c
    src/fft_octave.cc
    src/dream_error.c
    )

  add_library (oct_dreamrect_f MODULE
    ${oct_dreamrect_f_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamrect_f
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dreamrect_f PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamrect_f")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dreamrect_f_SOURCE_FILES
    mex/mex_dreamrect_f.cc
    src/att.cc
    src/affinity.cc
    src/dreamrect_f.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dreamrect_f MODULE
    ${mex_dreamrect_f_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamrect_f
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamrect_f PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamrect_f")

endif (Matlab_FOUND)


#
# dreamcirc
#

if (OCTAVE_FOUND)

  set (oct_dreamcirc_SOURCE_FILES
    oct/oct_dreamcirc.cc
    src/att.cc
    src/affinity.cc
    src/dreamcirc.c
    src/fft_octave.cc
    src/dream_error.c
    )

  if (OpenCL_FOUND)
    set (oct_dreamcirc_SOURCE_FILES ${oct_dreamcirc_SOURCE_FILES} "opencl/cl_dreamcirc.cc")
  endif (OpenCL_FOUND)

  add_library (oct_dreamcirc MODULE
    ${oct_dreamcirc_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamcirc
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${OpenCL_LIBRARIES}
    )

  set_target_properties (oct_dreamcirc PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamcirc")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dreamcirc_SOURCE_FILES
    mex/mex_dreamcirc.cc
    src/att.cc
    src/affinity.cc
    src/dreamcirc.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dreamcirc MODULE
    ${mex_dreamcirc_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamcirc
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamcirc PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamcirc")

endif (Matlab_FOUND)

#
# dreamcirc_f
#

if (OCTAVE_FOUND)

  set (oct_dreamcirc_f_SOURCE_FILES
    oct/oct_dreamcirc_f.cc
    src/att.cc
    src/affinity.cc
    src/dreamcirc_f.c
    src/fft_octave.cc
    src/dream_error.c
    )

  add_library (oct_dreamcirc_f MODULE
    ${oct_dreamcirc_f_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamcirc_f
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dreamcirc_f PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamcirc_f")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dreamcirc_f_SOURCE_FILES
    mex/mex_dreamcirc_f.cc
    src/att.cc
    src/affinity.cc
    src/dreamcirc_f.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dreamcirc_f MODULE
    ${mex_dreamcirc_f_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamcirc_f
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamcirc_f PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamcirc_f")

endif (Matlab_FOUND)

#
# dreamsphere_f
#

if (OCTAVE_FOUND)

  set (oct_dreamsphere_f_SOURCE_FILES
    oct/oct_dreamsphere_f.cc
    src/att.cc
    src/affinity.cc
    src/dreamsphere_f.c
    src/fft_octave.cc
    src/dream_error.c
    )

  add_library (oct_dreamsphere_f MODULE
    ${oct_dreamsphere_f_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamsphere_f
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dreamsphere_f PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamsphere_f")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dreamsphere_f_SOURCE_FILES
    mex/mex_dreamsphere_f.cc
    src/att.cc
    src/affinity.cc
    src/dreamsphere_f.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dreamsphere_f MODULE
    ${mex_dreamsphere_f_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamsphere_f
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamsphere_f PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamsphere_f")

endif (Matlab_FOUND)


#
# dreamsphere_d
#

if (OCTAVE_FOUND)

  set (oct_dreamsphere_d_SOURCE_FILES
    oct/oct_dreamsphere_d.cc
    src/att.cc
    src/affinity.cc
    src/dreamsphere_d.c
    src/fft_octave.cc
    src/dream_error.c
    )

  add_library (oct_dreamsphere_d MODULE
    ${oct_dreamsphere_d_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamsphere_d
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dreamsphere_d PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamsphere_d")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dreamsphere_d_SOURCE_FILES
    mex/mex_dreamsphere_d.cc
    src/att.cc
    src/affinity.cc
    src/dreamsphere_d.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dreamsphere_d MODULE
    ${mex_dreamsphere_d_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamsphere_d
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamsphere_d PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamsphere_d")

endif (Matlab_FOUND)

#
# dreamcylind_f
#

if (OCTAVE_FOUND)

  set (oct_dreamcylind_f_SOURCE_FILES
    oct/oct_dreamcylind_f.cc
    src/att.cc
    src/affinity.cc
    src/dreamcylind_f.c
    src/fft_octave.cc
    src/dream_error.c
    )

  add_library (oct_dreamcylind_f MODULE
    ${oct_dreamcylind_f_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamcylind_f
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dreamcylind_f PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamcylind_f")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dreamcylind_f_SOURCE_FILES
    mex/mex_dreamcylind_f.cc
    src/att.cc
    src/affinity.cc
    src/dreamcylind_f.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dreamcylind_f MODULE
    ${mex_dreamcylind_f_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamcylind_f
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamcylind_f PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamcylind_f")

endif (Matlab_FOUND)

#
# dreamcylind_d
#

if (OCTAVE_FOUND)

  set (oct_dreamcylind_d_SOURCE_FILES
    oct/oct_dreamcylind_d.cc
    src/att.cc
    src/affinity.cc
    src/dreamcylind_d.c
    src/fft_octave.cc
    src/dream_error.c
    )

  add_library (oct_dreamcylind_d MODULE
    ${oct_dreamcylind_d_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamcylind_d
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dreamcylind_d PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamcylind_d")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dreamcylind_d_SOURCE_FILES
    mex/mex_dreamcylind_d.cc
    src/att.cc
    src/affinity.cc
    src/dreamcylind_d.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dreamcylind_d MODULE
    ${mex_dreamcylind_d_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamcylind_d
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamcylind_d PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamcylind_d")


endif (Matlab_FOUND)

#
# dream_arr_rect
#

if (OCTAVE_FOUND)

  set (oct_dream_arr_rect_SOURCE_FILES
    oct/oct_dream_arr_rect.cc
    src/att.cc
    src/affinity.cc
    src/arr_functions.c
    src/dream_arr_rect.c
    src/fft_octave.cc
    src/dream_error.c
    )

  add_library (oct_dream_arr_rect MODULE
    ${oct_dream_arr_rect_SOURCE_FILES}
    )

  target_link_libraries (oct_dream_arr_rect
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dream_arr_rect PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dream_arr_rect")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dream_arr_circ_SOURCE_FILES
    mex/mex_dream_arr_circ.cc
    src/att.cc
    src/affinity.cc
    src/arr_functions.c
    src/dream_arr_circ.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dream_arr_circ MODULE
    ${mex_dream_arr_circ_SOURCE_FILES}
    )

  target_link_libraries (mex_dream_arr_circ
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dream_arr_circ PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dream_arr_circ")

endif (Matlab_FOUND)

#
# dream_arr_circ
#

if (OCTAVE_FOUND)

  set (oct_dream_arr_circ_SOURCE_FILES
    oct/oct_dream_arr_circ.cc
    src/att.cc
    src/affinity.cc
    src/arr_functions.c
    src/dream_arr_circ.c
    src/fft_octave.cc
    src/dream_error.c
    )

  add_library (oct_dream_arr_circ MODULE
    ${oct_dream_arr_circ_SOURCE_FILES}
    )

  target_link_libraries (oct_dream_arr_circ
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dream_arr_circ PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dream_arr_circ")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dream_arr_rect_SOURCE_FILES
    mex/mex_dream_arr_rect.cc
    src/att.cc
    src/affinity.cc
    src/arr_functions.c
    src/dream_arr_rect.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dream_arr_rect MODULE
    ${mex_dream_arr_rect_SOURCE_FILES}
    )

  target_link_libraries (mex_dream_arr_rect
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dream_arr_rect PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dream_arr_rect")

endif (Matlab_FOUND)

#
# dream_arr_cylind_f
#

if (OCTAVE_FOUND)

  set (oct_dream_arr_cylind_f_SOURCE_FILES
    oct/oct_dream_arr_cylind_f.cc
    src/att.cc
    src/affinity.cc
    src/arr_functions.c
    src/dream_arr_cylind_f.c
    src/fft_octave.cc
    src/dream_error.c
    )

  add_library (oct_dream_arr_cylind_f MODULE
    ${oct_dream_arr_cylind_f_SOURCE_FILES}
    )

  target_link_libraries (oct_dream_arr_cylind_f
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dream_arr_cylind_f PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dream_arr_cylind_f")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dream_arr_cylind_f_SOURCE_FILES
    mex/mex_dream_arr_cylind_f.cc
    src/att.cc
    src/affinity.cc
    src/arr_functions.c
    src/dream_arr_cylind_f.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dream_arr_cylind_f MODULE
    ${mex_dream_arr_cylind_f_SOURCE_FILES}
    )

  target_link_libraries (mex_dream_arr_cylind_f
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dream_arr_cylind_f PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dream_arr_cylind_f")

endif (Matlab_FOUND)

#
# dream_arr_cylind_d
#

if (OCTAVE_FOUND)

  set (oct_dream_arr_cylind_d_SOURCE_FILES
    oct/oct_dream_arr_cylind_d.cc
    src/att.cc
    src/affinity.cc
    src/arr_functions.c
    src/dream_arr_cylind_d.c
    src/fft_octave.cc
    src/dream_error.c
    )

  add_library (oct_dream_arr_cylind_d MODULE
    ${oct_dream_arr_cylind_d_SOURCE_FILES}
    )

  target_link_libraries (oct_dream_arr_cylind_d
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dream_arr_cylind_d PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dream_arr_cylind_d")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dream_arr_cylind_d_SOURCE_FILES
    mex/mex_dream_arr_cylind_d.cc
    src/att.cc
    src/affinity.cc
    src/arr_functions.c
    src/dream_arr_cylind_d.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dream_arr_cylind_d MODULE
    ${mex_dream_arr_cylind_d_SOURCE_FILES}
    )

  target_link_libraries (mex_dream_arr_cylind_d
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dream_arr_cylind_d PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dream_arr_cylind_d")

endif (Matlab_FOUND)

#
# dream_arr_annu
#

if (OCTAVE_FOUND)

  set (oct_dream_arr_annu_SOURCE_FILES
    oct/oct_dream_arr_annu.cc
    src/att.cc
    src/affinity.cc
    src/arr_functions.c
    src/dream_arr_annu.c
    src/fft_octave.cc
    src/dream_error.c
    )

  add_library (oct_dream_arr_annu MODULE
    ${oct_dream_arr_annu_SOURCE_FILES}
    )

  target_link_libraries (oct_dream_arr_annu
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dream_arr_annu PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dream_arr_annu")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dream_arr_annu_SOURCE_FILES
    mex/mex_dream_arr_annu.cc
    src/att.cc
    src/affinity.cc
    src/arr_functions.c
    src/dream_arr_annu.c
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dream_arr_annu MODULE
    ${mex_dream_arr_annu_SOURCE_FILES}
    )

  target_link_libraries (mex_dream_arr_annu
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dream_arr_annu PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dream_arr_annu")

endif (Matlab_FOUND)

#
# Misc routines.
#


#
# dream_apodwin
#

if (OCTAVE_FOUND)

  set (oct_dream_apodwin_SOURCE_FILES
    oct/oct_dream_apodwin.cc
    src/arr_functions.c
    src/dream_error.c
    )

  add_library (oct_dream_apodwin MODULE
    ${oct_dream_apodwin_SOURCE_FILES}
    )

  target_link_libraries (oct_dream_apodwin
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dream_apodwin PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dream_apodwin")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dream_apodwin_SOURCE_FILES
    mex/mex_dream_apodwin.cc
    src/arr_functions.c
    src/dream_error.c
    )

  add_library (mex_dream_apodwin MODULE
    ${mex_dream_apodwin_SOURCE_FILES}
    )

  target_link_libraries (mex_dream_apodwin
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dream_apodwin PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dream_apodwin")

endif (Matlab_FOUND)

#
# dream_att
#


if (OCTAVE_FOUND)

  set (oct_dream_att_SOURCE_FILES
    oct/oct_dream_att.cc
    src/att.cc
    src/affinity.cc
    src/fft_octave.cc
    src/dream_error.c
    )

  add_library (oct_dream_att MODULE
    ${oct_dream_att_SOURCE_FILES}
    )

  target_link_libraries (oct_dream_att
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dream_att PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dream_att")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_dream_att_SOURCE_FILES
    mex/mex_dream_att.c
    src/att.cc
    src/affinity.cc
    src/fft.c
    src/dream_error.c
    )

  add_library (mex_dream_att MODULE
    ${mex_dream_att_SOURCE_FILES}
    )

  target_link_libraries (mex_dream_att
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dream_att PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dream_att")

endif (Matlab_FOUND)


#
# Analytic SIRs.
#

# Circ

if (OCTAVE_FOUND)

  set (oct_circ_sir_SOURCE_FILES
    extra_src/oct_circ_sir.cc
    src/affinity.cc
    extra_src/circ_sir.c
    src/dream_error.c
    )

  add_library (oct_circ_sir MODULE
    ${oct_circ_sir_SOURCE_FILES}
    )

  target_link_libraries (oct_circ_sir
    ${OCTAVE_LIBRARIES}
    )

  set_target_properties (oct_circ_sir PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "circ_sir")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_circ_sir_SOURCE_FILES
    extra_src/mex_circ_sir.cc
    src/affinity.cc
    extra_src/circ_sir.c
    src/dream_error.c
    )

  add_library (mex_circ_sir MODULE
    ${mex_circ_sir_SOURCE_FILES}
    )

  target_link_libraries (mex_circ_sir
    ${MATLAB_LIBRARIES}
    )

  set_target_properties (mex_circ_sir PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "circ_sir")

endif (Matlab_FOUND)

# Sampled analytic circ

if (OCTAVE_FOUND)

  set (oct_scirc_sir_SOURCE_FILES
    extra_src/oct_scirc_sir.cc
    extra_src/scirc_sir.c
    src/dream_error.c
    )

  add_library (oct_scirc_sir MODULE
    ${oct_scirc_sir_SOURCE_FILES}
    )

  target_link_libraries (oct_scirc_sir
    ${OCTAVE_LIBRARIES}
    )

  set_target_properties (oct_scirc_sir PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "scirc_sir")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_scirc_sir_SOURCE_FILES
    extra_src/mex_scirc_sir.cc
    extra_src/scirc_sir.c
    src/dream_error.c
    )

  add_library (mex_scirc_sir MODULE
    ${mex_scirc_sir_SOURCE_FILES}
    )

  target_link_libraries (mex_scirc_sir
    ${MATLAB_LIBRARIES}
    )

  set_target_properties (mex_scirc_sir PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "scirc_sir")

endif (Matlab_FOUND)

# Rect

if (OCTAVE_FOUND)

  set (oct_rect_sir_SOURCE_FILES
    extra_src/oct_rect_sir.cc
    extra_src/rect_sir.c
    src/affinity.cc
    src/dream_error.c
    )

  if (OpenCL_FOUND)
    set (oct_rect_sir_SOURCE_FILES ${oct_rect_sir_SOURCE_FILES} "opencl/cl_rect_sir.cc")
  endif (OpenCL_FOUND)

  add_library (oct_rect_sir MODULE
    ${oct_rect_sir_SOURCE_FILES}
    )

  target_link_libraries (oct_rect_sir
    ${OCTAVE_LIBRARIES}
    ${OpenCL_LIBRARIES}
    )

  set_target_properties (oct_rect_sir PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "rect_sir")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_rect_sir_SOURCE_FILES
    extra_src/mex_rect_sir.cc
    extra_src/rect_sir.c
    src/affinity.cc
    src/dream_error.c
    )

  add_library (mex_rect_sir MODULE
    ${mex_rect_sir_SOURCE_FILES}
    )

  target_link_libraries (mex_rect_sir
    ${MATLAB_LIBRARIES}
    )

  set_target_properties (mex_rect_sir PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "rect_sir")

endif (Matlab_FOUND)

#
# Delay-and-sum (DAS)
#

# das

if (OCTAVE_FOUND)

  set (oct_das_SOURCE_FILES
    extra_src/oct_das.cc
    extra_src/das.c
    src/dream_error.c
    )

  add_library (oct_das MODULE
    ${oct_das_SOURCE_FILES}
    )

  target_link_libraries (oct_das
    ${OCTAVE_LIBRARIES}
    )

  set_target_properties (oct_das PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "das")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_das_SOURCE_FILES
    extra_src/mex_das.c
    extra_src/das.c
    src/dream_error.c
    )

  add_library (mex_das MODULE
    ${mex_das_SOURCE_FILES}
    )

  target_link_libraries (mex_das
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_das PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "das")

endif (Matlab_FOUND)

# das_err

if (OCTAVE_FOUND)

  set (oct_das_arr_SOURCE_FILES
    extra_src/oct_das_arr.cc
    extra_src/das_arr.c
    src/arr_functions.c
    src/dream_error.c
    )

  add_library (oct_das_arr MODULE
    ${oct_das_arr_SOURCE_FILES}
    )

  target_link_libraries (oct_das_arr
    ${OCTAVE_LIBRARIES}
    )

  set_target_properties (oct_das_arr PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "das_arr")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_das_arr_SOURCE_FILES
    extra_src/mex_das_arr.c
    extra_src/das_arr.c
    src/arr_functions.c
    src/dream_error.c
    )

  add_library (mex_das_arr MODULE
    ${mex_das_arr_SOURCE_FILES}
    )

  target_link_libraries (mex_das_arr
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_das_arr PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "das_arr")

endif (Matlab_FOUND)

#
# SAFT
#

if (OCTAVE_FOUND)

  set (oct_saft_SOURCE_FILES
    extra_src/oct_saft.cc
    src/affinity.cc
    src/dream_error.c
    )

  add_library (oct_saft MODULE
    ${oct_saft_SOURCE_FILES}
    )

  target_link_libraries (oct_saft
    ${OCTAVE_LIBRARIES}
    )

  set_target_properties (oct_saft PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "saft")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_saft_SOURCE_FILES
    extra_src/mex_saft.cc
    src/dream_error.c
    )

  add_library (mex_saft MODULE
    ${mex_saft_SOURCE_FILES}
    )

  target_link_libraries (mex_saft
    ${MATLAB_LIBRARIES}
    )

  set_target_properties (mex_saft PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "saft")

endif (Matlab_FOUND)

#
# Convolution algorithms
#

# conv_p

if (OCTAVE_FOUND)

  set (oct_conv_p_SOURCE_FILES
    extra_src/oct_conv_p.cc
    src/affinity.cc
    extra_src/conv.c
    src/dream_error.c
    )

  add_library (oct_conv_p MODULE
    ${oct_conv_p_SOURCE_FILES}
    )

  target_link_libraries (oct_conv_p
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_conv_p PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "conv_p")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_conv_p_SOURCE_FILES
    extra_src/mex_conv_p.cc
    src/affinity.cc
    extra_src/conv.c
    src/dream_error.c
    )

  add_library (mex_conv_p MODULE
    ${mex_conv_p_SOURCE_FILES}
    )

  target_link_libraries (mex_conv_p
    ${MATLAB_LIBRARIES}
    )

  set_target_properties (mex_conv_p PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "conv_p")

endif (Matlab_FOUND)

# fftconv_p

if (OCTAVE_FOUND)

  set (oct_fftconv_p_SOURCE_FILES
    extra_src/oct_fftconv_p.cc
    src/affinity.cc
    src/dream_error.c
    )

  add_library (oct_fftconv_p MODULE
    ${oct_fftconv_p_SOURCE_FILES}
    )

  target_link_libraries (oct_fftconv_p
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_fftconv_p PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "fftconv_p")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_fftconv_p_SOURCE_FILES
    extra_src/mex_fftconv_p.cc
    src/affinity.cc
    src/dream_error.c
    )

  add_library (mex_fftconv_p MODULE
    ${mex_fftconv_p_SOURCE_FILES}
    )

  target_link_libraries (mex_fftconv_p
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_fftconv_p PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "fftconv_p")

endif (Matlab_FOUND)


# sum_fftconv

if (OCTAVE_FOUND)

  set (oct_sum_fftconv_SOURCE_FILES
    extra_src/oct_sum_fftconv.cc
    src/affinity.cc
    src/dream_error.c
    )

  add_library (oct_sum_fftconv MODULE
    ${oct_sum_fftconv_SOURCE_FILES}
    )

  target_link_libraries (oct_sum_fftconv
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_sum_fftconv PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "sum_fftconv")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_sum_fftconv_SOURCE_FILES
    extra_src/mex_sum_fftconv.cc
    src/affinity.cc
    src/dream_error.c
    )

  add_library (mex_sum_fftconv MODULE
    ${mex_sum_fftconv_SOURCE_FILES}
    )

  target_link_libraries (mex_sum_fftconv
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_sum_fftconv PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "sum_fftconv")

endif (Matlab_FOUND)


# fftconv_ola

if (OCTAVE_FOUND)

  set (oct_fftconv_ola_SOURCE_FILES
    extra_src/oct_fftconv_ola.cc
    src/affinity.cc
    src/dream_error.c
    )

  add_library (oct_fftconv_ola MODULE
    ${oct_fftconv_ola_SOURCE_FILES}
    )

  target_link_libraries (oct_fftconv_ola
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${SFFTW_LIBRARIES}
    )

  set_target_properties (oct_fftconv_ola PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "fftconv_ola")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  # Not implemented yet.

endif (Matlab_FOUND)

#
# Copy
#

if (OCTAVE_FOUND)

  set (oct_copy_p_SOURCE_FILES
    extra_src/oct_copy_p.cc
    src/affinity.cc
    src/dream_error.c
    )

  add_library (oct_copy_p MODULE
    ${oct_copy_p_SOURCE_FILES}
    )

  target_link_libraries (oct_copy_p
    ${OCTAVE_LIBRARIES}
    )

  set_target_properties (oct_copy_p PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS "${OCT_FLAGS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "copy_p")

endif (OCTAVE_FOUND)

if (Matlab_FOUND)

  set (mex_copy_p_SOURCE_FILES
    extra_src/mex_copy_p.cc
    src/dream_error.c
    )

  add_library (mex_copy_p MODULE
    ${mex_copy_p_SOURCE_FILES}
    )

  target_link_libraries (mex_copy_p
    ${MATLAB_LIBRARIES}
    )

  set_target_properties (mex_copy_p PROPERTIES
    C_STANDARD 99
    CXX_STANDARD 11
    COMPILE_FLAGS ""
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "copy_p")

endif (Matlab_FOUND)

#
# User manual (N.B., we need CMake 3.x on Windows due to a MikTeX related bug).
#

if (BUILD_USERMAN)

  # Generate version files for the user manual.
  execute_process(COMMAND svnversion
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE SVNVER
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  file(WRITE  "${CMAKE_CURRENT_BINARY_DIR}/version.tex" "\\newcommand{\\version}{${DREAM_VERSION_MAJOR}.${DREAM_VERSION_MINOR}.${DREAM_VERSION_PATCH}\\xspace}\r")
  file(WRITE  "${CMAKE_CURRENT_BINARY_DIR}/revision.tex" "\\rhead{\\date Revision ${SVNVER} }\r")

  include (userman/UseLATEX.cmake)

  add_latex_document(userman/dream_userman.tex
    BIBFILES userman/dream_userman.bib
    IMAGE_DIRS userman/eps
    DEFAULT_SAFEPDF # Using dvips/ps2pdf
    #DEFAULT_PS
    #DEFAULT_PDF # Messes up page breakes and fig sizes.
    )

  if(WIN32 AND NOT UNIX)
    install (FILES ${CMAKE_CURRENT_BINARY_DIR}/dream_userman.pdf DESTINATION documents)
  endif(WIN32 AND NOT UNIX)

endif (BUILD_USERMAN)
