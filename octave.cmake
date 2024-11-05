#
# Copyright (C) 2022 Fredrik Lingvall
#

# From http://www.coolprop.org/coolprop/wrappers/Octave/index.html
find_package (Octave)

if (OCTAVE_FOUND)

  add_definitions(-DOCTAVE) # Needed to catch CTRL-C by using the OCTAVE_QUIT macro inside loops.

  # Octave oct flags.
  set (DREAM_OCT_FLAGS "-DDREAM_OCTAVE") # So that octave_idx_type is used for matrix/vector indexing.

  if (UNIX)
    set (OCT_LD_FLAGS "${CMAKE_LD_FLAGS} ${OCTAVE_LINK_FLAGS} -pthread") # We need phread to use pthread_setaffinity_np.
  else (UNIX)
    set (OCT_LD_FLAGS "${CMAKE_LD_FLAGS} ${OCTAVE_LINK_FLAGS}")
  endif (UNIX)

  #
  # Include paths
  #

  # Strip lead and trailing whitepasce
  string(STRIP "${OCTAVE_INCLUDE_DIRS}" OCTAVE_INCLUDE_DIRS)

  set (DREAM_OCT_INCLUDE_DIRS
    "${PROJECT_BINARY_DIR};"
    "${PROJECT_BINARY_DIR}/opencl;"
    "${PROJECT_SOURCE_DIR}/include;"
    "${FFTW_INCLUDE_DIR};"
    "${OCTAVE_INCLUDE_DIRS}"
    )
  message (STATUS "DREAM_OCT_INCLUDE_DIRS: ${DREAM_OCT_INCLUDE_DIRS}")

  #
  # dreamline
  #

  set (oct_dreamline_SOURCE_FILES
    oct/oct_dreamline.cc
    src/attenuation.cc
    src/affinity.cc
    src/dreamline.cc
    src/dream_error.cc
    )

  add_library (oct_dreamline MODULE
    ${oct_dreamline_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamline
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dreamline PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamline")

  #
  # dreamrect
  #

  set (oct_dreamrect_SOURCE_FILES
    oct/oct_dreamrect.cc
    src/dreamrect.cc
    src/attenuation.cc
    src/affinity.cc
    src/dream_error.cc
    )

  if (OpenCL_FOUND)
    set (oct_dreamrect_SOURCE_FILES ${oct_dreamrect_SOURCE_FILES} "opencl/cl_dreamrect.cc")
    file (COPY opencl/dreamrect.cl DESTINATION kernels)
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
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamrect")

  #
  # dreamrect_f
  #

  set (oct_dreamrect_f_SOURCE_FILES
    oct/oct_dreamrect_f.cc
    src/attenuation.cc
    src/affinity.cc
    src/dreamrect_f.cc
    src/dream_error.cc
    )

  add_library (oct_dreamrect_f MODULE
    ${oct_dreamrect_f_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamrect_f
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dreamrect_f PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamrect_f")

  #
  # dreamcirc
  #

  set (oct_dreamcirc_SOURCE_FILES
    oct/oct_dreamcirc.cc
    src/attenuation.cc
    src/affinity.cc
    src/dreamcirc.cc
    src/dream_error.cc
    )

  if (OpenCL_FOUND)
    set (oct_dreamcirc_SOURCE_FILES ${oct_dreamcirc_SOURCE_FILES} "opencl/cl_dreamcirc.cc")
    file (COPY opencl/dreamcirc.cl DESTINATION kernels)
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
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamcirc")

  #
  # dreamcirc_f
  #

  set (oct_dreamcirc_f_SOURCE_FILES
    oct/oct_dreamcirc_f.cc
    src/attenuation.cc
    src/affinity.cc
    src/dreamcirc_f.cc
    src/dream_error.cc
    )

  add_library (oct_dreamcirc_f MODULE
    ${oct_dreamcirc_f_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamcirc_f
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dreamcirc_f PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamcirc_f")

  #
  # dreamsphere
  #

  set (oct_dreamsphere_SOURCE_FILES
    oct/oct_dreamsphere.cc
    src/attenuation.cc
    src/affinity.cc
    src/dreamsphere.cc
    src/dream_error.cc
    )

  add_library (oct_dreamsphere MODULE
    ${oct_dreamsphere_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamsphere
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dreamsphere PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamsphere")


  #
  # dreamcylind
  #

  set (oct_dreamcylind_SOURCE_FILES
    oct/oct_dreamcylind.cc
    src/attenuation.cc
    src/affinity.cc
    src/dreamcylind.cc
    src/dream_error.cc
    )

  add_library (oct_dreamcylind MODULE
    ${oct_dreamcylind_SOURCE_FILES}
    )

  target_link_libraries (oct_dreamcylind
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dreamcylind PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dreamcylind")

  #
  # dream_arr_rect
  #

  set (oct_dream_arr_rect_SOURCE_FILES
    oct/oct_dream_arr_rect.cc
    src/attenuation.cc
    src/affinity.cc
    src/arr_functions.cc
    src/dream_arr_rect.cc
    src/dreamrect.cc
    src/dream_error.cc
    )

  add_library (oct_dream_arr_rect MODULE
    ${oct_dream_arr_rect_SOURCE_FILES}
    )

  target_link_libraries (oct_dream_arr_rect
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dream_arr_rect PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dream_arr_rect")
  #
  # dream_arr_circ
  #

  set (oct_dream_arr_circ_SOURCE_FILES
    oct/oct_dream_arr_circ.cc
    src/attenuation.cc
    src/affinity.cc
    src/arr_functions.cc
    src/dream_arr_circ.cc
    src/dreamcirc.cc
    src/dream_error.cc
    )

  add_library (oct_dream_arr_circ MODULE
    ${oct_dream_arr_circ_SOURCE_FILES}
    )

  target_link_libraries (oct_dream_arr_circ
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dream_arr_circ PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dream_arr_circ")

  #
  # dream_arr_cylind
  #

  set (oct_dream_arr_cylind_SOURCE_FILES
    oct/oct_dream_arr_cylind.cc
    src/attenuation.cc
    src/affinity.cc
    src/arr_functions.cc
    src/dream_arr_cylind.cc
    src/dreamcylind.cc
    src/dream_error.cc
    )

  add_library (oct_dream_arr_cylind MODULE
    ${oct_dream_arr_cylind_SOURCE_FILES}
    )

  target_link_libraries (oct_dream_arr_cylind
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dream_arr_cylind PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dream_arr_cylind")

  #
  # dream_arr_annu
  #

  set (oct_dream_arr_annu_SOURCE_FILES
    oct/oct_dream_arr_annu.cc
    src/attenuation.cc
    src/affinity.cc
    src/arr_functions.cc
    src/dream_arr_annu.cc
    src/dreamcirc.cc
    src/dream_error.cc
    )

  add_library (oct_dream_arr_annu MODULE
    ${oct_dream_arr_annu_SOURCE_FILES}
    )

  target_link_libraries (oct_dream_arr_annu
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dream_arr_annu PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dream_arr_annu")

  ############
  #
  # Misc routines.
  #
  ############

  #
  # dream_apodwin
  #

  set (oct_dream_apodwin_SOURCE_FILES
    oct/oct_dream_apodwin.cc
    src/arr_functions.cc
    src/dream_error.cc
    )

  add_library (oct_dream_apodwin MODULE
    ${oct_dream_apodwin_SOURCE_FILES}
    )

  target_link_libraries (oct_dream_apodwin
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dream_apodwin PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dream_apodwin")

  #
  # dream_att
  #

  set (oct_dream_att_SOURCE_FILES
    oct/oct_dream_att.cc
    src/attenuation.cc
    src/affinity.cc
    src/dream_error.cc
    )

  add_library (oct_dream_att MODULE
    ${oct_dream_att_SOURCE_FILES}
    )

  target_link_libraries (oct_dream_att
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_dream_att PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "dream_att")

  #
  # Analytic SIRs.
  #

  # Circ

  set (oct_circ_sir_SOURCE_FILES
    extra_src/oct_circ_sir.cc
    src/affinity.cc
    extra_src/circ_sir.cc
    src/dream_error.cc
    )

  add_library (oct_circ_sir MODULE
    ${oct_circ_sir_SOURCE_FILES}
    )

  target_link_libraries (oct_circ_sir
    ${OCTAVE_LIBRARIES}
    )

  set_target_properties (oct_circ_sir PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "circ_sir")


  # Sampled analytic circ

  set (oct_scirc_sir_SOURCE_FILES
    extra_src/oct_scirc_sir.cc
    src/affinity.cc
    extra_src/scirc_sir.cc
    src/dream_error.cc
    )

  add_library (oct_scirc_sir MODULE
    ${oct_scirc_sir_SOURCE_FILES}
    )

  target_link_libraries (oct_scirc_sir
    ${OCTAVE_LIBRARIES}
    )

  set_target_properties (oct_scirc_sir PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "scirc_sir")

  # Rect

  set (oct_rect_sir_SOURCE_FILES
    extra_src/oct_rect_sir.cc
    extra_src/rect_sir.cc
    src/affinity.cc
    src/dream_error.cc
    )

  if (OpenCL_FOUND)
    set (oct_rect_sir_SOURCE_FILES ${oct_rect_sir_SOURCE_FILES} "opencl/cl_rect_sir.cc")
    file (COPY opencl/rect_sir.cl DESTINATION kernels)
  endif (OpenCL_FOUND)

  add_library (oct_rect_sir MODULE
    ${oct_rect_sir_SOURCE_FILES}
    )

  target_link_libraries (oct_rect_sir
    ${OCTAVE_LIBRARIES}
    ${OpenCL_LIBRARIES}
    )

  set_target_properties (oct_rect_sir PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "rect_sir")

  #
  # Delay-and-sum (DAS)
  #

  # das

  set (oct_das_SOURCE_FILES
    extra_src/oct_das.cc
    extra_src/das.cc
    src/affinity.cc
    src/dream_error.cc
  )

  if (OpenCL_FOUND)
    set (oct_das_SOURCE_FILES ${oct_das_SOURCE_FILES} "opencl/cl_das.cc")
  endif (OpenCL_FOUND)

  add_library (oct_das MODULE
    ${oct_das_SOURCE_FILES}
  )

  target_link_libraries (oct_das
    ${OCTAVE_LIBRARIES}
    ${OpenCL_LIBRARIES}
    )

  set_target_properties (oct_das PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "das")

  #
  # Delay-and-sum (DAS) using uniform grids (OpenCL/GPU only)
  #

  if (OpenCL_FOUND)

    set (oct_das_uni_SOURCE_FILES
      extra_src/oct_das_uni.cc
      opencl/cl_das_uni.cc
    )

    file (COPY opencl/das_uni_float.cl DESTINATION kernels)
    file (COPY opencl/das_uni_double.cl DESTINATION kernels)

    add_library (oct_das_uni MODULE
      ${oct_das_uni_SOURCE_FILES}
    )

    target_link_libraries (oct_das_uni
      ${OCTAVE_LIBRARIES}
      ${OpenCL_LIBRARIES}
    )

    set_target_properties (oct_das_uni PROPERTIES
      CXX_STANDARD 14
      COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
      INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
      LINK_FLAGS ${OCT_LD_FLAGS}
      SUFFIX ".oct" PREFIX "" OUTPUT_NAME "das_uni")

  endif (OpenCL_FOUND)

  #
  # Convolution algorithms
  #

  # conv_p

  set (oct_conv_p_SOURCE_FILES
    extra_src/oct_conv_p.cc
    src/affinity.cc
    extra_src/conv.cc
    src/dream_error.cc
    )

  add_library (oct_conv_p MODULE
    ${oct_conv_p_SOURCE_FILES}
    )

  target_link_libraries (oct_conv_p
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_conv_p PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "conv_p")

  # fftconv_p

  set (oct_fftconv_p_SOURCE_FILES
    extra_src/oct_fftconv_p.cc
    extra_src/fftconv.cc
    src/affinity.cc
    src/dream_error.cc
    )

  add_library (oct_fftconv_p MODULE
    ${oct_fftconv_p_SOURCE_FILES}
    )

  target_link_libraries (oct_fftconv_p
    ${OCTAVE_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (oct_fftconv_p PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "fftconv_p")

  # sum_fftconv

  # This one fails on macOS
  if(NOT MACOSX)

    set (oct_sum_fftconv_SOURCE_FILES
      extra_src/oct_sum_fftconv.cc
      extra_src/fftconv.cc
      src/affinity.cc
      src/dream_error.cc
      )

    add_library (oct_sum_fftconv MODULE
      ${oct_sum_fftconv_SOURCE_FILES}
      )

    target_link_libraries (oct_sum_fftconv
      ${OCTAVE_LIBRARIES}
      ${FFTW_LIBRARIES}
      )

    set_target_properties (oct_sum_fftconv PROPERTIES
      CXX_STANDARD 14
      COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
      INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
      LINK_FLAGS ${OCT_LD_FLAGS}
      SUFFIX ".oct" PREFIX "" OUTPUT_NAME "sum_fftconv")

  endif(NOT MACOSX)

  # fftconv_ola

  set (oct_fftconv_ola_SOURCE_FILES
    extra_src/oct_fftconv_ola.cc
    extra_src/fftconv.cc
    src/affinity.cc
    src/dream_error.cc
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
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "fftconv_ola")

  #
  # Add
  #

  set (oct_add_p_SOURCE_FILES
    extra_src/oct_add_p.cc
    src/affinity.cc
    src/dream_error.cc
    )

  add_library (oct_add_p MODULE
    ${oct_add_p_SOURCE_FILES}
    )

  target_link_libraries (oct_add_p
    ${OCTAVE_LIBRARIES}
    )

  set_target_properties (oct_add_p PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "add_p")

  #
  # Copy
  #

  set (oct_copy_p_SOURCE_FILES
    extra_src/oct_copy_p.cc
    src/affinity.cc
    src/dream_error.cc
    )

  add_library (oct_copy_p MODULE
    ${oct_copy_p_SOURCE_FILES}
    )

  target_link_libraries (oct_copy_p
    ${OCTAVE_LIBRARIES}
    )

  set_target_properties (oct_copy_p PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_OCT_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_OCT_INCLUDE_DIRS}"
    LINK_FLAGS ${OCT_LD_FLAGS}
    SUFFIX ".oct" PREFIX "" OUTPUT_NAME "copy_p")

endif (OCTAVE_FOUND)
