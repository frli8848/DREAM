#
# Copyright (C) 2022,2023 Fredrik Lingvall
#

find_package (Matlab)
message (STATUS "Matlab_MEX_EXTENSION = ${Matlab_MEX_EXTENSION}")

if (Matlab_FOUND)

  # mex libs
  #set (MEX_LD_FLAGS "${CMAKE_LD_FLAGS} ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY}")
  message (STATUS "Matlab_ROOT_DIR ${Matlab_ROOT_DIR}")
  message (STATUS "Matlab_MAIN_PROGRAM ${Matlab_MAIN_PROGRAM}")
  message (STATUS "Matlab_INCLUDE_DIRS ${Matlab_INCLUDE_DIRS}")
  message (STATUS "Matlab_LIBRARIES  ${Matlab_LIBRARIES}")

  #find_package (MWFFTW)
  #if (MWFFTW_FOUND)
  #  message (STATUS "MWFFTW_LIB_PATH: ${MWFFTW_LIB_PATH}")
  #  message (STATUS "MWFFTW_LIBRARIES: ${MWFFTW_LIBRARIES}")
  #  #set (MEX_LD_FLAGS "${MEX_LD_FLAGS} -L${MWFFTW_LIB_PATH} -lmwfftw3")
  #  #set (MEX_LD_FLAGS "${MEX_LD_FLAGS} ${MWFFTW_LIBRARIES}")
  #  #set (FFTW_LIBRARIES "${MWFFTW_LIBRARIES}") # Override the system FFTW lib
  #endif (MWFFTW_FOUND)

  # Try static FFTW lib to avoid Matlab FFTW crashes
  #set (FFTW_LIBRARIES "${FFTW_STATIC_LIBRARIES}")

  #
  # Include paths
  #

  set (DREAM_MEX_INCLUDE_DIRS
    #${PROJECT_BINARY_DIR};"
    "${PROJECT_SOURCE_DIR}/include;"
    "${FFTW_INCLUDE_DIR};"
    "${Matlab_INCLUDE_DIRS}"
    )
  message (STATUS "DREAM_MEX_INCLUDE_DIRS: ${DREAM_MEX_INCLUDE_DIRS}")

  #
  # dreamline
  #

  set (mex_dreamline_SOURCE_FILES
    mex/mex_dreamline.cc
    src/attenuation.cc
    src/affinity.cc
    src/dreamline.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_dreamline MODULE
    SRC ${mex_dreamline_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "dreamline"
  )

  set_target_properties (mex_dreamline PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_dreamline PUBLIC DREAM_MATLAB)

  #
  # dreamrect
  #

  set (mex_dreamrect_SOURCE_FILES
    mex/mex_dreamrect.cc
    src/attenuation.cc
    src/affinity.cc
    src/dreamrect.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_dreamrect MODULE
    SRC ${mex_dreamrect_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "dreamrect"
  )

  set_target_properties (mex_dreamrect PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_dreamrect PUBLIC DREAM_MATLAB)

  #
  # dreamrect_f
  #

  set (mex_dreamrect_f_SOURCE_FILES
    mex/mex_dreamrect_f.cc
    src/attenuation.cc
    src/affinity.cc
    src/dreamrect_f.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_dreamrect_f MODULE
    SRC ${mex_dreamrect_f_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "dreamrect_f"
  )

  set_target_properties (mex_dreamrect_f PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_dreamrect_f PUBLIC DREAM_MATLAB)

  #
  # dreamcirc
  #

  set (mex_dreamcirc_SOURCE_FILES
    mex/mex_dreamcirc.cc
    src/attenuation.cc
    src/affinity.cc
    src/dreamcirc.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_dreamcirc MODULE
    SRC ${mex_dreamcirc_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "dreamcirc"
  )

  set_target_properties (mex_dreamcirc PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_dreamcirc PUBLIC DREAM_MATLAB)

  #
  # dreamcirc_f
  #

  set (mex_dreamcirc_f_SOURCE_FILES
    mex/mex_dreamcirc_f.cc
    src/attenuation.cc
    src/affinity.cc
    src/dreamcirc_f.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_dreamcirc_f MODULE
    SRC ${mex_dreamcirc_f_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "dreamcirc_f"
  )

  set_target_properties (mex_dreamcirc_f PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_dreamcirc_f PUBLIC DREAM_MATLAB)

  #
  # dreamsphere
  #

  set (mex_dreamsphere_SOURCE_FILES
    mex/mex_dreamsphere.cc
    src/attenuation.cc
    src/affinity.cc
    src/dreamsphere.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_dreamsphere MODULE
    SRC ${mex_dreamsphere_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "dreamsphere"
  )

  set_target_properties (mex_dreamsphere PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_dreamsphere PUBLIC DREAM_MATLAB)

  #
  # dreamcylind
  #

  set (mex_dreamcylind_SOURCE_FILES
    mex/mex_dreamcylind.cc
    src/attenuation.cc
    src/affinity.cc
    src/dreamcylind.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_dreamcylind MODULE
    SRC ${mex_dreamcylind_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "dreamcylind"
  )

  set_target_properties (mex_dreamcylind PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_dreamcylind PUBLIC DREAM_MATLAB)

  #
  # dream_arr_rect
  #

  set (mex_dream_arr_rect_SOURCE_FILES
    mex/mex_dream_arr_rect.cc
    src/attenuation.cc
    src/affinity.cc
    src/arr_functions.cc
    src/dream_arr_rect.cc
    src/dreamrect.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_dream_arr_rect MODULE
    SRC ${mex_dream_arr_rect_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "dream_arr_rect"
  )

  set_target_properties (mex_dream_arr_rect PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_dream_arr_rect PUBLIC DREAM_MATLAB)

  #
  # dream_arr_circ
  #

  set (mex_dream_arr_circ_SOURCE_FILES
    mex/mex_dream_arr_circ.cc
    src/attenuation.cc
    src/affinity.cc
    src/arr_functions.cc
    src/dream_arr_circ.cc
    src/dreamcirc.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_dream_arr_circ MODULE
    SRC ${mex_dream_arr_circ_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "dream_arr_circ"
  )

  set_target_properties (mex_dream_arr_circ PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_dream_arr_circ PUBLIC DREAM_MATLAB)

  #
  # dream_arr_cylind
  #

  set (mex_dream_arr_cylind_SOURCE_FILES
    mex/mex_dream_arr_cylind.cc
    src/attenuation.cc
    src/affinity.cc
    src/arr_functions.cc
    src/dream_arr_cylind.cc
    src/dreamcylind.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_dream_arr_cylind MODULE
    SRC ${mex_dream_arr_cylind_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "dream_arr_cylind"
  )

  set_target_properties (mex_dream_arr_cylind PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_dream_arr_cylind PUBLIC DREAM_MATLAB)

  #
  # dream_arr_annu
  #

  set (mex_dream_arr_annu_SOURCE_FILES
    mex/mex_dream_arr_annu.cc
    src/attenuation.cc
    src/affinity.cc
    src/arr_functions.cc
    src/dream_arr_annu.cc
    src/dreamcirc.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_dream_arr_annu MODULE
    SRC ${mex_dream_arr_annu_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "dream_arr_annu"
  )

  set_target_properties (mex_dream_arr_annu PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_dream_arr_annu PUBLIC DREAM_MATLAB)

  ############
  #
  # Misc routines.
  #
  ############

  #
  # dream_apodwin
  #

  set (mex_dream_apodwin_SOURCE_FILES
    mex/mex_dream_apodwin.cc
    src/arr_functions.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_dream_apodwin MODULE
    SRC ${mex_dream_apodwin_SOURCE_FILES}
    OUTPUT_NAME "dream_apodwin"
  )

  set_target_properties (mex_dream_apodwin PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_dream_apodwin PUBLIC DREAM_MATLAB)

  #
  # dream_att
  #

  set (mex_dream_att_SOURCE_FILES
    mex/mex_dream_att.cc
    src/attenuation.cc
    src/affinity.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_dream_att MODULE
    SRC ${mex_dream_att_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "dream_att"
  )

  set_target_properties (mex_dream_att PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_dream_att PUBLIC DREAM_MATLAB)

  #
  # Analytic SIRs.
  #

  # Circ

  set (mex_circ_sir_SOURCE_FILES
    extra_src/mex_circ_sir.cc
    src/affinity.cc
    extra_src/circ_sir.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_circ_sir MODULE
    ${mex_circ_sir_SOURCE_FILES}
    OUTPUT_NAME "circ_sir"
  )

  set_target_properties (mex_circ_sir PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_circ_sir PUBLIC DREAM_MATLAB)

  # Sampled analytic circ

  set (mex_scirc_sir_SOURCE_FILES
    extra_src/mex_scirc_sir.cc
    extra_src/scirc_sir.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_scirc_sir MODULE
    SRC ${mex_scirc_sir_SOURCE_FILES}
    OUTPUT_NAME "scirc_sir"
  )

  set_target_properties (mex_scirc_sir PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_scirc_sir PUBLIC DREAM_MATLAB)

  # Rect

  set (mex_rect_sir_SOURCE_FILES
    extra_src/mex_rect_sir.cc
    extra_src/rect_sir.cc
    src/affinity.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_rect_sir MODULE
    SRC ${mex_rect_sir_SOURCE_FILES}
    OUTPUT_NAME "rect_sir"
  )

  set_target_properties (mex_rect_sir PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_rect_sir PUBLIC DREAM_MATLAB)

  #
  # Delay-and-sum (DAS)
  #

  # das

  set (mex_das_SOURCE_FILES
    extra_src/mex_das.cc
    extra_src/das.cc
    src/affinity.cc
    src/dream_error.cc
  )

  if (OpenCL_FOUND)
    set (mex_das_SOURCE_FILES ${mex_das_SOURCE_FILES} "opencl/cl_das.cc")
    file (COPY opencl/das.cl DESTINATION kernels)
  endif (OpenCL_FOUND)

  matlab_add_mex (NAME mex_das MODULE
    SRC ${mex_das_SOURCE_FILES}
    LINK_TO ${OpenCL_LIBRARIES}
    OUTPUT_NAME "das"
  )

  set_target_properties (mex_das PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_das PUBLIC DREAM_MATLAB)

  #
  # Convolution algorithms
  #

  # conv_p

  set (mex_conv_p_SOURCE_FILES
    extra_src/mex_conv_p.cc
    src/affinity.cc
    extra_src/conv.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_conv_p MODULE
    SRC ${mex_conv_p_SOURCE_FILES}
    OUTPUT_NAME "conv_p"
  )

  set_target_properties (mex_conv_p PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_conv_p PUBLIC DREAM_MATLAB)

  # fftconv_p

  set (mex_fftconv_p_SOURCE_FILES
    extra_src/mex_fftconv_p.cc
    extra_src/fftconv.cc
    src/affinity.cc
    src/dream_error.cc
    )

  matlab_add_mex (NAME mex_fftconv_p MODULE
    SRC ${mex_fftconv_p_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "fftconv_p"
  )

  set_target_properties (mex_fftconv_p PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_fftconv_p PUBLIC DREAM_MATLAB)

  # sum_fftconv

  # This one fails on macOS and Windows
  if(NOT MACOSX AND NOT WIN32)

    set (mex_sum_fftconv_SOURCE_FILES
      extra_src/mex_sum_fftconv.cc
      extra_src/fftconv.cc
      src/affinity.cc
      src/dream_error.cc
    )

    matlab_add_mex (NAME mex_sum_fftconv MODULE
      SRC ${mex_sum_fftconv_SOURCE_FILES}
      LINK_TO ${FFTW_LIBRARIES}
      OUTPUT_NAME "sum_fftconv"
    )

    set_target_properties (mex_sum_fftconv PROPERTIES
      CXX_STANDARD 14
      INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
    )

    target_compile_definitions(mex_sum_fftconv PUBLIC DREAM_MATLAB)

  endif(NOT MACOSX AND NOT WIN32)

  # fftconv_ola

  set (mex_fftconv_ola_SOURCE_FILES
    extra_src/mex_fftconv_ola.cc
    extra_src/fftconv.cc
    src/affinity.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_fftconv_ola MODULE
    SRC ${mex_fftconv_ola_SOURCE_FILES}
    LINK_TO ${FFTW_LIBRARIES}
    OUTPUT_NAME "fftconv_ola"
  )

  set_target_properties (mex_fftconv_ola PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_fftconv_ola PUBLIC DREAM_MATLAB)

  #
  # Add
  #

  # No mex file yet

  #
  # Copy
  #

  set (mex_copy_p_SOURCE_FILES
    extra_src/mex_copy_p.cc
    src/affinity.cc
    src/dream_error.cc
  )

  matlab_add_mex (NAME mex_copy_p MODULE
    SRC ${mex_copy_p_SOURCE_FILES}
    LINK_TO ${Matlab_LIBRARIES}
    OUTPUT_NAME "copy_p"
  )

  set_target_properties (mex_copy_p PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
  )

  target_compile_definitions(mex_copy_p PUBLIC DREAM_MATLAB)

endif (Matlab_FOUND)
