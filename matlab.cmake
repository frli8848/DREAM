#
# Copyright (C) 2022 Fredrik Lingvall
#

find_package (Matlab)
message (STATUS "Matlab_MEX_EXTENSION = ${Matlab_MEX_EXTENSION}")
#message (STATUS "MATLAB_INCLUDE_DIR = ${MATLAB_INCLUDE_DIR}")
#message (STATUS "MATLAB_LIBARARIES = ${MATLAB_LIBRARIES}")

if (Matlab_FOUND)

  set (DREAM_MEX_FLAGS "-DDREAM_MATLAB") # Make Matlab use its internal FFT:s

  # mex libs
  set (MEX_LD_FLAGS "${CMAKE_LD_FLAGS} ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY} -lstdc++")
  message (STATUS "Matlab_INCLUDE_DIRS ${Matlab_INCLUDE_DIRS}")
  message (STATUS "Matlab_LIBRARIES  ${Matlab_LIBRARIES}")
  message (STATUS "MEX_LD_FLAGS ${MEX_LD_FLAGS}")
  if (FFTW_FOUND)
    set (MEX_LD_FLAGS "${MEX_LD_FLAGS} -lfftw3")
  endif (FFTW_FOUND)

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

  add_library (mex_dreamline MODULE
    ${mex_dreamline_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamline
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamline PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamline")


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

  add_library (mex_dreamrect MODULE
    ${mex_dreamrect_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamrect
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamrect PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamrect")

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

  add_library (mex_dreamrect_f MODULE
    ${mex_dreamrect_f_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamrect_f
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamrect_f PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamrect_f")

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

  add_library (mex_dreamcirc MODULE
    ${mex_dreamcirc_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamcirc
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamcirc PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamcirc")

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

  add_library (mex_dreamcirc_f MODULE
    ${mex_dreamcirc_f_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamcirc_f
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamcirc_f PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamcirc_f")

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

  add_library (mex_dreamsphere MODULE
    ${mex_dreamsphere_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamsphere
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamsphere PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamsphere")

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

  add_library (mex_dreamcylind MODULE
    ${mex_dreamcylind_SOURCE_FILES}
    )

  target_link_libraries (mex_dreamcylind
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dreamcylind PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dreamcylind")

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

  add_library (mex_dream_arr_rect MODULE
    ${mex_dream_arr_rect_SOURCE_FILES}
    )

  target_link_libraries (mex_dream_arr_rect
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dream_arr_rect PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dream_arr_rect")

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

  add_library (mex_dream_arr_circ MODULE
    ${mex_dream_arr_circ_SOURCE_FILES}
    )

  target_link_libraries (mex_dream_arr_circ
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dream_arr_circ PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dream_arr_circ")

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

  add_library (mex_dream_arr_cylind MODULE
    ${mex_dream_arr_cylind_SOURCE_FILES}
    )

  target_link_libraries (mex_dream_arr_cylind
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dream_arr_cylind PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dream_arr_cylind")

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

  add_library (mex_dream_arr_annu MODULE
    ${mex_dream_arr_annu_SOURCE_FILES}
    )

  target_link_libraries (mex_dream_arr_annu
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dream_arr_annu PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dream_arr_annu")

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

  add_library (mex_dream_apodwin MODULE
    ${mex_dream_apodwin_SOURCE_FILES}
    )

  target_link_libraries (mex_dream_apodwin
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dream_apodwin PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dream_apodwin")

  #
  # dream_att
  #

  set (mex_dream_att_SOURCE_FILES
    mex/mex_dream_att.cc
    src/attenuation.cc
    src/affinity.cc
    src/dream_error.cc
    )

  add_library (mex_dream_att MODULE
    ${mex_dream_att_SOURCE_FILES}
    )

  target_link_libraries (mex_dream_att
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_dream_att PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "dream_att")


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

  add_library (mex_circ_sir MODULE
    ${mex_circ_sir_SOURCE_FILES}
    )

  target_link_libraries (mex_circ_sir
    ${MATLAB_LIBRARIES}
    )

  set_target_properties (mex_circ_sir PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "circ_sir")

  # Sampled analytic circ

  set (mex_scirc_sir_SOURCE_FILES
    extra_src/mex_scirc_sir.cc
    extra_src/scirc_sir.cc
    src/dream_error.cc
    )

  add_library (mex_scirc_sir MODULE
    ${mex_scirc_sir_SOURCE_FILES}
    )

  target_link_libraries (mex_scirc_sir
    ${MATLAB_LIBRARIES}
    )

  set_target_properties (mex_scirc_sir PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "scirc_sir")

  # Rect

  set (mex_rect_sir_SOURCE_FILES
    extra_src/mex_rect_sir.cc
    extra_src/rect_sir.cc
    src/affinity.cc
    src/dream_error.cc
    )

  add_library (mex_rect_sir MODULE
    ${mex_rect_sir_SOURCE_FILES}
    )

  target_link_libraries (mex_rect_sir
    ${MATLAB_LIBRARIES}
    )

  set_target_properties (mex_rect_sir PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "rect_sir")

  #
  # Delay-and-sum (DAS)
  #

  # das

  set (mex_das_SOURCE_FILES
    extra_src/mex_das.cc
    extra_src/das.cc
    src/dream_error.cc
    )

  add_library (mex_das MODULE
    ${mex_das_SOURCE_FILES}
    )

  target_link_libraries (mex_das
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_das PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "das")

  # das_err

  set (mex_das_arr_SOURCE_FILES
    extra_src/mex_das_arr.cc
    extra_src/das_arr.cc
    src/arr_functions.cc
    src/dream_error.cc
    )

  add_library (mex_das_arr MODULE
    ${mex_das_arr_SOURCE_FILES}
    )

  target_link_libraries (mex_das_arr
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_das_arr PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "das_arr")

  #
  # SAFT
  #

  set (mex_saft_SOURCE_FILES
    extra_src/mex_saft.cc
    src/dream_error.cc
    )

  add_library (mex_saft MODULE
    ${mex_saft_SOURCE_FILES}
    )

  target_link_libraries (mex_saft
    ${MATLAB_LIBRARIES}
    )

  set_target_properties (mex_saft PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "saft")

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

  add_library (mex_conv_p MODULE
    ${mex_conv_p_SOURCE_FILES}
    )

  target_link_libraries (mex_conv_p
    ${MATLAB_LIBRARIES}
    )

  set_target_properties (mex_conv_p PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "conv_p")

  # fftconv_p

  set (mex_fftconv_p_SOURCE_FILES
    extra_src/mex_fftconv_p.cc
    src/affinity.cc
    src/dream_error.cc
    )

  add_library (mex_fftconv_p MODULE
    ${mex_fftconv_p_SOURCE_FILES}
    )

  target_link_libraries (mex_fftconv_p
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_fftconv_p PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "fftconv_p")

  # sum_fftconv

  set (mex_sum_fftconv_SOURCE_FILES
    extra_src/mex_sum_fftconv.cc
    src/affinity.cc
    src/dream_error.cc
    )

  add_library (mex_sum_fftconv MODULE
    ${mex_sum_fftconv_SOURCE_FILES}
    )

  target_link_libraries (mex_sum_fftconv
    ${MATLAB_LIBRARIES}
    ${FFTW_LIBRARIES}
    )

  set_target_properties (mex_sum_fftconv PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"

    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "sum_fftconv")

  # fftconv_ola

  # Not implemented yet.

  #
  # Add
  #

  # No mex file yet

  #
  # Copy
  #

  set (mex_copy_p_SOURCE_FILES
    extra_src/mex_copy_p.cc
    src/dream_error.cc
    )

  add_library (mex_copy_p MODULE
    ${mex_copy_p_SOURCE_FILES}
    )

  target_link_libraries (mex_copy_p
    ${MATLAB_LIBRARIES}
    )

  set_target_properties (mex_copy_p PROPERTIES
    CXX_STANDARD 14
    COMPILE_FLAGS "${DREAM_MEX_FLAGS}"
    INCLUDE_DIRECTORIES "${DREAM_MEX_INCLUDE_DIRS}"
    LINK_FLAGS ${MEX_LD_FLAGS}
    SUFFIX ".${Matlab_MEX_EXTENSION}" PREFIX "" OUTPUT_NAME "copy_p")

endif (Matlab_FOUND)
