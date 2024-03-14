#
# Copyright (C) 2023 Fredrik Lingvall
#

project(jl_dream)

# Get the CxxWrap prefix and Julia module paths
if (NOT DEFINED CXXWRAP_PREFIX)
  # NB. This does not work on Windows using julialauncher.exe
  # so one needs to add the prefix manually.
  execute_process(COMMAND julia -E "using CxxWrap; CxxWrap.prefix_path()" OUTPUT_VARIABLE CXXWRAP_PREFIX)
  string (STRIP ${CXXWRAP_PREFIX} CXXWRAP_PREFIX)
  string (REPLACE "\"" "" CXXWRAP_PREFIX ${CXXWRAP_PREFIX})
endif (NOT DEFINED CXXWRAP_PREFIX)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CXXWRAP_PREFIX}/lib/cmake/JlCxx/)
set (CXXWRAP_INCLUDE_DIR ${CXXWRAP_PREFIX}/include)

message (STATUS "Julia CxxWrap PREFIX ${CXXWRAP_PREFIX}")

# FIXME: The Julia CxxWrap macro @wrapmodule only seem to take hard coded paths
# so we cannot read the DREAM install path from a variable defined in a file
# like in ~/.julia/config/startup.jl, for example. For now we just let cmake
# generate the module files with hard coded paths on-the-fly.

set (DREAM_FUNC_LIST
  dreamline
  dreamrect
  dreamrect_f
  dreamcirc
  dreamcirc_f
  dreamsphere
  dreamcylind
  dream_arr_rect
  dream_arr_circ
  dream_arr_cylind
  dream_arr_annu
  conv_p
  fftconv_p
  das
  das_f
)

# FIXME: There must be a cleaner way than having a separate file for each Julia module!?
foreach (D_FUN IN LISTS DREAM_FUNC_LIST)
  #message (STATUS ${D_FUN})
  file (WRITE "${CMAKE_CURRENT_BINARY_DIR}/${D_FUN}_m.jl" "module ${D_FUN}_m\n")
  file (APPEND "${CMAKE_CURRENT_BINARY_DIR}/${D_FUN}_m.jl" "using CxxWrap\n")
  file (APPEND "${CMAKE_CURRENT_BINARY_DIR}/${D_FUN}_m.jl" "@wrapmodule(joinpath(\"${CMAKE_CURRENT_BINARY_DIR}\",\"${D_FUN}\"))\n")
  file (APPEND "${CMAKE_CURRENT_BINARY_DIR}/${D_FUN}_m.jl" "function __init__()\n")
  file (APPEND "${CMAKE_CURRENT_BINARY_DIR}/${D_FUN}_m.jl" "@initcxx\n")
  file (APPEND "${CMAKE_CURRENT_BINARY_DIR}/${D_FUN}_m.jl" "end\n")
  file (APPEND "${CMAKE_CURRENT_BINARY_DIR}/${D_FUN}_m.jl" "end\n")
endforeach()

# NB. CxxWrap uses the C++17 standard!

find_package(Julia REQUIRED)

set (DREAM_JULIA_INCLUDE_DIRS
  ${PROJECT_BINARY_DIR}
  "${PROJECT_BINARY_DIR}/../opencl;"
  "${PROJECT_SOURCE_DIR}/../include;"
  "${PROJECT_SOURCE_DIR}/include;"
  ${Julia_INCLUDE_DIRS}
  ${CXXWRAP_INCLUDE_DIR}
  ${FFTW_INCLUDE_DIR}
)

message (STATUS "DREAM_JULIA_INCLUDE_DIRS: ${DREAM_JULIA_INCLUDE_DIRS}")

set (DREAM_JULIA_FLAGS "-DDREAM_JULIA")
message (STATUS "JULIA_LD_FLAGS ${JULIA_LD_FLAGS}")

if (FFTW_FOUND)
  add_definitions( -DHAVE_FFTW )
  set (JULIA_LD_FLAGS "${JULIA_LD_FLAGS} -lfftw3")
endif (FFTW_FOUND)

#
# dreamline
#

set (jl_dreamline_SOURCE_FILES
  src/jl_dreamline.cc
  ../src/dreamline.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/dream_error.cc
)

add_library(jl_dreamline MODULE ${jl_dreamline_SOURCE_FILES})
target_link_libraries(jl_dreamline
  ${FFTW_LIBRARIES}
)

set_target_properties(jl_dreamline PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamline"
)

#
# dreamrect
#

set (jl_dreamrect_SOURCE_FILES
  src/jl_dreamrect.cc
  ../src/dreamrect.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/dream_error.cc
)

add_library(jl_dreamrect MODULE ${jl_dreamrect_SOURCE_FILES})
target_link_libraries(jl_dreamrect
  ${FFTW_LIBRARIES}
)

set_target_properties(jl_dreamrect PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamrect"
)


#
# dreamrect_f
#

set (jl_dreamrect_f_SOURCE_FILES
  src/jl_dreamrect_f.cc
  ../src/dreamrect_f.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(jl_dreamrect_f MODULE ${jl_dreamrect_f_SOURCE_FILES})
target_link_libraries(jl_dreamrect_f
  ${FFTW_LIBRARIES}
)

set_target_properties(jl_dreamrect_f PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamrect_f"
)

#
# dreamcirc
#

set (jl_dreamcirc_SOURCE_FILES
  src/jl_dreamcirc.cc
  ../src/dreamcirc.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/dream_error.cc
)

add_library(jl_dreamcirc MODULE ${jl_dreamcirc_SOURCE_FILES})
target_link_libraries(jl_dreamcirc
  ${FFTW_LIBRARIES}
)

set_target_properties(jl_dreamcirc PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamcirc"
)

#
# dreamcirc_f
#

set (jl_dreamcirc_f_SOURCE_FILES
  src/jl_dreamcirc_f.cc
  ../src/dreamcirc_f.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(jl_dreamcirc_f MODULE ${jl_dreamcirc_f_SOURCE_FILES})
target_link_libraries(jl_dreamcirc_f
  ${FFTW_LIBRARIES}
)

set_target_properties(jl_dreamcirc_f PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamcirc_f"
)

#
# dreamsphere
#

set (jl_dreamsphere_SOURCE_FILES
  src/jl_dreamsphere.cc
  ../src/dreamsphere.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/dream_error.cc
)

add_library(jl_dreamsphere MODULE ${jl_dreamsphere_SOURCE_FILES})
target_link_libraries(jl_dreamsphere
  ${FFTW_LIBRARIES}
)

set_target_properties(jl_dreamsphere PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamsphere"
)

#
# dreamcylind
#

set (jl_dreamcylind_SOURCE_FILES
  src/jl_dreamcylind.cc
  ../src/dreamcylind.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/dream_error.cc
)

add_library(jl_dreamcylind MODULE ${jl_dreamcylind_SOURCE_FILES})
target_link_libraries(jl_dreamcylind
  ${FFTW_LIBRARIES}
)

set_target_properties(jl_dreamcylind PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamcylind"
)

#
# dream_arr_rect
#

set (jl_dream_arr_rect_SOURCE_FILES
  src/jl_dream_arr_rect.cc
  ../src/dream_arr_rect.cc
  ../src/dreamrect.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(jl_dream_arr_rect MODULE ${jl_dream_arr_rect_SOURCE_FILES})
target_link_libraries(jl_dream_arr_rect
  ${FFTW_LIBRARIES}
)

set_target_properties(jl_dream_arr_rect PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dream_arr_rect"
)

#
# dream_arr_circ
#

set (jl_dream_arr_circ_SOURCE_FILES
  src/jl_dream_arr_circ.cc
  ../src/dream_arr_circ.cc
  ../src/dreamcirc.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(jl_dream_arr_circ MODULE ${jl_dream_arr_circ_SOURCE_FILES})
target_link_libraries(jl_dream_arr_circ
  ${FFTW_LIBRARIES}
)

set_target_properties(jl_dream_arr_circ PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dream_arr_circ"
)

#
# dream_arr_cylind
#

set (jl_dream_arr_cylind_SOURCE_FILES
  src/jl_dream_arr_cylind.cc
  ../src/dream_arr_cylind.cc
  ../src/dreamcylind.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(jl_dream_arr_cylind MODULE ${jl_dream_arr_cylind_SOURCE_FILES})
target_link_libraries(jl_dream_arr_cylind
  ${FFTW_LIBRARIES}
)

set_target_properties(jl_dream_arr_cylind PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dream_arr_cylind"
)

#
# dream_arr_annu
#

set (jl_dream_arr_annu_SOURCE_FILES
  src/jl_dream_arr_annu.cc
  ../src/dream_arr_annu.cc
  ../src/dreamcirc.cc
  ../src/attenuation.cc
  ../src/affinity.cc
  ../src/arr_functions.cc
  ../src/dream_error.cc
)

add_library(jl_dream_arr_annu MODULE ${jl_dream_arr_annu_SOURCE_FILES})
target_link_libraries(jl_dream_arr_annu
  ${FFTW_LIBRARIES}
)

set_target_properties(jl_dream_arr_annu PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dream_arr_annu"
)

#
# conv_p
#

set (jl_conv_p_SOURCE_FILES
  src/jl_conv_p.cc
  ../extra_src/conv.cc
  ../src/affinity.cc
)

add_library(jl_conv_p MODULE ${jl_conv_p_SOURCE_FILES})
target_link_libraries(jl_conv_p)

set_target_properties(jl_conv_p PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "conv_p"
)

#
# fftconv_p
#

set (jl_fftconv_p_SOURCE_FILES
  src/jl_fftconv_p.cc
  ../extra_src/fftconv.cc
  ../src/affinity.cc
)

add_library(jl_fftconv_p MODULE ${jl_fftconv_p_SOURCE_FILES})
target_link_libraries(jl_fftconv_p
  ${FFTW_LIBRARIES}
)

set_target_properties(jl_fftconv_p PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "fftconv_p"
)

#
# Delay-and-sum (DAS)
#

# single

set (jl_das_f_SOURCE_FILES
  src/jl_das.cc
  ../extra_src/das.cc
  ../src/affinity.cc
  ../src/dream_error.cc
)

if (OpenCL_FOUND)
  set (jl_das_f_SOURCE_FILES ${jl_das_f_SOURCE_FILES} "../opencl/cl_das.cc")
endif (OpenCL_FOUND)

add_library (jl_das_f MODULE ${jl_das_f_SOURCE_FILES})

target_link_libraries (jl_das_f
  ${OpenCL_LIBRARIES}
)

set_target_properties (jl_das_f PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS} -DUSE_FLOAT"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  SUFFIX "" PREFIX "" OUTPUT_NAME "das_f")

# Double

set (jl_das_SOURCE_FILES
  src/jl_das.cc
  ../extra_src/das.cc
  ../src/affinity.cc
  ../src/dream_error.cc
)

if (OpenCL_FOUND)
  set (jl_das_SOURCE_FILES ${jl_das_SOURCE_FILES} "../opencl/cl_das.cc")
endif (OpenCL_FOUND)

add_library (jl_das MODULE ${jl_das_SOURCE_FILES})

target_link_libraries (jl_das
  ${OpenCL_LIBRARIES}
)

set_target_properties (jl_das PROPERTIES
  CXX_STANDARD 17
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  SUFFIX "" PREFIX "" OUTPUT_NAME "das")
