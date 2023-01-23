#
# Copyright (C) 2023 Fredrik Lingvall
#

project(jl_dream)

find_package(FFTW)

# Get the CxxWrap prefix and Julia module paths
execute_process(COMMAND julia -E "using CxxWrap; CxxWrap.prefix_path()" OUTPUT_VARIABLE CXXWRAP_PREFIX)
string (STRIP ${CXXWRAP_PREFIX} CXXWRAP_PREFIX)
string (REPLACE "\"" "" CXXWRAP_PREFIX ${CXXWRAP_PREFIX})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CXXWRAP_PREFIX}/lib/cmake/JlCxx/)
set (CXXWRAP_INCLUDE_DIR ${CXXWRAP_PREFIX}/include)

message (STATUS "Julia CxxWrap PREFIX ${CXXWRAP_PREFIX}")

# FIXME: The Julia CxxWrap macro @@wrapmodule only seem to take hard coded paths
# so we cannot read the DREAM install path from a variable defined in a file
# like in ~/.julia/config/startup.jl, for example. For now we just let cmake
# generate the module files with hard coded paths on-the-fly.

set (DREAM_FUNC_LIST
  dreamline
  dreamrect
  dreamcirc
)

# FIXME: There must be a cleaner way than having a separate file for each Julia module!?
foreach (D_FUN IN LISTS DREAM_FUNC_LIST)
  message (STATUS ${D_FUN})
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
  "${PROJECT_SOURCE_DIR}/../include"
  ${Julia_INCLUDE_DIRS}
  ${CXXWRAP_INCLUDE_DIR}
  ${FFTW_INCLUDE_DIR}
)

message (STATUS "DREAM_JULIA_INCLUDE_DIRS: ${DREAM_JULIA_INCLUDE_DIRS}")

set (DREAM_JULIA_FLAGS "-DDREAM_JULIA")
message (STATUS "JULIA_LD_FLAGS ${JULIA_LD_FLAGS}")

if (FFTW_FOUND)
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
