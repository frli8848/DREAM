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


# The Julia CxxWrap macro @@wrapmodule only seem to take hard coded paths
# so we cannot read the DREAM install variable from variable definbed in in file
# like in ~/.julia/config/startup.jl, for example. For now we just let cmake
# generate the file with hard coded path on-the-fly.

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/DREAM.jl" "module DREAM\n")
file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/DREAM.jl" "using CxxWrap\n")
file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/DREAM.jl" "@wrapmodule(joinpath(\"${CMAKE_CURRENT_BINARY_DIR}\",\"dreamrect\"))\n")
file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/DREAM.jl" "function __init__()\n")
file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/DREAM.jl" "@initcxx\n")
file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/DREAM.jl" "end\n")
file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/DREAM.jl" "end\n")

#file (COPY jl_files/DREAM.jl DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
#message (STATUS "Julia DREAM MODULE PATH ${CMAKE_CURRENT_BINARY_DIR}")

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
