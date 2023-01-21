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

#set (JULIA_INCLUDE_DIR /usr/include/julia) # Gentoo

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
  CXX_STANDARD 14
  COMPILE_FLAGS "${DREAM_JULIA_FLAGS}"
  INCLUDE_DIRECTORIES "${DREAM_JULIA_INCLUDE_DIRS}"
  LINK_FLAGS ${JULIA_LD_FLAGS}
  PREFIX "" OUTPUT_NAME "dreamrect"
)
