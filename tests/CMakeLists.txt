cmake_minimum_required(VERSION 3.10)

message (STATUS "Current CTest bin dir: ${CMAKE_CURRENT_BINARY_DIR}")
message (STATUS "Current CTest PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

enable_testing()

#
# Here we have the C++ tests only. The tests for each binding is found in the
# corresponding directories.
#

#
# Platform specific settings - typically include directories and link libraries.
#
#include("${CMAKE_HOST_SYSTEM_NAME}.cmake")

#
# User specific settings - e.g. to override platform file locations or play with compiler options.
#
#include("user.cmake" OPTIONAL)

# Check that we have FFTW installed
find_package (FFTW)  # Will define FFTW_INCLUDE_DIR, FFTW_LIBRARIES, FFTW_STATIC_LIBRARIES, and FFTW_FOUND

set (UNIT_TEST_LIST
  test_fftconv_p
)

set (DREAM_SRC
  "${PROJECT_SOURCE_DIR}/src/affinity.cc"
  "${PROJECT_SOURCE_DIR}/extra_src/fftconv.cc"
)

# Use the list to make executables and what to link with (the dependency)
foreach(NAME IN LISTS UNIT_TEST_LIST)

  add_executable(${NAME} src/${NAME}.cc ${DREAM_SRC})

  set_target_properties (${NAME} PROPERTIES
    CXX_STANDARD 14
    INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/include;"
    )

  target_link_libraries(${NAME} PUBLIC "${FFTW_LIBRARIES}")

endforeach()

#
# C++ tests
#

# fftconv_p
add_test(test_fftconv_p_run ${CMAKE_CURRENT_BINARY_DIR}/test_fftconv_p)
