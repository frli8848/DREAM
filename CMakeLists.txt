#
# Copyright (C) 2016,2019,2021,2022,2023 Fredrik Lingvall
#

cmake_minimum_required(VERSION 3.14)

project (dream-toolbox)

# Add path to the various cmake files to detect Matlab, Octave, FTTW, etc.
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include(GetGitRevision)
git_version("DREAM")

message(STATUS "DREAM_VERSION_MAJOR: ${DREAM_VERSION_MAJOR}")
message(STATUS "DREAM_VERSION_MINOR: ${DREAM_VERSION_MINOR}")
message(STATUS "DREAM_VERSION_PATCH: ${DREAM_VERSION_PATCH}")
message(STATUS "DREAM_VERSION_COMMIT: ${DREAM_VERSION_COMMIT}")
message(STATUS "DREAM_VERSION_SHORT: ${DREAM_VERSION_SHORT}")
message(STATUS "DREAM_VERSION_LONG: ${DREAM_VERSION_LONG}")
message(STATUS "DREAM_VERSION_SHA1: ${DREAM_VERSION_SHA1}")
message(STATUS "DREAM_VERSION: ${DREAM_VERSION}")

# Version number.
#set (DREAM_VERSION_MAJOR 3)
#set (DREAM_VERSION_MINOR 0)
#set (DREAM_VERSION_PATCH 0)

message (STATUS "CMake system name: ${CMAKE_SYSTEM_NAME}")

option(BUILD_MEX "Enable building of the mex-files." OFF)
option(BUILD_OCT "Enable building of the oct-files." ON)
option(BUILD_PYTHON "Enable building of the Python bindings." OFF) # Work in progess
option(BUILD_JULIA "Enable building of the Julia bindings." OFF) # Work in progess
option(BUILD_USERMAN "Enable building of the user manual." OFF)
option(USE_OPENCL "Enable OpenCL support." OFF) # Experimental

#
# Testing
#

option(BUILD_TESTS "Determines whether to build tests." OFF)
if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif (BUILD_TESTS)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(MACOSX TRUE)
endif()

if (UNIX AND NOT MACOSX)

  # Detect which Linux distro were running on
  if(EXISTS /etc/gentoo-release)
    set(LINUX_DISTRO "GENTOO" CACHE STRING "Unix OS (options: UBUNTU/ubuntu or GENTOO/gentoo).")
    execute_process (COMMAND bash -c "cat /etc/gentoo-release | awk \'{ print $5 }\' | tr -d \'\n\'"
      OUTPUT_VARIABLE LINUX_DISTRO_VERSION
      )
  elseif(EXISTS /usr/bin/lsb_release)
    set(LINUX_DISTRO "UBUNTU" CACHE STRING "Unix OS (options: UBUNTU/ubuntu or GENTOO/gentoo).")
    execute_process (COMMAND bash -c "lsb_release -r | awk \'{ print $2 }\' | tr -d \'\n\'"
      OUTPUT_VARIABLE LINUX_DISTRO_VERSION
      )
  else()
    set(LINUX_DISTRO "Unknown" CACHE STRING "Unix OS (options: UBUNTU/ubuntu or GENTOO/gentoo).")
    message(WARNING "Unknown/non-supported Linux distrubution!")
  endif()

  message(STATUS " ")
  message(STATUS "****************************")
  message(STATUS " ${LINUX_DISTRO} ${LINUX_DISTRO_VERSION} Linux detected!")
  message(STATUS "****************************")
  message(STATUS " ")

endif (UNIX AND NOT MACOSX)

if (USE_OPENCL)

  find_package (OpenCL 1.2)

  if (OpenCL_FOUND)

    add_definitions(-DUSE_OPENCL)

    # Generate header files from the cl kernels so we can embed them
    # in the binaries directly.

    find_program(XXD xxd REQUIRED)

    # Put xxd:ed OpenCL kernel sources in the build folder
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/opencl")

    file(GLOB cl_files "opencl/*.cl")
    foreach(file ${cl_files})
      string(REGEX REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/opencl/" "" rel_path_file ${file})

      message (STATUS "${file} -> ${CMAKE_CURRENT_SOURCE_DIR}/include/${rel_path_file}.h")

      #file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${rel_path_file}.h" "const char* cl_kernel_string[] =\n")

      # This one eats all new lines...
      #file(STRINGS ${file} CL_SOURCES)
      #file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/opencl/${rel_path_file}.h" ${CL_SOURCES})

      # This one do not work with lines ending in "//"
      #file(READ ${file} CL_SOURCES)
      #foreach(line ${CL_SOURCES})
      #  file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/${rel_path_file}.h" "\"${line}\\n\"\n")
      #endforeach()
      #file(APPEND "${CMAKE_CURRENT_BINARY_DIR}/${rel_path_file}.h" ";\n")

      # Use external xxd program (part of app-editors/vim-core) # FIXME: Make cross-platform solution (preferable cmake only)
      set (TMP_VAR_NAME "__${rel_path_file}")
      string(REPLACE "." "_" C_VAR_NAME ${TMP_VAR_NAME})
      execute_process(COMMAND ${XXD} -i -n "${C_VAR_NAME}" ${file} "${CMAKE_CURRENT_BINARY_DIR}/opencl/${rel_path_file}.h")

    endforeach()

  endif (OpenCL_FOUND)

endif (USE_OPENCL)

# Check if we have FFTW installed
find_package (FFTW)  # Will define FFTW_INCLUDE_DIR, FFTW_LIBRARIES, FFTW_STATIC_LIBRARIES, and FFTW_FOUND

if (FFTW_FOUND)
  add_definitions( -DHAVE_FFTW )	# Build with FFTW support
endif (FFTW_FOUND)

# Check if we have pthreads (for setting the thread affinity)
find_package(Threads)
if (CMAKE_USE_PTHREADS_INIT)
  message (STATUS "We have pthread support -> enabing thread affinity control!")
  add_definitions( -DHAVE_PTHREAD )
endif(CMAKE_USE_PTHREADS_INIT)

#
# Octave
#

if (BUILD_OCT)
  include(octave.cmake)
endif (BUILD_OCT)

#
# Matlab
#

if (BUILD_MEX)
  include(matlab.cmake)
endif (BUILD_MEX)

#
# Python
#

if (BUILD_PYTHON)
  add_subdirectory(python)
endif (BUILD_PYTHON)

#
# Julia
#

if (BUILD_JULIA)
  add_subdirectory(julia)
endif (BUILD_JULIA)

#
# User manual (N.B., we need CMake 3.x on Windows due to a MikTeX related bug).
#

if (BUILD_USERMAN)

  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/version.tex" "\\newcommand{\\version}{\\texttt{${DREAM_VERSION_SHA1}}\\xspace}\r")
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/revision.tex" "\\rhead{\\date \\texttt{${DREAM_VERSION_SHA1}} }\r")

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
