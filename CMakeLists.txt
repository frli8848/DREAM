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
option(BUILD_PYTHON "Enable building of the Python binings." OFF) # Work in progess
option(BUILD_USERMAN "Enable building of the user manual." OFF)
option(USE_OPENCL "Enable OpenCL support." OFF) # Experimental

if (UNIX AND NOT DARWIN)

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
  message(STATUS "Octave bindings: ${LINUX_DISTRO} ${LINUX_DISTRO_VERSION} Linux detected!")
  message(STATUS "****************************")
  message(STATUS " ")

endif (UNIX AND NOT DARWIN)

if (USE_OPENCL)

  find_package (OpenCL 1.2)

  if (OpenCL_FOUND)
    add_definitions(-DUSE_OPENCL)
  endif (OpenCL_FOUND)

endif (USE_OPENCL)

# Check if we have FFTW installed
find_package (FFTW)  # Will define FFTW_INCLUDE_DIR, FFTW_LIBRARIES, and FFTW_FOUND
if (FFTW_FOUND)
  message (STATUS "${FFTW_LIBRARIES}")
  #add_definitions( -DHAVE_FFTW )	# Build with FFTW support
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
  include(python.cmake)
endif (BUILD_PYTHON)

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
