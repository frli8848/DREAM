# Found the FFTW lib bundeled with Matlab
#
# MWFFTW_INCLUDE_DIR = fftw3.h
# MWFFTW_LIBRARIES = libfftw3.a
# MWFFTW_FOUND = true if MWFFTW3 is found

if (WIN32)
  set(CMAKE_PREFIX_PATH "C:/MWFFTW64")
  #set(CMAKE_PREFIX_PATH "C:/MWFFTW")
endif(WIN32)

if (MACOSX)
  set(CMAKE_PREFIX_PATH "/opt/miniconda3")
  set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};/usr/local/miniconda")
endif(MACOSX)

message (STATUS "MWFFTW Matlab_ROOT_DIR ${Matlab_ROOT_DIR}")

# NB. Matlab do not distribute the header file so we
# so we use the upstream one.
find_path (MWFFTW_INCLUDE_DIR fftw3.h PATHS include)
if (NOT MWFFTW_INCLUDE_DIR)
  message (STATUS "Could not find fftw3.h")
endif (NOT MWFFTW_INCLUDE_DIR)

find_path (MWFFTW_LIB_PATH
  NAMES libmwfftw3 mwfftw3 libmwfftw3.so.3 mwfftw3.so.3 libmwfftw3.dylib mwfftw3.dylib
  PATHS "${Matlab_ROOT_DIR}/bin/"
  PATH_SUFFIXES glnxa64 maci64
  )

if (NOT MWFFTW_LIB_PATH)
  message (STATUS "Could not find Matlab library path")
endif (NOT MWFFTW_LIB_PATH)

find_library (MWFFTW_LIBRARIES
  NAMES libmwfftw3 mwfftw3 libmwfftw3.so.3 mwfftw3.so.3 libmwfftw3.dylib mwfftw3.dylib
  PATHS "${Matlab_ROOT_DIR}/bin/"
  PATH_SUFFIXES glnxa64 maci64
  NO_DEFAULT_PATH
  )

if (NOT MWFFTW_LIBRARIES)
  message (STATUS "Could not find the double precision Mathworks FFTW library")
endif (NOT MWFFTW_LIBRARIES)

find_library (SMWFFTW_LIBRARIES
  NAMES libmwfftw3f mwfftw3f libmwfftw3f.so.3 mwfftw3f.so.3 libmwfftw3f.dylib mwfftw3f.dylib
  PATHS "${Matlab_ROOT_DIR}/bin"
  PATH_SUFFIXES glnxa64
  NO_DEFAULT_PATH
  )

if (NOT SMWFFTW_LIBRARIES)
  message (STATUS "Could not find the single precision Mathworks FFTW library")
endif (NOT SMWFFTW_LIBRARIES)

if (MWFFTW_INCLUDE_DIR AND MWFFTW_LIBRARIES AND SMWFFTW_LIBRARIES)
  set (MWFFTW_FOUND TRUE)
endif (MWFFTW_INCLUDE_DIR AND MWFFTW_LIBRARIES AND SMWFFTW_LIBRARIES)

if (MWFFTW_FOUND)

  if (NOT MWFFTW_FIND_QUIETLY)
    message (STATUS "Found fftw3.h: ${MWFFTW_INCLUDE_DIR}")
  endif (NOT MWFFTW_FIND_QUIETLY)

  if (NOT MWFFTW_FIND_QUIETLY)
    message (STATUS "Found MWFFTW: ${MWFFTW_LIBRARIES}")
  endif (NOT MWFFTW_FIND_QUIETLY)

  if (NOT MWFFTW_FIND_QUIETLY)
    message (STATUS "Found MWFFTW: ${SMWFFTW_LIBRARIES}")
  endif (NOT MWFFTW_FIND_QUIETLY)

else (MWFFTW_FOUND)
  if (MWFFTW_FIND_REQUIRED)
    message (FATAL_ERROR "Could not find MWFFTW")
  endif (MWFFTW_FIND_REQUIRED)
endif (MWFFTW_FOUND)

mark_as_advanced (MWFFTW_INCLUDE_DIR
  MWFFTW_LIBRARIES
  MWFFTW_LIB_PATH
  SMWFFTW_LIBRARIES
  MWFFTW_FOUND)
