# FFTW_INCLUDE_DIR = fftw3.h
# FFTW_LIBRARIES = libfftw3.a
# FFTW_FOUND = true if FFTW3 is found


if (WIN32)
  set(CMAKE_PREFIX_PATH "C:/FFTW64")
  #set(CMAKE_PREFIX_PATH "C:/FFTW")
endif(WIN32)

find_path (FFTW_INCLUDE_DIR fftw3.h PATHS include)
if (NOT FFTW_INCLUDE_DIR)
  message (STATUS "Could not find fftw3.h")
endif (NOT FFTW_INCLUDE_DIR)

find_library (FFTW_LIBRARIES NAMES libfftw3-3 fftw3 PATHS lib)
if (NOT FFTW_LIBRARIES)
  message (STATUS "Could not find the double precision FFTW library")
endif (NOT FFTW_LIBRARIES)

find_library (SFFTW_LIBRARIES NAMES libfftw3f-3 fftw3f PATHS lib)
if (NOT SFFTW_LIBRARIES)
  message (STATUS "Could not find the single precision FFTW library")
endif (NOT SFFTW_LIBRARIES)

if (FFTW_INCLUDE_DIR AND FFTW_LIBRARIES AND SFFTW_LIBRARIES)
  set (FFTW_FOUND TRUE)
endif (FFTW_INCLUDE_DIR AND FFTW_LIBRARIES AND SFFTW_LIBRARIES)

if (FFTW_FOUND)

  if (NOT FFTW_FIND_QUIETLY)
    message (STATUS "Found fftw3.h: ${FFTW_INCLUDE_DIR}")
  endif (NOT FFTW_FIND_QUIETLY)

  if (NOT FFTW_FIND_QUIETLY)
    message (STATUS "Found FFTW: ${FFTW_LIBRARIES}")
  endif (NOT FFTW_FIND_QUIETLY)

  if (NOT FFTW_FIND_QUIETLY)
    message (STATUS "Found FFTW: ${SFFTW_LIBRARIES}")
  endif (NOT FFTW_FIND_QUIETLY)

else (FFTW_FOUND)
  if (FFTW_FIND_REQUIRED)
    message (FATAL_ERROR "Could not find FFTW")
  endif (FFTW_FIND_REQUIRED)
endif (FFTW_FOUND)

mark_as_advanced (FFTW_INCLUDE_DIR
  FFTW_LIBRARIES
  SFFTW_LIBRARIES
  FFTW_FOUND)
