# - Find LMFIT
# The following setings are defined
# LMFIT_ROOT_DIR, the root of the include and lib directory
# LMFIT_INCLUDE_DIR, the full path of the include dir (ADVANCED)
# LMFIT_LIBRARIES, the name of the LMFIT library (ADVANCED)
#####
# Find LMFIT

# Look for a root installation
FIND_PATH(LMFIT_ROOT_DIR include/
  /usr
  /usr/local/lib
  DOC "The root of an installed LMFIT installation"
)

# try to find the header
FIND_PATH(LMFIT_INCLUDE_DIR LMFIT/lmmin.h
  ${LMFIT_ROOT_DIR}/include
  /usr/include 
  /usr/local/include
)

# Find the library
FIND_LIBRARY(LMFIT_LIBRARY
   NAMES lmmin
   PATHS
     ${LMFIT_ROOT_DIR}/lib
     /usr/lib 
     /usr/local/lib
   DOC "The name of the LMFIT library"
)
IF (LMFIT_ROOT_DIR)
  IF (LMFIT_INCLUDE_DIR AND LMFIT_LIBRARY)
    SET (LMFIT_FOUND TRUE)
    SET (LMFIT_LIBRARIES "${LMFIT_LIBRARY}")
    # FIXME: There should be a better way of handling this?
    # FIXME: How can we test to see if the lib dir isn't 
    # FIXME: one of the default dirs?
    LINK_DIRECTORIES(${LMFIT_ROOT_DIR}/lib)
  ENDIF (LMFIT_INCLUDE_DIR AND LMFIT_LIBRARY)
ENDIF (LMFIT_ROOT_DIR)

IF (LMFIT_FOUND)
  IF (NOT LMFIT_FIND_QUIETLY)
    MESSAGE (STATUS "Found LMFIT: ${LMFIT_LIBRARY}")
  ENDIF (NOT LMFIT_FIND_QUIETLY)
ELSE (LMFIT_FOUND)
  IF (LMFIT_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find LMFIT")
  ENDIF (LMFIT_FIND_REQUIRED)
ENDIF (LMFIT_FOUND)

MARK_AS_ADVANCED(
  LMFIT_INCLUDE_DIR
  LMFIT_LIBRARY
)

