# - Find Geos
# The following setings are defined
# GEOS_ROOT_DIR, the root of the include and lib directory
# GEOS_INCLUDE_DIR, the full path of the include dir (ADVANCED)
# GEOS_LIBRARIES, the name of the Geos library (ADVANCED)
#####
# Find Geos

# Look for a root installation
FIND_PATH(GEOS_ROOT_DIR include/
  /usr
  DOC "The root of an installed Geos installation"
)

# try to find the header
FIND_PATH(GEOS_INCLUDE_DIR geos 
  ${GEOS_ROOT_DIR}/include
  /usr/include 
  /usr/local/include
  /usr/local/include/geos
)

# Find the library
FIND_LIBRARY(GEOS_LIBRARY
   NAMES geos 
   PATHS
     ${GEOS_ROOT_DIR}/lib
     /usr/lib 
     /usr/local/lib
   DOC "The name of the Geos library"
)
IF (GEOS_ROOT_DIR)
  IF (GEOS_INCLUDE_DIR AND GEOS_LIBRARY)
    SET (GEOS_FOUND TRUE)
    SET (GEOS_LIBRARIES "${GEOS_LIBRARY}")
    # FIXME: There should be a better way of handling this?
    # FIXME: How can we test to see if the lib dir isn't 
    # FIXME: one of the default dirs?
    LINK_DIRECTORIES(${GEOS_ROOT_DIR}/lib)
  ENDIF (GEOS_INCLUDE_DIR AND GEOS_LIBRARY)
ENDIF (GEOS_ROOT_DIR)

IF (GEOS_FOUND)
  IF (NOT GEOS_FIND_QUIETLY)
    MESSAGE (STATUS "Found Geos: ${GEOS_LIBRARY}")
  ENDIF (NOT GEOS_FIND_QUIETLY)
ELSE (GEOS_FOUND)
  IF (GEOS_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find Geos")
  ENDIF (GEOS_FIND_REQUIRED)
ENDIF (GEOS_FOUND)

MARK_AS_ADVANCED(
  GEOS_INCLUDE_DIR
  GEOS_LIBRARY
)

