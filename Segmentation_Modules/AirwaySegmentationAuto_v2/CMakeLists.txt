PROJECT(Tree_Auto_v2)
 
cmake_minimum_required(VERSION 2.6)

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


ADD_EXECUTABLE(Tree_Auto_v2 SegmentAirwayTree_FM.cxx)
TARGET_LINK_LIBRARIES(Tree_Auto_v2 ${ITK_LIBRARIES} )
