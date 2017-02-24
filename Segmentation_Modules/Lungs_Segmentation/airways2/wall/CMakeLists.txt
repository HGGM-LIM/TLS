PROJECT(AirwaySectioning)
 
     cmake_minimum_required(VERSION 2.4)

FIND_PACKAGE(XercesC)
FIND_PACKAGE(Boost COMPONENTS filesystem system REQUIRED)
FIND_PACKAGE(Geos REQUIRED)

# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

SET(Required_Packages
  ITK
)
execute_process(COMMAND g++ -c ${PROJECT_SOURCE_DIR}/tests/Measure-odb.cxx RESULT_VARIABLE rv)
message("rv='${rv}'")
execute_process(COMMAND g++ -c ${PROJECT_SOURCE_DIR}/tests/ImageEntity-odb.cxx RESULT_VARIABLE rv)
message("rv='${rv}'")


#ADD_EXECUTABLE(AirwayWall_Phantom_NoSQL AirwayWall_Phantom_NoSQL.cxx)
ADD_EXECUTABLE(AirwayWall_Phantom AirwayWall_Phantom.cxx Measure-odb.o ImageEntity-odb.o)
ADD_EXECUTABLE(AirwayWall_Human AirwayWall_Human.cxx Measure-odb.o ImageEntity-odb.o)
#ADD_EXECUTABLE(AirwayWall_Human_NoSQL AirwayWall_Human_NoSQL.cxx)


# TARGET_LINK_LIBRARIES(AirwayWall_Phantom_NoSQL
#ITKAlgorithms 
#ITKIO 
#ITKBasicFilters 
#ITKCommon
#ITKStatistics
#${GEOS_LIBRARY}
#${Boost_FILESYSTEM_LIBRARY} 
#${Boost_SYSTEM_LIBRARY}
#${LMFIT_LIBRARY}
#${XERCES_LIBRARY}
#${XERCESC_LIBRARY}
#)
 TARGET_LINK_LIBRARIES(AirwayWall_Phantom
ITKAlgorithms 
ITKIO 
ITKBasicFilters 
ITKCommon
ITKStatistics
${GEOS_LIBRARY}
${Boost_FILESYSTEM_LIBRARY} 
${Boost_SYSTEM_LIBRARY}
${LMFIT_LIBRARY}
#${XERCES_LIBRARY}
${XERCESC_LIBRARY}
odb
odb-mysql
)
# TARGET_LINK_LIBRARIES(AirwayWall_Human_NoSQL
#ITKAlgorithms 
#ITKIO 
#ITKBasicFilters 
#ITKCommon
#ITKStatistics
#${GEOS_LIBRARY}
#${Boost_FILESYSTEM_LIBRARY} 
#${Boost_SYSTEM_LIBRARY}
#${LMFIT_LIBRARY}
##${XERCES_LIBRARY}
#${XERCESC_LIBRARY}
#)
 TARGET_LINK_LIBRARIES(AirwayWall_Human
ITKAlgorithms 
ITKIO 
ITKBasicFilters 
ITKCommon
ITKStatistics
${GEOS_LIBRARY}
${Boost_FILESYSTEM_LIBRARY} 
${Boost_SYSTEM_LIBRARY}
${LMFIT_LIBRARY}
#${XERCES_LIBRARY}
${XERCESC_LIBRARY}
odb
odb-mysql
)

INSTALL(TARGETS AirwayWall_Human AirwayWall_Phantom
    RUNTIME DESTINATION ${RALPH_FILTERS_BIN_INSTALL_DIR}
)
