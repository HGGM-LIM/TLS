INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#FIND_PACKAGE( Boost 1.59 COMPONENTS program_options REQUIRED )
find_package(Boost COMPONENTS filesystem REQUIRED)
#INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/Segmentation_Modules/libs/)

#SET(OLD_CXX_FLAGS ${CMAKE_CXX_FLAGS})
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")


find_package(ITK REQUIRED)
include(${ITK_USE_FILE})


#ADD_EXECUTABLE(extract_air SingleThresholdSegmenter_main talk2manager segments manager propagator bifurcations globals_itk leakage log) 
#TARGET_LINK_LIBRARIES(extract_air ${ITK_LIBRARIES} ${LOG4CPP_LIBRARIES} ${Boost_SYSTEM_LIBRARY})
   
#ADD_EXECUTABLE(link_segments link_segments_main log) 
#TARGET_LINK_LIBRARIES(link_segments  ${ITK_LIBRARIES} ${LOG4CPP_LIBRARIES} ${Boost_SYSTEM_LIBRARY})

ADD_EXECUTABLE(evaluateAirways evaluateAirways.cxx)
TARGET_LINK_LIBRARIES(evaluateAirways ${ITK_LIBRARIES} ${Boost_LIBRARIES})

#INSTALL(TARGETS extract_air link_segments evaluateAirways  RUNTIME DESTINATION ${TLS_FILTERS_BIN_INSTALL_DIR})

#SET(CMAKE_CXX_FLAGS ${OLD_CXX_FLAGS})

#ADD_SUBDIRECTORY(wall)


