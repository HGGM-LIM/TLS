find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/Segmentation_Modules/libs/)

#set(CMAKE_CXX_FLAGS "-lboost_program_options")

#find_package(Boost 1.59 COMPONENTS program_options REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)

ADD_EXECUTABLE(enf_ext itkEmphysemaCalculationsImageFilter.txx enf_ext_main.cxx)
TARGET_LINK_LIBRARIES(enf_ext ${ITK_LIBRARIES} ${Boost_LIBRARIES})					

ADD_EXECUTABLE(quantify_enf_mask quantify_enf_mask.cxx)
TARGET_LINK_LIBRARIES(quantify_enf_mask ${ITK_LIBRARIES})

ADD_EXECUTABLE(quantify2d	 quantify2d.cxx)
TARGET_LINK_LIBRARIES(quantify2d ${ITK_LIBRARIES})

ADD_EXECUTABLE(nova NOVAFilter.cxx)
TARGET_LINK_LIBRARIES(nova ${ITK_LIBRARIES} ${Boost_LIBRARIES})
					
ADD_EXECUTABLE(nova2 NOVAFilter2.cxx)
TARGET_LINK_LIBRARIES(nova2 ${ITK_LIBRARIES} ${Boost_LIBRARIES})					

ADD_EXECUTABLE(nova3 NOVAFilter3.cxx)
TARGET_LINK_LIBRARIES(nova3 ${ITK_LIBRARIES} ${Boost_LIBRARIES})
					
ADD_EXECUTABLE(moving_aver MovingAverageFilter.cxx)
TARGET_LINK_LIBRARIES(moving_aver ${ITK_LIBRARIES} ${Boost_LIBRARIES})	
#ADD_EXECUTABLE(visual_score visual_score_main.cxx )
#TARGET_LINK_LIBRARIES(visual_score ${ITK_LIBRARIES} ${Boost_SYSTEM_LIBRARY})				

#INSTALL(TARGETS moving_aver nova3 enf_ext quantify_enf_mask quantify2d RUNTIME DESTINATION ${RALPH_FILTERS_BIN_INSTALL_DIR})
