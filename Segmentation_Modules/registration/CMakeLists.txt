
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

find_package(ITK REQUIRED)
include(${ITK_USE_FILE}  )

find_package(Boost COMPONENTS program_options REQUIRED)

ADD_EXECUTABLE(InvertDeformationField InvertDeformationField.cxx) 
TARGET_LINK_LIBRARIES(InvertDeformationField ${ITK_LIBRARIES} )

ADD_EXECUTABLE(TestFieldInversion TestFieldInversion_main.cxx) 
TARGET_LINK_LIBRARIES(TestFieldInversion ${ITK_LIBRARIES} ${Boost_LIBRARIES})

#INSTALL(TARGETS InvertDeformationField TestFieldInversion RUNTIME DESTINATION ${TLS_FILTERS_BIN_INSTALL_DIR})
