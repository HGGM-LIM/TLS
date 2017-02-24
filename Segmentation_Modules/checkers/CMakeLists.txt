#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
cmake_minimum_required(VERSION 2.8)
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

ADD_EXECUTABLE(quantify_volume quantify_volume_main.cxx)
TARGET_LINK_LIBRARIES(quantify_volume ${ITK_LIBRARIES})

ADD_EXECUTABLE(copy_metadata copy_metadata_main.cxx)
TARGET_LINK_LIBRARIES(copy_metadata ${ITK_LIBRARIES})

#INSTALL(TARGETS quantify_volume RUNTIME DESTINATION ${RALPH_FILTERS_BIN_INSTALL_DIR})
