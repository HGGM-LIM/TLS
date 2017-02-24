cmake_minimum_required(VERSION 2.6)
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/IJMacros.txt)

PROJECT(itkFastBilateralImageFilter)

#SET(Required_Packages ITK)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE}  )

find_package(Boost COMPONENTS program_options REQUIRED)

#FOREACH(Package ${Required_Packages})
#  LOADPACKAGE(${Package})
#ENDFOREACH(Package)

ADD_EXECUTABLE(FastBilateral2D_UC FastBilateral2D_UC_main.cxx)  
TARGET_LINK_LIBRARIES(FastBilateral2D_UC ${ITK_LIBRARIES})

ADD_EXECUTABLE(FastBilateral2D_SS FastBilateral2D_SS_main.cxx)  
TARGET_LINK_LIBRARIES(FastBilateral2D_SS ${ITK_LIBRARIES})

ADD_EXECUTABLE(FastBilateral3D_SS FastBilateral3D_SS_main.cxx)  
TARGET_LINK_LIBRARIES(FastBilateral3D_SS ${ITK_LIBRARIES} ${Boost_LIBRARIES})

#INSTALL(TARGETS FastBilateral3D_SS RUNTIME DESTINATION ${TLS_FILTERS_BIN_INSTALL_DIR})



