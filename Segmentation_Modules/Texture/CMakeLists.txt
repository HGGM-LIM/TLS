
cmake_minimum_required(VERSION 2.8)

project(TextureFeatures)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()



add_executable(TextureFeatures MACOSX_BUNDLE TextureFeatures.cxx)
add_executable(VesselnessMeasureImageFilter MACOSX_BUNDLE VesselnessMeasureImageFilter.cxx)
add_executable(SegmentBloodVesselsWithMultiScaleHessianBasedMeasure MACOSX_BUNDLE SegmentBloodVesselsWithMultiScaleHessianBasedMeasure.cxx)


if( "${ITK_VERSION_MAJOR}" LESS 4 ) 
  target_link_libraries(TextureFeatures ITKReview ${ITK_LIBRARIES})
  target_link_libraries(VesselnessMeasureImageFilter ITKReview ${ITK_LIBRARIES})
  target_link_libraries(SegmentBloodVesselsWithMultiScaleHessianBasedMeasure ITKReview ${ITK_LIBRARIES})

 else( "${ITK_VERSION_MAJOR}" LESS 4 ) 
  target_link_libraries(TextureFeatures ${ITK_LIBRARIES})
  target_link_libraries(VesselnessMeasureImageFilter ${ITK_LIBRARIES})
  target_link_libraries(SegmentBloodVesselsWithMultiScaleHessianBasedMeasure ${ITK_LIBRARIES})
  
  
endif( "${ITK_VERSION_MAJOR}" LESS 4 )

