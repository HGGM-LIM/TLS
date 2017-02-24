#Change PROJECT_NAME to the name of your project
PROJECT(propagator)

# required to avoid warnings while running cmake
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists


FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})


#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(Libraries
  ITKCommon
  ITKBasicFilters
  ITKIO
  ITKStatistics
)


# install devel files
#OPTION(INSTALL_DEVEL_FILES "Install C++ headers" ON)
#IF(INSTALL_DEVEL_FILES)
#	FILE(GLOB develFiles *.h *.txx) 
#	FOREACH(f ${develFiles})
#	  INSTALL_FILES(/include/InsightToolkit/BasicFilters FILES ${f})
#	ENDFOREACH(f)
#ENDIF(INSTALL_DEVEL_FILES)

# option for wrapping
OPTION(BUILD_WRAPPERS "Wrap restartable fast marching" OFF)
IF(BUILD_WRAPPERS)
   FIND_PACKAGE(WrapITK REQUIRED)



   WRAP_LIBRARY(RestartableFastMarching)

  SET(WRAPPER_LIBRARY_GROUPS
    itkFastMarchingRestartableImageFilter
  )

  SET(WRAPPER_LIBRARY_DEPENDS VXLNumerics Base Numerics)

  SET(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)

  AUTO_INCLUDE_MODULES()
  END_WRAP_LIBRARY()
ENDIF(BUILD_WRAPPERS)


#SET(CurrentExe "test_fm")
#ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx)
#TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries} ${ITK_LIBRARIES})

SET(CurrentExe "fast_segmenter")
ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

SET(CurrentExe "create_speed_image")
ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries} ${ITK_LIBRARIES})


#the following line is an example of how to add a test to your project.
#Testname is the title for this particular test.  ExecutableToRun is the
#program which will be running this test.  It can either be a part of this
#project or an external executable.  After that list any args that are needed
#for this test.  Include as many tests as you like.  If your project doesn't have
#any tests you can comment out or delete the following line.
# ADD_TEST(Testname ExecutableToRun arg1 arg2 arg3)

#ADD_TEST(EqualAttributesUnique1 ${TEST_COMMAND}
#  eq_attrib_unique
#  ${CMAKE_SOURCE_DIR}/images/simple_unique.png label_unique-1.png
#  1
#  --compare label_unique-1.png ${CMAKE_SOURCE_DIR}/images/label_unique-1.png
#)
