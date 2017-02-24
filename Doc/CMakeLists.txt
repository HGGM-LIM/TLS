try_compile(
  RESULT_VAR
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/CMakeDoxygenFilter.cxx
  COMPILE_DEFINITIONS "-DUSE_NAMESPACE=CMake"
  OUTPUT_VARIABLE COMPILATION_OUTPUT
  COPY_FILE ${CMAKE_CURRENT_BINARY_DIR}/CMakeDoxygenFilter${CMAKE_EXECUTABLE_SUFFIX})
if(NOT RESULT_VAR)
  message(FATAL_ERROR "Faild to compile Documentation/CMakeDoxygenFilter.cxx\n${COMPILATION_OUTPUT}")
endif()

set(M2CPP "echo")
if(USE_MATLAB)
  file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/m2cpp.pl
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  if(WIN32)
    set(M2CPP ${CMAKE_CURRENT_BINARY_DIR}/m2cpp.bat)
    file(WRITE ${M2CPP} "\"${CMA_MATLAB_DIR}/sys/perl/win32/bin/perl.exe\" m2cpp.pl %1")
  elseif(UNIX)
    set(M2CPP ${CMAKE_CURRENT_BINARY_DIR}/m2cpp.sh)
    set(M2CPP_ ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/m2cpp.sh)
    file(WRITE ${M2CPP_} "perl m2cpp.pl $1")
    file(
      COPY ${M2CPP_}
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
      FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
  else()
    message(FATAL_ERROR "Platform is not supported.")
  endif()
endif()

set(DOXYPYPY "echo")
if(USE_Python)
  file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/doxypypy.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  if(WIN32)
    set(DOXYPYPY ${CMAKE_CURRENT_BINARY_DIR}/doxypypy.bat)
    file(WRITE ${DOXYPYPY} "${CMA_PYTHON_EXECUTABLE} -m doxypypy -a -c %1")
  elseif(UNIX)
    set(DOXYPYPY ${CMAKE_CURRENT_BINARY_DIR}/doxypypy.sh)
    set(DOXYPYPY_ ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/doxypypy.sh)
    file(WRITE ${DOXYPYPY_} "${CMA_PYTHON_EXECUTABLE} -m doxypypy -a -c $1")
    file(
      COPY ${DOXYPYPY_}
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
      FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
  else()
    message(FATAL_ERROR "Platform is not supported.")
  endif()
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  @ONLY)
add_custom_target(doc
  COMMAND ${CMA_DOXYGEN_EXECUTABLE} Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
