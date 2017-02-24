INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SUBDIRS(lungs)
SUBDIRS(airways)
SUBDIRS(airways2)
SUBDIRS(emphysema)
SUBDIRS(managed_client_example1)


#ADD_SUBDIRECTORY(visitor)
ADD_SUBDIRECTORY(evaluators)

SET(modules lungs airways emphysema .)

SET(scripts lungs_analyze.sh)
FILE( COPY lungs_analyze.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
#INSTALL_FILES( ${TLS_FILTERS_UNPREFIXED_BIN_INSTALL_DIR} FILES ${scripts} )

#IF(INSTALL_DEVEL_FILES)
#    FOREACH(m ${modules})
#        FILE(GLOB develFiles ${m}/*.h ${m}/*.txx ${m}/*.cxx) 
#        FOREACH(f ${develFiles})
#          INSTALL_FILES(${TLS_FILTERS_DEVEL_INSTALL_DIR} FILES ${f})
#        ENDFOREACH(f)
#    ENDFOREACH(m)
#ENDIF(INSTALL_DEVEL_FILES)
