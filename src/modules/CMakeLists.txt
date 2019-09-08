MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

FOREACH(subdir ${SUBDIRS})    
    FILE(GLOB_RECURSE CPP_SRCS ${subdir}/*.cc)

    ADD_LIBRARY(as_${subdir} ${CPP_SRCS})

    TARGET_LINK_LIBRARIES(as_${subdir} ${OpenSeaSeis_LINKER_LIBS})

    INSTALL( TARGETS as_${subdir}
             LIBRARY DESTINATION lib
             ARCHIVE DESTINATION lib
             RUNTIME DESTINATION bin
    )
ENDFOREACH()

