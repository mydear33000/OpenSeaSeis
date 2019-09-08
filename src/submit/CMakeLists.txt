FILE(GLOB_RECURSE CPP_SRCS *.cc)

ADD_EXECUTABLE(as_submit ${CPP_SRCS})

TARGET_LINK_LIBRARIES(as_submit ${OpenSeaSeis_LINKER_LIBS}
                                ${CMAKE_DL_LIBS})

INSTALL( TARGETS as_submit
         LIBRARY DESTINATION lib
         ARCHIVE DESTINATION lib
         RUNTIME DESTINATION bin
)
