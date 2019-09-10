FILE(GLOB_RECURSE CPP_SRCS *.cc)

ADD_EXECUTABLE(as_cssubmit ${CPP_SRCS})

TARGET_LINK_LIBRARIES(as_cssubmit ${OpenSeaSeis_LINKER_LIBS}
                                ${CMAKE_DL_LIBS})

INSTALL( TARGETS as_cssubmit
         LIBRARY DESTINATION lib
         ARCHIVE DESTINATION lib
         RUNTIME DESTINATION bin
)
