SET( INCLUDES
    Dependencies
)

FOREACH( INC ${INCLUDES} )
    INCLUDE( ${CMKDIR}/${INC}.cmake )
ENDFOREACH()
