# In this directory, a target lokib_doxygen is configured.
# (This will not be built by default, but you can do e.g.
# 'make lokib_doxygen'

message(STATUS "Creating Doxygen configuration file...")
# set input and output files
set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/loki-b.doxy.in)
set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/loki-b.doxy)

# request to configure the file
configure_file(loki-b.doxy.in loki-b.doxy @ONLY)

# to always build the doxygen docs, comment in the ALL option
add_custom_target( lokib_doxygen
    # ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM )
