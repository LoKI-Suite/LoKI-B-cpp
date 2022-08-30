# Recursively install all ".in", ".json" and ".txt" files from this directory
# in directory "input/mixture" relative to the installation path.
# Exclude automake and cmake files.

install(
    DIRECTORY .
    DESTINATION Input
    FILES_MATCHING
    PATTERN "*.in"
    PATTERN "*.json"
    PATTERN "*.txt"
    PATTERN "Makefile.*" EXCLUDE
    PATTERN "CMakeLists.txt" EXCLUDE
  )

