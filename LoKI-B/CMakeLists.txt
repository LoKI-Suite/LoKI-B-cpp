# Install all headers in directory "LoKI-B/" relative to the installation path.
# Exclude automake and cmake files.

install(
    DIRECTORY .
    DESTINATION LoKI-B
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "Makefile.*" EXCLUDE
    PATTERN "CMakeLists.txt" EXCLUDE
  )

