# This function will prevent in-source builds
# Code based on:
# https://stackoverflow.com/questions/1208681/with-cmake-how-would-you-disable-in-source-builds 
# https://github.com/InsightSoftwareConsortium/ITK/blob/master/CMake/PreventInSourceBuilds.cmake
function(AssureOutOfSourceBuilds)
  # make sure the user doesn't play dirty with symlinks
  get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

  # disallow in-source builds
  if("${srcdir}" STREQUAL "${bindir}")
    message("######################################################")
    message("# Plasimo should not be configured & built in the Plasimo source directory")
    message("# You must run cmake in a build directory.")
    message("# For example:")
    message("# mkdir build-cmake ; cd build-cmake")
    message("# Then configure using:")
    message("# cmake .. # or ccmake, or cmake-gui ")
    message("# Finally, build using:")
    message("# make")
    message("#")
    message("# NOTE: Given that you already tried to make an in-source build")
    message("#       CMake have already created several files & directories")
    message("#       in your source tree. run 'git status' to find them and")
    message("#       remove them.")
    message("######################################################")
    message(FATAL_ERROR "Quitting configuration")
  endif()
endfunction()

AssureOutOfSourceBuilds()