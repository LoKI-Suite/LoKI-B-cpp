cmake_minimum_required(VERSION 3.13)
project(loki_cpp)

set(CMAKE_CXX_STANDARD 17)

# --- Compiler flags --- #

if(MSVC)
    message(STATUS "Using MSVC compiler")

    # #define _USE_MATH_DEFINES for M_PI and M_E
    # #define _CRT_SECURE_NO_DEPRECATE to surpress fopen warnings (fopen_s is not available on other platforms)
    # #define _MSVC to alter a regular expression in the Parse class (since MSVC treats the $ command differently)
    add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_DEPRECATE -D_MSVC)

    # warning flag: /W4
    # there is no alternative for -march=native in MSVC, therefore important relevant compiler flags have to be set manually for now (avx, avx2, fma etc.).
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG")

else() # clang or g++
    # warning flags: -Wall -Wextra -Werror
    message(STATUS "Using g++ or Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -march=native -DNDEBUG")

    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        # specific to clang
    else()
        # extra for g++ (and compilers different to msvc and clang)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    endif()
endif()

# --- loki-b --- #
set(src loki-b/source)
set(h   loki-b/headers)
set(pf  loki-b/property_functions)
set(inc include)

# Specify the libraries to link for different backends.
set(ob  openblas)
set(mkl mkl_intel_lp64 mkl_core)

# MSVC and Clang require to additionally link to the intel openmp library.
# g++ can also use intel openmp threading, however GNU openmp seemed to perform better.
# These settings might not yet work for Mac OS, as it is still untested.

# [TODO] Test compiling with MKL using Clang and alter libraries where necessary.
if(MSVC)
  set(mkl ${mkl} mkl_intel_thread libiomp5md)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(mkl ${mkl} mkl_intel_thread iomp5)
else() # g++
  set(mkl ${mkl} mkl_gnu_thread) #mkl_intel_thread iomp5 
endif()

# Specify the linear algebra backend to use (default is pure Eigen).
option(USE_MKL "Use Intel MKL" OFF)
option(USE_OPENBLAS "Use OpenBLAS" OFF)

if (USE_MKL AND USE_OPENBLAS)
  message(SEND_ERROR "Multiple backends defined, please just select one.")
endif()

if(USE_MKL)
  message(STATUS "Using Intel MKL")
  add_definitions(-DEIGEN_USE_MKL_ALL)
  set(backend ${mkl})
  unset(USE_MKL CACHE)
elseif(USE_OPENBLAS)
  message(STATUS "Using OpenBLAS")
  add_definitions(-DEIGEN_USE_LAPACKE -DEIGEN_USE_BLAS)
  set(backend ${ob})
  unset(USE_OPENBLAS CACHE)
else()
  message(STATUS "Using Eigen internal expressions")
  set(backend )
endif()

set(cpp ${src}/WorkingConditions.cpp
        ${src}/Grid.cpp ${src}/EedfGas.cpp ${src}/Enumeration.cpp
        ${src}/EedfState.cpp  ${src}/Setup.cpp
        ${src}/Simulation.cpp ${src}/ElectronKinetics.cpp
        ${src}/CrossSection.cpp ${src}/EedfGasMixture.cpp
        ${src}/EedfCollision.cpp ${src}/Power.cpp ${src}/Output.cpp
        ${src}/JobSystem.cpp ${src}/json.cpp)
set(headers ${h}/WorkingConditions.h ${h}/Grid.h
        ${h}/Gas.h ${h}/EedfGas.h ${h}/Constant.h
        ${h}/Collision.h ${h}/State.h ${h}/EedfState.h
        ${h}/Setup.h ${h}/Enumeration.h ${h}/Parse.h
        ${h}/Event.h ${h}/Simulation.h ${h}/ElectronKinetics.h
        ${h}/Log.h ${h}/GasMixture.h ${h}/CrossSection.h ${h}/Traits.h
        ${h}/EedfGasMixture.h ${h}/EedfCollision.h ${h}/InputStructures.h
        ${h}/Power.h ${h}/MacroscopicQuantities.h ${h}/Output.h
        ${h}/StandardPaths.h ${h}/JobSystem.h ${h}/json.h)
set(property_functions ${pf}/PropertyFunctions.h)

# Standard MKL directories on different operating systems. 
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
  set(mkl_lib "/opt/intel/mkl/lib/intel64")
  set(mkl_inc "/opt/intel/mkl/include")
  set(mkl_comp "/opt/intel/lib")
elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
  set(mkl_inc "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2020.0.166/windows/mkl/include")
  set(mkl_lib "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2020.0.166/windows/mkl/lib/intel64_win")
  set(mkl_comp "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2020.0.166/windows/compiler/lib/intel64_win")
elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
  # [TODO] Add Mac OS specific paths
  set(mkl_lib "")
  set(mkl_inc "")
  set(mkl_comp "")
endif()

# -- loki library -- #
add_library(loki_l ${headers} ${cpp} ${property_functions})

target_include_directories(loki_l PRIVATE ${h} ${pf} ${inc} ${mkl_inc})
target_link_directories(loki_l PRIVATE ${mkl_lib} ${mkl_comp})
target_link_libraries(loki_l PRIVATE ${backend})

# -- loki default front end -- #
add_executable(loki ${src}/main.cpp)

target_include_directories(loki PRIVATE ${h} ${pf} ${inc} ${mkl_inc})
target_link_directories(loki PRIVATE build ${mkl_lib} ${mkl_comp})
target_link_libraries(loki PRIVATE loki_l ${backend})

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # to always build the doxygen docs, add option ALL to the following list.
    add_custom_target( doc_doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
