cmake_minimum_required(VERSION 3.13)
project(loki_cpp)

set(CMAKE_CXX_STANDARD 17)

# --- Compiler flags --- #

if(MSVC)
    # warning flag: /W4
    # there is no alternative for -march=native in MSVC, therefore important relevant compiler flags have to be set manually for now (avx, avx2, fma etc.).
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
    set(CMAKE_CXX_FLAGS_RELEASE "{CMAKE_CXX_FLAGS_RELEASE} /O3 /DNDEBUG")

else() # clang or g++
    # warning flags: -Wall -Wextra -Werror
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")

    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        # specific to clang
    else()
        # extra for g++ (and compilers different to msvc and clang)
    endif()
endif()

# --- loki-b --- #
set(src loki-b/source)
set(h   loki-b/headers)
set(pf  loki-b/property_functions)
set(inc include)

# Specify the libraries to link for different backends.
set(mkl mkl_intel_lp64 mkl_gnu_thread mkl_core)
set(ob  openblas)

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
  # set(inc ${inc} include/MKL)
  unset(USE_MKL CACHE)
elseif(USE_OPENBLAS)
  message(STATUS "Using OpenBLAS")
  add_definitions(-DEIGEN_USE_LAPACKE -DEIGEN_USE_BLAS)
  set(backend ${ob})
  # set(inc ${inc} include/OpenBLAS)
  unset(USE_OPENBLAS CACHE)
else()
  message(STATUS "Using Eigen internal expressions")
  set(backend )
endif()

set(cpp ${src}/WorkingConditions.cpp
        ${src}/Grid.cpp ${src}/EedfGas.cpp
        ${src}/EedfState.cpp  ${src}/Setup.cpp
        ${src}/Simulation.cpp ${src}/ElectronKinetics.cpp
        ${src}/CrossSection.cpp ${src}/EedfGasMixture.cpp
        ${src}/EedfCollision.cpp ${src}/Power.cpp ${src}/Output.cpp
        ${src}/JobSystem.cpp)
set(headers ${h}/WorkingConditions.h ${h}/Grid.h
        ${h}/Gas.h ${h}/EedfGas.h ${h}/Constant.h
        ${h}/Collision.h ${h}/State.h ${h}/EedfState.h
        ${h}/Setup.h ${h}/Enumeration.h ${h}/Parse.h
        ${h}/Event.h ${h}/Simulation.h ${h}/ElectronKinetics.h
        ${h}/Log.h ${h}/GasMixture.h ${h}/CrossSection.h ${h}/Traits.h
        ${h}/EedfGasMixture.h ${h}/EedfCollision.h ${h}/InputStructures.h
        ${h}/Power.h ${h}/MacroscopicQuantities.h ${h}/Output.h
        ${h}/StandardPaths.h ${h}/JobSystem.h)
set(property_functions ${pf}/PropertyFunctions.h)

# -- loki library -- #
add_library(loki_l ${headers} ${cpp} ${property_functions})

target_include_directories(loki_l PRIVATE ${h} ${pf} ${inc})
target_link_directories(loki_l PRIVATE)
target_link_libraries(loki_l PRIVATE ${backend} stdc++fs)
target_link_options(loki_l PRIVATE -fopenmp)

# -- loki default front end -- #
add_executable(loki ${src}/main.cpp)

target_include_directories(loki PRIVATE ${h} ${pf} ${inc})
target_link_directories(loki PRIVATE build)
target_link_libraries(loki PRIVATE loki_l ${backend} stdc++fs)
target_link_options(loki PRIVATE -fopenmp)