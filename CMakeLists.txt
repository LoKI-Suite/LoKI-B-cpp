cmake_minimum_required(VERSION 3.13)
project(loki_cpp)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -fopenmp -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fopenmp")
set(CMAKE_CXX_STANDARD 17)

# --- loki-b --- #
set(src loki-b/source)
set(h   loki-b/headers)
set(pf  loki-b/property_functions)
set(inc include)
set(dep dependencies)

# Specify the libraries to link for different backends.
set(mkl mkl_intel_lp64 mkl_gnu_thread mkl_core)
set(ob  openblas)

# Specify the linear algebra backend to use.
option(USE_MKL "Use Intel MKL" OFF)
option(USE_OPENBLAS "Use OpenBLAS" OFF)

if (USE_MKL AND USE_OPENBLAS)
  message(SEND_ERROR "Multiple backends defined, please just select one.")
endif()

if(USE_MKL)
  message(STATUS "Using Intel MKL")
  add_definitions(-DEIGEN_USE_MKL_ALL)
  set(backend ${mkl})
  set(inc ${inc} include/MKL)
  unset(USE_MKL CACHE)
elseif(USE_OPENBLAS)
  message(STATUS "Using OpenBLAS")
  add_definitions(-DEIGEN_USE_LAPACKE)
  set(backend ${ob})
  set(inc ${inc} include/OpenBLAS)
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
target_link_directories(loki_l PRIVATE ${dep})
target_link_libraries(loki_l PRIVATE ${backend} stdc++fs)
target_link_options(loki_l PRIVATE -fopenmp)

# -- loki default front end -- #
add_executable(loki ${src}/main.cpp)

target_include_directories(loki PRIVATE ${h} ${pf} ${inc})
target_link_directories(loki PRIVATE build ${dep})
target_link_libraries(loki PRIVATE loki_l ${backend} stdc++fs)
target_link_options(loki PRIVATE -fopenmp)

# --- tests --- #

# -- loki_tests -- #
add_executable(loki_tests tests/loki_tests/main.cpp)
target_include_directories(loki_tests PRIVATE ${inc} ${h})

# -- regex_test -- #
add_executable(regex_test tests/regex_test/main.cpp)

# -- lua_test -- #
add_executable(lua_test tests/lua_test/main.cpp)

target_include_directories(lua_test PRIVATE tests/lua_test/include/)

target_link_directories(lua_test PRIVATE tests/lua_test/dependencies/lua)
target_link_libraries(lua_test liblua53.a dl)

# -- typedef_test -- #
add_executable(typedef_test tests/typedef_test/main.cpp)

# -- eigen_test -- #
add_executable(eigen_test tests/eigen_test/main.cpp tests/eigen_test/backend_abstraction.h loki-b/headers/LinearAlgebra.h)
target_include_directories(eigen_test PRIVATE include)

# -- event_listener -- #
add_executable(event_listener tests/event_listener/main.cpp)

# -- trait_test -- #
add_executable(trait_test tests/trait_test/main.cpp)

# -- trait_test_no_indexing -- #
add_executable(trait_test_no_indexing tests/trait_test_no_indexing/main.cpp tests/trait_test_no_indexing/Traits.h tests/trait_test_no_indexing/Gas.h tests/trait_test_no_indexing/Basic.h tests/trait_test_no_indexing/State.h tests/trait_test_no_indexing/GasMixture.h tests/trait_test_no_indexing/Working_v1.h)
target_include_directories(trait_test_no_indexing PRIVATE ${inc})

# -- lu_test -- #
add_executable(lu_test tests/lu_test/main.cpp)

# -- hess_solver -- #
add_executable(hess_solver tests/hess_solver/main.cpp)

# -- hess_reduction -- #
add_executable(hess_reduction tests/hess_reduction/main.cpp)

# -- reflection_test -- #
add_executable(reflection tests/reflection_test/main.cpp)

# -- tiled_lu_test -- #
add_executable(tiled_lu tests/tiled_lu_test/main.cpp tests/tiled_lu_test/LinearAlgebra.h tests/tiled_lu_test/TiledMatrix.h)

# OpenBLAS
target_include_directories(tiled_lu PRIVATE ${inc})
target_link_directories(tiled_lu PRIVATE ${dep})
target_link_libraries(tiled_lu PRIVATE openblas)

# Intel MKL
# target_include_directories(tiled_lu PRIVATE ${inc} ${mkl}/include)
# target_link_directories(tiled_lu PRIVATE ${mkl}/lib/intel64)
# target_link_libraries(tiled_lu PRIVATE mkl_intel_lp64 mkl_gnu_thread mkl_core)