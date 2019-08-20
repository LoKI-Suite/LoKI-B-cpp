cmake_minimum_required(VERSION 3.13)
project(loki_cpp)

# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fopenmp -mavx2 -mfma -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -fopenmp -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fopenmp")
set(CMAKE_CXX_STANDARD 17)

# --- loki-b --- #
set(src loki-b/source)
set(h   loki-b/headers)
set(pf  loki-b/property_functions)
set(inc include)
set(dep dependencies)
set(mkl /home/daan/intel/compilers_and_libraries/linux/mkl)

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

target_include_directories(loki_l PRIVATE ${h} ${pf} ${inc} ${mkl}/include)
target_link_directories(loki_l PRIVATE ${dep} ${mkl}/lib/intel64)
target_link_libraries(loki_l PRIVATE mkl_intel_lp64 mkl_gnu_thread mkl_core stdc++fs)

# -- loki default front end -- #
add_executable(loki ${src}/main.cpp)

target_include_directories(loki PRIVATE ${h} ${pf} ${inc} ${mkl}/include)
target_link_directories(loki PRIVATE build ${mkl}/lib/intel64)
target_link_libraries(loki PRIVATE loki_l mkl_intel_lp64 mkl_gnu_thread mkl_core stdc++fs)

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
add_executable(eigen_test tests/eigen_test/main.cpp loki-b/headers/LinearAlgebra.h)
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