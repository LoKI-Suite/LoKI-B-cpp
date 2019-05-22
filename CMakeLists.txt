cmake_minimum_required(VERSION 3.14)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
project(loki_cpp)

set(CMAKE_CXX_STANDARD 17)

# --- loki-b --- #
set(src loki-b/source)
set(h   loki-b/headers)
set(inc include)
set(dep dependencies)

set(cpp ${src}/WorkingConditions.cpp
        ${src}/Grid.cpp ${src}/EedfGas.cpp
        ${src}/EedfState.cpp  ${src}/Setup.cpp
        ${src}/Simulation.cpp ${src}/ElectronKinetics.cpp
        ${src}/CrossSection.cpp ${src}/EedfGasMixture.cpp
        ${src}/EedfCollision.cpp)
set(headers ${h}/WorkingConditions.h ${h}/Grid.h
        ${h}/Gas.h ${h}/EedfGas.h ${h}/Constant.h
        ${h}/Collision.h ${h}/State.h ${h}/EedfState.h
        ${h}/Setup.h ${h}/Enumeration.h ${h}/Parse.h
        ${h}/Event.h ${h}/Simulation.h ${h}/ElectronKinetics.h
        ${h}/Log.h ${h}/GasMixture.h ${h}/CrossSection.h ${h}/Traits.h
        ${h}/EedfGasMixture.h ${h}/EedfCollision.h)

add_executable(loki ${headers} ${cpp} ${src}/main.cpp)

target_include_directories(loki PRIVATE ${h} ${inc})
target_link_directories(loki PRIVATE ${dep})

# target_compile_options(loki PRIVATE -O3)

# -- library version -- #
add_library(loki_l ${headers} ${cpp})

target_include_directories(loki_l PRIVATE ${h} ${inc})
target_link_directories(loki_l PRIVATE ${dep})

# set_target_properties(loki_l PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
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