cmake_minimum_required(VERSION 3.14)
project(loki_cpp)

set(CMAKE_CXX_STANDARD 14)

# --- loki-b --- #
set(src loki-b/source)
set(h   loki-b/headers)
set(inc include)
set(dep dependencies)

set(cpp ${src}/main.cpp ${src}/WorkingConditions.cpp
        ${src}/Grid.cpp ${src}/Gas.cpp ${src}/EedfGas.cpp
        ${src}/Collision.cpp ${src}/State.cpp ${src}/EedfState.cpp
        ${src}/Setup.cpp ${src}/Simulation.cpp ${src}/Boltzmann.cpp)
set(headers ${h}/WorkingConditions.h ${h}/Grid.h
        ${h}/Gas.h ${h}/EedfGas.h ${h}/Constant.h
        ${h}/Collision.h ${h}/State.h ${h}/EedfState.h
        ${h}/Setup.h ${h}/Enumeration.h ${h}/Parse.h
        ${h}/Event.h ${h}/Simulation.h ${h}/Boltzmann.h)

add_executable(loki ${headers} ${cpp})

target_include_directories(loki PRIVATE ${h} ${inc})
target_link_directories(loki PRIVATE ${dep})


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
add_executable(eigen_test tests/eigen_test/main.cpp)
target_include_directories(eigen_test PRIVATE include)

# -- event_listener -- #
add_executable(event_listener tests/event_listener/main.cpp)