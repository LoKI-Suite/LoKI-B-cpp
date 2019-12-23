cmake_minimum_required(VERSION 3.13)
project(loki_cpp)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -fopenmp -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fopenmp")
set(CMAKE_CXX_STANDARD 17)

# --- loki-b --- #
set(src ../loki-b/source)
set(h   ../loki-b/headers)
set(pf  ../loki-b/property_functions)
set(inc ../include)
set(dep ../dependencies)

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
target_link_libraries(loki_l PRIVATE ${backend})# stdc++fs)
target_link_options(loki_l PRIVATE -fopenmp)

# -- loki default front end -- #
add_executable(loki main.cpp)#${src}/main.cpp)

target_include_directories(loki PRIVATE ${h} ${pf} ${inc})
target_link_directories(loki PRIVATE build ${dep})
target_link_libraries(loki PRIVATE loki_l ${backend})# stdc++fs)
target_link_options(loki PRIVATE -fopenmp)
set_target_properties(loki PROPERTIES LINK_FLAGS "-s EXPORTED_FUNCTIONS='[\"_run\", \"_main\"]' -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"cwrap\", \"lengthBytesUTF8\", \"stringToUTF8\"]' -s DISABLE_EXCEPTION_CATCHING=0 -s ALLOW_MEMORY_GROWTH=1")
