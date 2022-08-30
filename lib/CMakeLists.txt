set(cpp ../source/CrossSection.cpp
        ../source/Gas.cpp
        ../source/GasMixture.cpp
        ../source/Gnuplot.cpp
        ../source/Grid.cpp
        ../source/EedfCollisions.cpp
        ../source/EedfMixture.cpp
        ../source/ElectronKinetics.cpp
        ../source/Enumeration.cpp
        ../source/JobSystem.cpp
        ../source/json.cpp
        ../source/LinearAlgebra.cpp
        ../source/LookupTable.cpp
        ../source/Output.cpp
        ../source/Parse.cpp
        ../source/Power.cpp
        ../source/Setup.cpp
        ../source/Simulation.cpp
        ../source/StateEntry.cpp
        ../source/WorkingConditions.cpp)

add_library(loki-b SHARED ${cpp} ${property_functions})
target_include_directories(loki-b PUBLIC loki-b  ${PROJECT_SOURCE_DIR})
target_link_libraries(loki-b PUBLIC ${backend_libs} Eigen3::Eigen nlohmann_json::nlohmann_json)
target_compile_definitions(loki-b PUBLIC LOKIB_BUILDING_LIBLOKIB)
if(USE_OPENMP)
   target_compile_definitions(loki-b PUBLIC USE_OPENMP)
   target_link_libraries(loki-b PUBLIC OpenMP::OpenMP_CXX )
endif()

install(TARGETS loki-b)
