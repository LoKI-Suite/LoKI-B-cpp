set(inc ${CMAKE_CURRENT_SOURCE_DIR}/../loki-b ${inc_eigen} ${inc_nlohmann})
set(cpp ${src}/CrossSection.cpp
        ${src}/Gas.cpp
        ${src}/GasMixture.cpp
        ${src}/Gnuplot.cpp
        ${src}/Grid.cpp
        ${src}/EedfCollisions.cpp
        ${src}/EedfMixture.cpp
        ${src}/ElectronKinetics.cpp
        ${src}/Enumeration.cpp
        ${src}/JobSystem.cpp
        ${src}/json.cpp
        ${src}/LinearAlgebra.cpp
        ${src}/LookupTable.cpp
        ${src}/Output.cpp
        ${src}/Parse.cpp
        ${src}/Power.cpp
        ${src}/Setup.cpp
        ${src}/Simulation.cpp
        ${src}/StateEntry.cpp
        ${src}/WorkingConditions.cpp)

add_library(loki-b SHARED ${cpp} ${property_functions})
target_include_directories(loki-b PRIVATE loki-b ${inc} ${mkl_inc})
target_link_directories(loki-b PRIVATE ${mkl_lib} ${mkl_comp})
target_link_libraries(loki-b PRIVATE ${backend})

