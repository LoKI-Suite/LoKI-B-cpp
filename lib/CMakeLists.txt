add_library(loki-b SHARED)
target_sources(loki-b PRIVATE
   ../source/CrossSection.cpp
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
   ../source/WorkingConditions.cpp
)

target_include_directories(loki-b  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:.>  
)
target_link_libraries(loki-b PUBLIC ${backend_libs} Eigen3::Eigen nlohmann_json::nlohmann_json)
target_compile_definitions(loki-b PUBLIC LOKIB_BUILDING_LIBLOKIB)
if(USE_OPENMP)
   #OpenMP is only used in source files, so this dependency is private
   target_compile_definitions(loki-b PRIVATE USE_OPENMP)
   target_link_libraries(loki-b PRIVATE OpenMP::OpenMP_CXX )
endif()

install(
        TARGETS loki-b
        EXPORT loki_cpp-targets
)
