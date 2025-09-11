add_library(loki-b SHARED)
target_sources(loki-b PRIVATE
   ../source/CrossSection.cpp
   ../source/EedfCollisions.cpp
   ../source/EedfMixture.cpp
   ../source/EedfUtilities.cpp
   ../source/ElectronKinetics.cpp
   ../source/Enumeration.cpp
   ../source/Environment.cpp
   ../source/Gas.cpp
   ../source/GasMixture.cpp
   ../source/GasProperties.cpp
   ../source/Gnuplot.cpp
   ../source/Grid.cpp
   ../source/Iterators.cpp
   ../source/JobSystem.cpp
   ../source/json.cpp
   ../source/LegacyToJSON.cpp
   ../source/LinearAlgebra.cpp
   ../source/LookupTable.cpp
   ../source/OffSideToJSON.cpp
   ../source/Operators.cpp
   ../source/Output.cpp
   ../source/Parse.cpp
   ../source/Power.cpp
   ../source/Range.cpp
   ../source/Simulation.cpp
   ../source/StateEntry.cpp
   ../source/WorkingConditions.cpp
)

target_include_directories(loki-b  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:.>  
)

target_compile_definitions(loki-b PUBLIC LOKIB_BUILDING_LIBLOKIB)

if(LOKIB_USE_BUILTIN_NLOHMANN_JSON)
  target_include_directories(loki-b PUBLIC 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/nlohmann-3.7.0>
  )
else()
	target_link_libraries(loki-b PUBLIC nlohmann_json::nlohmann_json)
endif()

if(LOKIB_USE_BUILTIN_EIGEN)
  target_include_directories(loki-b PUBLIC 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/eigen-3.4.0>
  )
else()
	target_link_libraries(loki-b  PUBLIC Eigen3::Eigen)
endif()

if(LOKIB_USE_MKL)
  target_link_libraries(loki-b PUBLIC MKL::MKL)
  target_compile_definitions(loki-b PUBLIC EIGEN_USE_MKL_ALL)  
endif()

if(LOKIB_USE_OPENBLAS)
  #Once the minimum cmake version is bumped to 3.18 or higher, the following two lines
  #may be replaced with a more prefered:
  #target_link_libaries(loki-b INTERFACE BLAS::BLAS LAPACK::LAPACK)
  target_link_libraries(loki-b PUBLIC ${LAPACK_LINKER_FLAGS} ${BLAS_LIBRARIES}   )
  target_link_options(loki-b  PUBLIC ${LAPACK_LINKER_FLAGS} ${BLAS_LINKER_FLAGS})  
  target_compile_definitions(loki-b  PUBLIC EIGEN_USE_LAPACKE EIGEN_USE_BLAS)
endif()

if(LOKIB_USE_OPENMP)
   #OpenMP is only used in source files, so this dependency is private
   target_compile_definitions(loki-b PRIVATE LOKIB_USE_OPENMP)
   target_link_libraries(loki-b PRIVATE OpenMP::OpenMP_CXX )
endif()

if(LOKIB_USE_HDF5)
        target_link_libraries(loki-b PUBLIC HighFive)
endif()


install(
   TARGETS loki-b
   EXPORT loki_cpp-targets
)
