cmake_minimum_required(VERSION 3.14)
project(loki_cpp)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

# --- Compiler flags --- #

if(MSVC)
    message(STATUS "Using MSVC compiler")

    # #define _USE_MATH_DEFINES for M_PI and M_E
    # #define _CRT_SECURE_NO_DEPRECATE to surpress fopen warnings (fopen_s is not available on other platforms)
    # #define _MSVC to alter a regular expression in the Parse class (since MSVC treats the $ command differently)
    add_compile_definitions(_USE_MATH_DEFINES _CRT_SECURE_NO_DEPRECATE _MSVC)

    # warning flag: /W4
    # there is no alternative for -march=native in MSVC, therefore important relevant compiler flags have to be set manually for now (avx, avx2, fma etc.).
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG")

else() # clang or g++

	# warning flags: -Wall -Wextra -Werror
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -march=native -DNDEBUG")

endif()

# --- OpenMP --- #

#Search for OpenMP, if it is found an OpenMP::OpenMP_CXX target is exported.
#This target sets the appropriate compile/link flags (e.g. -fopenmp for GCC,
#-Qiopenmp for Intel,  -fopenmp=libomp for Clang etc.). Most compilers come with
#openmp support, but for some (e.g. Clang) you need to install additional
#packages (e.g. libomp9-devel for clang9, or libomp-dev) 
option(USE_OPENMP "Use OpenMP parallelization" ON)
if(USE_OPENMP)
	find_package(OpenMP)
	if(NOT OpenMP_CXX_FOUND)
		message(FATAL_ERROR "OpenMP could not be found. Install the missing OpenMP package(s), or disable OpenMP support with -DUSE_OPENMP=OFF")
	endif()
endif()

# --- Eigen3 --- #

option(USE_BUILTIN_EIGEN "Use builtin Eigen instead of system-installed" OFF)

if(USE_BUILTIN_EIGEN)
	#Use Eigen from include/eigen-3.3.7/eigen3
	#Do we need this?
	if(NOT TARGET Eigen3::Eigen)
		add_library(Eigen3::Eigen INTERFACE IMPORTED)
		target_include_directories(Eigen3::Eigen INTERFACE ${PROJECT_SOURCE_DIR}/include/eigen-3.3.7/eigen3)
	endif()
else()
	#Use Eigen found on the system, or download from git
	find_package(Eigen3 3.3 NO_MODULE QUIET)
	if(Eigen3_FOUND AND TARGET Eigen3::Eigen)
		#We found Eigen: done
		message(STATUS "System-installed Eigen found")
	else()
		#Eigen not found or bad configured: download source
		message(STATUS "Using FetchContent to download Eigen (avoid this by installing Eigen3 or build with -DUSE_BUILTIN_EIGEN=ON)")
		FetchContent_Declare(
			Eigen
			GIT_REPOSITORY    https://gitlab.com/libeigen/eigen.git
			GIT_TAG           3.4.0
		)
		FetchContent_MakeAvailable(Eigen)
	endif()
endif()


# --- nlohmann/json --- #

option(USE_BUILTIN_NLOHMANN_JSON "Use the version of nlohmann/json that is shipped with loki-b" OFF)

if(USE_BUILTIN_NLOHMANN_JSON)
	#Use nlohmann/json from include/nlohmann-3.7.0
	#Do we need this?
	if(NOT TARGET nlohmann_json::nlohmann_json)
		add_library(nlohmann_json::nlohmann_json INTERFACE IMPORTED)
		target_include_directories(nlohmann_json::nlohmann_json INTERFACE ${PROJECT_SOURCE_DIR}/include/nlohmann-3.7.0)
	endif()
else()
	#Use nlohmann/json found on the system, or download from git
	find_package(nlohmann_json 3.7.0 QUIET)
	if(nlohmann_json_FOUND AND TARGET nlohmann_json::nlohmann_json)
		#We found nlohmann/json: done
		message(STATUS "System-installed nlohmann/json found")
	else()
		#nlohmann/json not found or bad configured: download source (v 3.11.2)
		message(STATUS "Using FetchContent to download nlohmann/json (avoid this by installing nlohman/json or build with -DUSE_BUILTIN_NLOHMANN_JSON=ON)")
		FetchContent_Declare(
			nlohmann_json
			GIT_REPOSITORY    https://github.com/nlohmann/json.git
			GIT_TAG           bc889af
		)
		FetchContent_MakeAvailable(nlohmann_json)
	endif()
endif()

# Specify the linear algebra backend to use (default is pure Eigen).
option(USE_MKL "Use Intel MKL" OFF)
option(USE_OPENBLAS "Use OpenBLAS" OFF)

if (USE_MKL AND USE_OPENBLAS)
  message(SEND_ERROR "Multiple backends defined, please just select one.")
endif()

if(USE_MKL)
  message(STATUS "Using Intel MKL")
  target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_MKL_ALL)
  if(NOT BLA_VENDOR)
	 set(BLA_VENDOR Intel10_64lp)
  endif()
  find_package(BLAS REQUIRED)  
  if(${CMAKE_VERSION} VERSION_LESS "3.18")
 	#This can be simplified when cmake versions older than 3.18 get less widespread
	set(backend_libs  ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES})
  else()
  	set(backend_libs  BLAS::BLAS)
  endif()
elseif(USE_OPENBLAS)
  message(STATUS "Using OpenBLAS")
  if(NOT BLA_VENDOR)
	 set(BLA_VENDOR OpenBLAS)
  endif()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  if(${CMAKE_VERSION} VERSION_LESS "3.18")
 	#This can be simplified when cmake versions older than 3.18 get less widespread
	set(backend_libs  ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} 
					  ${BLAS_LINKER_FLAGS}   ${BLAS_LIBRARIES}  )
  else()
  	set(backend_libs  LAPACK::LAPACK BLAS::BLAS)
  endif()
  target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_LAPACKE EIGEN_USE_BLAS)
else()
  message(STATUS "Using Eigen internal expressions")
endif()


# --- loki-b --- #
set(src ${CMAKE_CURRENT_SOURCE_DIR}/loki-b/source)
set(h   ${CMAKE_CURRENT_SOURCE_DIR}/loki-b)
set(ideas loki-b/ideas)
set(tests loki-b/tests)

# set(headers ${h}/WorkingConditions.h ${h}/Grid.h
#         ${h}/Gas.h ${h}/GasMixture.h ${h}/EedfGas.h ${h}/Constant.h
#         ${h}/Collision.h
#         ${h}/Setup.h ${h}/Enumeration.h ${h}/Parse.h ${h}/Matrix2Picture.h
#         ${h}/Event.h ${h}/Simulation.h ${h}/ElectronKinetics.h
#         ${h}/Log.h ${h}/GasMixture.h ${h}/CrossSection.h
#         ${h}/EedfMixture.h ${h}/EedfCollision.h ${h}/StateEntry.h
#         ${h}/Power.h ${h}/MacroscopicQuantities.h ${h}/Output.h ${h}/PropertyFunctions.h
#         ${h}/StandardPaths.h ${h}/JobSystem.h ${h}/json.h)


# -- loki library -- #
add_subdirectory(lib)

# -- loki front ends -- #
if(NOT EMSCRIPTEN)
	# -- loki native default front end -- #
	add_executable(loki ${src}/main.cpp)
	target_link_libraries(loki PRIVATE loki-b)

	# -- loki tests -- #
	add_executable(test_parser ${tests}/test_parser.cpp)
	target_link_libraries(test_parser PRIVATE loki-b)

	add_executable(test_statebase ${tests}/test_statebase.cpp)
	target_link_libraries(test_statebase PRIVATE loki-b)

	add_executable(test_reaction_format ${tests}/test_reaction_format.cpp)
	target_link_libraries(test_reaction_format PRIVATE loki-b)

	add_executable(test_lut ${tests}/test_lut.cpp)
	target_link_libraries(test_lut PRIVATE loki-b)

	add_executable(bolsig_extract ${tests}/bolsig_extract.cpp)
	target_link_libraries(bolsig_extract PRIVATE loki-b)

	add_executable(compare_tables ${tests}/compare_tables.cpp)
	target_link_libraries(compare_tables PRIVATE loki-b)

	add_executable(test_druyvesteyn ${tests}/test_druyvesteyn.cpp)
	target_link_libraries(test_druyvesteyn PRIVATE loki-b)

	add_executable(test_statecontainer ${ideas}/test_statecontainer.cpp)
	target_link_libraries(test_statecontainer PRIVATE loki-b)

	add_executable(json2lxcat ${ideas}/json2lxcat.cpp)
	target_include_directories(json2lxcat PRIVATE ${inc})
	target_link_directories(json2lxcat PRIVATE build)
	target_link_libraries(json2lxcat PRIVATE nlohmann_json::nlohmann_json)

	enable_testing()
	add_test(lut test_lut)
	add_test(parser test_parser)
	add_test(reaction_format test_reaction_format)
	add_test(statebase test_statebase)
	add_test(statecontainer test_statecontainer)
	add_test(druyvesteyn test_druyvesteyn)
else()
	# -- WebAssembly -- #
	message(STATUS "Compiling to WebAssembly.")
	add_executable(loki_bindings loki-b/web/bindings.cpp)
	target_include_directories(loki_bindings PRIVATE loki-b ${inc} /usr/lib/emscripten/system/include)
	target_link_libraries(loki_bindings PRIVATE loki-b)
	set_target_properties(loki_bindings PROPERTIES LINK_FLAGS "\
		-s ASSERTIONS=1 \
		-s ALLOW_MEMORY_GROWTH=1 \
		-s MODULARIZE=1 \
		-s EXPORT_NAME=create_module \
		-s EXTRA_EXPORTED_RUNTIME_METHODS=[\"FS\"] \
		--bind \
		--preload-file Input/default_lokib_setup.json \
		--preload-file Input/N2_LXCat.json \
		--preload-file Input/N2_rot_LXCat.json \
		--preload-file Input/masses.txt \
		--preload-file Input/harmonicFrequencies.txt \
		--preload-file Input/rotationalConstants.txt \
		--preload-file Input/N2_vibpop.txt \
		--preload-file Input/OPBParameter.txt \
		--preload-file Input/quadrupoleMoment.txt \
		--preload-file Input/anharmonicFrequencies.txt")
endif()

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # to always build the doxygen docs, add option ALL to the following list.
    add_custom_target( doc_doxygen_loki_b
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen needs to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
