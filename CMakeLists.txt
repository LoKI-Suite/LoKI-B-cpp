cmake_minimum_required(VERSION 3.14)
project(loki_cpp)
file(STRINGS "VERSION.txt" PROJECT_VERSION)

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

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  #Put all dlls and executables in a single folder. This
  #enables all regtests and executables to find the dlls.
  #Probably there is a nicer way (a copy step at the end of
  #the compile that copies all dlls next to the executables?)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
else()
  #Set RPATH for linux systems, such that loki-b keeps working
  #when installed.
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
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

# --- OpenBLAS/MKL --- #

# Specify the linear algebra backend to use (default is pure Eigen).
option(USE_MKL "Use Intel MKL" OFF)
option(USE_OPENBLAS "Use (Open)BLAS" OFF)

if (USE_MKL AND USE_OPENBLAS)
  message(SEND_ERROR "Multiple backends defined, please just select one.")
endif()

if(USE_MKL)
	message(STATUS "Using Intel MKL")
	target_compile_definitions(Eigen3::Eigen INTERFACE EIGEN_USE_MKL_ALL)
	#This has only been tested with the Intel OneAPI toolkit, version 2021.3
	#According to documentation, version 2021.3 and higher should correctly 
	#export a MKLConfig.cmake file, which can be picked up by cmake:
	#https://www.intel.com/content/www/us/en/develop/documentation/onemkl-macos-developer-guide/top/getting-started/cmake-config-for-onemkl.html
	find_package(MKL CONFIG REQUIRED)
	message(STATUS "${MKL_IMPORTED_TARGETS}")
	set(backend_libs MKL::MKL)
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

# check if Doxygen is installed
find_package(Doxygen)

enable_testing()

add_subdirectory(lib)
if(EMSCRIPTEN)
    add_subdirectory(web)
else()
    add_subdirectory(app)
    add_subdirectory(tests)
    add_subdirectory(ideas)
endif(EMSCRIPTEN)
add_subdirectory(LoKI-B)
add_subdirectory(Input)
add_subdirectory(doc)


# --- Install --- #

## Installing / Exporting ##
#  Guided by: https://cmake.org/cmake/help/git-master/guide/tutorial/Adding%20Export%20Configuration.html
#  This SHOULD make the package findable by a simple find_package(loki_cpp) call

include(CMakePackageConfigHelpers)

# Export a file that lists the exported targets
install(EXPORT loki_cpp-targets
  FILE loki_cpp-targets.cmake
  DESTINATION lib/cmake/loki_cpp
  NAMESPACE loki_cpp::
)

# generate the config file itself (loaded and run by find_package)
configure_package_config_file(${PROJECT_SOURCE_DIR}/CMake/loki_cpp-config.cmake.in
  "${PROJECT_BINARY_DIR}/cmake/loki_cpp-config.cmake"
  INSTALL_DESTINATION "lib/cmake/loki_cpp"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

#Generate the version file (used by find_package to decide if the version is high enough)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/cmake/loki_cpp-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

#Install the version and config files
#typically to /usr/local/lib/cmake/BOLSIG
install(FILES
  ${PROJECT_BINARY_DIR}/cmake/loki_cpp-config.cmake
  ${PROJECT_BINARY_DIR}/cmake/loki_cpp-config-version.cmake
  DESTINATION lib/cmake/loki_cpp
)

#Make this project available in the build directory as well
export(EXPORT loki_cpp-targets
  FILE "${PROJECT_BINARY_DIR}/loki_cpp-targets.cmake"
  NAMESPACE loki_cpp::
)