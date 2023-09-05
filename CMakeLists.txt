cmake_minimum_required(VERSION 3.14)
project(loki_cpp)
file(STRINGS "VERSION.txt" PROJECT_VERSION)

#Some Utility macros and functions
include(CMake/LoKIUtilities.cmake)

set(CMAKE_CXX_STANDARD 17)

#Add protection against in-source builds
include(CMake/PreventInSourceBuilds.cmake)

#Enable fetching packages from git etc
include(FetchContent)

# --- Compiler flags --- #

# Adopt cmake policy 0092:
# do not ad /W3 by default for MSVC (and Intel on Windows?)
# since that gets overridden by /W4, which causes a warning
# for every file being compiled.
cmake_policy(SET CMP0092 NEW)

option(ALL_WARNINGS "-DALL_WARNINGS=ON|OFF : Show all warnings" OFF)
loki_constrain_configure_option(ALL_WARNINGS ON OFF)
if (MSVC)
    message(STATUS "Using MSVC compiler")

    # #define _USE_MATH_DEFINES for M_PI and M_E
    # #define _CRT_SECURE_NO_DEPRECATE to surpress fopen warnings (fopen_s is not available on other platforms)
    # #define _MSVC to alter a regular expression in the Parse class (since MSVC treats the $ command differently)
    add_compile_definitions(_USE_MATH_DEFINES _CRT_SECURE_NO_DEPRECATE _MSVC)

    # warning level 4
    loki_msvc_cxxflags(/W4)
    if(NOT ALL_WARNINGS)
        #disable some valid, but perhaps not so interesting warnings
        #disable warning C4100: "unreferenced formal parameter"
        loki_msvc_cxxflags(/wd4100)
        #disable warning C4244: "conversion from 'X' to 'Y', possible loss of data"
        loki_msvc_cxxflags(/wd4244)
        #disable warning C4267: "'return': conversion from 'size_t' to 'unsigned int', possible loss of data"
        loki_msvc_cxxflags(/wd4267)
        #disable warning C4996: 'XXX': This function or variable may be unsafe.
        loki_msvc_cxxflags(/wd4996)
        #disable warning C4389: '==': signed/unsigned mismatch
        loki_msvc_cxxflags(/wd4389)
        #disable warning C4127: conditional expression is constant
        loki_msvc_cxxflags(/wd4127)
        #disable warning C4459: declaration of 'attr' hides global declaration
        loki_msvc_cxxflags(/wd4459)
    endif()
else()
    # lots of warnings
    loki_ccxxflags(-Wall)
    if(ALL_WARNINGS)
        loki_ccxxflags(-Wextra)
    else()
        #disable some valid, but perhaps not so interesting warnings
        loki_ccxxflags(-Wno-unused-parameter)
        #loki_ccxxflags(-Wno-unused-parameter -Werror=maybe-uninitialized)
    endif()

    #Turn some warnings into errors
    loki_cxxflags(-Wall)
    # disabled because of sloppiness in Eigen-3.4.0 and SuperLU:
    # loki_cxxflags(-Werror=unused-but-set-parameter)
    loki_cxxflags(-Werror=return-type)
    loki_cxxflags(-Werror=ignored-attributes)
    loki_cxxflags(-Werror=overloaded-virtual)
    loki_cxxflags(-Werror=non-virtual-dtor)
    # Disabled because of Eigen/SuperLU
    # loki_cxxflags(-Werror=unused-variable)
    loki_cxxflags(-Werror=unused-result)
    loki_cxxflags(-Werror=ignored-qualifiers)
    loki_cxxflags(-Werror=zero-as-null-pointer-constant)
    loki_cxxflags(-Werror=misleading-indentation)
    loki_cxxflags(-Werror=reorder)
    # disabled, not supported by gcc 7.5.0:
    # loki_cxxflags(-Werror=cast-function-type)
    loki_cxxflags(-Werror=type-limits)
    loki_cxxflags(-Werror=missing-field-initializers)
    loki_cxxflags(-Werror=vla)
    loki_cxxflags(-Werror=invalid-pch)
    # TODO change into -Werror=shadow after fixing numerous local/member data
    # shadowing violations in the code. This should also be made compatible
    # with clang++, which has a slightly different syntax:
    #loki_cxxflags(-Werror=shadow=local)
    # We treat APPLE as a subset of non MSVC
    if(APPLE)
        # Disabled because boost generates a lot of these warnings
        # loki_cxxflags(-Werror=undef)
        add_link_options(LINKER:-undefined,error)
    else()
        loki_cxxflags(-Werror=undef)
        add_link_options(LINKER:-no-undefined)
    endif()
endif()


#CMake provides some default build types:
#                 |  build optimized?  |   debug info?    |   asserts?   |
# ========================================================================
#  Debug:         |        No          |      yes         |      yes     |
#  Release:       |   yes, for speed   |      no          |      no      |
#  RelWithDebinfo |   yes, for speed   |      yes         |      no      |
#  MinSizeRel     |   yes, for size    |      no          |      no      |


#In single-configuration generators, the build type can be set with
#   -D CMAKE_BUILD_TYPE = <build type>
#For multi-configuration generators, CMAKE_BUILD_TYPE is ignored, and you chose the
#configuration during compile time:
#    cmake --build build --config Debug

#The exact flags depend ofcourse on the compiler, but the defaults are somewhat
#sensible. Sensible enough at least, that these replace the autotools flags
# --enable-assert, --enable-debug and --enable-optimize

#For example with gcc:
# CMAKE_CXX_FLAGS_DEBUG is -g
# CMAKE_CXX_FLAGS_RELEASE is -O3 -DNDEBUG
# CMAKE_CXX_FLAGS_RELWITHDEBINFO is -O2 -g -DNDEBUG
# CMAKE_CXX_FLAGS_MINSIZEREL is -Os -DNDEBUG

#For MSVC:
#CMAKE_CXX_FLAGS_DEBUG is /Zi /Ob0 /Od /RTC1
#CMAKE_CXX_FLAGS_RELEASE is /O2 /Ob2 /DNDEBUG
#CMAKE_CXX_FLAGS_RELWITHDEBINFO is /Zi /O2 /Ob1 /DNDEBUG
#CMAKE_CXX_FLAGS_MINSIZEREL is /O1 /Ob1 /DNDEBUG

###############################################################################
###### build configuration ####################################################
###############################################################################

# Set the default build type to Release if none was specified, only for
# single-configuration generators
# https://www.kitware.com/cmake-and-the-default-build-type/
set(default_build_type "Release")
if(NOT CMAKE_CONFIGURATION_TYPES)
 #On single configuration generators CMAKE_CONFIGURATION_TYPES will be empty
 if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
  endif()
  # Set the possible values of build type for cmake-gui
  loki_constrain_configure_option(CMAKE_BUILD_TYPE Release Debug RelWithDebInfo MinSizeRel)
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
loki_constrain_configure_option(USE_OPENMP ON OFF)
if(USE_OPENMP)
	find_package(OpenMP)
	if(NOT OpenMP_CXX_FOUND)
		message(FATAL_ERROR "OpenMP could not be found. Install the missing OpenMP package(s), or disable OpenMP support with -DUSE_OPENMP=OFF")
	endif()
endif()


# --- nlohmann/json --- #

option(USE_BUILTIN_NLOHMANN_JSON "Use the version of nlohmann/json that is shipped with loki-b" OFF)

if(USE_BUILTIN_NLOHMANN_JSON)
	#Use built-in nlohmann/json from include/nlohmann-3.7.0
	message(STATUS "Using built-in nlohmann_json")
else()
	#Use nlohmann/json found on the system, or download from git
	message(STATUS "Using system-installed nlohmann/json.")
	find_package(nlohmann_json 3.7.0 REQUIRED)
endif()


# --- Eigen3 --- #

option(USE_BUILTIN_EIGEN "Use builtin Eigen instead of system-installed" OFF)

if(USE_BUILTIN_EIGEN)
	#Use built-in Eigen from include/eigen-3.3.7/eigen3
	message(STATUS "Using built-in Eigen")
else()
	#Use Eigen found on the system, or download from git
	message(STATUS "Using System-installed Eigen")
	find_package(Eigen3 3.3 NO_MODULE REQUIRED)
endif()

# --- OpenBLAS/MKL --- #

# Specify the linear algebra backend to use (default is pure Eigen).
option(USE_MKL "Use Intel MKL as backend for Eigen" OFF)
option(USE_OPENBLAS "Use (Open)BLAS as backend for Eigen" OFF)

if (USE_MKL AND USE_OPENBLAS)
  message(SEND_ERROR "Multiple backends defined, please just select one.")
endif()

if(USE_MKL)
  message(STATUS "Using Intel MKL")
  #This has only been tested with the Intel OneAPI toolkit, version 2021.3
  #According to documentation, version 2021.3 and higher should correctly 
  #export a MKLConfig.cmake file, which can be picked up by cmake:
  #https://www.intel.com/content/www/us/en/develop/documentation/onemkl-macos-developer-guide/top/getting-started/cmake-config-for-onemkl.html
  find_package(MKL CONFIG REQUIRED)
  message(STATUS "${MKL_IMPORTED_TARGETS}")
elseif(USE_OPENBLAS)
  message(STATUS "Using (Open)BLAS")
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
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
add_subdirectory(input)
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
