cmake_minimum_required(VERSION 3.13)
project(loki_cpp)

set(CMAKE_CXX_STANDARD 17)

# --- Compiler flags --- #

if(MSVC)
    message(STATUS "Using MSVC compiler")

    # #define _USE_MATH_DEFINES for M_PI and M_E
    # #define _CRT_SECURE_NO_DEPRECATE to surpress fopen warnings (fopen_s is not available on other platforms)
    # #define _MSVC to alter a regular expression in the Parse class (since MSVC treats the $ command differently)
    add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_DEPRECATE -D_MSVC)

    # warning flag: /W4
    # there is no alternative for -march=native in MSVC, therefore important relevant compiler flags have to be set manually for now (avx, avx2, fma etc.).
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG")

else() # clang or g++
    # warning flags: -Wall -Wextra -Werror
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -march=native -DNDEBUG")

    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        # specific to clang
				message(STATUS "Using Clang")
    else()
        # extra for g++ (and compilers different to msvc and clang)
    		message(STATUS "Using g++")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    endif()
endif()

# the following two must be enabled to make the versions of Eigen / nlohmann
# that are shipped with LoKI-B available. When disabled, eigen3/Eigen/...
# and nlohmann/json.hpp will be resolved using the default include paths.

#todo see how/when to set these two. For now, we rely on a 'system' eigen3/,
#except for builds with MSVC: then we use the built-in Eigen.
if(MSVC OR EMSCRIPTEN)
  set(inc_eigen include/eigen-3.3.7/eigen3)
else()
  if ("x${EIGEN_PATH}" STREQUAL "x")
    set (EIGEN_PATH "/usr/include/eigen3")
  endif()
  set(inc_eigen ${EIGEN_PATH})
endif()
message(STATUS "Using Eigen library at ${inc_eigen}.")

option(USE_BUILTIN_NLOHMANN_JSON "Use the version of nlohmann/json that is shipped with loki-b" OFF)
if (USE_BUILTIN_NLOHMANN_JSON)
  message("Using the nlohmann/json distribution that is shipped with LoKI-B")
  get_filename_component(NLOHMANN_ABSOLUTE_PATH include/nlohmann-3.7.0 ABSOLUTE)
  set(inc_nlohmann ${NLOHMANN_ABSOLUTE_PATH})
else()
  message("Using the system's nlohmann/json distribution")
  set(inc_nlohmann )
endif()

# --- loki-b --- #
set(src ${CMAKE_CURRENT_SOURCE_DIR}/loki-b/source)
set(h   ${CMAKE_CURRENT_SOURCE_DIR}/loki-b)
set(inc ${inc_eigen} ${inc_nlohmann})
set(ideas loki-b/ideas)
set(tests loki-b/tests)

# Specify the libraries to link for different backends.
set(ob  openblas)
set(mkl mkl_intel_lp64 mkl_core)

# MSVC and Clang require to additionally link to the intel openmp library.
# g++ can also use intel openmp threading, however GNU openmp seemed to perform better.
# These settings might not yet work for Mac OS, as it is still untested.

# [TODO] Test compiling with MKL using Clang and alter libraries where necessary.
if(MSVC)
  set(mkl ${mkl} mkl_intel_thread libiomp5md)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(mkl ${mkl} mkl_intel_thread iomp5)
else() # g++
  set(mkl ${mkl} mkl_gnu_thread) #mkl_intel_thread iomp5
endif()

# Specify the linear algebra backend to use (default is pure Eigen).
option(USE_MKL "Use Intel MKL" OFF)
option(USE_OPENBLAS "Use OpenBLAS" OFF)

if (USE_MKL AND USE_OPENBLAS)
  message(SEND_ERROR "Multiple backends defined, please just select one.")
endif()

if(USE_MKL)
  message(STATUS "Using Intel MKL")
  add_definitions(-DEIGEN_USE_MKL_ALL)
  set(backend ${mkl})
  unset(USE_MKL CACHE)
elseif(USE_OPENBLAS)
  message(STATUS "Using OpenBLAS")
  add_definitions(-DEIGEN_USE_LAPACKE -DEIGEN_USE_BLAS)
  set(backend ${ob})
  unset(USE_OPENBLAS CACHE)
else()
  message(STATUS "Using Eigen internal expressions")
  set(backend )
endif()

# set(headers ${h}/WorkingConditions.h ${h}/Grid.h
#         ${h}/Gas.h ${h}/GasMixture.h ${h}/EedfGas.h ${h}/Constant.h
#         ${h}/Collision.h
#         ${h}/Setup.h ${h}/Enumeration.h ${h}/Parse.h ${h}/Matrix2Picture.h
#         ${h}/Event.h ${h}/Simulation.h ${h}/ElectronKinetics.h
#         ${h}/Log.h ${h}/GasMixture.h ${h}/CrossSection.h
#         ${h}/EedfMixture.h ${h}/EedfCollision.h ${h}/StateEntry.h
#         ${h}/Power.h ${h}/MacroscopicQuantities.h ${h}/Output.h ${h}/PropertyFunctions.h
#         ${h}/StandardPaths.h ${h}/JobSystem.h ${h}/json.h)

# Standard MKL directories on different operating systems.
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
  set(mkl_lib "/opt/intel/mkl/lib/intel64")
  set(mkl_inc "/opt/intel/mkl/include")
  set(mkl_comp "/opt/intel/lib")
elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
  set(mkl_inc "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2020.0.166/windows/mkl/include")
  set(mkl_lib "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2020.0.166/windows/mkl/lib/intel64_win")
  set(mkl_comp "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2020.0.166/windows/compiler/lib/intel64_win")
elseif("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
  # [TODO] Add Mac OS specific paths
  set(mkl_lib "")
  set(mkl_inc "")
  set(mkl_comp "")
endif()

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
	target_link_libraries(json2lxcat)

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
