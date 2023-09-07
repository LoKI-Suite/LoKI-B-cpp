# LoKI-B C++

[![action](https://github.com/DAANBOER/luxurious-loki/workflows/CMake/badge.svg)](../../actions?query=workflow%3A%22CMake%22)
[![codecov](https://codecov.io/gh/LoKI-Suite/LoKI-B/graph/badge.svg?token=GRKRVNPAAZ)](https://codecov.io/gh/LoKI-Suite/LoKI-B)

Numerical Boltzmann solver for the two-term approximation of the Boltzmann equation for electrons.

The clean branch only contains the files needed to compile the loki-b project (apart from the basic input to get you started). Furthermore, the CMakeLists.txt file has been updated to support compilers different from g++. It should compile out of the box on g++, MSVC and Clang (untested). Additional, compiler specific, flags can always be added in the CMakeLists.txt file.

Note that when compiling loki with a non-default backend, currently MKL or OpenBLAS, the relevant files should be available to the compiler. The standard paths to the MKL in Linux and Windows have been added to the CMakeLists.txt file, you should edit them if the library is found in a different location. Alternatively they can be added to the PATH.

## Clang format

A Clang format style file is present at the root of this repository (`.clang-format`). To format all the source files, you can issue the following command from the repository root.

 - `clang-format -i -style=file source/*.cpp LoKI-B/*.h app/*.cpp tests/*.cpp ideas/*.cpp web/bindings.cpp`

## Compilation instructions:

1. Make sure Git and CMake are installed on your system as well as a suitable C/C++ compiler (i.e. gcc/g++ on Linux, MSVC on Windows and Clang on Mac OS).
1. Open a terminal and navigate to the folder where you want to save the LoKI-B C++ project folder.
1. Clone the `luxurious-loki` repository: `git clone git@github.com:DAANBOER/luxurious-loki.git` and `cd` into the newly created `luxurious-loki` folder. The download might take a while, since the master branch contains a lot of extra data.
1. Switch to the `clean` branch, using: `git switch clean`.
1. Create the build directory (`mkdir build`) and `cd` into it.

### Linux
1. run: `cmake -DCMAKE_BUILD_TYPE=Release -D<BACKEND_FLAG>=ON ..`
    - where `<BACKEND_FLAG>=USE_MKL/USE_OPENBLAS`, specifying the backend to supply to Eigen
    - this flag can also be omitted to build with pure Eigen
    - LoKI-B assumes that Eigen is available in the directory /usr/include/eigen3. To specify
      another path, run cmake with an additional option like -DEIGEN_PATH=/opt/include/eigen3
1. run: `make -j <NUM_JOBS>`
    - where `<NUM_JOBS>` is the maximum number of jobs to run simultaneously when compiling; just use the number of physical cores in your system. Omit this flag to use the default settings.
    
### Windows
1. run: `cmake -D<BACKEND_FLAG>=ON ..`
    - where `<BACKEND_FLAG>=USE_MKL/USE_OPENBLAS`, specifying the backend to supply to Eigen
    - this flag can also be omitted to build with pure Eigen
1. run: `cmake --build . --config Release -j <NUM_JOBS>`
    - where `<NUM_JOBS>` is the maximum number of jobs to run simultaneously when compiling; just use the number of physical cores in your system. Omit this flag to use the default settings.
1. The loki exec executable should always run from the build directory (or any directory that is a direct subdirectory of the main folder). When compiling with MSVC the executable can be placed in the `build\Release` directory by default, therefore it should be moved to `build` before executing. Additionally, when compiling on Windows with the MKL, copy the `libiomp5md.dll` file from the main folder or, preferably from its default location in your MKL installation (see `${mkl_comp}` in CMakeLists.txt), to the build directory.

## Running and Plotting

As of yet, LoKI-B C++ is run from the command line. To plot the computed eedf you can install gnuplot. The execution on different operating systems is very similar.

### Linux
1. Open a terminal and navigate to the `luxurious-loki/build directory.
1. Run the following command: `./loki <INPUT_FILE> | gnuplot --persist`.
    - Where `<INPUT_FILE>` is the full input file name to use, including the `.in` extension, e.g. `default_lokib_setup.in`. The input files need to be present in the `loki/input` folder.
    - The output from LoKI-B is then piped into gnuplot, where the `--persist` flag avoids gnuplot from immediately closing after plotting.

### Windows
1. Open a terminal and navigate to the `luxurious-loki\build` directory.
1. Run the following command: `loki.exe <INPUT_FILE> | gnuplot --persist`.
    - Where `<INPUT_FILE>` is the full input file name to use, including the `.in` extension, e.g. `default_lokib_setup.in`. The input files need to be present in the `loki\input` folder.
    - The output from LoKI-B is then piped into gnuplot, where the `--persist` flag avoids gnuplot from immediately closing after plotting.
