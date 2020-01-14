# LoKI-B C++
Numerical Boltzmann solver for the two-term approximation of the Boltzmann equation for electrons.

The clean branch only contains the files needed to compile the loki-b project (apart from the basic input to get you started). Furthermore, the CMakeLists.txt file has been updated to support compilers different from g++. It should compile out of the box on g++, MSVC and Clang (untested). Additional, compiler specific, flags can always be added in the CMakeLists.txt file.

Note that when compiling loki with a non-default backend, currently MKL or OpenBLAS, the relevant files should be available to the compiler. The standard paths to the MKL in Linux and Windows have been added to the CMakeLists.txt file, you should edit them if the library is found in a different location. Alternatively they can be added to the PATH.

## Compilation instructions:

1. Make sure Git and CMake are installed on your system as well as a suitable C/C++ compiler (i.e. gcc/g++ on Linux, MSVC on Windows and Clang on Mac OS).
2. Open a terminal and cd into luxurious-loki
3. Switch to the 'clean' branch, using: 'git switch clean'.
4. Create the build directory ('mkdir build') and cd into it 

### Linux
5. run: cmake -DCMAKE\_BUILD\_TYPE=Release -D<BACKEND\_FLAG>=ON ..
    - where BACKEND\_FLAG=USE\_MKL/USE\_OPENBLAS, specifying the backend to supply to Eigen
    - this flag can also be omitted to build with pure Eigen
6. run: make -j <NUM\_JOBS>
    - where <NUM\_JOBS> is the maximum number of jobs to run simultaneously when compiling; just use the number of physical cores in your system. Omit this flag to use the default settings.
    
### Windows
5. run: cmake -D<BACKEND\_FLAG>=ON ..
    - where <BACKEND\_FLAG>=USE\_MKL/USE\_OPENBLAS, specifying the backend to supply to Eigen
    - this flag can also be omitted to build with pure Eigen
6. run: cmake --build . --config Release -j <NUM\_JOBS>
    - where <NUM\_JOBS> is the maximum number of jobs to run simultaneously when compiling; just use the number of physical cores in your system. Omit this flag to use the default settings.
7. The loki exec executable should always run from the build directory (or any directory that is a direct subdirectory of the main folder). When compiling with MSVC the executable can be placed in the 'build\\Release directory' by default, therefore it should be moved to 'build' before executing. Additionally, when compiling on Windows with the MKL, copy the 'libiomp5md.dll' file from the main folder or, preferably from its default location in your MKL installation (see ${mkl\_comp} in CMakeLists.txt), to the build directory.

## Running and Plotting

As of yet, LoKI-B C++ is run from the command line. To plot the data you can install gnuplot. The execution on different operating systems is very similar.

### Linux
1. Open a terminal and navigate to the 'luxurious-loki/build' directory.
2. Run the following command: './loki <INPUT\_FILE> | gnuplot --persist'.
    - Where <INPUT\_FILE> is the full input file name to use, including the '.in' extension, e.g. 'default_lokib_setup.in'. The input files need to be present in the 'luxurious-loki/Input' folder.
    - The output from LoKI-B is then piped into gnuplot, where the '--persist' flag avoids gnuplot from immediately closing after plotting.

### Windows
1. Open a terminal and navigate to the 'luxurious-loki\\build' directory.
2. Run the following command: 'loki.exe <INPUT\_FILE> | gnuplot --persist'.
    - Where <INPUT\_FILE> is the full input file name to use, including the '.in' extension, e.g. 'default_lokib_setup.in'. The input files need to be present in the 'luxurious-loki\\Input' folder.
    - The output from LoKI-B is then piped into gnuplot, where the '--persist' flag avoids gnuplot from immediately closing after plotting.
