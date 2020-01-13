# luxurious-loki
Numerical Boltzmann solver for the two-term approximation of the Boltzmann equation for electrons.

The clean branch only contains the files needed to compile the loki-b project (apart from the basic input to get you started). Furthermore, the CMakeLists.txt file has been updated to support compilers different from g++. It should compile out of the box on GCC, MSVC and Clang (untested). Additional, compiler specific, flags can always be added in the CMakeLists.txt file.

To compile using cmake:

1. Open a terminal and cd into luxurious-loki
2. Create the build directory ('mkdir build') and cd into it 
3. run: cmake -DCMAKE_BUILD_TYPE=Release -D<BACKEND\_FLAG>=ON ../
    - where BACKEND\_FLAG=USE\_MKL/USE\_OPENBLAS, specifying the backend to supply to Eigen
    - this flag can also be omitted to build with pure Eigen
4. run: make -j <NUM\_JOBS>
    - where <NUM\_JOBS> is the number of jobs to run simultaneously when compiling; just use the number of physical cores in your system.

Note that when compiling loki with a non-default backend, currently MKL or OpenBLAS, the relevant files should be available to the compiler. The standard paths to the MKL in Linux and Windows have been added to the CMakeLists.txt file, you should edit them if the library is found in a different location. Alternatively they can be added to the PATH.
