# luxurious-loki
Numerical Boltzmann solver for the two-term approximation of the Boltzmann equation for electrons.

The clean branch only contains the files needed to compile the loki-b project (apart from the basic input to get you started). Furthermore, the CMakeLists.txt file has been updated to support compilers different from g++. It should compile out of the box on g++, MSVC and clang. Additional, compiler specific, flags can always be added in the CMakeLists.txt file.

To compile using cmake:

1. Open a terminal and cd into luxurious-loki
2. Create the build directory ('mkdir build') and cd into it 
3. run: cmake -D<BACKEND\_FLAG>=ON ../
    - where BACKEND\_FLAG=USE\_MKL/USEi\_OPENBLAS, specifying the backend to supply to Eigen
    - this flag can also be omitted to build with pure Eigen
4. run: cmake --build . --target loki --config Release -- -j <NUM\_JOBS>
    - where <NUM\_JOBS> is the number of jobs to run simultaneously when compiling; just use the number of physical cores in your system.

Note that when compiling loki with a non-default backend, currently MKL or OpenBLAS, the relevant files should be available to the compiler (i.e. added to the PATH). Alternatively the CMakeLists.txt file can be edited for the files to be found in a custom location.
