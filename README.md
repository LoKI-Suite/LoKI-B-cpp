# luxurious-loki
Numerical Boltzmann solver for the two-term approximation of the Boltzmann equation for electrons.

To compile using cmake:

1. Open a terminal and cd into luxurious-loki/build
2. run: cmake -D<BACKEND_FLAG>=ON ../
    - where BACKEND_FLAG=USE_MKL/USE_OPENBLAS, specifying the backend to supply to Eigen
    - this flag can also be omitted to build with pure Eigen
3. run: cmake --build . --target loki --config Release -- -j <NUM_JOBS>
    - where <NUM_JOBS> is the number of jobs to run simultaneously when compiling; just use the number of physical cores in your system.
    
NOTE that by default loki will probably only compile under linux when supplying a backend (MKL/OpenBLAS). When compiling with a backend on different operating systems, make sure the required libraries are either installed in the default location on your system or present in the "dependencies" folder.
