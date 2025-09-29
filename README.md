# LoKI-B++

[![action](https://github.com/DAANBOER/luxurious-loki/workflows/CMake/badge.svg)](../../actions?query=workflow%3A%22CMake%22)
[![codecov](https://codecov.io/gh/LoKI-Suite/LoKI-B/graph/badge.svg?token=GRKRVNPAAZ)](https://codecov.io/gh/LoKI-Suite/LoKI-B)

LoKI-B++ is a numerical Boltzmann solver for the two-term approximation of the
Boltzmann equation for electrons. It is based on the equivalent matlab version
[LoKI-B MATLAB](https://github.com/LoKI-Suite/LoKI-B).

## Clang format

A Clang format style file is present at the root of this repository
(`.clang-format`). To format all the source files, you can issue the following
command from the repository root.

- `clang-format -i -style=file source/*.cpp LoKI-B/*.h app/*.cpp tests/*.cpp ideas/*.cpp web/bindings.cpp`

## Compilation instructions:

1. Make sure [Git](https://git-scm.com/) and [CMake](https://cmake.org/) are
   installed on your system as well as a suitable C/C++ compiler (i.e.
   [gcc/g++](https://gcc.gnu.org/) on Linux,
   [MSVC](https://visualstudio.microsoft.com/vs/features/cplusplus/) on Windows
   and [Apple Clang](https://developer.apple.com/xcode/cpp/) on Mac OS).
1. LoKI-B relies on two dependencies, [Eigen](https://eigen.tuxfamily.org/) for
   linear algebra, and [nlohmann-json](https://github.com/nlohmann/json) for
   JSON handling.
1. Open a terminal and navigate to the folder where you want to save the
   LoKI-B++ project folder.
1. Clone the `LoKI-B-cpp` repository:
   `git clone git@github.com:LoKI-Suite/LoKI-B-cpp.git` and `cd` into the newly
   created `LoKI-B-cpp` folder.
1. Create the build directory (`mkdir build`).

### Linux

1. run: `cmake -DCMAKE_BUILD_TYPE=Release -D<BACKEND_FLAG>=ON -B build`
   - Where `<BACKEND_FLAG>=LOKIB_USE_MKL/LOKIB_USE_OPENBLAS`, is a flag to set
     the backend to supply to Eigen,
     [OpenBLAS](http://www.openmathlib.org/OpenBLAS/) or
     [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html).
     This flag can also be omitted to build with pure Eigen.
1. run: `cmake --build build -j <NUM_JOBS>`
   - Where `<NUM_JOBS>` is the maximum number of jobs to run simultaneously when
     compiling; just use the number of physical cores in your system. Omit this
     flag to use the default settings.

### Linux using [Nix](https://nixos.org/)

1. Setting up the LoKI-B development environment.

   To launch a development shell use
   ```bash
   nix develop
   ```
   in the root of the repository. This shell will contain all necessary
   dependencies to build LoKI-B. Additionally, all binaries can be built using
   ```bash
   nix build
   ```
   Alternatively, a Cobertura coverage report can be generated using
   ```bash
   nix build .#coverage
   ```
   The resulting output will be available in `<repo_root>/result`.

2. Installing LoKI-B without manually cloning the repository.

   The LoKI-B binary can be built and run anywhere by issuing
   ```bash
   nix run github:loki-suite/loki-b <input_file>
   ```

   Similarly, a shell can be launched with access to the `loki` binary using
   ```bash
   nix shell github:loki-suite/loki-b
   ```

   Additionally, the web version (using WebAssembly) can be built using the
   `loki-web` package.
   ```bash
   nix build github:loki-suite/loki-b#loki-web
   ```

   To serve the pages you can then e.g. use the python http server.
   ```bash
   python -m http.server -d result/share/loki-web
   ```

### Windows

1. run: `cmake -D<BACKEND_FLAG>=ON -B build`
   - where `<BACKEND_FLAG>=LOKIB_USE_MKL/LOKIB_USE_OPENBLAS`, specifying the
     backend to supply to Eigen
   - this flag can also be omitted to build with pure Eigen
1. run: `cmake --build build --config Release -j <NUM_JOBS>`
   - where `<NUM_JOBS>` is the maximum number of jobs to run simultaneously when
     compiling; just use the number of physical cores in your system. Omit this
     flag to use the default settings.

## Running and Plotting

As of yet, LoKI-B++ is run from the command line. To plot the computed eedf you
can install gnuplot. The execution on different operating systems is very
similar.

### Linux

1. Make sure [gnuplot](http://gnuplot.info/) is available on your system.
1. After compiling the `loki` executable, run the following command in the root
   of the repository: `./build/app/loki <INPUT_FILE> | gnuplot --persist`.
   - Where `<INPUT_FILE>` is the path to the input file relative to the current
     working directory.
   - The output from LoKI-B++ is then piped into gnuplot, where the `--persist`
     flag avoids gnuplot from immediately closing after plotting.

### Windows

1. Make sure [gnuplot](http://gnuplot.info/) is available on your system.
1. Run the following command:
   `.\build\app\loki.exe <INPUT_FILE> | gnuplot --persist`.
   - Where `<INPUT_FILE>` is the path to the input file relative to the current
     working directory.
   - The output from LoKI-B++ is then piped into gnuplot, where the `--persist`
     flag avoids gnuplot from immediately closing after plotting.
