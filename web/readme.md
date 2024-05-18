## How to compile LoKI-B to WebAssembly

Prerequisites:
 - Make sure the Emscripten SDK is installed. See
   <https://emscripten.org/docs/getting_started/downloads.html>

To compile the project using emscripten and cmake, run the following commands
in the project root.
```
emcmake cmake -DUSE_OPENMP=OFF -DUSE_BUILTIN_EIGEN=ON -DUSE_BUILTIN_NLOHMANN_JSON=ON -DCMAKE_BUILD_TYPE=Release -B build .
cmake --build build
```
The web version of LoKI-B can be launched by e.g. starting a local Python
server in `web`: `python3 -m http.server 8000`. Visit <http://localhost:8000>
to view the demo user interface.

### Using Nix

LoKI-B provides a Nix derivation for the web build called `loki-web`. Therefore,
you can also build the web build using
```bash
nix build .#
```

or from anywhere, without the need to clone the repository using
```bash
nix build github:loki-suite/loki-b#loki-web
```

The resulting files will be available in `./result/share/loki-web`.
