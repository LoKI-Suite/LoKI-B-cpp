## How to compile LoKI-B to WebAssembly

Prerequisites:
 - Make sure the Emscripten SDK is installed. See https://emscripten.org/docs/getting_started/downloads.html

Run the following commands in the current folder `/web`.
 - Compile the project using emscripten and cmake. First, do `emcmake cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release ..`. Next, run `make`.
 - Setup a local server using Python: `python3 -m http.server 8000`.

Then visit <http://localhost:8000> to view the demo user interface.
