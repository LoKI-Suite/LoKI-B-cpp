# -- WebAssembly -- #
message(STATUS "Compiling to WebAssembly.")
add_executable(loki_bindings bindings.cpp)

target_include_directories(loki_bindings PRIVATE loki-b /usr/lib/emscripten/system/include)
target_link_libraries(loki_bindings PRIVATE loki-b)

set_target_properties(loki_bindings PROPERTIES LINK_FLAGS "\
  -s ASSERTIONS=1 \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s MODULARIZE=1 \
  -s EXPORT_NAME=create_module \
  -s EXPORTED_RUNTIME_METHODS=FS,run \
  --bind \
  --preload-file ${PROJECT_SOURCE_DIR}/Input/JSON/N2/N2_LXCat.json@JSON/N2/ \
  --preload-file ${PROJECT_SOURCE_DIR}/Input/JSON/N2/N2_rot_LXCat.json@JSON/N2/ \
  --preload-file ${PROJECT_SOURCE_DIR}/Input/Databases/masses.txt@Input/Databases/ \
  --preload-file ${PROJECT_SOURCE_DIR}/Input/Databases/harmonicFrequencies.txt@Input/Databases/ \
  --preload-file ${PROJECT_SOURCE_DIR}/Input/Databases/rotationalConstants.txt@Input/Databases/ \
  --preload-file ${PROJECT_SOURCE_DIR}/Input/Databases/OPBParameter.txt@Input/Databases/ \
  --preload-file ${PROJECT_SOURCE_DIR}/Input/Databases/quadrupoleMoment.txt@Input/Databases/ \
  --preload-file ${PROJECT_SOURCE_DIR}/Input/Databases/anharmonicFrequencies.txt@Input/Databases/")

# FIXME: The generated files should be placed in `web/wasm`. However, this is
# currently not possible due to a bug in emscripten
# (https://github.com/emscripten-core/emscripten/issues/16240).
set_target_properties(loki_bindings 
  PROPERTIES 
  ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/web"
  LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/web"
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/web"
)
