# -- WebAssembly -- #
message(STATUS "Compiling to WebAssembly.")
add_executable(loki_bindings bindings.cpp)

target_include_directories(loki_bindings PRIVATE loki-b /usr/lib/emscripten/system/include)
target_link_libraries(loki_bindings PRIVATE loki_cpp::loki-b)

set_target_properties(loki_bindings PROPERTIES LINK_FLAGS "\
  -s ASSERTIONS=1 \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s MODULARIZE=1 \
  -s EXPORT_NAME=create_module \
  -s EXPORTED_RUNTIME_METHODS=[FS,run] \
  -s ENVIRONMENT=web,worker \
  -s NO_DISABLE_EXCEPTION_CATCHING \
  -lembind \
  --emit-tsd loki_bindings.d.ts \
  --preload-file ${PROJECT_SOURCE_DIR}/input/JSON/N2/N2_LXCat.json@JSON/N2/ \
  --preload-file ${PROJECT_SOURCE_DIR}/input/JSON/N2/N2_rot_LXCat.json@JSON/N2/ \
  --preload-file ${PROJECT_SOURCE_DIR}/input/Databases/masses.txt@Databases/ \
  --preload-file ${PROJECT_SOURCE_DIR}/input/Databases/harmonicFrequencies.txt@Databases/ \
  --preload-file ${PROJECT_SOURCE_DIR}/input/Databases/rotationalConstants.txt@Databases/ \
  --preload-file ${PROJECT_SOURCE_DIR}/input/Databases/OPBParameter.txt@Databases/ \
  --preload-file ${PROJECT_SOURCE_DIR}/input/Databases/quadrupoleMoment.txt@Databases/ \
  --preload-file ${PROJECT_SOURCE_DIR}/input/Databases/anharmonicFrequencies.txt@Databases/")

# FIXME: The generated files should be placed in `web/wasm`. However, this is
# currently not possible due to a bug in emscripten
# (https://github.com/emscripten-core/emscripten/issues/16240).
set_target_properties(loki_bindings
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/web"
  LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/web"
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/web"
)
