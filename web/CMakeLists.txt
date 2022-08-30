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
	-s EXTRA_EXPORTED_RUNTIME_METHODS=[\"FS\"] \
	--bind \
	--preload-file ${PROJECT_SOURCE_DIR}/web/Input/default_lokib_setup.json@Input/ \
	--preload-file ${PROJECT_SOURCE_DIR}/web/Input/N2_LXCat.json@Input/ \
	--preload-file ${PROJECT_SOURCE_DIR}/web/Input/N2_rot_LXCat.json@Input/ \
	--preload-file ${PROJECT_SOURCE_DIR}/web/Input/masses.txt@Input/ \
	--preload-file ${PROJECT_SOURCE_DIR}/web/Input/harmonicFrequencies.txt@Input/ \
	--preload-file ${PROJECT_SOURCE_DIR}/web/Input/rotationalConstants.txt@Input/ \
	--preload-file ${PROJECT_SOURCE_DIR}/web/Input/N2_vibpop.txt@Input/ \
	--preload-file ${PROJECT_SOURCE_DIR}/web/Input/OPBParameter.txt@Input/ \
	--preload-file ${PROJECT_SOURCE_DIR}/web/Input/quadrupoleMoment.txt@Input/ \
	--preload-file ${PROJECT_SOURCE_DIR}/web/Input/anharmonicFrequencies.txt@Input/")

