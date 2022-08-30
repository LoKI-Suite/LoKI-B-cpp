# -- WebAssembly -- #
message(STATUS "Compiling to WebAssembly.")
add_executable(loki_bindings bindings.cpp)
target_include_directories(loki_bindings PRIVATE loki-b ${inc} /usr/lib/emscripten/system/include)
target_link_libraries(loki_bindings PRIVATE loki-b)
set_target_properties(loki_bindings PROPERTIES LINK_FLAGS "\
	-s ASSERTIONS=1 \
	-s ALLOW_MEMORY_GROWTH=1 \
	-s MODULARIZE=1 \
	-s EXPORT_NAME=create_module \
	-s EXTRA_EXPORTED_RUNTIME_METHODS=[\"FS\"] \
	--bind \
	--preload-file Input/default_lokib_setup.json \
	--preload-file Input/N2_LXCat.json \
	--preload-file Input/N2_rot_LXCat.json \
	--preload-file Input/masses.txt \
	--preload-file Input/harmonicFrequencies.txt \
	--preload-file Input/rotationalConstants.txt \
	--preload-file Input/N2_vibpop.txt \
	--preload-file Input/OPBParameter.txt \
	--preload-file Input/quadrupoleMoment.txt \
	--preload-file Input/anharmonicFrequencies.txt")

