add_executable(loki_gasproperties loki_gasproperties.cpp)
target_link_libraries(loki_gasproperties PRIVATE loki_cpp::loki-b)

add_executable(loki_legacytojson loki_legacytojson.cpp)
target_link_libraries(loki_legacytojson PRIVATE loki_cpp::loki-b)

add_executable(loki loki_main.cpp)
target_link_libraries(loki PRIVATE loki_cpp::loki-b)

add_executable(loki_offsidetojson loki_offsidetojson.cpp)
target_link_libraries(loki_offsidetojson PRIVATE loki_cpp::loki-b)

install(TARGETS loki loki_legacytojson loki_offsidetojson)
