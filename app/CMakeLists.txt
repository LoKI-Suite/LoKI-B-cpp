add_executable(loki loki_main.cpp)
target_link_libraries(loki PRIVATE loki-b)

add_executable(loki_legacytojson loki_legacytojson.cpp)
target_link_libraries(loki_legacytojson PRIVATE loki-b)

add_executable(loki_offsidetojson loki_offsidetojson.cpp)
target_link_libraries(loki_offsidetojson PRIVATE loki-b)

install(TARGETS loki loki_legacytojson loki_offsidetojson)
