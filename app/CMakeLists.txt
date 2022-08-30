# -- loki native default front end -- #
add_executable(loki loki_main.cpp)
target_link_libraries(loki PRIVATE loki-b)

install(TARGETS loki)
