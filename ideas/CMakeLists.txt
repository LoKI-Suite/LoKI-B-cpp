add_executable(json2lxcat json2lxcat.cpp)
target_link_libraries(json2lxcat PRIVATE nlohmann_json::nlohmann_json)

add_executable(test_statecontainer test_statecontainer.cpp)
target_link_libraries(test_statecontainer PRIVATE loki-b)
add_test(statecontainer test_statecontainer)


