#This custom target_add_test function also handles target relocation 
#by the CMAKE_RUNTIME_OUTPUT_DIRECTORY directive. see also: 
#https://stackoverflow.com/questions/21394768/cmake-cannot-find-test-if-cmake-runtime-output-directory-is-changed
function (target_add_test target)
  add_test (NAME ${target} COMMAND $<TARGET_FILE:${target}>)
endfunction (target_add_test)

add_executable(json2lxcat json2lxcat.cpp)
if(USE_BUILTIN_NLOHMANN_JSON)
	target_include_directories(json2lxcat PRIVATE ${PROJECT_SOURCE_DIR}/include/nlohmann-3.7.0)
else()
	target_link_libraries(json2lxcat PRIVATE nlohmann_json::nlohmann_json)
endif()

add_executable(loki_steps loki_steps.cpp)
target_link_libraries(loki_steps PRIVATE loki-b)

add_executable(test_statecontainer test_statecontainer.cpp)
target_link_libraries(test_statecontainer PRIVATE loki-b)
target_add_test(test_statecontainer)
