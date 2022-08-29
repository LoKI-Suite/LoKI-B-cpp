add_executable(bolsig_extract bolsig_extract.cpp)
target_link_libraries(bolsig_extract PRIVATE loki-b)

add_executable(compare_tables compare_tables.cpp)
target_link_libraries(compare_tables PRIVATE loki-b)

add_executable(test_lut test_lut.cpp)
target_link_libraries(test_lut PRIVATE loki-b)

add_executable(test_parser test_parser.cpp)
target_link_libraries(test_parser PRIVATE loki-b)

add_executable(test_reaction_format test_reaction_format.cpp)
target_link_libraries(test_reaction_format PRIVATE loki-b)

add_executable(test_statebase test_statebase.cpp)
target_link_libraries(test_statebase PRIVATE loki-b)

add_executable(test_druyvesteyn test_druyvesteyn.cpp)
target_link_libraries(test_druyvesteyn PRIVATE loki-b)

add_test(lut test_lut)
add_test(parser test_parser)
add_test(reaction_format test_reaction_format)
add_test(statebase test_statebase)
add_test(druyvesteyn test_druyvesteyn)
