#This custom target_add_test function also handles target relocation 
#by the CMAKE_RUNTIME_OUTPUT_DIRECTORY directive. see also: 
#https://stackoverflow.com/questions/21394768/cmake-cannot-find-test-if-cmake-runtime-output-directory-is-changed
function (target_add_test target)
  add_test (NAME ${target} COMMAND $<TARGET_FILE:${target}>)

  if (LOKIB_COVERAGE)
    set(LOKIB_COVERAGE_TESTS ${LOKIB_COVERAGE_TESTS} ${target} PARENT_SCOPE)
  endif()
endfunction (target_add_test)

enable_testing()

# Option -Dcoverage: configure whether coverage information is collected.
option(coverage "-Dcoverage(=OFF|ON): Collect test coverage statistics." OFF)
set(LOKIB_COVERAGE ${coverage})
if(LOKIB_COVERAGE)
  include(../CMake/CodeCoverage.cmake)
  append_coverage_compiler_flags_to_target(loki-b)
  set(LOKIB_COVERAGE_TESTS loki-b)
endif()

add_executable(bolsig_extract bolsig_extract.cpp)
target_link_libraries(bolsig_extract PRIVATE loki-b)

add_executable(compare_tables compare_tables.cpp)
target_link_libraries(compare_tables PRIVATE loki-b)

add_executable(test_druyvesteyn test_druyvesteyn.cpp)
target_link_libraries(test_druyvesteyn PRIVATE loki-b)

add_executable(test_ee test_ee.cpp)
target_link_libraries(test_ee PRIVATE loki-b)

add_executable(test_ee_relaxation test_ee_relaxation.cpp)
target_link_libraries(test_ee_relaxation PRIVATE loki-b)

add_executable(test_events test_events.cpp)
target_link_libraries(test_events PRIVATE loki-b)

add_executable(test_lut test_lut.cpp)
target_link_libraries(test_lut PRIVATE loki-b)

add_executable(test_parser test_parser.cpp)
target_link_libraries(test_parser PRIVATE loki-b)

add_executable(test_reaction_format test_reaction_format.cpp)
target_link_libraries(test_reaction_format PRIVATE loki-b)

add_executable(test_statebase test_statebase.cpp)
target_link_libraries(test_statebase PRIVATE loki-b)

add_executable(test_trans_coefs test_trans_coefs.cpp)
target_link_libraries(test_trans_coefs PRIVATE loki-b)

target_add_test(test_druyvesteyn)
target_add_test(test_ee)
target_add_test(test_ee_relaxation)
target_add_test(test_events)
target_add_test(test_lut)
target_add_test(test_parser)
target_add_test(test_reaction_format)
target_add_test(test_statebase)
target_add_test(test_trans_coefs)

if (LOKIB_COVERAGE)
  setup_target_for_coverage_gcovr_xml(
    NAME coverage 
    EXECUTABLE ctest --test-dir tests 
    BASE_DIRECTORY "${PROJECT_SOURCE_DIR}"
    DEPENDENCIES ${LOKIB_COVERAGE_TESTS}
  )
endif()
