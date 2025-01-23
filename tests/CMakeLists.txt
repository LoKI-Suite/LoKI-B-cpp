#This custom target_add_test function also handles target relocation 
#by the CMAKE_RUNTIME_OUTPUT_DIRECTORY directive. see also: 
#https://stackoverflow.com/questions/21394768/cmake-cannot-find-test-if-cmake-runtime-output-directory-is-changed
function (target_add_test target)
  add_test (NAME ${target} COMMAND $<TARGET_FILE:${target}>)
  set_property(TEST ${target} PROPERTY ENVIRONMENT "LOKI_TEST_INPUT_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
  if (LOKIB_COVERAGE)
    set(LOKIB_COVERAGE_TESTS ${LOKIB_COVERAGE_TESTS} ${target} PARENT_SCOPE)
  endif()
endfunction (target_add_test)

enable_testing()

# Option -Dcoverage: configure whether coverage information is collected.
option(coverage "-Dcoverage(=OFF|ON): Collect test coverage statistics." OFF)
loki_constrain_configure_option(coverage ON OFF)
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

add_executable(loki_test_derivatives test_derivatives.cpp)
target_link_libraries(loki_test_derivatives PRIVATE loki-b)

add_executable(loki_test_eedf_utilities test_eedf_utilities.cpp)
target_link_libraries(loki_test_eedf_utilities PRIVATE loki-b)

add_executable(loki_test_gas_properties test_gas_properties.cpp)
target_link_libraries(loki_test_gas_properties PRIVATE loki-b)

add_executable(loki_test_interpolation test_interpolation.cpp)
target_link_libraries(loki_test_interpolation PRIVATE loki-b)

add_executable(loki_test_offside_to_json test_offside_to_json.cpp)
target_link_libraries(loki_test_offside_to_json PRIVATE loki-b)

add_executable(test_druyvesteyn test_druyvesteyn.cpp)
target_link_libraries(test_druyvesteyn PRIVATE loki-b)

add_executable(test_ee test_ee.cpp)
target_link_libraries(test_ee PRIVATE loki-b)

add_executable(test_ee_relaxation test_ee_relaxation.cpp)
target_link_libraries(test_ee_relaxation PRIVATE loki-b)

add_executable(test_events test_events.cpp)
target_link_libraries(test_events PRIVATE loki-b)

add_executable(test_grid test_grid.cpp)
target_link_libraries(test_grid PRIVATE loki-b)

add_executable(loki_test_integrals test_integrals.cpp)
target_link_libraries(loki_test_integrals PRIVATE loki-b)

add_executable(test_linalg test_linalg.cpp)
target_link_libraries(test_linalg PRIVATE loki-b)

add_executable(test_lut test_lut.cpp)
target_link_libraries(test_lut PRIVATE loki-b)

add_executable(test_parser test_parser.cpp)
target_link_libraries(test_parser PRIVATE loki-b)

add_executable(loki_test_range test_range.cpp)
target_link_libraries(loki_test_range PRIVATE loki-b)

add_executable(test_reaction_format test_reaction_format.cpp)
target_link_libraries(test_reaction_format PRIVATE loki-b)

add_executable(test_statebase test_statebase.cpp)
target_link_libraries(test_statebase PRIVATE loki-b)

add_executable(test_trans_coefs test_trans_coefs.cpp)
target_link_libraries(test_trans_coefs PRIVATE loki-b)

add_executable(test_nonuniform_elastic test_nonuniform_elastic.cpp)
target_link_libraries(test_nonuniform_elastic PRIVATE loki-b)

add_executable(test_maxwell test_maxwell.cpp)
target_link_libraries(test_maxwell PRIVATE loki-b)

add_executable(test_nonuniform_field test_nonuniform_field.cpp)
target_link_libraries(test_nonuniform_field PRIVATE loki-b)

add_executable(test_nonuniform_druyvesteyn test_nonuniform_druyvesteyn.cpp)
target_link_libraries(test_nonuniform_druyvesteyn PRIVATE loki-b)

add_executable(test_nonuniform_CAR test_nonuniform_CAR.cpp)
target_link_libraries(test_nonuniform_CAR PRIVATE loki-b)

add_executable(test_nonuniform_inelastic test_nonuniform_inelastic.cpp)
target_link_libraries(test_nonuniform_inelastic PRIVATE loki-b)

add_executable(test_nonuniform_operator_distribution test_nonuniform_operator_distribution.cpp)
target_link_libraries(test_nonuniform_operator_distribution PRIVATE loki-b)

add_executable(test_invalid_inputs test_invalid_inputs.cpp)
target_link_libraries(test_invalid_inputs PRIVATE loki-b)

add_executable(test_valid_inputs test_valid_inputs.cpp)
target_link_libraries(test_valid_inputs PRIVATE loki-b)

target_add_test(loki_test_offside_to_json)
target_add_test(loki_test_derivatives)
target_add_test(loki_test_eedf_utilities)
target_add_test(test_druyvesteyn)
target_add_test(test_ee)
target_add_test(test_ee_relaxation)
target_add_test(test_events)
target_add_test(test_grid)
target_add_test(loki_test_integrals)
target_add_test(loki_test_interpolation)
target_add_test(test_linalg)
target_add_test(test_lut)
target_add_test(test_parser)
target_add_test(loki_test_range)
target_add_test(test_reaction_format)
target_add_test(test_statebase)
target_add_test(test_trans_coefs)
target_add_test(test_nonuniform_elastic)
target_add_test(test_maxwell)
target_add_test(test_nonuniform_field)
target_add_test(test_nonuniform_druyvesteyn)
target_add_test(test_nonuniform_CAR)
target_add_test(test_nonuniform_inelastic)
target_add_test(test_nonuniform_operator_distribution)
target_add_test(test_invalid_inputs)
target_add_test(test_valid_inputs)

if (LOKIB_COVERAGE)
  setup_target_for_coverage_gcovr_xml(
    NAME coverage 
    EXECUTABLE ctest --test-dir tests 
    BASE_DIRECTORY "${PROJECT_SOURCE_DIR}"
    DEPENDENCIES ${LOKIB_COVERAGE_TESTS}
  )
endif()
