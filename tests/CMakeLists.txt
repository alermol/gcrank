add_executable(gcrank_tests
  test_main.cpp
  test_bitvector.cpp
  test_fasta_reader.cpp
  test_manifest.cpp
  test_mmap.cpp
  test_parallel.cpp
  test_roundtrip.cpp)
target_link_libraries(gcrank_tests PRIVATE gcrank)
target_compile_options(gcrank_tests PRIVATE ${GCRANK_WARNING_FLAGS})

# doctest.h lives in this directory and is pulled in with <doctest.h>, so the
# SYSTEM include keeps its warnings out of our build.
target_include_directories(gcrank_tests SYSTEM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(gcrank_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(gcrank_tests PRIVATE
  GCRANK_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data"
  GCRANK_TEST_TMP_DIR="${CMAKE_CURRENT_BINARY_DIR}/tmp")

add_test(NAME unit COMMAND gcrank_tests)
set_tests_properties(unit PROPERTIES TIMEOUT 900)

if(GCRANK_BUILD_CLI)
  add_test(NAME cli
    COMMAND ${CMAKE_COMMAND}
      -DCLI=$<TARGET_FILE:gcrank-build>
      -DDATA=${CMAKE_CURRENT_SOURCE_DIR}/data
      -DTMP=${CMAKE_CURRENT_BINARY_DIR}/cli-tmp
      -P ${CMAKE_CURRENT_SOURCE_DIR}/cli_test.cmake)
  set_tests_properties(cli PROPERTIES TIMEOUT 300)
endif()
