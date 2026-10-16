add_executable(gcrank_tests
  test_main.cpp
  test_bitvector.cpp
  test_fasta_reader.cpp
  test_manifest.cpp
  test_mmap.cpp
  test_parallel.cpp
  test_roundtrip.cpp)
target_link_libraries(gcrank_tests PRIVATE gcrank)
target_compile_options(gcrank_tests PRIVATE
  ${GCRANK_WARNING_FLAGS} ${GCRANK_TSAN_COMPILE_OPTIONS})
target_link_options(gcrank_tests PRIVATE ${GCRANK_TSAN_LINK_OPTIONS})

# doctest.h lives in this directory and is pulled in with <doctest.h>, so the
# SYSTEM include keeps its warnings out of our build.
target_include_directories(gcrank_tests SYSTEM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(gcrank_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(gcrank_tests PRIVATE
  GCRANK_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data"
  GCRANK_TEST_TMP_DIR="${CMAKE_CURRENT_BINARY_DIR}/tmp")

add_test(NAME unit COMMAND gcrank_tests)
set_tests_properties(unit PROPERTIES TIMEOUT 900)

# Doc-test for the on-disk bit layout: an independent Python reader decodes the
# golden index using only the rule documented in format.h/README.md and checks
# it against the reference values. Skipped when Python is not available.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
  add_test(NAME spec
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/data-raw/check_spec.py)
  set_tests_properties(spec PROPERTIES TIMEOUT 120)
endif()

if(GCRANK_BUILD_CLI)
  add_test(NAME cli
    COMMAND ${CMAKE_COMMAND}
      -DCLI=$<TARGET_FILE:gcrank-build>
      -DQUERY=$<TARGET_FILE:gcrank-query>
      -DDATA=${CMAKE_CURRENT_SOURCE_DIR}/data
      -DTMP=${CMAKE_CURRENT_BINARY_DIR}/cli-tmp
      -P ${CMAKE_CURRENT_SOURCE_DIR}/cli_test.cmake)
  set_tests_properties(cli PROPERTIES TIMEOUT 300)

  # The headline claim, checked against the real bedtools when it is installed:
  # pct_gc identity on an exhaustive window set plus the adversarial rounding
  # cases. Skips (exit 77) when bedtools is not available.
  if(Python3_Interpreter_FOUND)
    find_program(GCRANK_BEDTOOLS bedtools)
    add_test(NAME bedtools
      COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/data-raw/check_bedtools_parity.py
        --query $<TARGET_FILE:gcrank-query>
        --build $<TARGET_FILE:gcrank-build>
        --bedtools ${GCRANK_BEDTOOLS}
        --data ${CMAKE_CURRENT_SOURCE_DIR}/data
        --tmp ${CMAKE_CURRENT_BINARY_DIR}/bedtools-tmp)
    set_tests_properties(bedtools PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 900)
  endif()
endif()

# Installs into a scratch prefix and builds a C++20 consumer with
# find_package(gcrank), exercising the export and the installed headers.
# Skipped for sanitizer builds: a consumer of an instrumented static library
# would have to link the sanitizer runtime, which is not a supported scenario.
set(GCRANK_SANITIZED OFF)
if(GCRANK_TSAN OR CMAKE_CXX_FLAGS MATCHES "fsanitize" OR CMAKE_EXE_LINKER_FLAGS MATCHES "fsanitize")
  set(GCRANK_SANITIZED ON)
endif()
if(GCRANK_SANITIZED)
  message(STATUS "sanitizer build: skipping the install/consumer test")
else()
  add_test(NAME consumer
    COMMAND ${CMAKE_COMMAND}
      -DBUILD_DIR=${CMAKE_BINARY_DIR}
      -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
      -DCONFIG=$<CONFIG>
      -DTMP=${CMAKE_CURRENT_BINARY_DIR}/consumer-tmp
      -P ${CMAKE_CURRENT_SOURCE_DIR}/consumer_test.cmake)
  set_tests_properties(consumer PROPERTIES TIMEOUT 600)
endif()
