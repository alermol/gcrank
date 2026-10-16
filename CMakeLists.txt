cmake_minimum_required(VERSION 3.20)

project(gcrank
  VERSION 0.1.0
  DESCRIPTION "Standalone C++23 library and CLI for O(1) GC-content queries"
  LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(GCRANK_BUILD_TESTS "Build the gcrank test suite" ON)
option(GCRANK_BUILD_CLI "Build the gcrank-build command line tool" ON)
option(GCRANK_TSAN "Build gcrank with ThreadSanitizer" OFF)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(GCRANK_IS_TOP_LEVEL ON)
else()
  set(GCRANK_IS_TOP_LEVEL OFF)
endif()

# --- dependencies -----------------------------------------------------------

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# --- targets ----------------------------------------------------------------

add_library(gcrank
  src/bitvector.cpp
  src/crc32.cpp
  src/fasta_reader.cpp
  src/gc_rank.cpp
  src/input_stream.cpp
  src/line_reader.cpp
  src/manifest.cpp
  src/sha256.cpp
  src/thread_pool.cpp)
add_library(gcrank::gcrank ALIAS gcrank)

target_compile_features(gcrank PUBLIC cxx_std_20)
target_compile_features(gcrank PRIVATE cxx_std_23)
target_include_directories(gcrank PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(gcrank PRIVATE ZLIB::ZLIB Threads::Threads)

set_target_properties(gcrank PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})

# --- sanitizer plumbing ------------------------------------------------------
# ThreadSanitizer flags are applied PRIVATE to every target of this project, so
# they can never leak into the installed CMake export (a consumer must not be
# forced to link the sanitizer runtime).
set(GCRANK_TSAN_COMPILE_OPTIONS "")
set(GCRANK_TSAN_LINK_OPTIONS "")
if(GCRANK_TSAN)
  set(GCRANK_TSAN_COMPILE_OPTIONS -fsanitize=thread -g -O1)
  set(GCRANK_TSAN_LINK_OPTIONS -fsanitize=thread)
endif()

# Warning flags shared by the library, the CLI and the tests.
if(MSVC)
  set(GCRANK_WARNING_FLAGS /W4)
else()
  set(GCRANK_WARNING_FLAGS -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

target_compile_options(gcrank PRIVATE ${GCRANK_WARNING_FLAGS} ${GCRANK_TSAN_COMPILE_OPTIONS})
if(GCRANK_TSAN)
  target_link_options(gcrank PRIVATE ${GCRANK_TSAN_LINK_OPTIONS})
endif()

if(GCRANK_BUILD_CLI)
  add_executable(gcrank-build cli/gcrank_build.cpp)
  target_compile_options(gcrank-build PRIVATE
    ${GCRANK_WARNING_FLAGS} ${GCRANK_TSAN_COMPILE_OPTIONS})
  target_link_options(gcrank-build PRIVATE ${GCRANK_TSAN_LINK_OPTIONS})
  target_link_libraries(gcrank-build PRIVATE gcrank)

  add_executable(gcrank-query cli/gcrank_query.cpp)
  target_compile_options(gcrank-query PRIVATE
    ${GCRANK_WARNING_FLAGS} ${GCRANK_TSAN_COMPILE_OPTIONS})
  target_link_options(gcrank-query PRIVATE ${GCRANK_TSAN_LINK_OPTIONS})
  target_link_libraries(gcrank-query PRIVATE gcrank)
endif()

if(GCRANK_BUILD_TESTS)
  # Compiles the library example from README.md so the documentation cannot
  # drift away from the API (see examples/readme.cpp).
  add_executable(gcrank-example-readme examples/readme.cpp)
  target_compile_options(gcrank-example-readme PRIVATE
    ${GCRANK_WARNING_FLAGS} ${GCRANK_TSAN_COMPILE_OPTIONS})
  target_link_options(gcrank-example-readme PRIVATE ${GCRANK_TSAN_LINK_OPTIONS})
  target_link_libraries(gcrank-example-readme PRIVATE gcrank)
endif()

if(GCRANK_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# --- installation -----------------------------------------------------------

install(TARGETS gcrank EXPORT gcrankTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(GCRANK_BUILD_CLI)
  install(TARGETS gcrank-build gcrank-query RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# Only the public headers are installed; everything under detail/ is internal.
install(FILES
  include/gcrank/error.h
  include/gcrank/format.h
  include/gcrank/gc_rank.h
  include/gcrank/manifest.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gcrank)

set(GCRANK_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/gcrank)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/gcrankConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/gcrankConfig.cmake
  INSTALL_DESTINATION ${GCRANK_CMAKE_DIR})
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/gcrankConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/gcrankConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/gcrankConfigVersion.cmake
  DESTINATION ${GCRANK_CMAKE_DIR})
install(EXPORT gcrankTargets
  FILE gcrankTargets.cmake
  NAMESPACE gcrank::
  DESTINATION ${GCRANK_CMAKE_DIR})
