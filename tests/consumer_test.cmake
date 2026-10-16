# Install-and-consume test for the exported CMake package.
#
# Invoked by CTest as:
#   cmake -DBUILD_DIR=<build> -DSOURCE_DIR=<src> -DCONFIG=<cfg> -DTMP=<scratch> -P consumer_test.cmake
#
# Installs the build into a scratch prefix and builds a C++20 consumer project
# with find_package(gcrank), so the exported target and the installed headers
# are exercised the way a downstream user would use them.
foreach(required BUILD_DIR SOURCE_DIR TMP)
  if(NOT DEFINED ${required})
    message(FATAL_ERROR "consumer_test.cmake: -D${required}=... is required")
  endif()
endforeach()

file(REMOVE_RECURSE "${TMP}")
file(MAKE_DIRECTORY "${TMP}/prefix")

set(config_args)
if(DEFINED CONFIG AND NOT "${CONFIG}" STREQUAL "")
  list(APPEND config_args --config "${CONFIG}")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}" --install "${BUILD_DIR}" --prefix "${TMP}/prefix" ${config_args}
  RESULT_VARIABLE install_code
  OUTPUT_VARIABLE install_out
  ERROR_VARIABLE install_err
  WORKING_DIRECTORY "${TMP}")
if(NOT install_code EQUAL 0)
  message(FATAL_ERROR "install failed (${install_code}):\n${install_out}\n${install_err}")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}" -S "${SOURCE_DIR}/tests/consumer" -B "${TMP}/build"
          "-DCMAKE_PREFIX_PATH=${TMP}/prefix"
  RESULT_VARIABLE configure_code
  OUTPUT_VARIABLE configure_out
  ERROR_VARIABLE configure_err
  WORKING_DIRECTORY "${TMP}")
if(NOT configure_code EQUAL 0)
  message(FATAL_ERROR "consumer configure failed (${configure_code}):\n${configure_out}\n${configure_err}")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}" --build "${TMP}/build" ${config_args}
  RESULT_VARIABLE build_code
  OUTPUT_VARIABLE build_out
  ERROR_VARIABLE build_err
  WORKING_DIRECTORY "${TMP}")
if(NOT build_code EQUAL 0)
  message(FATAL_ERROR "consumer build failed (${build_code}):\n${build_out}\n${build_err}")
endif()

if(WIN32)
  set(consumer_exe "${TMP}/build/consumer.exe")
else()
  set(consumer_exe "${TMP}/build/consumer")
endif()
execute_process(
  COMMAND "${consumer_exe}"
  RESULT_VARIABLE run_code
  OUTPUT_VARIABLE run_out
  ERROR_VARIABLE run_err
  WORKING_DIRECTORY "${TMP}")
if(NOT run_code EQUAL 0)
  message(FATAL_ERROR "consumer run failed (${run_code}):\n${run_out}\n${run_err}")
endif()

message(STATUS "gcrank install/consumer test passed: ${run_out}")
