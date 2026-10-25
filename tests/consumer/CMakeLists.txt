cmake_minimum_required(VERSION 3.20)

project(gcrank_consumer LANGUAGES CXX)

# Deliberately C++20: the installed headers and the exported target must work
# without the C++23 the library itself is built with.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(gcrank 0.3 REQUIRED)

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE gcrank::gcrank)
