cmake_minimum_required(VERSION 3.28)

project(
    linalg_mini_library
    VERSION 0.1.0
    DESCRIPTION "Toy linalg library"
    LANGUAGES CXX
)


set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


add_executable(
        hello
        tests/hello.cpp
)

add_executable(
        test
        tests/test.cpp
)
    
add_library(lml
    src/matrix.cpp
)

target_link_libraries(
    test
    PRIVATE
        lml
)

target_include_directories(
    lml
    PUBLIC
        include
)
