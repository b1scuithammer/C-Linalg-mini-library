cmake_minimum_required(VERSION 3.28)

project(
    linalg_mini_library
    VERSION 0.1.0
    DESCRIPTION "Finite Element Method Solver"
    LANGUAGES CXX
)


set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(${PROJECT_NAME}
        tests/hello.cpp
)
    
add_library(linalg_mini_lib
    src/matrix.cpp
)

target_include_directories(
    linalg_mini_lib
    PUBLIC
        include
)
