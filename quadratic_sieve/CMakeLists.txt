cmake_minimum_required(VERSION 3.12)

project(QuadraticSieve VERSION 1.0 LANGUAGES CXX)

include_directories(.)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-pthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -flto -march=native -mpopcnt")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(
    main
    src/main.cpp
    src/quadratic_sieve.cpp
    src/auxiliary.cpp
)

add_executable(
    test
    src/test.cpp
    src/quadratic_sieve.cpp
    src/auxiliary.cpp
)
