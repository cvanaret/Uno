# Copyright (c) 2025 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

cmake_minimum_required(VERSION 3.7)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
   cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(tests_hs015 LANGUAGES C)

# compile options
if (MSVC)
   add_compile_options("$<$<COMPILE_LANGUAGE:C>:/utf-8>")
else()
   set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wnon-virtual-dtor -pedantic -Wunused-value -Wconversion")
   set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG") # disable asserts
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmaybe-uninitialized")
endif()

find_library(UNO REQUIRED)
add_executable(tests_hs015 tests_hs015.c)
# link the math library and Uno
target_link_libraries(tests_hs015 PUBLIC m ${UNO})