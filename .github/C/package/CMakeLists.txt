# Copyright (c) 2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# Consumer of the installed CMake package of Uno (find_package + Uno::uno), used by CI.
# The C++ language is enabled because the static library of Uno needs the C++ linker and runtime.

cmake_minimum_required(VERSION 3.12)
project(uno_package_consumer LANGUAGES C CXX)

set(UNO_TARGET "Uno::uno" CACHE STRING "Uno target to link against (Uno::uno, Uno::uno_static or Uno::uno_shared)")

find_package(Uno 2 REQUIRED CONFIG)
message(STATUS "Found Uno ${Uno_VERSION} in ${Uno_DIR}")

add_executable(tests_hs015 ${CMAKE_CURRENT_SOURCE_DIR}/../tests_hs015.c)
target_link_libraries(tests_hs015 PRIVATE ${UNO_TARGET})
if(NOT WIN32)
   target_link_libraries(tests_hs015 PRIVATE m)
endif()
