# Copyright (c) 2025 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

cmake_minimum_required(VERSION 3.7)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
   cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(example_hs015 LANGUAGES C)

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

# external dependencies
set(LIBRARIES "")
# additional libraries to link against
set(AUXILIARY_LIBRARIES "" CACHE STRING "Additional libraries to link against")

# Uno
find_library(UNO REQUIRED)
list(APPEND LIBRARIES m ${UNO}) # math library and Uno

# LAPACK
# convert relative paths of LAPACK to absolute paths
if(LAPACK_LIBRARIES)
   set(ABSOLUTE_LAPACK_LIBRARIES)
   foreach(relative_path IN LISTS LAPACK_LIBRARIES)
      get_filename_component(absolute_path "${relative_path}" ABSOLUTE)
      list(APPEND ABSOLUTE_LAPACK_LIBRARIES "${absolute_path}")
   endforeach()
   list(APPEND LIBRARIES ${ABSOLUTE_LAPACK_LIBRARIES})
else()
   find_package(LAPACK REQUIRED)
   list(APPEND LIBRARIES LAPACK::LAPACK)
endif()

# BLAS
# convert relative paths of BLAS to absolute paths
if(BLAS_LIBRARIES)
   set(ABSOLUTE_BLAS_LIBRARIES)
   foreach(relative_path IN LISTS BLAS_LIBRARIES)
      get_filename_component(absolute_path "${relative_path}" ABSOLUTE)
      list(APPEND ABSOLUTE_BLAS_LIBRARIES "${absolute_path}")
   endforeach()
   list(APPEND LIBRARIES ${ABSOLUTE_BLAS_LIBRARIES})
else()
   find_package(BLAS REQUIRED)
   list(APPEND LIBRARIES BLAS::BLAS)
endif()

# add the auxiliary libraries last
if(AUXILIARY_LIBRARIES)
   set(ABSOLUTE_AUXILIARY_LIBRARIES)
   foreach(relative_path IN LISTS AUXILIARY_LIBRARIES)
      get_filename_component(absolute_path "${relative_path}" ABSOLUTE)
      list(APPEND ABSOLUTE_AUXILIARY_LIBRARIES "${absolute_path}")
   endforeach()
   list(APPEND LIBRARIES ${ABSOLUTE_AUXILIARY_LIBRARIES})
endif()

# example
add_executable(example_hs015 example_hs015.c)
target_link_libraries(example_hs015 PUBLIC ${LIBRARIES})