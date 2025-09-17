# Copyright (c) 2018-2025 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

cmake_minimum_required(VERSION 3.7)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
   cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

######################
# project definition #
######################

# define the project
project(Uno VERSION 2.1.0
        DESCRIPTION "Uno (Unifying Nonlinear Optimization)"
        LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)

# compile options
if (MSVC)
   add_compile_options("$<$<COMPILE_LANGUAGE:C>:/utf-8>")
   add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:/utf-8>")
   # https://github.com/google/googletest/tree/main/googletest#visual-studio-dynamic-vs-static-runtimes
   set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
else()
   set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wnon-virtual-dtor -pedantic -Wunused-value -Wconversion")
   set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG") # disable asserts
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmaybe-uninitialized")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake-library/finders)

# options
option(BUILD_STATIC_LIBS "Build using static libraries" ON)
option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)
if(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
   message(FATAL_ERROR "At least one of BUILD_SHARED_LIBS or BUILD_STATIC_LIBS must be ON.")
endif()

# determine whether a Fortran compiler is required, based on the available optional dependencies
find_library(HSL hsl)
find_library(MA57 ma57)
find_library(MA27 ma27)
find_library(BQPD bqpd)
if(HSL OR MA57 OR MA27 OR BQPD)
   message(STATUS "Fortran compiler required")
   enable_language(Fortran)
   include(FortranCInterface)
   FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/include/fortran_interface.h
           MACRO_NAMESPACE "FC_"
           SYMBOL_NAMESPACE "FC_")
endif()

# directories
set(DIRECTORIES uno ${CMAKE_BINARY_DIR}/include)

# source files
file(GLOB UNO_SOURCE_FILES
   uno/Uno.cpp
   uno/ingredients/constraint_relaxation_strategies/*.cpp
   uno/ingredients/globalization_mechanisms/*.cpp
   uno/ingredients/globalization_strategies/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/filter_methods/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/filter_methods/filters/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/funnel_methods/*.cpp
   uno/ingredients/hessian_models/*.cpp
   uno/ingredients/inequality_handling_methods/*.cpp
   uno/ingredients/inequality_handling_methods/inequality_constrained_methods/*.cpp
   uno/ingredients/inequality_handling_methods/interior_point_methods/*.cpp
   uno/ingredients/regularization_strategies/*.cpp
   uno/ingredients/subproblem/*.cpp
   uno/ingredients/subproblem_solvers/*.cpp
   uno/model/*.cpp
   uno/optimization/*.cpp
   uno/options/*.cpp
   uno/tools/*.cpp
   bindings/C/Uno_C_API.cpp
)

# unit test source files
file(GLOB TESTS_UNO_SOURCE_FILES
   unotest/unotest.cpp
   unotest/unit_tests/CollectionAdapterTests.cpp
   unotest/unit_tests/ConcatenationTests.cpp
   unotest/unit_tests/COOSparseStorageTests.cpp
   unotest/unit_tests/CSCSparseStorageTests.cpp
   unotest/unit_tests/RangeTests.cpp
   unotest/unit_tests/ScalarMultipleTests.cpp
   unotest/unit_tests/SparseVectorTests.cpp
   unotest/unit_tests/SumTests.cpp
   unotest/unit_tests/VectorTests.cpp
   unotest/unit_tests/VectorViewTests.cpp
)

#########################
# external dependencies #
#########################
set(LIBRARIES "")

# function that links an existing library to Uno
function(link_to_uno library_name library_path)
   # add the library
   set(LIBRARIES ${LIBRARIES} ${library_path} PARENT_SCOPE)
   # add a preprocessor definition
   string(TOUPPER ${library_name} library_name_upper)
   add_definitions("-D HAS_${library_name_upper}")
   # include the corresponding directory
   get_filename_component(directory ${library_path} DIRECTORY)
   set(DIRECTORIES ${DIRECTORIES} ${directory} PARENT_SCOPE)
   message(STATUS "Found ${library_name}")
endfunction()

# HSL or MA57
find_library(HSL hsl)
if(HSL)
   link_to_uno(hsl ${HSL})
else()
   find_library(MA57 ma57)
   if(MA57)
      link_to_uno(ma57 ${MA57})
   endif()
   find_library(MA27 ma27)
   if(MA27)
      link_to_uno(ma27 ${MA27})
   endif()
endif()
if(HSL OR MA57)
   list(APPEND UNO_SOURCE_FILES uno/ingredients/subproblem_solvers/MA57/MA57Solver.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/MA57SolverTests.cpp)

   find_package(BLAS REQUIRED)
   list(APPEND LIBRARIES ${BLAS_LIBRARIES})
endif()
if(HSL OR MA27)
   list(APPEND UNO_SOURCE_FILES uno/ingredients/subproblem_solvers/MA27/MA27Solver.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/MA27SolverTests.cpp)

   find_package(BLAS REQUIRED)
   list(APPEND LIBRARIES ${BLAS_LIBRARIES})
endif()

# METIS
find_library(METIS metis)
if(METIS)
   link_to_uno(metis ${METIS})
endif()

# BQPD
find_library(BQPD bqpd)
if(BQPD)
   list(APPEND UNO_SOURCE_FILES
           uno/ingredients/subproblem_solvers/BQPD/BQPDSolver.cpp
           uno/ingredients/subproblem_solvers/BQPD/BQPDEvaluationSpace.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/BQPDSolverTests.cpp)
   link_to_uno(bqpd ${BQPD})
endif()

# HiGHS
find_package(HIGHS)
if(HIGHS_FOUND)
   list(APPEND UNO_SOURCE_FILES
           uno/ingredients/subproblem_solvers/HiGHS/HiGHSSolver.cpp
           uno/ingredients/subproblem_solvers/HiGHS/HiGHSEvaluationSpace.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/HiGHSSolverTests.cpp)
   list(APPEND LIBRARIES highs::highs)
   add_definitions("-D HAS_HIGHS")
   message(STATUS "Found HiGHS")
endif()

# MUMPS
find_package(MUMPS)
if(MUMPS_LIBRARY)
   list(APPEND UNO_SOURCE_FILES uno/ingredients/subproblem_solvers/MUMPS/MUMPSSolver.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/functional_tests/MUMPSSolverTests.cpp)
   list(APPEND LIBRARIES ${MUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${MUMPS_PORD_LIBRARY})

   list(APPEND DIRECTORIES ${MUMPS_INCLUDE_DIR})

   if(NOT MUMPS_MPISEQ_LIBRARY)
      # parallel
      add_definitions("-D MUMPS_PARALLEL")
      find_package(MPI REQUIRED)
      list(APPEND LIBRARIES MPI::MPI_CXX MPI::MPI_Fortran)
      add_definitions("-D HAS_MPI")

      find_package(BLACS REQUIRED)
      list(APPEND LIBRARIES ${BLACS_LIBRARY})
      list(APPEND DIRECTORIES ${BLACS_INCLUDE_DIRS})

      find_package(ScaLAPACK REQUIRED)
      list(APPEND LIBRARIES ${ScaLAPACK_LIBRARY})
      list(APPEND DIRECTORIES ${ScaLAPACK_INCLUDE_DIRS})
   else()
      # link dummy parallel library (provided by MUMPS distribution)
      link_to_uno(MUMPS_MPISEQ_LIBRARY ${MUMPS_MPISEQ_LIBRARY})
   endif()

   find_package(METIS REQUIRED)
   list(APPEND LIBRARIES ${METIS_LIBRARY})
   list(APPEND DIRECTORIES ${METIS_INCLUDE_DIRS})

   find_package(LAPACK REQUIRED)
   list(APPEND LIBRARIES ${LAPACK_LIBRARIES})

   find_package(BLAS REQUIRED)
   list(APPEND LIBRARIES ${BLAS_LIBRARIES})

   find_package(OpenMP REQUIRED)
   list(APPEND LIBRARIES OpenMP::OpenMP_CXX)

   add_definitions("-D HAS_MUMPS")
   message(STATUS "Found MUMPS")
endif()

###############
# Uno library #
###############
# link with Fortran libraries (filtered for Apple architectures) even when static libraries are linked
set(FORTRAN_LIBS ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
if(APPLE)
   list(REMOVE_ITEM FORTRAN_LIBS gcc)
   list(REMOVE_ITEM FORTRAN_LIBS emutls_w)
   list(REMOVE_ITEM FORTRAN_LIBS heapt_w)
endif()

# compile the Uno source files once, whatever the compilation type is (static or shared)
add_library(compiled_uno_files OBJECT ${UNO_SOURCE_FILES})
set_property(TARGET compiled_uno_files PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(compiled_uno_files PUBLIC ${DIRECTORIES})
target_link_libraries(compiled_uno_files PUBLIC ${LIBRARIES} ${FORTRAN_LIBS})

if(BUILD_STATIC_LIBS)
   add_library(uno_static STATIC $<TARGET_OBJECTS:compiled_uno_files>)
   target_include_directories(uno_static PUBLIC ${DIRECTORIES})
   target_link_libraries(uno_static PUBLIC ${LIBRARIES} ${FORTRAN_LIBS})
   set_target_properties(uno_static PROPERTIES OUTPUT_NAME "uno" POSITION_INDEPENDENT_CODE ON)
   set_target_properties(uno_static PROPERTIES PUBLIC_HEADER "bindings/C/Uno_C_API.h")
endif()

if(BUILD_SHARED_LIBS)
   add_library(uno_shared SHARED $<TARGET_OBJECTS:compiled_uno_files>)
   target_include_directories(uno_shared PUBLIC ${DIRECTORIES})
   target_link_libraries(uno_shared PUBLIC ${LIBRARIES} ${FORTRAN_LIBS})
   set_target_properties(uno_shared PROPERTIES OUTPUT_NAME "uno")
   set_target_properties(uno_shared PROPERTIES PUBLIC_HEADER "bindings/C/Uno_C_API.h")
endif()

if(BUILD_STATIC_LIBS)
   # by default, the target to install is the static library
   set(DEFAULT_UNO_LIB uno_static)
   if (MSVC)
      add_compile_options(/MT)
   endif()
else()
   set(DEFAULT_UNO_LIB uno_shared)
endif()

######################
# optional AMPL main #
######################
find_library(AMPLSOLVER amplsolver)
if(AMPLSOLVER)
   message(STATUS "Found amplsolver")
   add_executable(uno_ampl EXCLUDE_FROM_ALL bindings/AMPL/AMPLModel.cpp bindings/AMPL/uno_ampl.cpp)
   target_include_directories(uno_ampl PUBLIC ${DIRECTORIES})
   target_link_libraries(uno_ampl PUBLIC ${DEFAULT_UNO_LIB} ${AMPLSOLVER} ${LIBRARIES} ${CMAKE_DL_LIBS} ${FORTRAN_LIBS})
   add_definitions("-D HAS_AMPLSOLVER")
   # include the corresponding directory
   get_filename_component(directory ${AMPLSOLVER} DIRECTORY)
   include_directories(${directory})
endif()

############################
# optional Python bindings #
############################

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11)
if(pybind11_DIR)
   # link-time optimization
   if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto=auto")
   endif()
   
   file(GLOB PYTHON_MODULE_FILES
      bindings/Python/cpp_classes/*.cpp
      bindings/Python/python_modules/*.cpp
      bindings/Python/unopy.cpp
   )

   add_definitions("-D PYBIND11_DETAILED_ERROR_MESSAGES")
   pybind11_add_module(unopy SHARED EXCLUDE_FROM_ALL $<TARGET_OBJECTS:compiled_uno_files> ${PYTHON_MODULE_FILES})
   target_link_libraries(unopy PUBLIC ${LIBRARIES} ${FORTRAN_LIBS})
   target_include_directories(unopy PUBLIC ${DIRECTORIES})
endif()

##################################
# optional GoogleTest unit tests #
##################################
find_package(GTest CONFIG)
if(GTest_DIR)
   message(STATUS "Found GoogleTest")
   add_executable(run_unotest EXCLUDE_FROM_ALL ${TESTS_UNO_SOURCE_FILES})
   target_include_directories(run_unotest PUBLIC ${DIRECTORIES} ${GTEST_INCLUDE_DIR})
   target_link_libraries(run_unotest PUBLIC GTest::gtest ${DEFAULT_UNO_LIB} ${LIBRARIES})
   if (MSVC)
      set_target_properties(run_unotest PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
   endif()
endif()

#########################################
# install library (and AMPL executable) #
#########################################
set(INSTALL_TARGETS "")
if(TARGET uno_static)
   list(APPEND INSTALL_TARGETS uno_static)
endif()
if(TARGET uno_shared)
   list(APPEND INSTALL_TARGETS uno_shared)
endif()
if(TARGET uno_ampl)
   list(APPEND INSTALL_TARGETS uno_ampl)
endif()
install(TARGETS ${INSTALL_TARGETS}
   COMPONENT libuno
   EXPORT UnoTargets
   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}  # DLLs on windows
   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}  # static libraries
   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}  # shared libraries
   PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # headers
)