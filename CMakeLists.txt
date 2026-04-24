# Copyright (c) 2018-2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

cmake_minimum_required(VERSION 3.12)

# define the project
project(Uno VERSION 2.7.1
        DESCRIPTION "Uno (Unifying Nonlinear Optimization)"
        LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)

# ======================
# Options
# ======================

option(BUILD_STATIC_LIBS "Build using static libraries" ON)
option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)
option(ENABLE_TESTS "Build unit tests" OFF)

if(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
   message(FATAL_ERROR "At least one of BUILD_SHARED_LIBS or BUILD_STATIC_LIBS must be ON.")
endif()

# use Fortran compiler for name mangling
enable_language(Fortran)
include(FortranCInterface)
FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/include/fortran_interface.h
        MACRO_NAMESPACE "FC_"
        SYMBOL_NAMESPACE "FC_")

# compile options
if (MSVC)
   add_compile_options("$<$<COMPILE_LANGUAGE:C>:/utf-8>")
   add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:/utf-8>")
   # https://github.com/google/googletest/tree/main/googletest#visual-studio-dynamic-vs-static-runtimes
else()
   add_compile_options(-Wall -Wextra -Wnon-virtual-dtor -pedantic -Wunused-value -Wconversion)
   set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG") # disable asserts
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmaybe-uninitialized")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake-library/finders)

# ======================
# Core targets
# ======================

# compile the Uno source files once, whatever the compilation type is (static or shared)
add_library(compiled_uno_files OBJECT)
set_property(TARGET compiled_uno_files PROPERTY POSITION_INDEPENDENT_CODE ON)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   target_compile_definitions(compiled_uno_files PUBLIC _GLIBCXX_USE_CXX11_ABI=1) # C++ ABI
endif()

add_library(uno_dependencies INTERFACE)
target_include_directories(uno_dependencies INTERFACE
   ${CMAKE_CURRENT_SOURCE_DIR}/uno
   ${CMAKE_BINARY_DIR}/include
)

if(ENABLE_TESTS)
   add_executable(run_unotest EXCLUDE_FROM_ALL)
endif()

# ======================
# Sources
# ======================

file(GLOB UNO_SOURCE_FILES CONFIGURE_DEPENDS
   interfaces/C/Uno_C_API.cpp
   uno/Uno.cpp
   uno/ingredients/constraint_relaxation_strategies/*.cpp
   uno/ingredients/constraint_relaxation_strategies/relaxed_problems/*.cpp
   uno/ingredients/globalization_mechanisms/*.cpp
   uno/ingredients/globalization_strategies/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/filter_methods/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/filter_methods/filters/*.cpp
   uno/ingredients/globalization_strategies/switching_methods/funnel_methods/*.cpp
   uno/ingredients/hessian_models/*.cpp
   uno/ingredients/hessian_models/quasi_newton/*.cpp
   uno/ingredients/hessian_models/quasi_newton/direct/*.cpp
   uno/ingredients/hessian_models/quasi_newton/inverse/*.cpp
   uno/ingredients/inequality_handling_methods/*.cpp
   uno/ingredients/inequality_handling_methods/interior_point_methods/*.cpp
   uno/ingredients/inequality_handling_methods/interior_point_methods/barrier_problems/*.cpp
   uno/ingredients/inertia_correction_strategies/*.cpp
   uno/ingredients/subproblem/*.cpp
   uno/ingredients/subproblem_solvers/*.cpp
   uno/model/*.cpp
   uno/optimization/*.cpp
   uno/options/*.cpp
   uno/tools/*.cpp
)
target_sources(compiled_uno_files PRIVATE ${UNO_SOURCE_FILES})

file(GLOB TESTS CONFIGURE_DEPENDS
   unotest/unotest.cpp
   unotest/functional_tests/BLASTests.cpp
   unotest/functional_tests/LAPACKTests.cpp
   unotest/unit_tests/*.cpp
)
if(ENABLE_TESTS)
   target_sources(run_unotest PRIVATE ${TESTS})
endif()

# ======================
# Helper function
# ======================

# function that links an existing library to Uno
function(link_to_uno library_name library_path)
   target_link_libraries(uno_dependencies INTERFACE ${library_path})

   # add a preprocessor definition
   string(TOUPPER ${library_name} library_name_upper)
   target_compile_definitions(uno_dependencies INTERFACE HAS_${library_name_upper})

   # include the corresponding directory
   get_filename_component(directory ${library_path} DIRECTORY)
   target_include_directories(uno_dependencies SYSTEM INTERFACE ${directory})

   message(STATUS "Found ${library_name}")
endfunction()

# ======================
# Dependencies
# ======================

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
   target_sources(compiled_uno_files PRIVATE uno/ingredients/subproblem_solvers/MA57/MA57Solver.cpp)
   if(ENABLE_TESTS)
      target_sources(run_unotest PRIVATE unotest/functional_tests/MA57SolverTests.cpp)
   endif()
endif()
if(HSL OR MA27)
   target_sources(compiled_uno_files PRIVATE uno/ingredients/subproblem_solvers/MA27/MA27Solver.cpp)
   if(ENABLE_TESTS)
      target_sources(run_unotest PRIVATE unotest/functional_tests/MA27SolverTests.cpp)
   endif()
endif()

# METIS
find_library(METIS metis)
if(METIS)
   link_to_uno(metis ${METIS})
endif()

# BQPD
find_library(BQPD bqpd)
if(BQPD)
   link_to_uno(bqpd ${BQPD})
   target_sources(compiled_uno_files PRIVATE
           uno/ingredients/subproblem_solvers/BQPD/BQPDSolver.cpp
           uno/ingredients/subproblem_solvers/BQPD/BQPDWorkspace.cpp)
   if(ENABLE_TESTS)
      target_sources(run_unotest PRIVATE unotest/functional_tests/BQPDSolverTests.cpp)
   endif()
endif()

# HiGHS
find_library(HIGHS highs)
if(HIGHS)
   target_link_libraries(uno_dependencies INTERFACE ${HIGHS})
   target_compile_definitions(uno_dependencies INTERFACE HAS_HIGHS)

   get_filename_component(HIGHS_DIRECTORY ${HIGHS} DIRECTORY)
   target_include_directories(uno_dependencies SYSTEM INTERFACE ${HIGHS_DIRECTORY}/../include/highs)

   target_sources(compiled_uno_files PRIVATE
           uno/ingredients/subproblem_solvers/HiGHS/HiGHSSolver.cpp
           uno/ingredients/subproblem_solvers/HiGHS/HiGHSWorkspace.cpp)
   if(ENABLE_TESTS)
      target_sources(run_unotest PRIVATE unotest/functional_tests/HiGHSSolverTests.cpp)
   endif()

   message(STATUS "Found HIGHS")
endif()

# MUMPS
find_package(MUMPS)
if(MUMPS_LIBRARY)
   target_link_libraries(uno_dependencies INTERFACE ${MUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${MUMPS_PORD_LIBRARY})
   target_compile_definitions(uno_dependencies INTERFACE HAS_MUMPS)

   target_sources(compiled_uno_files PRIVATE uno/ingredients/subproblem_solvers/MUMPS/MUMPSSolver.cpp)
   target_include_directories(uno_dependencies SYSTEM INTERFACE ${MUMPS_INCLUDE_DIR})
   if(ENABLE_TESTS)
      target_sources(run_unotest PRIVATE unotest/functional_tests/MUMPSSolverTests.cpp)
   endif()

   if(NOT MUMPS_MPISEQ_LIBRARY) # parallel
      find_package(MPI REQUIRED)
      target_link_libraries(uno_dependencies INTERFACE MPI::MPI_CXX MPI::MPI_Fortran)
      target_compile_definitions(uno_dependencies INTERFACE HAS_MPI MUMPS_PARALLEL)

      find_package(BLACS REQUIRED)
      target_link_libraries(uno_dependencies INTERFACE ${BLACS_LIBRARY})
      target_include_directories(uno_dependencies SYSTEM INTERFACE ${BLACS_INCLUDE_DIRS})

      find_package(ScaLAPACK REQUIRED)
      target_link_libraries(uno_dependencies INTERFACE ${ScaLAPACK_LIBRARY})
      target_include_directories(uno_dependencies SYSTEM INTERFACE ${ScaLAPACK_INCLUDE_DIRS})
   else()
      # link dummy parallel library (provided by MUMPS distribution)
      link_to_uno(MUMPS_MPISEQ_LIBRARY ${MUMPS_MPISEQ_LIBRARY})
   endif()
   
   find_package(METIS REQUIRED)
   target_link_libraries(uno_dependencies INTERFACE ${METIS_LIBRARY})
   target_include_directories(uno_dependencies SYSTEM INTERFACE ${METIS_INCLUDE_DIRS})

   find_package(OpenMP REQUIRED)
   target_link_libraries(uno_dependencies INTERFACE OpenMP::OpenMP_CXX)

   message(STATUS "Found MUMPS")
endif()

# SPRAL SSIDS
find_library(SPRAL NAMES spral)
if(SPRAL)
   target_sources(compiled_uno_files PRIVATE uno/ingredients/subproblem_solvers/SSIDS/SSIDSSolver.cpp)
   if(ENABLE_TESTS)
      target_sources(run_unotest PRIVATE unotest/functional_tests/SSIDSSolverTests.cpp)
   endif()
   link_to_uno(spral ${SPRAL})
endif()

# helper function: resolve a mixed list of -lflags and relative paths to absolute paths
function(resolve_library_paths input_list output_list)
   set(result)
   foreach(entry IN LISTS ${input_list})
      if(entry MATCHES "^-l")
         list(APPEND result "${entry}")
      else()
         get_filename_component(absolute_entry "${entry}" ABSOLUTE)
         list(APPEND result "${absolute_entry}")
      endif()
   endforeach()
   set(${output_list} "${result}" PARENT_SCOPE)
endfunction()

# LAPACK
if(LAPACK_LIBRARIES)
   resolve_library_paths(LAPACK_LIBRARIES ABSOLUTE_LAPACK_LIBRARIES)
   target_link_libraries(uno_dependencies INTERFACE ${ABSOLUTE_LAPACK_LIBRARIES})
else()
   find_package(LAPACK REQUIRED)
   target_link_libraries(uno_dependencies INTERFACE LAPACK::LAPACK)
endif()

# BLAS
if(BLAS_LIBRARIES)
   resolve_library_paths(BLAS_LIBRARIES ABSOLUTE_BLAS_LIBRARIES)
   target_link_libraries(uno_dependencies INTERFACE ${ABSOLUTE_BLAS_LIBRARIES})
else()
   find_package(BLAS REQUIRED)
   target_link_libraries(uno_dependencies INTERFACE BLAS::BLAS)
endif()

# additional libraries to link against
set(AUXILIARY_LIBRARIES "" CACHE STRING "Additional libraries to link against")
if(AUXILIARY_LIBRARIES)
   resolve_library_paths(AUXILIARY_LIBRARIES ABSOLUTE_AUXILIARY_LIBRARIES)
   target_link_libraries(uno_dependencies INTERFACE ${ABSOLUTE_AUXILIARY_LIBRARIES})
endif()

# ======================
# Link core
# ======================

# Fortran runtime (filtered for Apple architectures) even when static libraries are linked
set(FORTRAN_LIBS ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
if(APPLE)
   list(REMOVE_ITEM FORTRAN_LIBS gcc)
   list(REMOVE_ITEM FORTRAN_LIBS emutls_w)
   list(REMOVE_ITEM FORTRAN_LIBS heapt_w)
endif()
target_link_libraries(uno_dependencies INTERFACE ${FORTRAN_LIBS})

target_link_libraries(compiled_uno_files PUBLIC uno_dependencies)

# ======================
# Final libraries
# ======================

set(UNO_PUBLIC_HEADERS
   interfaces/C/Uno_C_API.h
   interfaces/C/uno_int.h
   interfaces/Fortran/uno_c.f90
   interfaces/Fortran/uno_fortran.f90
)

if(BUILD_STATIC_LIBS)
   add_library(uno_static STATIC $<TARGET_OBJECTS:compiled_uno_files>)
   target_link_libraries(uno_static PUBLIC uno_dependencies)
   set_target_properties(uno_static PROPERTIES OUTPUT_NAME "uno" POSITION_INDEPENDENT_CODE ON)
   set_target_properties(uno_static PROPERTIES PUBLIC_HEADER "${UNO_PUBLIC_HEADERS}")
endif()

if(BUILD_SHARED_LIBS)
   add_library(uno_shared SHARED $<TARGET_OBJECTS:compiled_uno_files>)
   target_link_libraries(uno_shared PUBLIC uno_dependencies)
   set_target_properties(uno_shared PROPERTIES OUTPUT_NAME "uno")
   set_target_properties(uno_shared PROPERTIES PUBLIC_HEADER "${UNO_PUBLIC_HEADERS}")
   if (MSVC)
      set_target_properties(uno_shared PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
   endif()
endif()

if(BUILD_STATIC_LIBS)
   set(DEFAULT_UNO_LIB uno_static)
else()
   set(DEFAULT_UNO_LIB uno_shared)
endif()

# ======================
# AMPL executable
# ======================

find_library(AMPLSOLVER amplsolver)
if(AMPLSOLVER)
   message(STATUS "Found amplsolver")

   add_executable(uno_ampl EXCLUDE_FROM_ALL interfaces/AMPL/AMPLModel.cpp interfaces/AMPL/uno_ampl.cpp)
   target_link_libraries(uno_ampl PUBLIC ${DEFAULT_UNO_LIB} ${AMPLSOLVER} uno_dependencies ${CMAKE_DL_LIBS})
   target_compile_definitions(uno_ampl PRIVATE HAS_AMPLSOLVER)

   # include the corresponding directory
   get_filename_component(directory ${AMPLSOLVER} DIRECTORY)
   target_include_directories(uno_ampl PRIVATE ${directory})
endif()

# =========================
# unopy (pybind11 bindings)
# =========================

if(SKBUILD)
   find_package(pybind11 CONFIG REQUIRED)
   message(STATUS "Building Python extension (SKBUILD detected).")

   file(GLOB PYTHON_MODULE_FILES
      interfaces/Python/cpp_classes/*.cpp
      interfaces/Python/python_modules/*.cpp
      interfaces/Python/unopy.cpp
   )

   add_compile_definitions(PYBIND11_DETAILED_ERROR_MESSAGES)
   pybind11_add_module(unopy SHARED $<TARGET_OBJECTS:compiled_uno_files> ${PYTHON_MODULE_FILES})

   target_link_libraries(unopy PUBLIC uno_dependencies)

   if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      # link-time optimization
      target_compile_options(unopy PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-flto=auto>)
      target_link_options(unopy PRIVATE -flto=auto)
   endif()

   install(TARGETS unopy
           LIBRARY DESTINATION unopy    # Linux/macOS: unopy/unopy.so
           RUNTIME DESTINATION unopy    # Windows: unopy/unopy.pyd
   )
endif()

# ======================
# GoogleTest unit tests
# ======================

if(ENABLE_TESTS)
   include(FetchContent)

   FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
   )
   # do not install gtest in install()
   set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
   set(INSTALL_GMOCK OFF CACHE BOOL "" FORCE)

   set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
   FetchContent_MakeAvailable(googletest)
   message(STATUS "Fetched GoogleTest")

   enable_testing()

   target_link_libraries(run_unotest PUBLIC GTest::gtest ${DEFAULT_UNO_LIB} uno_dependencies ${CMAKE_DL_LIBS})
   if (MSVC)
      set_target_properties(run_unotest PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
   endif()

   include(GoogleTest)
   gtest_discover_tests(run_unotest DISCOVERY_MODE PRE_TEST)
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
include(GNUInstallDirs)  # ensures CMAKE_INSTALL_INCLUDEDIR is set to "include"
install(TARGETS ${INSTALL_TARGETS}
   COMPONENT libuno
   EXPORT UnoTargets
   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}  # DLLs on windows
   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}  # static libraries
   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}  # shared libraries
   PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/uno # headers
)