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
option(ENABLE_LTO "Enable link-time optimization" ON)

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

file(GLOB TESTS CONFIGURE_DEPENDS
   unotest/unotest.cpp
   unotest/functional_tests/BLASTests.cpp
   unotest/functional_tests/LAPACKTests.cpp
   unotest/unit_tests/CollectionAdapterTests.cpp
   unotest/unit_tests/ConcatenationTests.cpp
   unotest/unit_tests/COOSparseStorageTests.cpp
   unotest/unit_tests/CSCSparseStorageTests.cpp
   unotest/unit_tests/DenseMatrixTests.cpp
   unotest/unit_tests/RangeTests.cpp
   unotest/unit_tests/ScalarMultipleTests.cpp
   unotest/unit_tests/SparseVectorTests.cpp
   unotest/unit_tests/SumTests.cpp
   unotest/unit_tests/VectorTests.cpp
   unotest/unit_tests/VectorViewTests.cpp
)

target_sources(compiled_uno_files PRIVATE ${UNO_SOURCE_FILES})

# ======================
# Dependencies
# ======================

# directories
set(DIRECTORIES uno ${CMAKE_BINARY_DIR}/include)

# list of all libraries
set(LIBRARIES "")
# additional libraries to link against
set(AUXILIARY_LIBRARIES "" CACHE STRING "Additional libraries to link against")

# function that links an existing library to Uno
function(link_to_uno library_name library_path)
   # add the library
   set(LIBRARIES ${LIBRARIES} ${library_path} PARENT_SCOPE)
   # add a preprocessor definition
   string(TOUPPER ${library_name} library_name_upper)
   add_compile_definitions(HAS_${library_name_upper})
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
   target_sources(compiled_uno_files PRIVATE uno/ingredients/subproblem_solvers/MA57/MA57Solver.cpp)
   list(APPEND TESTS unotest/functional_tests/MA57SolverTests.cpp)
endif()
if(HSL OR MA27)
   target_sources(compiled_uno_files PRIVATE uno/ingredients/subproblem_solvers/MA27/MA27Solver.cpp)
   list(APPEND TESTS unotest/functional_tests/MA27SolverTests.cpp)
endif()

# METIS
find_library(METIS metis)
if(METIS)
   link_to_uno(metis ${METIS})
endif()

# BQPD
find_library(BQPD bqpd)
if(BQPD)
   target_sources(compiled_uno_files PRIVATE
           uno/ingredients/subproblem_solvers/BQPD/BQPDSolver.cpp
           uno/ingredients/subproblem_solvers/BQPD/BQPDWorkspace.cpp)
   list(APPEND TESTS unotest/functional_tests/BQPDSolverTests.cpp)
   link_to_uno(bqpd ${BQPD})
endif()

# HiGHS
find_library(HIGHS highs)
if(HIGHS)
   target_sources(compiled_uno_files PRIVATE
           uno/ingredients/subproblem_solvers/HiGHS/HiGHSSolver.cpp
           uno/ingredients/subproblem_solvers/HiGHS/HiGHSWorkspace.cpp)
   list(APPEND TESTS unotest/functional_tests/HiGHSSolverTests.cpp)
   # add the library
   set(LIBRARIES ${LIBRARIES} ${HIGHS})
   # add a preprocessor definition
   add_compile_definitions(HAS_HIGHS)
   # include the corresponding directory
   get_filename_component(HIGHS_DIRECTORY ${HIGHS} DIRECTORY)
   set(DIRECTORIES ${DIRECTORIES} ${HIGHS_DIRECTORY}/../include/highs)
   message(STATUS "Found HIGHS")
endif()

# MUMPS
find_package(MUMPS)
if(MUMPS_LIBRARY)
   target_sources(compiled_uno_files PRIVATE uno/ingredients/subproblem_solvers/MUMPS/MUMPSSolver.cpp)
   list(APPEND TESTS unotest/functional_tests/MUMPSSolverTests.cpp)
   list(APPEND LIBRARIES ${MUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${MUMPS_PORD_LIBRARY})

   list(APPEND DIRECTORIES ${MUMPS_INCLUDE_DIR})

   if(NOT MUMPS_MPISEQ_LIBRARY)
      # parallel
      add_compile_definitions(MUMPS_PARALLEL)
      find_package(MPI REQUIRED)
      list(APPEND LIBRARIES MPI::MPI_CXX MPI::MPI_Fortran)
      add_compile_definitions(HAS_MPI)

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

   find_package(OpenMP REQUIRED)
   list(APPEND LIBRARIES OpenMP::OpenMP_CXX)

   add_compile_definitions(HAS_MUMPS)
   message(STATUS "Found MUMPS")
endif()

# SPRAL SSIDS
find_library(SPRAL NAMES spral)
if(SPRAL)
   target_sources(compiled_uno_files PRIVATE uno/ingredients/subproblem_solvers/SSIDS/SSIDSSolver.cpp)
   list(APPEND TESTS unotest/functional_tests/SSIDSSolverTests.cpp)
   link_to_uno(spral ${SPRAL})
endif()

# LAPACK
# convert relative paths of LAPACK to absolute paths
if(LAPACK_LIBRARIES)
   set(ABSOLUTE_LAPACK_LIBRARIES)
   foreach(relative_path IN LISTS LAPACK_LIBRARIES)
      if(relative_path MATCHES "^-l")
         # pass raw linker flags through unchanged
         list(APPEND ABSOLUTE_LAPACK_LIBRARIES "${relative_path}")
      else()
         get_filename_component(absolute_path "${relative_path}" ABSOLUTE)
         list(APPEND ABSOLUTE_LAPACK_LIBRARIES "${absolute_path}")
      endif()
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
      if(relative_path MATCHES "^-l")
         # pass raw linker flags through unchanged
         list(APPEND ABSOLUTE_BLAS_LIBRARIES "${relative_path}")
      else()
         get_filename_component(absolute_path "${relative_path}" ABSOLUTE)
         list(APPEND ABSOLUTE_BLAS_LIBRARIES "${absolute_path}")
      endif()
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
      if(relative_path MATCHES "^-l")
         # pass raw linker flags through unchanged
         list(APPEND ABSOLUTE_AUXILIARY_LIBRARIES "${relative_path}")
      else()
         get_filename_component(absolute_path "${relative_path}" ABSOLUTE)
         list(APPEND ABSOLUTE_AUXILIARY_LIBRARIES "${absolute_path}")
      endif()
   endforeach()
   list(APPEND LIBRARIES ${ABSOLUTE_AUXILIARY_LIBRARIES})
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

target_include_directories(compiled_uno_files SYSTEM PUBLIC ${DIRECTORIES})
target_link_libraries(compiled_uno_files PUBLIC ${LIBRARIES} ${FORTRAN_LIBS})

if(BUILD_STATIC_LIBS)
   add_library(uno_static STATIC $<TARGET_OBJECTS:compiled_uno_files>)
   target_include_directories(uno_static SYSTEM PUBLIC ${DIRECTORIES})
   target_link_libraries(uno_static PUBLIC uno_dependencies ${LIBRARIES} ${FORTRAN_LIBS})
   set_target_properties(uno_static PROPERTIES OUTPUT_NAME "uno" POSITION_INDEPENDENT_CODE ON)
   set_target_properties(uno_static PROPERTIES PUBLIC_HEADER "interfaces/C/Uno_C_API.h;interfaces/C/uno_int.h;interfaces/Fortran/uno_c.f90;interfaces/Fortran/uno_fortran.f90")
endif()

if(BUILD_SHARED_LIBS)
   add_library(uno_shared SHARED $<TARGET_OBJECTS:compiled_uno_files>)
   target_include_directories(uno_shared SYSTEM PUBLIC ${DIRECTORIES})
   target_link_libraries(uno_shared PUBLIC uno_dependencies ${LIBRARIES} ${FORTRAN_LIBS})
   set_target_properties(uno_shared PROPERTIES OUTPUT_NAME "uno")
   set_target_properties(uno_shared PROPERTIES
      PUBLIC_HEADER "interfaces/C/Uno_C_API.h;interfaces/C/uno_int.h;interfaces/Fortran/uno_c.f90;interfaces/Fortran/uno_fortran.f90"
   )
   if (MSVC)
      set_target_properties(uno_shared PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
   endif()
endif()

if(BUILD_STATIC_LIBS)
   set(DEFAULT_UNO_LIB uno_static)
else()
   set(DEFAULT_UNO_LIB uno_shared)
endif()

######################
# optional AMPL main #
######################
find_library(AMPLSOLVER amplsolver)
if(AMPLSOLVER)
   message(STATUS "Found amplsolver")
   add_executable(uno_ampl EXCLUDE_FROM_ALL interfaces/AMPL/AMPLModel.cpp interfaces/AMPL/uno_ampl.cpp)
   target_include_directories(uno_ampl SYSTEM PUBLIC ${DIRECTORIES})
   target_link_libraries(uno_ampl PUBLIC ${DEFAULT_UNO_LIB} ${AMPLSOLVER} ${LIBRARIES} ${CMAKE_DL_LIBS} ${FORTRAN_LIBS})
   add_compile_definitions(HAS_AMPLSOLVER)
   # include the corresponding directory
   get_filename_component(directory ${AMPLSOLVER} DIRECTORY)
   target_include_directories(uno_ampl PRIVATE ${directory})
endif()

#############################
# unopy (pybind11 bindings) #
#############################
if(SKBUILD)
   message(STATUS "Building Python extension (SKBUILD detected).")

   # set(PYBIND11_FINDPYTHON ON)
   find_package(pybind11 CONFIG REQUIRED)

   file(GLOB PYTHON_MODULE_FILES
      interfaces/Python/cpp_classes/*.cpp
      interfaces/Python/python_modules/*.cpp
      interfaces/Python/unopy.cpp
   )

   add_compile_definitions(PYBIND11_DETAILED_ERROR_MESSAGES)
   pybind11_add_module(unopy SHARED $<TARGET_OBJECTS:compiled_uno_files> ${PYTHON_MODULE_FILES})

   target_link_libraries(unopy PUBLIC ${LIBRARIES} ${FORTRAN_LIBS})
   target_include_directories(unopy SYSTEM PUBLIC ${DIRECTORIES})

   if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      # link-time optimization
      target_compile_options(unopy PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-flto=auto>)
      target_link_options(unopy PRIVATE -flto=auto)
   endif()

   install(TARGETS unopy
           LIBRARY DESTINATION unopy    # Linux/macOS: unopy/unopy.so
           RUNTIME DESTINATION unopy    # Windows: unopy/unopy.pyd
   )
else()
   message(STATUS "SKBUILD not detected, skipping Python extension.")
endif()

##################################
# optional GoogleTest unit tests #
##################################
if(ENABLE_TESTS)
   include(FetchContent)

   FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
   )
   set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
   FetchContent_MakeAvailable(googletest)
   message(STATUS "Fetched GoogleTest")

   enable_testing()

   add_executable(run_unotest EXCLUDE_FROM_ALL ${TESTS})
   target_include_directories(run_unotest SYSTEM PUBLIC ${DIRECTORIES} ${GTEST_INCLUDE_DIR})
   target_link_libraries(run_unotest PUBLIC GTest::gtest ${DEFAULT_UNO_LIB} ${LIBRARIES} ${CMAKE_DL_LIBS} ${FORTRAN_LIBS})
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