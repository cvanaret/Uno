# Copyright (c) 2018-2024 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

cmake_minimum_required(VERSION 3.7)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# define the project name
project(Uno VERSION 1.0
        DESCRIPTION "Uno (Unifying Nonconvex Optimization)"
        LANGUAGES CXX C Fortran)
set(CMAKE_CXX_STANDARD 17)

include(GNUInstallDirs)
set(MKL_LINK static)
set(MKL_THREADING sequential)
set(MKL_INTERFACE lp64)
find_package(MKL REQUIRED)

include($ENV{CMAKE_INSTALL_PREFIX}/cmake/GKlibTargets.cmake)
include($ENV{CMAKE_INSTALL_PREFIX}/cmake/metisTargets.cmake)
include($ENV{CMAKE_INSTALL_PREFIX}/cmake/AMPLTargets.cmake)
include($ENV{CMAKE_INSTALL_PREFIX}/cmake/bqpdTargets.cmake)
include($ENV{CMAKE_INSTALL_PREFIX}/cmake/mumpsTargets.cmake)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)


FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/include/fortran_interface.h
                            MACRO_NAMESPACE "FC_"
                            SYMBOL_NAMESPACE "FC_")

# directories
set(DIRECTORIES uno ${CMAKE_BINARY_DIR}/include)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake-library/finders)

# optional Gtest
option(WITH_GTEST "Enable GoogleTest" OFF)
message(STATUS "GoogleTest: WITH_GTEST=${WITH_GTEST}")


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
   uno/ingredients/subproblems/*.cpp
   uno/ingredients/subproblems/inequality_constrained_methods/*.cpp
   uno/ingredients/subproblems/interior_point_methods/*.cpp
   uno/model/*.cpp
   uno/optimization/*.cpp
   uno/options/*.cpp
   uno/preprocessing/*.cpp
   uno/tools/*.cpp
)

# unit test source files
file(GLOB TESTS_UNO_SOURCE_FILES
   unotest/unotest.cpp
   unotest/CollectionAdapterTests.cpp
   unotest/ConcatenationTests.cpp
   unotest/COOSparseStorageTests.cpp
   unotest/CSCSparseStorageTests.cpp
   unotest/MatrixVectorProductTests.cpp
   unotest/RangeTests.cpp
   unotest/ScalarMultipleTests.cpp
   unotest/SparseVectorTests.cpp
   unotest/SumTests.cpp
   unotest/VectorTests.cpp
   unotest/VectorViewTests.cpp
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
   message(STATUS "Library ${library_name} was found.")
endfunction()

# HSL or MA57
find_library(HSL hsl)
if(HSL)
   link_to_uno(hsl ${HSL})
else()
   message(WARNING "Optional library HSL was not found.")
   find_library(MA57 ma57)
   if(MA57)
      link_to_uno(ma57 ${MA57})
   else()
      message(WARNING "Optional library MA57 was not found.")
   endif()
endif()
if(HSL OR MA57)
   list(APPEND UNO_SOURCE_FILES uno/solvers/MA57/MA57Solver.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/MA57SolverTests.cpp)
   
   find_package(BLAS REQUIRED)
   list(APPEND LIBRARIES ${BLAS_LIBRARIES})
endif()

# find_library(metis metis)
# if(${metis} STREQUAL "metis-NOTFOUND")
#    message(WARNING "Optional library METIS was not found.")
# else()
   # link_to_uno(metis ${metis})
   set(LIBRARIES ${LIBRARIES} metis::metis)
   add_definitions("-D HAS_METIS")
# endif()

# find_library(bqpd bqpd)
# if(${bqpd} STREQUAL "bqpd-NOTFOUND")
#    message(WARNING "Optional library BQPD was not found.")
# else()
   list(APPEND UNO_SOURCE_FILES uno/solvers/BQPD/BQPDSolver.cpp uno/solvers/BQPD/wdotd.f)
   add_definitions("-D HAS_BQPD")
   set(LIBRARIES ${LIBRARIES} bqpd::bqpd_sparse)
#    link_to_uno(bqpd ${bqpd})
# endif()

# find_package(MUMPS)
# if(${MUMPS_LIBRARY} STREQUAL "MUMPS_LIBRARY-NOTFOUND")
#    message(WARNING "Optional library MUMPS was not found.")
# else()
   list(APPEND UNO_SOURCE_FILES uno/solvers/MUMPS/MUMPSSolver.cpp)
   list(APPEND TESTS_UNO_SOURCE_FILES unotest/MUMPSSolverTests.cpp)
   list(APPEND LIBRARIES mumps::libdmumps)
   
   # list(APPEND DIRECTORIES ${MUMPS_INCLUDE_DIR})
   

   # if(${MUMPS_MPISEQ_LIBRARY} STREQUAL "MUMPS_MPISEQ_LIBRARY-NOTFOUND")
   #    # parallel
   #    add_definitions("-D MUMPS_PARALLEL")
   #    find_package(MPI REQUIRED)
   #    list(APPEND LIBRARIES MPI::MPI_CXX MPI::MPI_Fortran)
   #    add_definitions("-D HAS_MPI")

   #    find_package(BLACS REQUIRED)
   #    list(APPEND LIBRARIES ${BLACS_LIBRARY})
   #    list(APPEND DIRECTORIES ${BLACS_INCLUDE_DIRS})

   #    find_package(ScaLAPACK REQUIRED)
   #    list(APPEND LIBRARIES ${ScaLAPACK_LIBRARY})
   #    list(APPEND DIRECTORIES ${ScaLAPACK_INCLUDE_DIRS})
   # else()
      # sequential
      add_definitions("-D MUMPS_SEQUENTIAL")
      # link dummy parallel library (provided by MUMPS distribution)
      # link_to_uno(MUMPS_MPISEQ_LIBRARY ${MUMPS_MPISEQ_LIBRARY})
   # endif()
   
   # find_package(METIS REQUIRED)
   list(APPEND LIBRARIES metis::metis)
   # list(APPEND DIRECTORIES ${METIS_INCLUDE_DIRS})
   
   # find_package(BLAS REQUIRED)
   # list(APPEND LIBRARIES ${BLAS_LIBRARIES})
   
   # include(FindOpenMP)
   # if (WIN32)
   # find_package(OpenMP REQUIRED)
   # else()
   # find_package(OpenMP REQUIRED C)
   # list(APPEND LIBRARIES OpenMP::OpenMP_C)
   # endif()

   add_definitions("-D HAS_MUMPS")
# endif()

###############
# Uno library #
###############
add_library(uno STATIC ${UNO_SOURCE_FILES})
set_target_properties(uno PROPERTIES 
      POSITION_INDEPENDENT_CODE ON
      MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
target_include_directories(uno PUBLIC ${DIRECTORIES})
target_link_libraries(uno PUBLIC ${LIBRARIES})

if (MSVC)
   target_compile_options(uno PUBLIC /utf-8)
endif(MSVC)

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(uno PUBLIC 
                              $<$<COMPILE_LANGUAGE:CXX>:/W4>
                              $<$<COMPILE_LANGUAGE:Fortran>:-Wall -Wextra -Wnon-virtual-dtor -pedantic -Wunused-value -Wconversion>)
else()
   target_compile_options(uno PUBLIC -Wall -Wextra -Wnon-virtual-dtor -pedantic -Wunused-value -Wconversion)
   if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmaybe-uninitialized")
   endif()
   # set(CMAKE_CXX_FLAGS_DEBUG "-pg")
   # set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG") # disable asserts
endif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
   add_compile_options(-fsanitize=address)
   add_link_options(-fsanitize=address -static-libasan)
endif()


# copy the option file
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/uno.options DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

######################
# optional AMPL main #
######################
# find_library(amplsolver amplsolver)
# if(${amplsolver} STREQUAL "amplsolver-NOTFOUND")
#    message(WARNING "Optional library amplsolver (ASL) was not found.")
# else()
   # message(STATUS "Library amplsolver was found.")
   add_executable(uno_ampl bindings/AMPL/AMPLModel.cpp bindings/AMPL/uno_ampl.cpp)
   
   target_link_libraries(uno_ampl PUBLIC uno ampl::ampl)
   add_definitions("-D HAS_AMPLSOLVER")
   set_target_properties(uno_ampl PROPERTIES 
      POSITION_INDEPENDENT_CODE ON
      MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
   # include the corresponding directory
   # get_filename_component(directory ${amplsolver} DIRECTORY)
   # include_directories(${directory})
# endif()

##################################
# optional GoogleTest unit tests #
##################################
if(WITH_GTEST)
   find_package(GTest CONFIG REQUIRED)
   if(NOT GTest_DIR)
      message(WARNING "Optional library GTest was not found.")
   else()
      add_executable(run_unotest ${TESTS_UNO_SOURCE_FILES})
      target_link_libraries(run_unotest PUBLIC GTest::gtest uno)
   endif()
endif()

add_subdirectory(localtest)