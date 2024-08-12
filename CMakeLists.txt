# Copyright (c) 2018-2024 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

cmake_minimum_required(VERSION 3.7)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# define the project name
project(Uno VERSION 1.0
			DESCRIPTION "Uno (Unifying Nonlinear Optimization)"
            LANGUAGES CXX C Fortran)

# set C++17 and enable other languages
set(CMAKE_CXX_STANDARD 17)
enable_language(CXX C Fortran)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wnon-virtual-dtor -pedantic -Wunused-value -Wconversion -Wmaybe-uninitialized")
set(CMAKE_CXX_FLAGS_DEBUG "-pg")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG") # disable asserts

# optional Gtest
option(WITH_GTEST "Enable GoogleTest" OFF)
message(STATUS "GoogleTest: WITH_GTEST=${WITH_GTEST}")

# directories
set(DIRECTORIES uno)

# source files
file(GLOB UNO_SOURCE_FILES
    uno/Uno.cpp
    uno/ingredients/globalization_mechanism/*.cpp
    uno/ingredients/globalization_strategy/*.cpp
    uno/ingredients/globalization_strategy/filter_method/*.cpp
    uno/ingredients/globalization_strategy/filter_method/filter/*.cpp
    uno/ingredients/constraint_relaxation_strategy/*.cpp
    uno/ingredients/subproblem/*.cpp
    uno/ingredients/subproblem/inequality_constrained_methods/*.cpp
    uno/ingredients/subproblem/interior_point_methods/*.cpp
    uno/model/*.cpp
    uno/optimization/*.cpp
    uno/preprocessing/*.cpp
    uno/tools/*.cpp
)

# find libraries
set(LIBRARIES "")
set(OPTIONAL_LIBRARIES amplsolver ma57 metis bqpd CACHE STRING "Optional libraries")
set(REQUIRED_LIBRARIES dl blas lapack)

# automatic detection of optional libraries
foreach(library_name IN LISTS OPTIONAL_LIBRARIES)
    if(${library_name} STREQUAL ma57)
        # Try to find either libma57 or libhsl
        find_library(ma57 NAMES ma57 hsl)
    else()
        find_library(${library_name} ${library_name})
    endif()
	# ${library_name} is the name of the library (e.g. ma57)
	# ${${library_name}} is the path of the library if found, otherwise ${library_name}-NOTFOUND
	if(${${library_name}} STREQUAL "${library_name}-NOTFOUND")
		message(WARNING "Optional library ${library_name} was not found. Use ccmake to configure its path.")
	else()
		# add the library
		list(APPEND LIBRARIES ${${library_name}})
		# add a preprocessor definition
		string(TOUPPER ${library_name} library_name_upper)
		add_definitions("-D HAS_${library_name_upper}")
		# include the corresponding directory
		get_filename_component(directory ${${library_name}} DIRECTORY)
		include_directories(${directory})
        message(STATUS "Library ${library_name} was found")

        # add the corresponding sources
        if(${library_name} STREQUAL amplsolver)
            list(APPEND UNO_SOURCE_FILES uno/interfaces/AMPL/AMPLModel.cpp)
        elseif(${library_name} STREQUAL bqpd)
            list(APPEND UNO_SOURCE_FILES uno/solvers/QP/BQPDSolver.cpp uno/solvers/QP/wdotd.f)
        elseif(${library_name} STREQUAL ma57)
            list(APPEND UNO_SOURCE_FILES uno/solvers/linear/MA57Solver.cpp)
        endif()
	endif()
endforeach(library_name)

# automatic detection of required libraries
foreach(library_name IN LISTS REQUIRED_LIBRARIES)
	find_library(${library_name} ${library_name})
	# ${library_name} is the name of the library(e.g. blas)
	# ${${library_name}} is the path of the library if found, otherwise ${library_name}-NOTFOUND
	if(${${library_name}} STREQUAL "${library_name}-NOTFOUND")
		message(FATAL_ERROR "Required library ${library_name} was not found. Use ccmake to configure its path.")
	else()
		# add the library
		list(APPEND LIBRARIES ${${library_name}})
	endif()
endforeach(library_name)
find_package(OpenMP REQUIRED)
list(APPEND LIBRARIES OpenMP::OpenMP_CXX)

#######
# Uno #
#######
add_library(uno STATIC ${UNO_SOURCE_FILES})
target_include_directories(uno PUBLIC ${DIRECTORIES})

# link the libraries
target_link_libraries(uno PUBLIC ${LIBRARIES})

# copy the option file
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/uno.options DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#############
# AMPL main #
#############
add_executable(uno_ampl uno/main.cpp)
target_link_libraries(uno_ampl PUBLIC uno)

#########################
# GoogleTest unit tests #
#########################
if(WITH_GTEST)
    find_package(GTest CONFIG REQUIRED)
    if(NOT ${GTest}_DIR STREQUAL "${GTest}-NOTFOUND")
        file(GLOB TESTS_UNO_SOURCE_FILES
            unotest/*.cpp
        )
        add_executable(run_unotest ${TESTS_UNO_SOURCE_FILES})
        target_link_libraries(run_unotest PUBLIC GTest::gtest uno)
    endif()
endif()
