cmake_minimum_required(VERSION 3.7)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# define the project name
project(Uno VERSION 1.0
			DESCRIPTION "Uno: Unifying Framework for Optimization" 
            LANGUAGES CXX C Fortran)

# set C++17 and enable other languages
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
enable_language(CXX C Fortran)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wnon-virtual-dtor -pedantic -Wunused-value -Wconversion")
set(CMAKE_CXX_FLAGS_DEBUG "-pg")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# optional Gtest
option(WITH_GTEST "Enable GoogleTest" OFF)
message(STATUS "GoogleTest: WITH_GTEST=${WITH_GTEST}")

# directories
set(DIRECTORIES uno)

# source files
set(SRC_FILES
    uno/Uno.cpp

    uno/ingredients/mechanism/GlobalizationMechanism.cpp
    uno/ingredients/mechanism/GlobalizationMechanismFactory.cpp
    uno/ingredients/mechanism/BacktrackingLineSearch.cpp
    uno/ingredients/mechanism/TrustRegion.cpp

    uno/ingredients/strategy/GlobalizationStrategyFactory.cpp
    uno/ingredients/strategy/FilterStrategy.cpp
    uno/ingredients/strategy/Filter.cpp
    uno/ingredients/strategy/l1MeritFunction.cpp

    uno/ingredients/constraint_relaxation/ConstraintRelaxationStrategy.cpp
    uno/ingredients/constraint_relaxation/ConstraintRelaxationStrategyFactory.cpp
    uno/ingredients/constraint_relaxation/l1Relaxation.cpp
    uno/ingredients/constraint_relaxation/FeasibilityRestoration.cpp
    
    uno/ingredients/subproblem/Subproblem.cpp
    uno/ingredients/subproblem/SubproblemFactory.cpp
    uno/ingredients/subproblem/Direction.cpp
    uno/ingredients/subproblem/SLP.hpp
    uno/ingredients/subproblem/SQP.hpp
    uno/ingredients/subproblem/InteriorPoint.hpp
    uno/ingredients/subproblem/HessianModel.hpp
    uno/ingredients/subproblem/PredictedReductionModel.cpp

    uno/linear_algebra/SymmetricMatrix.cpp
    uno/linear_algebra/COOSymmetricMatrix.cpp
    uno/linear_algebra/CSCSymmetricMatrix.cpp
    uno/linear_algebra/SparseVector.cpp
    uno/linear_algebra/Vector.cpp
    
    uno/optimization/Problem.cpp
    uno/optimization/Constraint.cpp
    uno/optimization/Preprocessing.cpp
    uno/optimization/Iterate.cpp

    uno/solvers/linear/LinearSolverFactory.hpp
	uno/solvers/QP/QPSolverFactory.hpp
    uno/solvers/QP/LPSolver.hpp
    uno/solvers/QP/QPSolver.hpp

    uno/tools/Logger.hpp
    uno/tools/Statistics.cpp
    uno/tools/Timer.cpp
    uno/tools/Options.cpp
)

# find libraries
set(LIBS "")
set(OPTIONAL_LIBS amplsolver ma57 pardiso metis bqpd CACHE STRING "Optional libraries")
set(REQUIRED_LIBS dl blas lapack gomp)

# automatic detection of optional libraries
foreach(library_name IN LISTS OPTIONAL_LIBS)
	find_library(${library_name} ${library_name})
	# ${library_name} is the name of the library(e.g. ma57)
	# ${${library_name}} is the path of the library if found, otherwise ${library_name}-NOTFOUND
	if(${${library_name}} STREQUAL "${library_name}-NOTFOUND")
		message(WARNING "Optional library ${library_name} was not found. Use ccmake to configure its path.")
	else()
		# add the library
		list(APPEND LIBS ${${library_name}})
		# add a preprocessor definition
		string(TOUPPER ${library_name} library_name_upper)
		add_definitions("-D HAS_${library_name_upper}")
		# include the corresponding directory
		get_filename_component(directory ${${library_name}} DIRECTORY)
		include_directories(${directory})
        message(STATUS "Library ${library_name} was found")
        
        # add the corresponding sources
        if(${library_name} STREQUAL amplsolver)
            list(APPEND SRC_FILES uno/interfaces/AMPL/AMPLModel.cpp)
        elseif(${library_name} STREQUAL bqpd)
            list(APPEND SRC_FILES uno/solvers/QP/BQPDSolver.cpp uno/solvers/QP/wdotd.f)
        elseif(${library_name} STREQUAL ma57)
            list(APPEND SRC_FILES uno/solvers/linear/MA57Solver.cpp)
        elseif(${library_name} STREQUAL pardiso)
            list(APPEND SRC_FILES uno/solvers/linear/PardisoSolver.cpp)
        endif()
	endif()
endforeach(library_name)

# automatic detection of required libraries
foreach(library_name IN LISTS REQUIRED_LIBS)
	find_library(${library_name} ${library_name})
	# ${library_name} is the name of the library(e.g. blas)
	# ${${library_name}} is the path of the library if found, otherwise ${library_name}-NOTFOUND
	if(${${library_name}} STREQUAL "${library_name}-NOTFOUND")
		message(FATAL_ERROR "Required library ${library_name} was not found. Use ccmake to configure its path.")
	else()
		# add the library
		list(APPEND LIBS ${${library_name}})
	endif()
endforeach(library_name)

#######
# Uno #
#######
add_library(uno STATIC ${SRC_FILES})
target_include_directories(uno PUBLIC ${DIRECTORIES})

# link the libraries
target_link_libraries(uno PUBLIC ${LIBS})

# copy the configuration file
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/uno.cfg DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

########
# main #
########
add_executable(uno_ampl uno/main.cpp)
target_link_libraries(uno_ampl PUBLIC uno)

#########################
# GoogleTest unit tests #
#########################
if(WITH_GTEST)
    find_package(GTest CONFIG REQUIRED)
    if(NOT ${GTest}_DIR STREQUAL "${GTest}-NOTFOUND")
        set(TESTS_SRC_FILES
            unotest/unotest.cpp
            unotest/MA57SolverTests.cpp
            unotest/CSCSymmetricMatrixTests.cpp
            unotest/PredictedReductionModelTests.cpp
            unotest/SparseVectorTests.cpp
        )
        add_executable(run_unotest ${TESTS_SRC_FILES})
        target_link_libraries(run_unotest PUBLIC GTest::gtest uno)
    endif()
endif()

#############
# filter-al #
#############
#add_executable(filter-al filter-al/filter-al.cpp)
#target_include_directories(filter-al PUBLIC ${DIRECTORIES})
#target_link_libraries(filter-al PUBLIC ${LIBS})
