cmake_minimum_required(VERSION 3.7)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# define the project name
project(UNO VERSION 1.0
			DESCRIPTION "Unifying Framework for Optimization" 
            LANGUAGES CXX C Fortran)

# set C++17 and enable languages
set(CMAKE_CXX_STANDARD 17)
enable_language(CXX C Fortran)

#if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
#endif()

#set(CMAKE_Fortran_FLAGS "-cpp")
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-pg")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# directories
set(DIRECTORIES base base/optimization_problem base/constraint_relaxation base/mechanism base/strategy base/subproblem interfaces interfaces/AMPL linear_algebra logger solvers solvers/linear solvers/QP solvers/NLP)

# source files
set(SRC_FILES
    base/Uno.cpp
    base/Iterate.cpp base/Preprocessing.cpp
    base/Statistics.cpp base/Timer.cpp
    
    base/optimization_problem/Problem.cpp
    base/optimization_problem/Constraint.cpp
    
    base/mechanism/GlobalizationMechanism.cpp
    base/mechanism/GlobalizationMechanismFactory.cpp
    base/mechanism/BacktrackingLineSearch.cpp
    base/mechanism/TrustRegion.cpp

    base/strategy/GlobalizationStrategyFactory.cpp
    base/strategy/FilterStrategy.cpp base/strategy/Filter.cpp
    base/strategy/l1MeritFunction.cpp

    base/constraint_relaxation/ConstraintRelaxationStrategy.cpp
    base/constraint_relaxation/ConstraintRelaxationStrategyFactory.cpp
    base/constraint_relaxation/l1Relaxation.cpp
    base/constraint_relaxation/FeasibilityRestoration.cpp
    
    base/subproblem/Subproblem.cpp
    base/subproblem/SubproblemFactory.cpp
    base/subproblem/Direction.cpp
    base/subproblem/InteriorPoint.cpp
    base/subproblem/SLP.cpp
    base/subproblem/SQP.cpp
    base/subproblem/HessianEvaluation.cpp
    
    linear_algebra/Matrix.cpp linear_algebra/Vector.cpp
    
    logger/Logger.hpp

    solvers/linear/LinearSolver.hpp
    solvers/linear/LinearSolverFactory.cpp
	solvers/QP/QPSolverFactory.cpp
    solvers/QP/LPSolver.hpp solvers/QP/QPSolver.hpp
)

# find libraries
set(LIBS "")
set(OPTIONAL_LIBS amplsolver ma57 metis bqpd CACHE STRING "Optional libraries")
set(REQUIRED_LIBS blas lapack dl)

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
        
        # add the corresponding sources
        if(${library_name} STREQUAL amplsolver)
            list(APPEND SRC_FILES interfaces/AMPL/AMPLModel.cpp)
        elseif(${library_name} STREQUAL bqpd)
            list(APPEND SRC_FILES solvers/QP/BQPDSolver.cpp solvers/QP/wdotd.f)
        elseif(${library_name} STREQUAL ma57)
            list(APPEND SRC_FILES solvers/linear/MA57Solver.cpp)
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
# UNO #
#######
add_library(uno STATIC ${SRC_FILES})
target_include_directories(uno PUBLIC ${DIRECTORIES})
# TODO: create a library with ${SRC_FILES} and an executable with main.cpp
#install(TARGETS uno DESTINATION bin)
#add_executable(uno main.cpp)

if(MSVC)
  target_compile_options(uno PRIVATE /W4 /WX)
else()
  target_compile_options(uno PRIVATE -Wall -pedantic)
endif()

# link the libraries
target_link_libraries(uno PUBLIC ${LIBS})

########
# main #
########
add_executable(uno_ampl main.cpp)
target_link_libraries(uno_ampl PUBLIC uno)

####################
# gtest unit tests #
####################
find_package(GTest CONFIG REQUIRED)
if(NOT ${GTest}_DIR STREQUAL "${GTest}-NOTFOUND")
	add_executable(run_unit_tests unit_tests/tests.cpp)
	target_include_directories(run_unit_tests PUBLIC ${DIRECTORIES})
	target_link_libraries(run_unit_tests PUBLIC GTest::gtest ${LIBS})
endif()

#############
# filter-al #
#############
#add_executable(filter-al filter-al/filter-al.cpp)
#target_include_directories(filter-al PUBLIC ${DIRECTORIES})
#target_link_libraries(filter-al PUBLIC ${LIBS})
