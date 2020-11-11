cmake_minimum_required(VERSION 3.7)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
	cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# define the project name
project(Argonot VERSION 1.0
				DESCRIPTION "Unifying Framework for Optimization" 
				LANGUAGES CXX C Fortran)

# set C++17 and enable languages
set (CMAKE_CXX_STANDARD 17)
enable_language(CXX C Fortran)

# directories
set(DIRECTORIES base base/mechanism base/strategy base/subproblem interfaces interfaces/AMPL logger solvers solvers/linear solvers/QP solvers/NLP)

# source files
set (SRC_FILES
    base/Argonot.cpp base/Constraint.cpp base/HessianEvaluation.cpp
    base/Iterate.cpp base/Matrix.cpp base/Phase.hpp base/Problem.cpp
    base/Statistics.cpp base/Utils.cpp
    
    base/mechanism/GlobalizationMechanism.cpp
    base/mechanism/GlobalizationMechanismFactory.cpp
    base/mechanism/LineSearch.cpp base/mechanism/TrustRegion.cpp
    
    base/strategy/Filter.cpp base/strategy/FilterStrategy.cpp
    base/strategy/GlobalizationStrategy.cpp 
    base/strategy/GlobalizationStrategyFactory.cpp
    base/strategy/PenaltyMeritFunction.cpp 
    
    base/subproblem/ActiveSetMethod.cpp base/subproblem/Direction.cpp
    base/subproblem/InteriorPoint.cpp base/subproblem/Sl1QP.cpp
    base/subproblem/SLP.cpp base/subproblem/SLPEQP.cpp
    base/subproblem/SQP.cpp base/subproblem/Subproblem.cpp
    base/subproblem/SubproblemFactory.cpp

    interfaces/AMPL/AMPLModel.cpp

    logger/Logger.hpp

    solvers/linear/LinearSolver.hpp
    solvers/linear/LinearSolverFactory.cpp
    solvers/linear/MA57Solver.cpp
    solvers/QP/BQPDSolver.cpp
	solvers/QP/QPSolverFactory.cpp
    solvers/QP/LPSolver.hpp solvers/QP/QPSolver.hpp
    solvers/QP/wdotd.f
)

###########
# Argonot #
###########
add_executable(argonot ${SRC_FILES} main.cpp)
target_include_directories(argonot PUBLIC ${DIRECTORIES})
# TODO: create a library with ${SRC_FILES} and an executable with main.cpp
#install(TARGETS argonot DESTINATION bin)
#add_executable (argonot_ampl main.cpp)

if(MSVC)
  target_compile_options(argonot PRIVATE /W4 /WX)
else()
  target_compile_options(argonot PRIVATE -Wall -pedantic)
endif()

# find libraries
set(LIBS "")
set(OPTIONAL_LIBS amplsolver ma57 metis bqpd CACHE STRING "Optional libraries")
set(REQUIRED_LIBS blas lapack dl)

# automatic detection of optional libraries
foreach(library_name IN LISTS OPTIONAL_LIBS)
	find_library(${library_name} ${library_name})
	# ${library_name} is the name of the library (e.g. ma57)
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
	endif()
endforeach(library_name)

# automatic detection of required libraries
foreach(library_name IN LISTS REQUIRED_LIBS)
	find_library(${library_name} ${library_name})
	# ${library_name} is the name of the library (e.g. blas)
	# ${${library_name}} is the path of the library if found, otherwise ${library_name}-NOTFOUND
	if(${${library_name}} STREQUAL "${library_name}-NOTFOUND")
		message(FATAL_ERROR "Required library ${library_name} was not found. Use ccmake to configure its path.")
	else()
		# add the library
		list(APPEND LIBS ${${library_name}})
	endif()
endforeach(library_name)

# link the libraries
target_link_libraries(argonot PUBLIC ${LIBS})

####################
# gtest unit tests #
####################
find_package(GTest CONFIG REQUIRED)
if(NOT ${GTest}_DIR STREQUAL "${GTest}-NOTFOUND")
	add_executable(run_unit_tests unit_tests/tests.cpp)
	target_include_directories(run_unit_tests PUBLIC ${DIRECTORIES})
	target_link_libraries(run_unit_tests PUBLIC GTest::gtest ${LIBS})
endif()
