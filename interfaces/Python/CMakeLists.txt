cmake_minimum_required(VERSION 3.18)
project(unopy LANGUAGES CXX C)

# set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

# Fail early if BQPD is not provided or invalid
if(NOT DEFINED BQPD)
    message(FATAL_ERROR "BQPD must be provided to CMake.")
endif()

if(NOT EXISTS "${BQPD}")
    message(FATAL_ERROR "Static library not found: ${BQPD}")
endif()

add_definitions("-D HAS_BQPD")
message(STATUS "Fortran compiler required")
enable_language(Fortran)
include(FortranCInterface)
FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/include/fortran_interface.h
		  MACRO_NAMESPACE "FC_"
		  SYMBOL_NAMESPACE "FC_")

# link-time optimization
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto=auto")
endif()

# link with Fortran libraries (filtered for Apple architectures) even when static libraries are linked
set(FORTRAN_LIBS ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
if(APPLE)
   list(REMOVE_ITEM FORTRAN_LIBS gcc)
   list(REMOVE_ITEM FORTRAN_LIBS emutls_w)
   list(REMOVE_ITEM FORTRAN_LIBS heapt_w)
endif()

# source files
file(GLOB UNO_SOURCE_FILES
   ../interfaces/C/Uno_C_API.cpp
   ../../uno/Uno.cpp
   ../../uno/ingredients/constraint_relaxation_strategies/*.cpp
   ../../uno/ingredients/constraint_relaxation_strategies/relaxed_problems/*.cpp
   ../../uno/ingredients/globalization_mechanisms/*.cpp
   ../../uno/ingredients/globalization_strategies/*.cpp
   ../../uno/ingredients/globalization_strategies/switching_methods/*.cpp
   ../../uno/ingredients/globalization_strategies/switching_methods/filter_methods/*.cpp
   ../../uno/ingredients/globalization_strategies/switching_methods/filter_methods/filters/*.cpp
   ../../uno/ingredients/globalization_strategies/switching_methods/funnel_methods/*.cpp
   ../../uno/ingredients/hessian_models/*.cpp
   ../../uno/ingredients/inequality_handling_methods/*.cpp
   ../../uno/ingredients/inequality_handling_methods/inequality_constrained_methods/*.cpp
   ../../uno/ingredients/inequality_handling_methods/interior_point_methods/*.cpp
   ../../uno/ingredients/inequality_handling_methods/interior_point_methods/barrier_problems/*.cpp
   ../../uno/ingredients/inertia_correction_strategies/*.cpp
   ../../uno/ingredients/subproblem/*.cpp
   ../../uno/ingredients/subproblem_solvers/*.cpp
   ../../uno/ingredients/subproblem_solvers/BQPD/*.cpp
   ../../uno/model/*.cpp
   ../../uno/optimization/*.cpp
   ../../uno/options/*.cpp
   ../../uno/tools/*.cpp
)
add_library(compiled_uno_files OBJECT ${UNO_SOURCE_FILES})
set_property(TARGET compiled_uno_files PROPERTY POSITION_INDEPENDENT_CODE ON)
set(DIRECTORIES ../../uno ${CMAKE_BINARY_DIR}/include)
target_include_directories(compiled_uno_files PUBLIC ${DIRECTORIES})
target_link_libraries(compiled_uno_files PUBLIC ${BQPD} ${FORTRAN_LIBS})

file(GLOB PYTHON_MODULE_FILES
	cpp_classes/*.cpp
	python_modules/*.cpp
	unopy.cpp
)

add_definitions("-D PYBIND11_DETAILED_ERROR_MESSAGES")
pybind11_add_module(unopy SHARED $<TARGET_OBJECTS:compiled_uno_files> ${PYTHON_MODULE_FILES})

target_link_libraries(unopy PUBLIC ${LIBRARIES} ${FORTRAN_LIBS})
target_include_directories(unopy PUBLIC ${DIRECTORIES})

if(SKBUILD)
	install(TARGETS unopy LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()