cmake_minimum_required(VERSION 3.18)
project(unopy LANGUAGES CXX C)

# set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

set(LIBRARIES "")
set(DIRECTORIES ../../uno ${CMAKE_BINARY_DIR}/include)

# Fail early if BQPD is not provided
if(NOT DEFINED BQPD)
    message(FATAL_ERROR "BQPD must be provided to CMake.")
endif()
find_library(BQPD REQUIRED)
message(STATUS "Found BQPD: ${BQPD}")
list(APPEND LIBRARIES ${BQPD})
add_definitions("-D HAS_BQPD")

# Fail early if MUMPS is not provided
find_package(MUMPS REQUIRED)
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
	list(APPEND LIBRARIES ${MUMPS_MPISEQ_LIBRARY})
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
message(STATUS "Found MUMPS: ${MUMPS_LIBRARY}")
add_definitions("-D HAS_MUMPS")
list(APPEND LIBRARIES ${BQPD} ${MUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${MUMPS_PORD_LIBRARY})

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
   ../../uno/ingredients/subproblem_solvers/MUMPS/*.cpp
   ../../uno/model/*.cpp
   ../../uno/optimization/*.cpp
   ../../uno/options/*.cpp
   ../../uno/tools/*.cpp
)
add_library(compiled_uno_files OBJECT ${UNO_SOURCE_FILES})
set_property(TARGET compiled_uno_files PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(compiled_uno_files PUBLIC ${DIRECTORIES})
target_link_libraries(compiled_uno_files PUBLIC ${LIBRARIES} ${FORTRAN_LIBS})

file(GLOB PYTHON_MODULE_FILES
	cpp_classes/*.cpp
	python_modules/*.cpp
	unopy.cpp
)

add_definitions("-D PYBIND11_DETAILED_ERROR_MESSAGES")
pybind11_add_module(unopy SHARED $<TARGET_OBJECTS:compiled_uno_files> ${PYTHON_MODULE_FILES})

target_link_libraries(unopy PUBLIC ${LIBRARIES} ${FORTRAN_LIBS})
target_include_directories(unopy PUBLIC ${DIRECTORIES})

install(TARGETS unopy LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})