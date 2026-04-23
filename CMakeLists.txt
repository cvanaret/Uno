cmake_minimum_required(VERSION 3.16)

project(Uno VERSION 2.7.1
        DESCRIPTION "Uno (Unifying Nonlinear Optimization)"
        LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

# ======================
# Options
# ======================

option(BUILD_STATIC_LIBS "Build static libraries" ON)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(ENABLE_TESTS "Build unit tests" OFF)
option(ENABLE_LTO "Enable link-time optimization" ON)

if(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
    message(FATAL_ERROR "At least one of BUILD_STATIC_LIBS or BUILD_SHARED_LIBS must be ON.")
endif()

# ======================
# Fortran interface
# ======================

enable_language(Fortran)
include(FortranCInterface)

FortranCInterface_HEADER(${CMAKE_BINARY_DIR}/include/fortran_interface.h
    MACRO_NAMESPACE "FC_"
    SYMBOL_NAMESPACE "FC_"
)

# ======================
# Core targets
# ======================

add_library(compiled_uno_files OBJECT)

add_library(uno_dependencies INTERFACE)

target_include_directories(uno_dependencies INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/uno
    ${CMAKE_BINARY_DIR}/include
)


# ======================
# Helper
# ======================

function(link_to_uno library_name library_path)
    string(TOUPPER ${library_name} library_name_upper)

    target_link_libraries(uno_dependencies INTERFACE ${library_path})
    target_compile_definitions(uno_dependencies INTERFACE HAS_${library_name_upper})

    get_filename_component(dir ${library_path} DIRECTORY)
    target_include_directories(uno_dependencies INTERFACE ${dir})

    message(STATUS "Found ${library_name}")
endfunction()

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
        unotest/unit_tests/*.cpp
)

target_sources(compiled_uno_files PRIVATE ${UNO_SOURCE_FILES})

# ======================
# Dependencies
# ======================

# HSL / MA57 / MA27
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
endif()
if(HSL OR MA27)
    target_sources(compiled_uno_files PRIVATE uno/ingredients/subproblem_solvers/MA27/MA27Solver.cpp)
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
        uno/ingredients/subproblem_solvers/BQPD/BQPDWorkspace.cpp
    )
endif()

# HiGHS
find_library(HIGHS highs)
if(HIGHS)
    target_link_libraries(uno_dependencies INTERFACE ${HIGHS})
    target_compile_definitions(uno_dependencies INTERFACE HAS_HIGHS)

    get_filename_component(HIGHS_DIR ${HIGHS} DIRECTORY)
    target_include_directories(uno_dependencies INTERFACE ${HIGHS_DIR}/../include/highs)

    target_sources(compiled_uno_files PRIVATE
        uno/ingredients/subproblem_solvers/HiGHS/HiGHSSolver.cpp
        uno/ingredients/subproblem_solvers/HiGHS/HiGHSWorkspace.cpp
    )

    message(STATUS "Found HIGHS")
endif()

# MUMPS
find_package(MUMPS)
if(MUMPS_LIBRARY)
    target_compile_definitions(uno_dependencies INTERFACE HAS_MUMPS)

    target_link_libraries(uno_dependencies INTERFACE
        ${MUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${MUMPS_PORD_LIBRARY}
    )
    target_include_directories(uno_dependencies INTERFACE ${MUMPS_INCLUDE_DIR})

    target_sources(compiled_uno_files PRIVATE
        uno/ingredients/subproblem_solvers/MUMPS/MUMPSSolver.cpp
    )

    if(NOT MUMPS_MPISEQ_LIBRARY)
        find_package(MPI REQUIRED)
        target_link_libraries(uno_dependencies INTERFACE MPI::MPI_CXX MPI::MPI_Fortran)
        target_compile_definitions(uno_dependencies INTERFACE HAS_MPI MUMPS_PARALLEL)

        find_package(BLACS REQUIRED)
        target_link_libraries(uno_dependencies INTERFACE ${BLACS_LIBRARY})
        target_include_directories(uno_dependencies INTERFACE ${BLACS_INCLUDE_DIRS})

        find_package(ScaLAPACK REQUIRED)
        target_link_libraries(uno_dependencies INTERFACE ${ScaLAPACK_LIBRARY})
        target_include_directories(uno_dependencies INTERFACE ${ScaLAPACK_INCLUDE_DIRS})
    else()
        link_to_uno(MUMPS_MPISEQ_LIBRARY ${MUMPS_MPISEQ_LIBRARY})
    endif()

    find_package(OpenMP REQUIRED)
    target_link_libraries(uno_dependencies INTERFACE OpenMP::OpenMP_CXX)

    message(STATUS "Found MUMPS")
endif()

# LAPACK / BLAS
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

target_link_libraries(uno_dependencies INTERFACE LAPACK::LAPACK BLAS::BLAS)

# ======================
# Link core
# ======================

target_link_libraries(compiled_uno_files PUBLIC uno_dependencies)

# Fortran runtime
set(FORTRAN_LIBS ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
target_link_libraries(compiled_uno_files PUBLIC ${FORTRAN_LIBS})

# ======================
# LTO (clean)
# ======================

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)

if(ipo_supported AND ENABLE_LTO)
    set_property(TARGET compiled_uno_files PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

# ======================
# Final libraries
# ======================

if(BUILD_STATIC_LIBS)
    add_library(uno_static STATIC $<TARGET_OBJECTS:compiled_uno_files>)
    target_link_libraries(uno_static PUBLIC uno_dependencies ${FORTRAN_LIBS})
    set_target_properties(uno_static PROPERTIES OUTPUT_NAME uno)
endif()

if(BUILD_SHARED_LIBS)
    add_library(uno_shared SHARED $<TARGET_OBJECTS:compiled_uno_files>)
    target_link_libraries(uno_shared PUBLIC uno_dependencies ${FORTRAN_LIBS})
    set_target_properties(uno_shared PROPERTIES OUTPUT_NAME uno)
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
    add_executable(uno_ampl EXCLUDE_FROM_ALL
        interfaces/AMPL/AMPLModel.cpp
        interfaces/AMPL/uno_ampl.cpp
    )

    target_link_libraries(uno_ampl
        PUBLIC ${DEFAULT_UNO_LIB} uno_dependencies ${AMPLSOLVER} ${CMAKE_DL_LIBS}
    )

    target_compile_definitions(uno_ampl PRIVATE HAS_AMPLSOLVER)

    # include the corresponding directory
    get_filename_component(directory ${AMPLSOLVER} DIRECTORY)
    target_include_directories(uno_ampl PRIVATE ${directory})
endif()

# ======================
# Python module
# ======================

if(SKBUILD)
    find_package(pybind11 CONFIG REQUIRED)

    file(GLOB PYTHON_MODULE_FILES CONFIGURE_DEPENDS
        interfaces/Python/**/*.cpp
    )

    pybind11_add_module(unopy SHARED
        $<TARGET_OBJECTS:compiled_uno_files>
        ${PYTHON_MODULE_FILES}
    )

    target_link_libraries(unopy PUBLIC uno_dependencies ${FORTRAN_LIBS})
endif()

# ======================
# Tests
# ======================

if(ENABLE_TESTS)
    include(FetchContent)

    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    )
    FetchContent_MakeAvailable(googletest)

    enable_testing()

    add_executable(run_unotest ${TESTS})

    target_link_libraries(run_unotest
        PUBLIC GTest::gtest ${DEFAULT_UNO_LIB} uno_dependencies ${CMAKE_DL_LIBS}
    )

    include(GoogleTest)
    gtest_discover_tests(run_unotest)
endif()