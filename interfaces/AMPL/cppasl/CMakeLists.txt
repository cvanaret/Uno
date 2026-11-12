cmake_minimum_required(VERSION 3.16)
project(cppasl VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

option(CPPASL_USE_LAPACK "Use LAPACK (dpotrf) for the convexity checks" ON)
option(CPPASL_BUILD_TESTS "Build the unit tests" ON)
set(CPPASL_ASL2_ROOT "" CACHE PATH "ASL source tree with a built libasl2.a (build/lib), for the benchmarks")

add_library(cppasl
   src/nl_reader.cpp
   src/nl_model.cpp
   src/model_builder.cpp
   src/function_decomposition.cpp
   src/tape.cpp
   src/problem_classification.cpp)
target_include_directories(cppasl PUBLIC include PRIVATE src)
target_compile_options(cppasl PRIVATE -Wall -Wextra)
if(CPPASL_USE_LAPACK)
   find_package(LAPACK REQUIRED)
   target_compile_definitions(cppasl PRIVATE CPPASL_USE_LAPACK)
   target_link_libraries(cppasl PRIVATE LAPACK::LAPACK)
endif()

# tools: .nl writer and the Mittelmann qcqp generator
add_executable(generate_qcqp_nl tools/generate_qcqp_nl.cpp)
target_include_directories(generate_qcqp_nl PRIVATE tools)
add_executable(generate_nlp_nl tools/generate_nlp_nl.cpp)
target_include_directories(generate_nlp_nl PRIVATE tools)

if(CPPASL_BUILD_TESTS)
   enable_testing()
   add_executable(cppasl_tests tests/test_main.cpp tests/test_reader.cpp tests/test_quadratic.cpp
      tests/test_derivatives.cpp tests/test_classification.cpp tests/test_large.cpp)
   target_include_directories(cppasl_tests PRIVATE tests tools src)
   target_link_libraries(cppasl_tests PRIVATE cppasl)
   add_test(NAME cppasl_tests COMMAND cppasl_tests)
endif()

if(CPPASL_ASL2_ROOT)
   add_library(asl2_session STATIC benchmarks/asl2_session.c)
   target_include_directories(asl2_session PRIVATE ${CPPASL_ASL2_ROOT}/src/solvers2 ${CPPASL_ASL2_ROOT}/build/include)
   target_link_libraries(asl2_session PUBLIC ${CPPASL_ASL2_ROOT}/build/lib/libasl2.a m dl)
   add_executable(compare_with_asl2 benchmarks/compare_with_asl2.cpp)
   target_link_libraries(compare_with_asl2 PRIVATE cppasl asl2_session)
endif()
