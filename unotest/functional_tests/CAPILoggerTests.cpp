// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// The logger stream and the logger level of the C API belong to a solver: they apply to the optimizations of this
// solver only, and do not leak to other solvers or to the rest of the program

#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <type_traits>
#include "CAPITestProblem.hpp"
#include "tools/Logger.hpp"

using namespace uno;
using namespace uno::capi_test;

static_assert(std::is_same_v<decltype(uno_set_logger_stream_callback(static_cast<uno_solver*>(nullptr), nullptr, nullptr)),
   uno_int>);
static_assert(std::is_same_v<decltype(uno_reset_logger_stream(static_cast<uno_solver*>(nullptr))), uno_int>);

namespace {
   uno_int capture(const char* buffer, uno_int length, void* user_data) {
      static_cast<std::string*>(user_data)->append(buffer, static_cast<size_t>(length));
      return length;
   }

   class CAPILoggerTests: public ::testing::Test {
   protected:
      void SetUp() override {
         if (available_preset() == nullptr) {
            GTEST_SKIP() << "no subproblem solver available";
         }
         this->model = create_hs015_model();
         ASSERT_NE(this->model, nullptr);
      }

      void TearDown() override {
         uno_destroy_model(this->model);
      }

      [[nodiscard]] uno_solver* create_solver(const char* logger) const {
         uno_solver* solver = uno_create_solver();
         EXPECT_NE(solver, nullptr);
         EXPECT_EQ(uno_set_solver_preset(solver, available_preset()), UNO_OPTIONS_SUCCESS);
         EXPECT_EQ(uno_set_solver_string_option(solver, "logger", logger), UNO_OPTIONS_SUCCESS);
         return solver;
      }

      uno_model* model{nullptr};
   };
} // namespace

TEST_F(CAPILoggerTests, NullSolverIsInvalidArgument) {
   std::string output;
   EXPECT_EQ(uno_set_logger_stream_callback(nullptr, capture, &output), UNO_SOLVER_INVALID_ARGUMENT);
   EXPECT_EQ(uno_reset_logger_stream(nullptr), UNO_SOLVER_INVALID_ARGUMENT);
}

// the stream of a solver receives the output of its optimizations only
TEST_F(CAPILoggerTests, StreamBelongsToTheSolver) {
   std::string output_a;
   uno_solver* solver_a = this->create_solver("INFO");
   uno_solver* solver_b = this->create_solver("SILENT");
   ASSERT_EQ(uno_set_logger_stream_callback(solver_a, capture, &output_a), UNO_SOLVER_SUCCESS);

   // solver B does not write to the stream of solver A
   ASSERT_EQ(uno_optimize(solver_b, this->model), UNO_OPTIMIZE_SUCCESS);
   EXPECT_TRUE(output_a.empty());

   // solver A writes to its stream
   ASSERT_EQ(uno_optimize(solver_a, this->model), UNO_OPTIMIZE_SUCCESS);
   EXPECT_NE(output_a.find("Optimization status"), std::string::npos);

   // after the reset, solver A no longer writes to the callback
   ASSERT_EQ(uno_reset_logger_stream(solver_a), UNO_SOLVER_SUCCESS);
   ASSERT_EQ(uno_set_solver_string_option(solver_a, "logger", "SILENT"), UNO_OPTIONS_SUCCESS);
   output_a.clear();
   ASSERT_EQ(uno_optimize(solver_a, this->model), UNO_OPTIMIZE_SUCCESS);
   EXPECT_TRUE(output_a.empty());

   uno_destroy_solver(solver_b);
   uno_destroy_solver(solver_a);
}

// destroying a solver with a stream does not leave a dangling global stream
TEST_F(CAPILoggerTests, DestroyedSolverLeavesNoDanglingStream) {
   std::string output_a;
   uno_solver* solver_a = this->create_solver("INFO");
   ASSERT_EQ(uno_set_logger_stream_callback(solver_a, capture, &output_a), UNO_SOLVER_SUCCESS);
   uno_destroy_solver(solver_a);
   EXPECT_EQ(Logger::stream, &std::cout);
}

// the logger level and the stream are restored after the optimization
TEST_F(CAPILoggerTests, LoggerStateIsRestoredAfterTheOptimization) {
   const Level previous_level = Logger::level;
   std::ostream* previous_stream = Logger::stream;
   Logger::level = WARNING;
   std::string output;
   uno_solver* solver = this->create_solver("SILENT");
   ASSERT_EQ(uno_set_logger_stream_callback(solver, capture, &output), UNO_SOLVER_SUCCESS);
   EXPECT_EQ(Logger::stream, previous_stream);
   ASSERT_EQ(uno_optimize(solver, this->model), UNO_OPTIMIZE_SUCCESS);
   EXPECT_EQ(Logger::level, WARNING);
   EXPECT_EQ(Logger::stream, previous_stream);
   uno_destroy_solver(solver);
   Logger::level = previous_level;
}
