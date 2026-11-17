// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Tests of the functions of the C API that deal with options: they report errors with status codes and never throw

#include <gtest/gtest.h>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include "CAPITestProblem.hpp"
#include "options/Options.hpp"

using namespace uno;
using namespace uno::capi_test;

// the declarations are noexcept for C++ callers
static_assert(noexcept(uno_set_solver_integer_option(nullptr, nullptr, 0)));
static_assert(noexcept(uno_set_solver_double_option(nullptr, nullptr, 0.)));
static_assert(noexcept(uno_set_solver_bool_option(nullptr, nullptr, false)));
static_assert(noexcept(uno_set_solver_string_option(nullptr, nullptr, nullptr)));
static_assert(noexcept(uno_get_solver_option_type(nullptr, nullptr)));
static_assert(noexcept(uno_load_solver_option_file(nullptr, nullptr)));
static_assert(noexcept(uno_set_solver_preset(nullptr, nullptr)));
static_assert(noexcept(uno_options_begin_iterator()));
static_assert(noexcept(uno_optimize(nullptr, nullptr)));

namespace {
   // owns a C API solver
   class CAPIOptionsTests: public ::testing::Test {
   protected:
      void SetUp() override {
         this->solver = uno_create_solver();
         ASSERT_NE(this->solver, nullptr);
      }

      void TearDown() override {
         uno_destroy_solver(this->solver);
      }

      void* solver{nullptr};
   };

   // writes a temporary option file and removes it at destruction
   class TemporaryFile {
   public:
      TemporaryFile(std::string name, const std::string& content): name(std::move(name)) {
         std::ofstream file(this->name);
         file << content;
      }
      ~TemporaryFile() {
         std::remove(this->name.c_str());
      }
      TemporaryFile(const TemporaryFile&) = delete;
      TemporaryFile& operator=(const TemporaryFile&) = delete;

      const std::string name;
   };
} // namespace

TEST_F(CAPIOptionsTests, SettersSucceed) {
   EXPECT_EQ(uno_set_solver_integer_option(this->solver, "max_iterations", 50), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(uno_set_solver_double_option(this->solver, "primal_tolerance", 1e-6), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(uno_set_solver_bool_option(this->solver, "print_solution", true), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(uno_set_solver_string_option(this->solver, "hessian_model", "LBFGS"), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(uno_get_solver_integer_option(this->solver, "max_iterations"), 50);
   EXPECT_EQ(uno_get_solver_double_option(this->solver, "primal_tolerance"), 1e-6);
   EXPECT_TRUE(uno_get_solver_bool_option(this->solver, "print_solution"));
   EXPECT_STREQ(uno_get_solver_string_option(this->solver, "hessian_model"), "LBFGS");
}

TEST_F(CAPIOptionsTests, NullArgumentsAreInvalid) {
   EXPECT_EQ(uno_set_solver_integer_option(nullptr, "max_iterations", 50), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_solver_integer_option(this->solver, nullptr, 50), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_solver_double_option(nullptr, "primal_tolerance", 1e-6), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_solver_double_option(this->solver, nullptr, 1e-6), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_solver_bool_option(nullptr, "print_solution", true), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_solver_bool_option(this->solver, nullptr, true), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_solver_string_option(nullptr, "logger", "INFO"), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_solver_string_option(this->solver, nullptr, "INFO"), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_solver_string_option(this->solver, "logger", nullptr), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_solver_preset(nullptr, "ipopt"), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_set_solver_preset(this->solver, nullptr), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_load_solver_option_file(nullptr, "uno.opt"), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_load_solver_option_file(this->solver, nullptr), UNO_OPTIONS_INVALID_ARGUMENT);
}

TEST_F(CAPIOptionsTests, UnknownOptionIsNotFound) {
   EXPECT_EQ(uno_set_solver_integer_option(this->solver, "max_iteration", 50), UNO_OPTIONS_NOT_FOUND);
   EXPECT_EQ(uno_set_solver_double_option(this->solver, "primal_tol", 1e-6), UNO_OPTIONS_NOT_FOUND);
   EXPECT_EQ(uno_set_solver_bool_option(this->solver, "print_solutions", true), UNO_OPTIONS_NOT_FOUND);
   EXPECT_EQ(uno_set_solver_string_option(this->solver, "loger", "INFO"), UNO_OPTIONS_NOT_FOUND);
}

TEST_F(CAPIOptionsTests, WrongTypeIsReported) {
   EXPECT_EQ(uno_set_solver_integer_option(this->solver, "primal_tolerance", 1), UNO_OPTIONS_WRONG_TYPE);
   EXPECT_EQ(uno_set_solver_double_option(this->solver, "max_iterations", 50.), UNO_OPTIONS_WRONG_TYPE);
   EXPECT_EQ(uno_set_solver_bool_option(this->solver, "logger", true), UNO_OPTIONS_WRONG_TYPE);
   EXPECT_EQ(uno_set_solver_string_option(this->solver, "print_solution", "yes"), UNO_OPTIONS_WRONG_TYPE);
}

// regression test for #1070: the misspelled value is rejected by the setter, and the solver can still be used
TEST_F(CAPIOptionsTests, InvalidStringValueIsRejected) {
   EXPECT_EQ(uno_set_solver_string_option(this->solver, "logger", "SOLENT"), UNO_OPTIONS_INVALID_VALUE);
   EXPECT_EQ(uno_set_solver_string_option(this->solver, "hessian_model", "BFGS"), UNO_OPTIONS_INVALID_VALUE);
   EXPECT_EQ(uno_set_solver_string_option(this->solver, "linear_solver", "MA97"), UNO_OPTIONS_INVALID_VALUE);
}

TEST_F(CAPIOptionsTests, FailedSetterLeavesTheOptionsUnchanged) {
   ASSERT_EQ(uno_set_solver_string_option(this->solver, "logger", "DISCRETE"), UNO_OPTIONS_SUCCESS);
   ASSERT_EQ(uno_set_solver_string_option(this->solver, "logger", "SOLENT"), UNO_OPTIONS_INVALID_VALUE);
   EXPECT_STREQ(uno_get_solver_string_option(this->solver, "logger"), "DISCRETE");
   ASSERT_EQ(uno_set_solver_preset(this->solver, "ipopt"), UNO_OPTIONS_SUCCESS);
   ASSERT_EQ(uno_set_solver_preset(this->solver, "snopt"), UNO_OPTIONS_INVALID_VALUE);
   EXPECT_STREQ(uno_get_solver_string_option(this->solver, "preset"), "ipopt");
   EXPECT_STREQ(uno_get_solver_string_option(this->solver, "inequality_handling_method"), "interior_point");
}

TEST_F(CAPIOptionsTests, PresetCanBeSetAsAnOption) {
   EXPECT_EQ(uno_set_solver_string_option(this->solver, "preset", "filtersqp"), UNO_OPTIONS_SUCCESS);
   EXPECT_STREQ(uno_get_solver_string_option(this->solver, "globalization_mechanism"), "TR");
   EXPECT_EQ(uno_set_solver_string_option(this->solver, "preset", "ipopt"), UNO_OPTIONS_SUCCESS);
   EXPECT_STREQ(uno_get_solver_string_option(this->solver, "globalization_mechanism"), "LS");
   EXPECT_EQ(uno_set_solver_string_option(this->solver, "preset", "snopt"), UNO_OPTIONS_INVALID_VALUE);
   EXPECT_EQ(uno_set_solver_preset(this->solver, "auto"), UNO_OPTIONS_SUCCESS);
}

TEST_F(CAPIOptionsTests, OptionTypes) {
   EXPECT_EQ(uno_get_solver_option_type(this->solver, "max_iterations"), UNO_OPTION_TYPE_INTEGER);
   EXPECT_EQ(uno_get_solver_option_type(this->solver, "primal_tolerance"), UNO_OPTION_TYPE_DOUBLE);
   EXPECT_EQ(uno_get_solver_option_type(this->solver, "print_solution"), UNO_OPTION_TYPE_BOOL);
   EXPECT_EQ(uno_get_solver_option_type(this->solver, "logger"), UNO_OPTION_TYPE_STRING);
   EXPECT_EQ(uno_get_solver_option_type(this->solver, "loger"), UNO_OPTION_TYPE_NOT_FOUND);
   EXPECT_EQ(uno_get_solver_option_type(this->solver, nullptr), UNO_OPTION_TYPE_NOT_FOUND);
   // used to return false, that is UNO_OPTION_TYPE_INTEGER
   EXPECT_EQ(uno_get_solver_option_type(nullptr, "logger"), UNO_OPTION_TYPE_NOT_FOUND);
}

TEST_F(CAPIOptionsTests, OptionFileIsLoaded) {
   const TemporaryFile file("unotest_capi_options.opt", "# comment\nmax_iterations 42\nlogger DISCRETE\n");
   EXPECT_EQ(uno_load_solver_option_file(this->solver, file.name.c_str()), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(uno_get_solver_integer_option(this->solver, "max_iterations"), 42);
   EXPECT_STREQ(uno_get_solver_string_option(this->solver, "logger"), "DISCRETE");
}

TEST_F(CAPIOptionsTests, MissingOptionFileIsInvalid) {
   EXPECT_EQ(uno_load_solver_option_file(this->solver, "this_option_file_does_not_exist.opt"), UNO_OPTIONS_INVALID_FILE);
}

TEST_F(CAPIOptionsTests, MalformedOptionFileIsInvalid) {
   const TemporaryFile file("unotest_capi_malformed_options.opt", "max_iterations 10 20\n");
   EXPECT_EQ(uno_load_solver_option_file(this->solver, file.name.c_str()), UNO_OPTIONS_INVALID_FILE);
}

TEST_F(CAPIOptionsTests, OptionFileErrorsLeaveTheOptionsUnchanged) {
   ASSERT_EQ(uno_set_solver_integer_option(this->solver, "max_iterations", 7), UNO_OPTIONS_SUCCESS);
   {
      const TemporaryFile file("unotest_capi_unknown_option.opt", "max_iterations 42\nloger INFO\n");
      EXPECT_EQ(uno_load_solver_option_file(this->solver, file.name.c_str()), UNO_OPTIONS_NOT_FOUND);
   }
   {
      const TemporaryFile file("unotest_capi_invalid_value.opt", "max_iterations 42\nlogger SOLENT\n");
      EXPECT_EQ(uno_load_solver_option_file(this->solver, file.name.c_str()), UNO_OPTIONS_INVALID_VALUE);
   }
   {
      const TemporaryFile file("unotest_capi_invalid_preset.opt", "max_iterations 42\npreset snopt\n");
      EXPECT_EQ(uno_load_solver_option_file(this->solver, file.name.c_str()), UNO_OPTIONS_INVALID_VALUE);
   }
   EXPECT_EQ(uno_get_solver_integer_option(this->solver, "max_iterations"), 7);
}

TEST_F(CAPIOptionsTests, IteratorVisitsAllOptions) {
   size_t number_options = 0;
   for (uno_options_iterator iterator = uno_options_begin_iterator(); iterator != uno_options_end_iterator();
         uno_options_next_iterator(&iterator)) {
      const char* name = uno_options_iterator_name(iterator);
      ASSERT_NE(name, nullptr);
      EXPECT_EQ(uno_options_iterator_type(iterator), uno_get_solver_option_type(this->solver, name)) << name;
      ++number_options;
   }
   EXPECT_EQ(number_options, Options::option_types.size());
   // null iterators are handled
   uno_options_next_iterator(nullptr);
   uno_options_destroy_iterator(nullptr);
   EXPECT_EQ(uno_options_iterator_name(nullptr), nullptr);
   EXPECT_EQ(uno_options_iterator_type(nullptr), UNO_OPTION_TYPE_NOT_FOUND);
}

TEST(CAPIOptimize, NullArgumentsAreRejected) {
   void* solver = uno_create_solver();
   void* model = create_hs015_model();
   EXPECT_FALSE(uno_optimize(nullptr, model));
   EXPECT_FALSE(uno_optimize(solver, nullptr));
   EXPECT_FALSE(uno_optimize(nullptr, nullptr));
   uno_destroy_model(model);
   uno_destroy_solver(solver);
}

TEST(CAPIOptimize, ModelWithoutObjectiveAndConstraintsIsRejected) {
   void* solver = uno_create_solver();
   void* model = uno_create_unconstrained_model(UNO_PROBLEM_NONLINEAR, 2, UNO_ZERO_BASED_INDEXING);
   EXPECT_FALSE(uno_optimize(solver, model));
   uno_destroy_model(model);
   uno_destroy_solver(solver);
}

// #1070: an invalid logger used to throw from uno_optimize. It is now rejected by the setter and the solve goes ahead
TEST(CAPIOptimize, SolveAfterRejectedOption) {
   if (available_preset() == nullptr) {
      GTEST_SKIP() << "no subproblem solver available";
   }
   void* solver = uno_create_solver();
   void* model = create_hs015_model();
   ASSERT_EQ(uno_set_solver_preset(solver, available_preset()), UNO_OPTIONS_SUCCESS);
   ASSERT_EQ(uno_set_solver_string_option(solver, "logger", "SILENT"), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(uno_set_solver_string_option(solver, "logger", "SOLENT"), UNO_OPTIONS_INVALID_VALUE);
   EXPECT_TRUE(uno_optimize(solver, model));
   EXPECT_EQ(uno_get_optimization_status(solver), UNO_SUCCESS);
   EXPECT_EQ(uno_get_solution_status(solver), UNO_FEASIBLE_KKT_POINT);
   uno_destroy_model(model);
   uno_destroy_solver(solver);
}
