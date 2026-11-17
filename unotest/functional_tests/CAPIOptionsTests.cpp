// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Tests of the functions of the C API that deal with options: they report errors with status codes and never throw

#include <gtest/gtest.h>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <set>
#include <string>
#include <vector>
#include "CAPITestProblem.hpp"
#include "options/DefaultOptions.hpp"
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
static_assert(noexcept(uno_get_solver_integer_option(nullptr, nullptr, nullptr)));
static_assert(noexcept(uno_get_solver_double_option(nullptr, nullptr, nullptr)));
static_assert(noexcept(uno_get_solver_bool_option(nullptr, nullptr, nullptr)));
static_assert(noexcept(uno_get_solver_string_option(nullptr, nullptr, nullptr)));

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

      // values of options that must be available
      uno_int integer_option(const char* option_name) const {
         uno_int value{-1};
         EXPECT_EQ(uno_get_solver_integer_option(this->solver, option_name, &value), UNO_OPTIONS_SUCCESS) << option_name;
         return value;
      }

      double double_option(const char* option_name) const {
         double value{-1.};
         EXPECT_EQ(uno_get_solver_double_option(this->solver, option_name, &value), UNO_OPTIONS_SUCCESS) << option_name;
         return value;
      }

      bool bool_option(const char* option_name) const {
         bool value{false};
         EXPECT_EQ(uno_get_solver_bool_option(this->solver, option_name, &value), UNO_OPTIONS_SUCCESS) << option_name;
         return value;
      }

      std::string string_option(const char* option_name) const {
         const char* value{nullptr};
         EXPECT_EQ(uno_get_solver_string_option(this->solver, option_name, &value), UNO_OPTIONS_SUCCESS) << option_name;
         return (value != nullptr) ? std::string{value} : std::string{"<null>"};
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
   EXPECT_EQ(this->integer_option("max_iterations"), 50);
   EXPECT_EQ(this->double_option("primal_tolerance"), 1e-6);
   EXPECT_TRUE(this->bool_option("print_solution"));
   EXPECT_EQ(this->string_option("hessian_model"), "LBFGS");
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
   EXPECT_EQ(this->string_option("logger"), "DISCRETE");
   ASSERT_EQ(uno_set_solver_preset(this->solver, "ipopt"), UNO_OPTIONS_SUCCESS);
   ASSERT_EQ(uno_set_solver_preset(this->solver, "snopt"), UNO_OPTIONS_INVALID_VALUE);
   EXPECT_EQ(this->string_option("preset"), "ipopt");
   EXPECT_EQ(this->string_option("inequality_handling_method"), "interior_point");
}

TEST_F(CAPIOptionsTests, PresetCanBeSetAsAnOption) {
   EXPECT_EQ(uno_set_solver_string_option(this->solver, "preset", "filtersqp"), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(this->string_option("globalization_mechanism"), "TR");
   EXPECT_EQ(uno_set_solver_string_option(this->solver, "preset", "ipopt"), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(this->string_option("globalization_mechanism"), "LS");
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
   EXPECT_EQ(this->integer_option("max_iterations"), 42);
   EXPECT_EQ(this->string_option("logger"), "DISCRETE");
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
   EXPECT_EQ(this->integer_option("max_iterations"), 7);
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

TEST_F(CAPIOptionsTests, GettersReportUnknownOptions) {
   uno_int integer_value{-1};
   double double_value{-1.};
   bool bool_value{false};
   const char* string_value{nullptr};
   EXPECT_EQ(uno_get_solver_integer_option(this->solver, "max_iteration", &integer_value), UNO_OPTIONS_NOT_FOUND);
   EXPECT_EQ(uno_get_solver_double_option(this->solver, "primal_tol", &double_value), UNO_OPTIONS_NOT_FOUND);
   EXPECT_EQ(uno_get_solver_bool_option(this->solver, "print_solutions", &bool_value), UNO_OPTIONS_NOT_FOUND);
   EXPECT_EQ(uno_get_solver_string_option(this->solver, "loger", &string_value), UNO_OPTIONS_NOT_FOUND);
   // the values are untouched
   EXPECT_EQ(integer_value, -1);
   EXPECT_EQ(double_value, -1.);
   EXPECT_FALSE(bool_value);
   EXPECT_EQ(string_value, nullptr);
}

TEST_F(CAPIOptionsTests, GettersReportWrongTypes) {
   uno_int integer_value{-1};
   double double_value{-1.};
   bool bool_value{false};
   const char* string_value{nullptr};
   EXPECT_EQ(uno_get_solver_integer_option(this->solver, "primal_tolerance", &integer_value), UNO_OPTIONS_WRONG_TYPE);
   EXPECT_EQ(uno_get_solver_double_option(this->solver, "max_iterations", &double_value), UNO_OPTIONS_WRONG_TYPE);
   EXPECT_EQ(uno_get_solver_bool_option(this->solver, "logger", &bool_value), UNO_OPTIONS_WRONG_TYPE);
   EXPECT_EQ(uno_get_solver_string_option(this->solver, "print_solution", &string_value), UNO_OPTIONS_WRONG_TYPE);
   EXPECT_EQ(integer_value, -1);
   EXPECT_EQ(double_value, -1.);
   EXPECT_FALSE(bool_value);
   EXPECT_EQ(string_value, nullptr);
}

TEST_F(CAPIOptionsTests, GettersRejectNullArguments) {
   uno_int integer_value{-1};
   const char* string_value{nullptr};
   EXPECT_EQ(uno_get_solver_integer_option(nullptr, "max_iterations", &integer_value), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_solver_integer_option(this->solver, nullptr, &integer_value), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_solver_integer_option(this->solver, "max_iterations", nullptr), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_solver_double_option(this->solver, "primal_tolerance", nullptr), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_solver_bool_option(this->solver, "print_solution", nullptr), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_solver_string_option(nullptr, "logger", &string_value), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_solver_string_option(this->solver, nullptr, &string_value), UNO_OPTIONS_INVALID_ARGUMENT);
   EXPECT_EQ(uno_get_solver_string_option(this->solver, "logger", nullptr), UNO_OPTIONS_INVALID_ARGUMENT);
}

TEST_F(CAPIOptionsTests, GettersReturnTheDefaultValues) {
   Options default_options;
   DefaultOptions::load(default_options);
   EXPECT_EQ(this->integer_option("max_iterations"), default_options.get_int("max_iterations"));
   EXPECT_EQ(this->double_option("primal_tolerance"), default_options.get_double("primal_tolerance"));
   EXPECT_EQ(this->bool_option("print_solution"), default_options.get_bool("print_solution"));
   EXPECT_EQ(this->string_option("logger"), default_options.get_string("logger"));
   EXPECT_EQ(this->string_option("preset"), default_options.get_string("preset"));
}

// #925: the options that have no default value (preset-dependent, option_file) are reported as not set
TEST_F(CAPIOptionsTests, OptionsWithoutDefaultAreNotSet) {
   const char* string_value{nullptr};
   EXPECT_EQ(this->string_option("preset"), "auto");
   EXPECT_EQ(uno_get_solver_string_option(this->solver, "globalization_mechanism", &string_value), UNO_OPTIONS_NOT_SET);
   EXPECT_EQ(uno_get_solver_string_option(this->solver, "hessian_model", &string_value), UNO_OPTIONS_NOT_SET);
   EXPECT_EQ(uno_get_solver_string_option(this->solver, "option_file", &string_value), UNO_OPTIONS_NOT_SET);
   EXPECT_EQ(string_value, nullptr);
   // a preset sets them
   ASSERT_EQ(uno_set_solver_preset(this->solver, "filtersqp"), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(this->string_option("globalization_mechanism"), "TR");
   EXPECT_EQ(this->string_option("hessian_model"), "exact");
   // back to auto: not set again
   ASSERT_EQ(uno_set_solver_preset(this->solver, "auto"), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(uno_get_solver_string_option(this->solver, "globalization_mechanism", &string_value), UNO_OPTIONS_NOT_SET);
   // the user can set them
   ASSERT_EQ(uno_set_solver_string_option(this->solver, "globalization_mechanism", "LS"), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(this->string_option("globalization_mechanism"), "LS");
}

// #925: iterating over the options and getting their values never throws. Exactly the options without a default value
// are not set: the preset-dependent options, the options that are paths (option_file, libhsl_path, dump_kkt_path) and,
// in a build without the corresponding solvers, QP_solver, LP_solver and linear_solver
TEST_F(CAPIOptionsTests, AllOptionsCanBeQueried) {
   Options default_options;
   DefaultOptions::load(default_options);
   std::set<std::string> options_not_set{};
   for (uno_options_iterator iterator = uno_options_begin_iterator(); iterator != uno_options_end_iterator();
         uno_options_next_iterator(&iterator)) {
      const char* name = uno_options_iterator_name(iterator);
      const uno_int type = uno_options_iterator_type(iterator);
      uno_int status{UNO_OPTIONS_ERROR};
      if (type == UNO_OPTION_TYPE_INTEGER) {
         uno_int value{};
         status = uno_get_solver_integer_option(this->solver, name, &value);
      }
      else if (type == UNO_OPTION_TYPE_DOUBLE) {
         double value{};
         status = uno_get_solver_double_option(this->solver, name, &value);
      }
      else if (type == UNO_OPTION_TYPE_BOOL) {
         bool value{};
         status = uno_get_solver_bool_option(this->solver, name, &value);
      }
      else if (type == UNO_OPTION_TYPE_STRING) {
         const char* value{nullptr};
         status = uno_get_solver_string_option(this->solver, name, &value);
         if (status == UNO_OPTIONS_SUCCESS) {
            EXPECT_NE(value, nullptr) << name;
         }
      }
      else {
         ADD_FAILURE() << "unknown type for " << name;
      }
      const uno_int expected_status = default_options.contains(name) ? UNO_OPTIONS_SUCCESS : UNO_OPTIONS_NOT_SET;
      EXPECT_EQ(status, expected_status) << name;
      if (status == UNO_OPTIONS_NOT_SET) {
         options_not_set.insert(name);
      }
   }
   for (const char* preset_dependent_option: {"constraint_relaxation_strategy", "inequality_handling_method",
         "globalization_mechanism", "globalization_strategy", "hessian_model", "inertia_correction_strategy"}) {
      EXPECT_EQ(options_not_set.count(preset_dependent_option), 1) << preset_dependent_option;
   }
   EXPECT_EQ(options_not_set.count("option_file"), 1);
}

TEST_F(CAPIOptionsTests, StringOverriddenBySolverIsReportedAfterTheSolve) {
   if (available_preset() == nullptr) {
      GTEST_SKIP() << "no subproblem solver available";
   }
   // without a Hessian, the solver falls back to L-BFGS
   void* model = uno_create_unconstrained_model(UNO_PROBLEM_NONLINEAR, number_variables, UNO_ZERO_BASED_INDEXING);
   double constraints_lower_bounds[] = {1., 0.};
   double constraints_upper_bounds[] = {Inf, Inf};
   uno_int jacobian_row_indices[] = {0, 1, 0, 1};
   uno_int jacobian_column_indices[] = {0, 0, 1, 1};
   double x0[] = {-2., 1.};
   uno_set_objective(model, UNO_MINIMIZE, objective_function, objective_gradient);
   uno_set_constraints(model, number_constraints, constraint_functions, constraints_lower_bounds,
      constraints_upper_bounds, 4, jacobian_row_indices, jacobian_column_indices, jacobian);
   uno_set_variable_upper_bound(model, 0, 0.5);
   uno_set_initial_primal_iterate(model, x0);
   ASSERT_EQ(uno_set_solver_preset(this->solver, available_preset()), UNO_OPTIONS_SUCCESS);
   ASSERT_EQ(uno_set_solver_string_option(this->solver, "logger", "SILENT"), UNO_OPTIONS_SUCCESS);
   ASSERT_EQ(this->string_option("hessian_model"), "exact");
   ASSERT_TRUE(uno_optimize(this->solver, model));
   EXPECT_EQ(this->string_option("hessian_model"), "LBFGS");
   // after a change of options, the value of the user/preset/default options is reported again
   ASSERT_EQ(uno_set_solver_integer_option(this->solver, "max_iterations", 100), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(this->string_option("hessian_model"), "exact");
   uno_destroy_model(model);
}

// after a successful solve, the getters return the options used by the solve, including the values picked by the
// "auto" preset, until the options are modified or the next optimization
TEST(CAPIOptimize, AutoPresetValuesCanBeReadAfterTheSolve) {
   if (available_preset() == nullptr) {
      GTEST_SKIP() << "no subproblem solver available";
   }
   void* solver = uno_create_solver();
   void* model = create_hs015_model();
   ASSERT_EQ(uno_set_solver_string_option(solver, "logger", "SILENT"), UNO_OPTIONS_SUCCESS);
   ASSERT_EQ(uno_set_solver_integer_option(solver, "max_iterations", 123), UNO_OPTIONS_SUCCESS);
   const char* value{nullptr};
   ASSERT_EQ(uno_get_solver_string_option(solver, "preset", &value), UNO_OPTIONS_SUCCESS);
   ASSERT_STREQ(value, "auto");
   value = nullptr;
   for (const char* option_name: {"globalization_mechanism", "globalization_strategy", "inequality_handling_method"}) {
      EXPECT_EQ(uno_get_solver_string_option(solver, option_name, &value), UNO_OPTIONS_NOT_SET) << option_name;
   }
   EXPECT_EQ(value, nullptr);

   ASSERT_TRUE(uno_optimize(solver, model));
   // the values picked by the "auto" preset
   for (const char* option_name: {"globalization_mechanism", "globalization_strategy", "inequality_handling_method",
         "constraint_relaxation_strategy", "inertia_correction_strategy"}) {
      value = nullptr;
      ASSERT_EQ(uno_get_solver_string_option(solver, option_name, &value), UNO_OPTIONS_SUCCESS) << option_name;
      const std::vector<std::string> possible_values = Options::get_possible_values(option_name);
      EXPECT_NE(std::find(possible_values.begin(), possible_values.end(), std::string{value}), possible_values.end())
         << option_name << " = " << value;
   }
   // the model has a Hessian
   ASSERT_EQ(uno_get_solver_string_option(solver, "hessian_model", &value), UNO_OPTIONS_SUCCESS);
   EXPECT_STREQ(value, "exact");
   // the user and default options used by the solve
   uno_int max_iterations{};
   EXPECT_EQ(uno_get_solver_integer_option(solver, "max_iterations", &max_iterations), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(max_iterations, 123);
   Options default_options;
   DefaultOptions::load(default_options);
   bool print_solution{};
   EXPECT_EQ(uno_get_solver_bool_option(solver, "print_solution", &print_solution), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(print_solution, default_options.get_bool("print_solution"));
   // a numerical option set by the preset picked by "auto"
   double primal_tolerance{};
   EXPECT_EQ(uno_get_solver_double_option(solver, "primal_tolerance", &primal_tolerance), UNO_OPTIONS_SUCCESS);
   EXPECT_GT(primal_tolerance, 0.);
   // options that are not set by the solve either
   EXPECT_EQ(uno_get_solver_string_option(solver, "option_file", &value), UNO_OPTIONS_NOT_SET);

   // modifying an option discards the options of the solve
   ASSERT_EQ(uno_set_solver_integer_option(solver, "max_iterations", 50), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(uno_get_solver_string_option(solver, "globalization_mechanism", &value), UNO_OPTIONS_NOT_SET);
   EXPECT_EQ(uno_get_solver_integer_option(solver, "max_iterations", &max_iterations), UNO_OPTIONS_SUCCESS);
   EXPECT_EQ(max_iterations, 50);
   // so does a failed optimization
   ASSERT_TRUE(uno_optimize(solver, model));
   EXPECT_EQ(uno_get_solver_string_option(solver, "globalization_mechanism", &value), UNO_OPTIONS_SUCCESS);
   EXPECT_FALSE(uno_optimize(solver, nullptr));
   EXPECT_EQ(uno_get_solver_string_option(solver, "globalization_mechanism", &value), UNO_OPTIONS_NOT_SET);
   uno_destroy_model(model);
   uno_destroy_solver(solver);
}

// an option overridden by the solver is reported after the solve, even if it had no value (here with the "auto" preset)
TEST(CAPIOptimize, OverrideOfAnOptionWithoutValueIsReportedAfterTheSolve) {
   if (available_preset() == nullptr) {
      GTEST_SKIP() << "no subproblem solver available";
   }
   // hs015 without a Hessian: the solver falls back to L-BFGS
   void* model = uno_create_unconstrained_model(UNO_PROBLEM_NONLINEAR, number_variables, UNO_ZERO_BASED_INDEXING);
   double constraints_lower_bounds[] = {1., 0.};
   double constraints_upper_bounds[] = {Inf, Inf};
   uno_int jacobian_row_indices[] = {0, 1, 0, 1};
   uno_int jacobian_column_indices[] = {0, 0, 1, 1};
   double x0[] = {-2., 1.};
   uno_set_objective(model, UNO_MINIMIZE, objective_function, objective_gradient);
   uno_set_constraints(model, number_constraints, constraint_functions, constraints_lower_bounds,
      constraints_upper_bounds, 4, jacobian_row_indices, jacobian_column_indices, jacobian);
   uno_set_variable_upper_bound(model, 0, 0.5);
   uno_set_initial_primal_iterate(model, x0);
   void* solver = uno_create_solver();
   ASSERT_EQ(uno_set_solver_string_option(solver, "logger", "SILENT"), UNO_OPTIONS_SUCCESS);
   const char* value{nullptr};
   ASSERT_EQ(uno_get_solver_string_option(solver, "hessian_model", &value), UNO_OPTIONS_NOT_SET);
   ASSERT_TRUE(uno_optimize(solver, model));
   ASSERT_EQ(uno_get_solver_string_option(solver, "hessian_model", &value), UNO_OPTIONS_SUCCESS);
   EXPECT_STREQ(value, "LBFGS");
   uno_destroy_solver(solver);
   uno_destroy_model(model);
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
