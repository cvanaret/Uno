// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <optional>
#include <string>
#include <vector>
#include "ingredients/constraint_relaxation_strategies/ConstraintRelaxationStrategyFactory.hpp"
#include "ingredients/globalization_mechanisms/GlobalizationMechanismFactory.hpp"
#include "ingredients/globalization_strategies/GlobalizationStrategyFactory.hpp"
#include "ingredients/hessian_models/HessianSubproblemSolverJointFactory.hpp"
#include "ingredients/inequality_handling_methods/InequalityHandlingMethodFactory.hpp"
#include "ingredients/inertia_correction_strategies/InertiaCorrectionStrategyFactory.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"
#include "options/Presets.hpp"

using namespace uno;

namespace {
   // the string options that accept any string
   const std::vector<std::string> free_string_options{"option_file", "libhsl_path", "dump_kkt_path"};

   bool is_possible_value(const std::string& option_name, const std::string& value) {
      const std::vector<std::string> possible_values = Options::get_possible_values(option_name);
      return std::find(possible_values.begin(), possible_values.end(), value) != possible_values.end();
   }

   template <typename Container>
   void expect_possible_values(const std::string& option_name, const Container& values) {
      for (const auto& value: values) {
         EXPECT_TRUE(is_possible_value(option_name, value)) << value << " is not a possible value of " << option_name;
      }
   }

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

// regression test for #1070: a misspelled value is rejected when it is set, not when the solver runs
TEST(Options, InvalidLoggerValueIsRejectedWhenSet) {
   Options options;
   EXPECT_THROW(options.set_string("logger", "SOLENT"), InvalidOptionValue);
   EXPECT_FALSE(options.get_string_optional("logger").has_value());
   EXPECT_NO_THROW(options.set_string("logger", "SILENT"));
   EXPECT_EQ(options.get_string("logger"), "SILENT");
}

TEST(Options, InvalidValueLeavesThePreviousValue) {
   Options options;
   options.set_string("hessian_model", "LBFGS");
   EXPECT_THROW(options.set_string("hessian_model", "BFGS"), InvalidOptionValue);
   EXPECT_EQ(options.get_string("hessian_model"), "LBFGS");
}

TEST(Options, InvalidValueErrorListsThePossibleValues) {
   Options options;
   try {
      options.set_string("globalization_mechanism", "trust_region");
      FAIL() << "the value should have been rejected";
   }
   catch (const InvalidOptionValue& exception) {
      const std::string message = exception.what();
      EXPECT_NE(message.find("trust_region"), std::string::npos);
      EXPECT_NE(message.find("TR, LS"), std::string::npos);
   }
}

TEST(Options, AllPossibleValuesAreAccepted) {
   for (const auto& [option_name, option_type]: Options::option_types) {
      for (const std::string& value: Options::get_possible_values(option_name)) {
         Options options;
         EXPECT_NO_THROW(options.set_string(option_name, value)) << option_name << " = " << value;
      }
   }
}

TEST(Options, EnumeratedStringOptionsAreCaseSensitiveUnlessTheFactoryIsNot) {
   Options options;
   EXPECT_THROW(options.set_string("hessian_model", "lbfgs"), InvalidOptionValue);
   EXPECT_THROW(options.set_string("linear_solver", "ma57"), InvalidOptionValue);
   EXPECT_THROW(options.set_string("QP_solver", "bqpd"), InvalidOptionValue);
   // the factories compare HIGHS and HIPO case-insensitively
   EXPECT_NO_THROW(options.set_string("QP_solver", "highs"));
   EXPECT_NO_THROW(options.set_string("LP_solver", "HiGHS"));
   EXPECT_NO_THROW(options.set_string("linear_solver", "hipo"));
}

TEST(Options, FreeStringOptionsAcceptAnyValue) {
   for (const std::string& option_name: free_string_options) {
      Options options;
      EXPECT_TRUE(Options::get_possible_values(option_name).empty()) << option_name;
      EXPECT_NO_THROW(options.set_string(option_name, "any value/at all.txt")) << option_name;
   }
}

// a new string option must either declare its possible values or be explicitly free
TEST(Options, EveryStringOptionHasPossibleValuesOrIsFree) {
   for (const auto& [option_name, option_type]: Options::option_types) {
      if (option_type == OptionType::STRING) {
         const bool is_free = std::find(free_string_options.begin(), free_string_options.end(), option_name) !=
            free_string_options.end();
         EXPECT_NE(is_free, !Options::get_possible_values(option_name).empty()) << option_name;
      }
      else {
         EXPECT_TRUE(Options::get_possible_values(option_name).empty()) << option_name;
      }
   }
}

TEST(Options, DefaultValuesArePossibleValues) {
   Options default_options;
   EXPECT_NO_THROW(DefaultOptions::load(default_options));
   for (const auto& [option_name, option_type]: Options::option_types) {
      if (option_type == OptionType::STRING && !Options::get_possible_values(option_name).empty()) {
         const std::optional<std::string> value = default_options.get_string_optional(option_name);
         if (value.has_value()) {
            EXPECT_TRUE(is_possible_value(option_name, *value)) << option_name << " = " << *value;
         }
      }
   }
}

TEST(Options, PresetValuesArePossibleValues) {
   for (const std::string& preset: Options::get_possible_values("preset")) {
      if (preset != "auto") {
         Options options;
         EXPECT_NO_THROW(Presets::set(options, preset)) << preset;
      }
   }
}

// the possible values must contain everything the factories can create (in any build)
TEST(Options, PossibleValuesContainTheFactoryStrategies) {
   expect_possible_values("constraint_relaxation_strategy", ConstraintRelaxationStrategyFactory::available_strategies);
   expect_possible_values("globalization_mechanism", GlobalizationMechanismFactory::available_strategies);
   expect_possible_values("globalization_strategy", GlobalizationStrategyFactory::available_strategies);
   expect_possible_values("hessian_model", HessianSubproblemSolverJointFactory::available_strategies);
   expect_possible_values("inertia_correction_strategy", InertiaCorrectionStrategyFactory::available_strategies);
   expect_possible_values("inequality_handling_method", InequalityHandlingMethodFactory::available_strategies());
   expect_possible_values("QP_solver", QPSolverFactory::available_solvers);
   expect_possible_values("LP_solver", LPSolverFactory::available_solvers);
   expect_possible_values("linear_solver", SymmetricIndefiniteLinearSolverFactory::available_solvers());
}

TEST(Options, UnknownOptionIsNotFound) {
   Options options;
   EXPECT_THROW(options.set_string("loger", "SILENT"), OptionNotFound);
   EXPECT_THROW(options.set_integer("max_iteration", 10), OptionNotFound);
   EXPECT_THROW(options.set("max_iteration", "10"), OptionNotFound);
   EXPECT_THROW(static_cast<void>(options.get_option_type("max_iteration")), OptionNotFound);
}

TEST(Options, WrongTypeIsReported) {
   Options options;
   EXPECT_THROW(options.set_integer("primal_tolerance", 1), OptionWrongType);
   EXPECT_THROW(options.set_double("max_iterations", 10.), OptionWrongType);
   EXPECT_THROW(options.set_bool("logger", true), OptionWrongType);
   EXPECT_THROW(options.set_string("print_solution", "yes"), OptionWrongType);
}

TEST(Options, UnparsableValueIsInvalid) {
   Options options;
   EXPECT_THROW(options.set("max_iterations", "ten"), InvalidOptionValue);
   EXPECT_THROW(options.set("max_iterations", "99999999999999999999"), InvalidOptionValue);
   EXPECT_THROW(options.set("primal_tolerance", "small"), InvalidOptionValue);
   EXPECT_THROW(options.set("print_solution", "maybe"), InvalidOptionValue);
   EXPECT_THROW(options.set("logger", "SOLENT"), InvalidOptionValue);
   EXPECT_NO_THROW(options.set("logger", "INFO"));
}

TEST(Options, MissingOptionFileIsInvalid) {
   Options options;
   EXPECT_THROW(Options::load_option_file(options, "this_option_file_does_not_exist.opt"), InvalidOptionFile);
}

TEST(Options, MalformedOptionFileIsInvalid) {
   const TemporaryFile file("unotest_malformed_options.opt", "max_iterations 10 20\n");
   Options options;
   EXPECT_THROW(Options::load_option_file(options, file.name), InvalidOptionFile);
}

TEST(Options, OptionFileWithInvalidValueIsRejected) {
   const TemporaryFile file("unotest_invalid_value_options.opt", "max_iterations 10\nlogger SOLENT\n");
   Options options;
   EXPECT_THROW(Options::load_option_file(options, file.name), InvalidOptionValue);
}

TEST(Options, ExceptionsDeriveFromTheStandardExceptions) {
   Options options;
   EXPECT_THROW(options.set_string("loger", "SILENT"), std::out_of_range);
   EXPECT_THROW(options.set_integer("primal_tolerance", 1), std::invalid_argument);
   EXPECT_THROW(options.set_string("logger", "SOLENT"), std::invalid_argument);
   EXPECT_THROW(Options::load_option_file(options, "this_option_file_does_not_exist.opt"), std::invalid_argument);
}
