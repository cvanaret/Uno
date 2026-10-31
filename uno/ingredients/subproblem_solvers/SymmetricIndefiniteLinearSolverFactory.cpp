// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include <string>
#include "SymmetricIndefiniteLinearSolverFactory.hpp"
#include "DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ScalingStrategy.hpp"
#include "linear_algebra/Vector.hpp"
#include "options/Options.hpp"

#if defined(HAS_HSL) || defined(HAS_MA57)
#include "ingredients/subproblem_solvers/MA57/MA57Solver.hpp"
#endif

#if defined(HAS_HSL) || defined(HAS_MA27)
#include "ingredients/subproblem_solvers/MA27/MA27Solver.hpp"
#endif

#ifdef HAS_MA86
#include "ingredients/subproblem_solvers/MA86/MA86Solver.hpp"
#endif

#ifdef HAS_HSL
namespace uno {
   extern "C" {
      bool LIBHSL_isfunctional();
   }
}
#endif

#ifdef HSL_RUNTIME_LOADING
#include "ingredients/subproblem_solvers/HSL/HSLLoader.hpp"
#endif

#ifdef HAS_MUMPS
#include "ingredients/subproblem_solvers/MUMPS/MUMPSSolver.hpp"
#endif

#ifdef HAS_SPRAL
#include "ingredients/subproblem_solvers/SSIDS/SSIDSSolver.hpp"
#endif

namespace uno {
   namespace {
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> create_solver(const Options& options, bool native_scaling_allowed);

      // options without the defaults (e.g. a standalone linear solver): no scaling
      LinearSystemScalingParameters read_scaling_parameters(const Options& options) {
         LinearSystemScalingParameters parameters;
         parameters.method = linear_system_scaling_method_from_string(options.get_string_optional("linear_system_scaling").value_or("none"));
         parameters.on_demand = options.get_bool_optional("linear_system_scaling_on_demand").value_or(parameters.on_demand);
         parameters.delayed_pivots_threshold = options.get_double_optional("linear_system_scaling_delayed_pivots_threshold")
            .value_or(parameters.delayed_pivots_threshold);
         return parameters;
      }
   } // namespace

   std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> SymmetricIndefiniteLinearSolverFactory::create(
         const Options& options, bool allow_scaling) {
      const LinearSystemScalingParameters parameters = read_scaling_parameters(options);
      const bool linear_system_scaling = (parameters.method != LinearSystemScalingMethod::NONE);
      // the solvers' own scalings (MA57_use_scaling, MUMPS_use_scaling) are disabled when Uno scales the linear systems
      auto linear_solver = create_solver(options, /* native_scaling_allowed = */ !linear_system_scaling);
      if (allow_scaling && linear_system_scaling) {
         linear_solver->set_scaling_strategy(std::make_unique<ScalingStrategy>(parameters));
      }
      return linear_solver;
   }

   namespace {
   std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> create_solver(const Options& options, bool native_scaling_allowed) {
      const auto& linear_solver = options.get_string("linear_solver");
#ifdef HSL_RUNTIME_LOADING
      // honor the requested HSL library name before LIBHSL_isfunctional() triggers the (cached) load
      load_hsl_library(options.get_string_optional("libhsl_path").value_or(""));
#endif
#if defined(HAS_HSL) || defined(HAS_MA57)
      if (linear_solver == "MA57"
   #if defined(HSL_RUNTIME_LOADING)
         && ma57_symbols_available()
   #elif defined(HAS_HSL)
         && LIBHSL_isfunctional()
   #endif
            ) {
         return std::make_unique<MA57Solver>(native_scaling_allowed && options.get_bool_optional("MA57_use_scaling").value_or(false));
      }
#endif

#if defined(HAS_HSL) || defined(HAS_MA27)
      if (linear_solver == "MA27"
   #if defined(HSL_RUNTIME_LOADING)
         && ma27_symbols_available()
   #elif defined(HAS_HSL)
         && LIBHSL_isfunctional()
   #endif
      ) {
         return std::make_unique<MA27Solver>();
      }
#endif // HAS_HSL || HAS_MA27

#ifdef HAS_MA86
      if (linear_solver == "MA86"
   #if defined(HSL_RUNTIME_LOADING)
         && ma86_symbols_available()
   #elif defined(HAS_HSL)
         && LIBHSL_isfunctional()
   #endif
      ) {
         return std::make_unique<MA86Solver>(Indexing::C_indexing);
      }
#endif

#ifdef HAS_MUMPS
      if (linear_solver == "MUMPS") {
         return std::make_unique<MUMPSSolver>(native_scaling_allowed && options.get_bool_optional("MUMPS_use_scaling").value_or(true));
      }
#endif
#ifdef HAS_SPRAL
      if (linear_solver == "SSIDS") {
         return std::make_unique<SSIDSSolver>();
      }
#endif
      // the HSL solver is supported by this build, but the HSL library is not functional
      bool supported_hsl_solver = false;
#if defined(HAS_HSL) || defined(HAS_MA57)
      supported_hsl_solver = supported_hsl_solver || (linear_solver == "MA57");
#endif
#if defined(HAS_HSL) || defined(HAS_MA27)
      supported_hsl_solver = supported_hsl_solver || (linear_solver == "MA27");
#endif
#ifdef HAS_MA86
      supported_hsl_solver = supported_hsl_solver || (linear_solver == "MA86");
#endif
      if (supported_hsl_solver) {
         throw std::invalid_argument("The linear solver " + linear_solver + " was requested, but the HSL library is not "
            "functional or does not provide it (with runtime loading, check libhsl_path or the UNO_HSL_LIBRARY environment variable)");
      }
      std::string message = "The linear solver ";
      message.append(linear_solver).append(" is unknown").append("\n").append("The following values are available: ")
            .append(join(SymmetricIndefiniteLinearSolverFactory::available_solvers(), ", "));
      throw std::invalid_argument(message);
   }
   } // namespace

   // return the list of available solvers
   std::vector<std::string> SymmetricIndefiniteLinearSolverFactory::available_solvers() {
      std::vector<std::string> solvers{};
#if defined(HSL_RUNTIME_LOADING)
      // the runtime library may provide only some of the solvers (e.g. MA27 only, from the HSL Archive)
      if (ma57_symbols_available()) {
         solvers.emplace_back("MA57");
      }
      if (ma27_symbols_available()) {
         solvers.emplace_back("MA27");
      }
#elif defined(HAS_HSL)
      if (LIBHSL_isfunctional()) {
         solvers.emplace_back("MA57");
         solvers.emplace_back("MA27");
      }
#else
   #ifdef HAS_MA57
      solvers.emplace_back("MA57");
   #endif
   #ifdef HAS_MA27
      solvers.emplace_back("MA27");
   #endif
#endif

#ifdef HAS_MA86
   #if defined(HSL_RUNTIME_LOADING)
      if (ma86_symbols_available()) {
         solvers.emplace_back("MA86");
      }
   #elif defined(HAS_HSL)
      if (LIBHSL_isfunctional()) {
         solvers.emplace_back("MA86");
      }
   #else
      solvers.emplace_back("MA86");
   #endif
#endif

#ifdef HAS_MUMPS
      solvers.emplace_back("MUMPS");
#endif

#ifdef HAS_SPRAL
      solvers.emplace_back("SSIDS");
#endif
      return solvers;
   }
} // namespace