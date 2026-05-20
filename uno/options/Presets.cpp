// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "Presets.hpp"
#include "Options.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"

namespace uno {
   void Presets::set_default(Options& options) {
      // default preset
      const auto linear_solvers = SymmetricIndefiniteLinearSolverFactory::available_solvers();
      if constexpr (0 < QPSolverFactory::available_solvers.size()) {
         Presets::set(options, "filtersqp");
      }
      else if (!linear_solvers.empty()) {
         Presets::set(options, "ipopt");
      }
      else if constexpr (0 < LPSolverFactory::available_solvers.size()) {
         Presets::set(options, "filterslp");
      }
   }

   void Presets::set(Options& options, const std::string& preset_name) {
      // shortcuts for state-of-the-art combinations
      if (preset_name == "ipopt") {
         options.set_string("constraint_relaxation_strategy", "feasibility_restoration");
         options.set_string("inequality_handling_method", "interior_point");
         options.set_string("barrier_function", "log");
         options.set_string("hessian_model", "exact");
         options.set_string("inertia_correction_strategy", "primal_dual");
         options.set_string("globalization_mechanism", "LS");
         options.set_string("globalization_strategy", "waechter_filter_method");
         options.set_string("filter_type", "standard");
         options.set_double("filter_beta", 0.99999);
         options.set_double("filter_gamma", 1e-8);
         options.set_double("switching_delta", 1);
         options.set_double("filter_ubd", 1e4);
         options.set_double("filter_fact", 1e4);
         options.set_double("filter_switching_infeasibility_exponent", 1.1);
         options.set_double("armijo_decrease_fraction", 1e-8);
         options.set_double("LS_backtracking_ratio", 0.5);
         options.set_double("LS_min_step_length", 5e-7);
         options.set_double("barrier_tau_min", 0.99);
         options.set_double("barrier_damping_factor", 1e-5);
         options.set_double("l1_constraint_violation_coefficient", 1000.);
         options.set_string("progress_norm", "L1");
         options.set_string("residual_norm", "INF");
         options.set_double("primal_tolerance", 1e-8);
         options.set_double("dual_tolerance", 1e-8);
         options.set_double("loose_primal_tolerance", 1e-6);
         options.set_double("loose_dual_tolerance", 1e-6);
         options.set_integer("loose_tolerance_iteration_threshold", 15);
         options.set_bool("switch_to_optimality_requires_linearized_feasibility", false);
         options.set_bool("LS_scale_duals_with_step_length", true);
         options.set_bool("protect_actual_reduction_against_roundoff", true);
      }
      else if (preset_name == "filtersqp") {
         options.set_string("constraint_relaxation_strategy", "feasibility_restoration");
         options.set_string("inequality_handling_method", "inequality_constrained");
         options.set_string("hessian_model", "exact");
         options.set_string("inertia_correction_strategy", "none");
         options.set_string("globalization_mechanism", "TR");
         options.set_string("globalization_strategy", "fletcher_filter_method");
         options.set_string("filter_type", "standard");
         options.set_string("progress_norm", "L1");
         options.set_string("residual_norm", "L2");
         options.set_double("TR_radius", 10.);
         options.set_double("l1_constraint_violation_coefficient", 1.);
         options.set_double("primal_tolerance", 1e-6);
         options.set_double("dual_tolerance", 1e-6);
         options.set_bool("switch_to_optimality_requires_linearized_feasibility", true);
         options.set_bool("protect_actual_reduction_against_roundoff", false);
      }
      else if (preset_name == "funnelsqp") {
         options.set_string("constraint_relaxation_strategy", "feasibility_restoration");
         options.set_string("inequality_handling_method", "inequality_constrained");
         options.set_string("hessian_model", "exact");
         options.set_string("inertia_correction_strategy", "none");
         options.set_string("globalization_mechanism", "TR");
         options.set_string("globalization_strategy", "funnel_method");
         options.set_string("progress_norm", "L1");
         options.set_string("residual_norm", "L2");
         options.set_double("TR_radius", 10);
         options.set_double("l1_constraint_violation_coefficient", 1.);
         options.set_double("primal_tolerance", 1e-6);
         options.set_double("dual_tolerance", 1e-6);
         options.set_bool("switch_to_optimality_requires_acceptance", false);
         options.set_bool("switch_to_optimality_requires_linearized_feasibility", true);

         options.set_double("funnel_beta", 0.9999);
         options.set_double("funnel_gamma", 0.001);
         options.set_double("switching_delta", 0.999);
         options.set_double("funnel_kappa", 0.5);
         options.set_double("funnel_ubd", 1.0);
         options.set_double("funnel_fact", 1.5);
         options.set_double("funnel_switching_infeasibility_exponent", 2);
         options.set_integer("funnel_update_strategy", 2);
      }
      else if (preset_name == "filterslp") {
         options.set_string("constraint_relaxation_strategy", "feasibility_restoration");
         options.set_string("inequality_handling_method", "inequality_constrained");
         options.set_string("hessian_model", "zero");
         options.set_string("inertia_correction_strategy", "none");
         options.set_string("globalization_mechanism", "TR");
         options.set_string("globalization_strategy", "fletcher_filter_method");
         options.set_string("filter_type", "standard");
         options.set_string("progress_norm", "L1");
         options.set_string("residual_norm", "L2");
         options.set_double("TR_radius", 10);
         options.set_double("l1_constraint_violation_coefficient", 1.);
         options.set_double("primal_tolerance", 1e-6);
         options.set_double("dual_tolerance", 1e-5);
         options.set_bool("switch_to_optimality_requires_linearized_feasibility", true);
         options.set_bool("protect_actual_reduction_against_roundoff", false);
      }
      else {
         throw std::runtime_error("The preset " + preset_name + " is not known");
      }
   }
} // namespace