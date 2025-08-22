// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "Presets.hpp"
#include "Options.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"

namespace uno {
   Options Presets::get_preset_options(const std::optional<std::string>& optional_preset) {
      Options options;

      /** optional user preset **/
      if (optional_preset.has_value()) {
         Presets::set(options, *optional_preset);
      }
      else {
         /** default preset **/
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
      return options;
   }

   void Presets::set(Options& options, const std::string& preset_name) {
      // shortcuts for state-of-the-art combinations
      if (preset_name == "ipopt") {
         options.set("constraint_relaxation_strategy", "feasibility_restoration");
         options.set("inequality_handling_method", "primal_dual_interior_point");
         options.set("hessian_model", "exact");
         options.set("regularization_strategy", "primal_dual");
         options.set("globalization_mechanism", "LS");
         options.set("globalization_strategy", "waechter_filter_method");
         options.set("filter_type", "standard");
         options.set("filter_beta", "0.99999");
         options.set("filter_gamma", "1e-8");
         options.set("switching_delta", "1");
         options.set("filter_ubd", "1e4");
         options.set("filter_fact", "1e4");
         options.set("filter_switching_infeasibility_exponent", "1.1");
         options.set("armijo_decrease_fraction", "1e-8");
         options.set("LS_backtracking_ratio", "0.5");
         options.set("LS_min_step_length", "5e-7");
         options.set("barrier_tau_min", "0.99");
         options.set("barrier_damping_factor", "1e-5");
         options.set("l1_constraint_violation_coefficient", "1000.");
         options.set("progress_norm", "L1");
         options.set("residual_norm", "INF");
         options.set("scale_functions", "yes");
         options.set("primal_tolerance", "1e-8");
         options.set("dual_tolerance", "1e-8");
         options.set("loose_dual_tolerance", "1e-6");
         options.set("loose_tolerance_consecutive_iteration_threshold", "15");
         options.set("switch_to_optimality_requires_linearized_feasibility", "no");
         options.set("LS_scale_duals_with_step_length", "yes");
         options.set("protect_actual_reduction_against_roundoff", "yes");
      }
      else if (preset_name == "filtersqp") {
         options.set("constraint_relaxation_strategy", "feasibility_restoration");
         options.set("inequality_handling_method", "inequality_constrained");
         options.set("hessian_model", "exact");
         options.set("regularization_strategy", "none");
         options.set("globalization_mechanism", "TR");
         options.set("globalization_strategy", "fletcher_filter_method");
         options.set("filter_type", "standard");
         options.set("progress_norm", "L1");
         options.set("residual_norm", "L2");
         options.set("TR_radius", "10");
         options.set("l1_constraint_violation_coefficient", "1.");
         options.set("enforce_linear_constraints", "yes");
         options.set("primal_tolerance", "1e-6");
         options.set("dual_tolerance", "1e-6");
         options.set("loose_dual_tolerance", "1e-6");
         options.set("TR_min_radius", "1e-8");
         options.set("switch_to_optimality_requires_linearized_feasibility", "yes");
         options.set("protect_actual_reduction_against_roundoff", "no");
      }
      else if (preset_name == "funnelsqp") {
         options.set("constraint_relaxation_strategy", "feasibility_restoration");
         options.set("inequality_handling_method", "inequality_constrained");
         options.set("hessian_model", "exact");
         options.set("regularization_strategy", "none");
         options.set("globalization_mechanism", "TR");
         options.set("globalization_strategy", "funnel_method");
         options.set("progress_norm", "L1");
         options.set("residual_norm", "L2");
         options.set("TR_radius", "10");
         options.set("l1_constraint_violation_coefficient", "1.");
         options.set("enforce_linear_constraints", "yes");
         options.set("primal_tolerance", "1e-6");
         options.set("dual_tolerance", "1e-6");
         options.set("loose_dual_tolerance", "1e-6");
         options.set("TR_min_radius", "1e-8");
         options.set("switch_to_optimality_requires_acceptance", "no");
         options.set("switch_to_optimality_requires_linearized_feasibility", "yes");

         options.set("funnel_beta", "0.9999");
         options.set("funnel_gamma", "0.001");
         options.set("switching_delta", "0.999");
         options.set("funnel_kappa", "0.5");
         options.set("funnel_ubd", "1.0");
         options.set("funnel_fact", "1.5");
         options.set("funnel_switching_infeasibility_exponent", "2");
         options.set("funnel_update_strategy", "2");
      }
      else if (preset_name == "filterslp") {
         options.set("constraint_relaxation_strategy", "feasibility_restoration");
         options.set("inequality_handling_method", "inequality_constrained");
         options.set("hessian_model", "zero");
         options.set("regularization_strategy", "none");
         options.set("globalization_mechanism", "TR");
         options.set("globalization_strategy", "fletcher_filter_method");
         options.set("filter_type", "standard");
         options.set("progress_norm", "L1");
         options.set("residual_norm", "L2");
         options.set("TR_radius", "10");
         options.set("l1_constraint_violation_coefficient", "1.");
         options.set("enforce_linear_constraints", "yes");
         options.set("primal_tolerance", "1e-6");
         options.set("dual_tolerance", "1e-5");
         options.set("loose_dual_tolerance", "1e-4");
         options.set("TR_min_radius", "1e-8");
         options.set("switch_to_optimality_requires_linearized_feasibility", "yes");
         options.set("protect_actual_reduction_against_roundoff", "no");
      }
      else {
         throw std::runtime_error("The preset " + preset_name + " is not known");
      }
   }
} // namespace