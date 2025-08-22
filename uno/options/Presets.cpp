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
      Options preset_options;
      // shortcuts for state-of-the-art combinations
      if (preset_name == "ipopt") {
         preset_options.set("constraint_relaxation_strategy", "feasibility_restoration");
         preset_options.set("inequality_handling_method", "primal_dual_interior_point");
         preset_options.set("hessian_model", "exact");
         preset_options.set("regularization_strategy", "primal_dual");
         preset_options.set("globalization_mechanism", "LS");
         preset_options.set("globalization_strategy", "waechter_filter_method");
         preset_options.set("filter_type", "standard");
         preset_options.set("filter_beta", "0.99999");
         preset_options.set("filter_gamma", "1e-8");
         preset_options.set("switching_delta", "1");
         preset_options.set("filter_ubd", "1e4");
         preset_options.set("filter_fact", "1e4");
         preset_options.set("filter_switching_infeasibility_exponent", "1.1");
         preset_options.set("armijo_decrease_fraction", "1e-8");
         preset_options.set("LS_backtracking_ratio", "0.5");
         preset_options.set("LS_min_step_length", "5e-7");
         preset_options.set("barrier_tau_min", "0.99");
         preset_options.set("barrier_damping_factor", "1e-5");
         preset_options.set("l1_constraint_violation_coefficient", "1000.");
         preset_options.set("progress_norm", "L1");
         preset_options.set("residual_norm", "INF");
         preset_options.set("scale_functions", "yes");
         preset_options.set("primal_tolerance", "1e-8");
         preset_options.set("dual_tolerance", "1e-8");
         preset_options.set("loose_dual_tolerance", "1e-6");
         preset_options.set("loose_tolerance_consecutive_iteration_threshold", "15");
         preset_options.set("switch_to_optimality_requires_linearized_feasibility", "no");
         preset_options.set("LS_scale_duals_with_step_length", "yes");
         preset_options.set("protect_actual_reduction_against_roundoff", "yes");
      }
      else if (preset_name == "filtersqp") {
         preset_options.set("constraint_relaxation_strategy", "feasibility_restoration");
         preset_options.set("inequality_handling_method", "inequality_constrained");
         preset_options.set("hessian_model", "exact");
         preset_options.set("regularization_strategy", "none");
         preset_options.set("globalization_mechanism", "TR");
         preset_options.set("globalization_strategy", "fletcher_filter_method");
         preset_options.set("filter_type", "standard");
         preset_options.set("progress_norm", "L1");
         preset_options.set("residual_norm", "L2");
         preset_options.set("TR_radius", "10");
         preset_options.set("l1_constraint_violation_coefficient", "1.");
         preset_options.set("enforce_linear_constraints", "yes");
         preset_options.set("primal_tolerance", "1e-6");
         preset_options.set("dual_tolerance", "1e-6");
         preset_options.set("loose_dual_tolerance", "1e-6");
         preset_options.set("TR_min_radius", "1e-8");
         preset_options.set("switch_to_optimality_requires_linearized_feasibility", "yes");
         preset_options.set("protect_actual_reduction_against_roundoff", "no");
      }
      else if (preset_name == "funnelsqp") {
         preset_options.set("constraint_relaxation_strategy", "feasibility_restoration");
         preset_options.set("inequality_handling_method", "inequality_constrained");
         preset_options.set("hessian_model", "exact");
         preset_options.set("regularization_strategy", "none");
         preset_options.set("globalization_mechanism", "TR");
         preset_options.set("globalization_strategy", "funnel_method");
         preset_options.set("progress_norm", "L1");
         preset_options.set("residual_norm", "L2");
         preset_options.set("TR_radius", "10");
         preset_options.set("l1_constraint_violation_coefficient", "1.");
         preset_options.set("enforce_linear_constraints", "yes");
         preset_options.set("primal_tolerance", "1e-6");
         preset_options.set("dual_tolerance", "1e-6");
         preset_options.set("loose_dual_tolerance", "1e-6");
         preset_options.set("TR_min_radius", "1e-8");
         preset_options.set("switch_to_optimality_requires_acceptance", "no");
         preset_options.set("switch_to_optimality_requires_linearized_feasibility", "yes");

         preset_options.set("funnel_beta", "0.9999");
         preset_options.set("funnel_gamma", "0.001");
         preset_options.set("switching_delta", "0.999");
         preset_options.set("funnel_kappa", "0.5");
         preset_options.set("funnel_ubd", "1.0");
         preset_options.set("funnel_fact", "1.5");
         preset_options.set("funnel_switching_infeasibility_exponent", "2");
         preset_options.set("funnel_update_strategy", "2");
      }
      else if (preset_name == "filterslp") {
         preset_options.set("constraint_relaxation_strategy", "feasibility_restoration");
         preset_options.set("inequality_handling_method", "inequality_constrained");
         preset_options.set("hessian_model", "zero");
         preset_options.set("regularization_strategy", "none");
         preset_options.set("globalization_mechanism", "TR");
         preset_options.set("globalization_strategy", "fletcher_filter_method");
         preset_options.set("filter_type", "standard");
         preset_options.set("progress_norm", "L1");
         preset_options.set("residual_norm", "L2");
         preset_options.set("TR_radius", "10");
         preset_options.set("l1_constraint_violation_coefficient", "1.");
         preset_options.set("enforce_linear_constraints", "yes");
         preset_options.set("primal_tolerance", "1e-6");
         preset_options.set("dual_tolerance", "1e-5");
         preset_options.set("loose_dual_tolerance", "1e-4");
         preset_options.set("TR_min_radius", "1e-8");
         preset_options.set("switch_to_optimality_requires_linearized_feasibility", "yes");
         preset_options.set("protect_actual_reduction_against_roundoff", "no");
      }
      else {
         throw std::runtime_error("The preset " + preset_name + " is not known");
      }
      options.overwrite_with(preset_options);
   }
} // namespace