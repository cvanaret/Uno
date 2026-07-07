// Copyright (c) 2024-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "Presets.hpp"
#include "Options.hpp"
#include "model/Model.hpp"
#include "tools/Logger.hpp"

namespace uno {
   Preset Presets::from_string(const std::string& preset) {
      if (preset == "filtersqp") {
         return Preset::FILTERSQP;
      }
      if (preset == "ipopt") {
         return Preset::IPOPT;
      }
      if (preset == "funnelsqp") {
         return Preset::FUNNELSQP;
      }
      if (preset == "filterslp") {
         return Preset::FILTERSLP;
      }
      throw std::runtime_error("Unknown preset");
   }

   Preset Presets::pick_auto_preset(const Model& model) {
      // thresholds encode "dense active-set QP stops scaling around here"
      // TODO: calibrate on a representative benchmark set (CUTEst, MINLPTests)
      constexpr size_t large_problem_dimension = 2000;       // n + m
      constexpr size_t large_problem_nonzeros = 50'000;      // Jacobian + Hessian nnz
      // constexpr size_t many_inequalities_floor = 500;        // absolute floor for the ratio rule

      const size_t total_dimension = model.number_variables + model.number_constraints;
      size_t total_nonzeros = model.number_jacobian_nonzeros();
      if (model.has_hessian_matrix()) {
         total_nonzeros += model.number_hessian_nonzeros();
      }

      // 1. large/sparse => barrier
      if (large_problem_dimension <= total_dimension || large_problem_nonzeros <= total_nonzeros) {
         INFO << "Automatically picked the ipopt preset\n";
         return Preset::IPOPT;
      }
      /*
      // 2. many potential active constraints, and more than the DOF => barrier
      if (many_inequalities_floor < number_potential_active_constraints &&
            model.number_variables < number_potential_active_constraints) {
         INFO << "Automatically picked the ipopt preset\n";
         return Preset::IPOPT;
      }
      */
      // 3. small/medium with modest inequalities (and the specialized cases) => active-set SQP
      else {
         INFO << "Automatically picked the filtersqp preset\n";
         return Preset::FILTERSQP;
      }
   }

   void Presets::set(Options& options, const std::string& preset) {
      if (preset == "auto") {
         throw std::runtime_error("The auto preset requires the model");
      }
      else {
         set(options, from_string(preset));
      }
   }

   void Presets::set(const Model& model, Options& options, const std::string& preset) {
      if (preset == "auto") {
         set(options, pick_auto_preset(model));
      }
      else {
         set(options, from_string(preset));
      }
   }

   void Presets::set(Options& options, Preset preset) {
      // shortcuts for state-of-the-art combinations
      if (preset == Preset::IPOPT) {
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
      else if (preset == Preset::FILTERSQP) {
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
      else if (preset == Preset::FUNNELSQP) {
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
      else if (preset == Preset::FILTERSLP) {
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
         throw std::runtime_error("The preset is not known");
      }
   }
} // namespace