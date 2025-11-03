// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_OPTIONS_H
#define UNO_OPTIONS_H

#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>

namespace uno {
   enum class OptionType {INTEGER, DOUBLE, BOOL, STRING};

   class Options {
   public:
      Options() = default;

      void set_integer(const std::string& option_name, int32_t option_value, bool flag_as_overwritten = false);
      void set_double(const std::string& option_name, double option_value, bool flag_as_overwritten = false);
      void set_bool(const std::string& option_name, bool option_value, bool flag_as_overwritten = false);
      void set_string(const std::string& option_name, const std::string& option_value, bool flag_as_overwritten = false);
      // setter for option with unknown type
      void set(const std::string& option_name, const std::string& option_value, bool flag_as_overwritten = false);

      [[nodiscard]] int get_int(const std::string& option_name) const;
      [[nodiscard]] size_t get_unsigned_int(const std::string& option_name) const;
      [[nodiscard]] double get_double(const std::string& option_name) const;
      [[nodiscard]] bool get_bool(const std::string& option_name) const;
      [[nodiscard]] const std::string& get_string(const std::string& option_name) const;
      [[nodiscard]] std::optional<std::string> get_string_optional(const std::string& option_name) const;
      [[nodiscard]] OptionType get_option_type(const std::string& option_name) const;
      [[nodiscard]] std::map<std::string, int32_t> get_integer_options() const;
      [[nodiscard]] std::map<std::string, double> get_double_options() const;
      [[nodiscard]] std::map<std::string, bool> get_bool_options() const;
      [[nodiscard]] std::map<std::string, std::string> get_string_options() const;

      [[nodiscard]] static Options get_command_line_options(int argc, char* argv[], size_t offset);
      [[nodiscard]] static Options load_option_file(const std::string& file_name);

      void overwrite_with(const Options& overwriting_options);
      void print_used_overwritten() const;

   private:
      std::map<std::string, int32_t> integer_options{};
      std::map<std::string, double> double_options{};
      std::map<std::string, bool> bool_options{};
      std::map<std::string, std::string> string_options{};

      mutable std::map<std::string, bool> used{};
      mutable std::map<std::string, bool> overwritten_options{};

      const std::unordered_map<std::string, OptionType> option_types = {
         {"primal_tolerance", OptionType::DOUBLE},
         {"dual_tolerance", OptionType::DOUBLE},
         {"loose_primal_tolerance", OptionType::DOUBLE},
         {"loose_dual_tolerance", OptionType::DOUBLE},
         {"loose_tolerance_consecutive_iteration_threshold", OptionType::INTEGER},
         {"max_iterations", OptionType::INTEGER},
         {"time_limit", OptionType::DOUBLE},
         {"print_solution", OptionType::BOOL},
         {"unbounded_objective_threshold", OptionType::DOUBLE},
         {"enforce_linear_constraints", OptionType::BOOL},
         {"statistics_major_column_order", OptionType::INTEGER},
         {"statistics_minor_column_order", OptionType::INTEGER},
         {"statistics_penalty_parameter_column_order", OptionType::INTEGER},
         {"statistics_barrier_parameter_column_order", OptionType::INTEGER},
         {"statistics_SOC_column_order", OptionType::INTEGER},
         {"statistics_TR_radius_column_order", OptionType::INTEGER},
         {"statistics_LS_step_length_column_order", OptionType::INTEGER},
         {"statistics_restoration_phase_column_order", OptionType::INTEGER},
         {"statistics_primal_regularization_column_order", OptionType::INTEGER},
         {"statistics_funnel_width_column_order", OptionType::INTEGER},
         {"statistics_step_norm_column_order", OptionType::INTEGER},
         {"statistics_objective_column_order", OptionType::INTEGER},
         {"statistics_primal_feasibility_column_order", OptionType::INTEGER},
         {"statistics_dual_feasibility_column_order", OptionType::INTEGER},
         {"statistics_stationarity_column_order", OptionType::INTEGER},
         {"statistics_complementarity_column_order", OptionType::INTEGER},
         {"statistics_status_column_order", OptionType::INTEGER},
         {"logger", OptionType::STRING},
         {"constraint_relaxation_strategy", OptionType::STRING},
         {"inequality_handling_method", OptionType::STRING},
         {"globalization_mechanism",OptionType::STRING},
         {"globalization_strategy", OptionType::STRING},
         {"hessian_model", OptionType::STRING},
         {"inertia_correction_strategy", OptionType::STRING},
         {"scale_functions", OptionType::BOOL},
         {"function_scaling_threshold", OptionType::DOUBLE},
         {"function_scaling_factor", OptionType::DOUBLE},
         {"scale_residuals", OptionType::BOOL},
         {"progress_norm", OptionType::STRING},
         {"residual_norm", OptionType::STRING},
         {"residual_scaling_threshold", OptionType::DOUBLE},
         {"protect_actual_reduction_against_roundoff", OptionType::BOOL},
         {"print_subproblem", OptionType::BOOL},
         {"armijo_decrease_fraction", OptionType::DOUBLE},
         {"armijo_tolerance", OptionType::DOUBLE},
         {"switching_delta", OptionType::DOUBLE},
         {"switching_infeasibility_exponent", OptionType::DOUBLE},
         {"filter_type", OptionType::STRING},
         {"filter_beta", OptionType::DOUBLE},
         {"filter_gamma", OptionType::DOUBLE},
         {"filter_ubd", OptionType::DOUBLE},
         {"filter_fact", OptionType::DOUBLE},
         {"filter_capacity", OptionType::INTEGER},
         {"filter_sufficient_infeasibility_decrease_factor", OptionType::DOUBLE},
         {"nonmonotone_filter_number_dominated_entries", OptionType::INTEGER},
         {"funnel_kappa", OptionType::DOUBLE},
         {"funnel_beta", OptionType::DOUBLE},
         {"funnel_gamma", OptionType::DOUBLE},
         {"funnel_ubd", OptionType::DOUBLE},
         {"funnel_fact", OptionType::DOUBLE},
         {"funnel_update_strategy", OptionType::INTEGER},
         {"funnel_require_acceptance_wrt_current_iterate", OptionType::BOOL},
         {"LS_backtracking_ratio", OptionType::DOUBLE},
         {"LS_min_step_length", OptionType::DOUBLE},
         {"LS_scale_duals_with_step_length", OptionType::BOOL},
         {"regularization_failure_threshold", OptionType::DOUBLE},
         {"regularization_initial_value", OptionType::DOUBLE},
         {"regularization_increase_factor", OptionType::DOUBLE},
         {"primal_regularization_initial_factor", OptionType::DOUBLE},
         {"dual_regularization_fraction", OptionType::DOUBLE},
         {"primal_regularization_lb", OptionType::DOUBLE},
         {"primal_regularization_decrease_factor", OptionType::DOUBLE},
         {"primal_regularization_fast_increase_factor", OptionType::DOUBLE},
         {"primal_regularization_slow_increase_factor", OptionType::DOUBLE},
         {"threshold_unsuccessful_attempts", OptionType::INTEGER},
         {"TR_radius", OptionType::DOUBLE},
         {"TR_increase_factor", OptionType::DOUBLE},
         {"TR_decrease_factor", OptionType::DOUBLE},
         {"TR_aggressive_decrease_factor", OptionType::DOUBLE},
         {"TR_activity_tolerance", OptionType::DOUBLE},
         {"TR_min_radius", OptionType::DOUBLE},
         {"TR_radius_reset_threshold", OptionType::DOUBLE},
         {"switch_to_optimality_requires_linearized_feasibility", OptionType::BOOL},
         {"l1_constraint_violation_coefficient", OptionType::DOUBLE},
         {"barrier_initial_parameter", OptionType::DOUBLE},
         {"barrier_default_multiplier", OptionType::DOUBLE},
         {"barrier_tau_min", OptionType::DOUBLE},
         {"barrier_k_sigma", OptionType::DOUBLE},
         {"barrier_smax", OptionType::DOUBLE},
         {"barrier_k_mu", OptionType::DOUBLE},
         {"barrier_theta_mu", OptionType::DOUBLE},
         {"barrier_k_epsilon", OptionType::DOUBLE},
         {"barrier_update_fraction", OptionType::DOUBLE},
         {"barrier_regularization_exponent", OptionType::DOUBLE},
         {"barrier_small_direction_factor", OptionType::DOUBLE},
         {"barrier_push_variable_to_interior_k1", OptionType::DOUBLE},
         {"barrier_push_variable_to_interior_k2", OptionType::DOUBLE},
         {"barrier_damping_factor", OptionType::DOUBLE},
         {"least_square_multiplier_max_norm", OptionType::DOUBLE},
         {"BQPD_kmax", OptionType::INTEGER},
         {"QP_solver", OptionType::STRING},
         {"LP_solver", OptionType::STRING},
         {"linear_solver", OptionType::STRING},
         {"preset", OptionType::STRING},
      };
   };
} // namespace

#endif // UNO_OPTIONS_H