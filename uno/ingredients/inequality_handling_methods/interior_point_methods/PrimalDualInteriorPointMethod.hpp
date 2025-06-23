// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALINTERIORPOINTMETHOD_H
#define UNO_PRIMALDUALINTERIORPOINTMETHOD_H

#include <memory>
#include "../InequalityHandlingMethod.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "BarrierParameterUpdateStrategy.hpp"

namespace uno {
   // forward references
   class DualResiduals;
   class Subproblem;

   struct InteriorPointParameters {
      double tau_min;
      double k_sigma;
      double regularization_exponent;
      double small_direction_factor;
      double push_variable_to_interior_k1;
      double push_variable_to_interior_k2;
   };

   class PrimalDualInteriorPointMethod : public InequalityHandlingMethod {
   public:
      explicit PrimalDualInteriorPointMethod(const Options& options);

      void initialize(const OptimizationProblem& problem, const HessianModel& hessian_model,
         RegularizationStrategy<double>& regularization_strategy) override;
      void initialize_statistics(Statistics& statistics, const Options& options) override;
      void generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) override;
      void set_initial_point(const Vector<double>& point) override;

      void initialize_feasibility_problem(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
      void exit_feasibility_problem(const OptimizationProblem& problem, Iterate& trial_iterate) override;
      void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& constraint_index) override;
      [[nodiscard]] double proximal_coefficient() const override;

      void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,  const Multipliers& current_multipliers,
         Direction& direction, HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy,
         double trust_region_radius, WarmstartInformation& warmstart_information) override;
      [[nodiscard]] double hessian_quadratic_product(const Vector<double>& vector) const override;

      void set_auxiliary_measure(const OptimizationProblem& problem, Iterate& iterate) override;
      [[nodiscard]] double compute_predicted_auxiliary_reduction_model(const Model& model, const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const override;

      void postprocess_iterate(const OptimizationProblem& problem, Vector<double>& primals, Multipliers& multipliers) override;

      [[nodiscard]] std::string get_name() const override;

   protected:
      SparseVector<double> objective_gradient; /*!< Sparse Jacobian of the objective */
      std::vector<double> constraints; /*!< Constraint values (size \f$m)\f$ */
      Vector<double> solution{};
      const std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<size_t, double>> linear_solver;

      BarrierParameterUpdateStrategy barrier_parameter_update_strategy;
      double previous_barrier_parameter;
      const double default_multiplier;
      const InteriorPointParameters parameters;
      const double least_square_multiplier_max_norm;
      const double damping_factor; // (Section 3.7 in IPOPT paper)
      const double l1_constraint_violation_coefficient; // (rho in Section 3.3.1 in IPOPT paper)

      bool solving_feasibility_problem{false};
      bool first_feasibility_iteration{false};

      [[nodiscard]] double barrier_parameter() const;
      [[nodiscard]] double push_variable_to_interior(double variable_value, double lower_bound, double upper_bound) const;
      void update_barrier_parameter(const OptimizationProblem& problem, const Iterate& current_iterate, const Multipliers& current_multipliers,
         const DualResiduals& residuals);
      [[nodiscard]] bool is_small_step(const OptimizationProblem& problem, const Vector<double>& current_primals, const Vector<double>& direction_primals) const;
      [[nodiscard]] double evaluate_subproblem_objective(const Direction& direction) const;
      [[nodiscard]] double compute_barrier_term_directional_derivative(const Model& model, const Iterate& current_iterate,
         const Vector<double>& primal_direction) const;
      [[nodiscard]] static double primal_fraction_to_boundary(const OptimizationProblem& problem, const Vector<double>& current_primals,
         const Vector<double>& primal_direction, double tau);
      [[nodiscard]] static double dual_fraction_to_boundary(const OptimizationProblem& problem, const Multipliers& current_multipliers,
         Multipliers& direction_multipliers, double tau);
      void assemble_primal_dual_direction(const OptimizationProblem& problem, const Vector<double>& current_primals, const Multipliers& current_multipliers,
         Vector<double>& direction_primals, Multipliers& direction_multipliers);
      void compute_bound_dual_direction(const OptimizationProblem& problem, const Vector<double>& current_primals, const Multipliers& current_multipliers,
         const Vector<double>& primal_direction, Multipliers& direction_multipliers);
   };
} // namespace

#endif // UNO_PRIMALDUALINTERIORPOINTMETHOD_H