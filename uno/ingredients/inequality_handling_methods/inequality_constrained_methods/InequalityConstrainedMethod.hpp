// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYCONSTRAINEDMETHOD_H
#define UNO_INEQUALITYCONSTRAINEDMETHOD_H

#include "../InequalityHandlingMethod.hpp"
#include "ingredients/subproblem_solvers/InequalityQPSolver.hpp"

namespace uno {
   class InequalityConstrainedMethod : public InequalityHandlingMethod {
   public:
      InequalityConstrainedMethod(size_t number_variables, size_t number_constraints,
         size_t number_objective_gradient_nonzeros, size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros,
         const Options& options);

      ~InequalityConstrainedMethod() override;

      // number of nonzeros
      [[nodiscard]] size_t compute_number_hessian_nonzeros(const OptimizationProblem& problem,
         const HessianModel& hessian_model) const override;

      void initialize_statistics(Statistics& statistics, const Options& options) override;
      void generate_initial_iterate(Statistics& /*statistics*/, const OptimizationProblem& /*problem*/,
         Iterate& /*initial_iterate*/) override;
      [[nodiscard]] SubproblemStatus solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Vector<double>& direction_primals, Multipliers& direction_multipliers, double& subproblem_objective,
         double trust_region_radius, WarmstartInformation& warmstart_information) override;

      void initialize_feasibility_problem(const l1RelaxedProblem& /*problem*/, Iterate& /*current_iterate*/) override;
      void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
      [[nodiscard]] double proximal_coefficient(const Iterate& /*current_iterate*/) const override;
      void exit_feasibility_problem(Statistics& /*statistics*/, const OptimizationProblem& /*problem*/,
         Iterate& /*trial_iterate*/) override;

      [[nodiscard]] double hessian_quadratic_product(const Vector<double>& primal_direction) const override;
      // auxiliary measure is 0 in inequality-constrained methods
      void set_auxiliary_measure(const Model& /*model*/, Iterate& iterate) override;
      [[nodiscard]] double compute_predicted_auxiliary_reduction_model(const Model& /*model*/, const Iterate& /*current_iterate*/,
         const Vector<double>& /*primal_direction*/, double /*step_length*/) const override;

      void postprocess_iterate(const OptimizationProblem& /*problem*/, Iterate& /*iterate*/) override;

      void set_initial_point(const Vector<double>& point) override;

   protected:
      const bool enforce_linear_constraints_at_initial_iterate;
      // pointer to allow polymorphism
      const std::unique_ptr<InequalityQPSolver> subproblem_solver;

      [[nodiscard]] SubproblemStatus solve_inequality_subproblem(Statistics& statistics, LagrangeNewtonSubproblem& subproblem,
         Vector<double>& direction_primals, Multipliers& direction_multipliers, double& subproblem_objective,
         WarmstartInformation& warmstart_information);
      static void compute_dual_displacements(const Multipliers& current_multipliers, Multipliers& direction_multipliers);
   };
} // namespace

#endif // UNO_INEQUALITYCONSTRAINEDMETHOD_H