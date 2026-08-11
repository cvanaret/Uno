// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INVERSENEWTONSOLVER_H
#define UNO_INVERSENEWTONSOLVER_H

#include "SubproblemSolver.hpp"
#include "SolverWorkspace.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/Direction.hpp"

namespace uno {
   class NewtonWorkspace: public SolverWorkspace {
   public:
      NewtonWorkspace() = default;

      double compute_hessian_quadratic_form(const Subproblem& /*subproblem*/, const Iterate& /*current_iterate*/,
            const Vector<double>& /*vector*/) const override {
         // no explicit Hessian. Since the (inverse) Hessian model is positive definite, the predicted reduction can
         // be kept first order
         return 0.;
      }
   };

   // forward declaration
   class InverseLBFGSHessian;

   class InverseNewtonSolver: public SubproblemSolver {
   public:
      explicit InverseNewtonSolver(InverseLBFGSHessian& hessian_model);
      ~InverseNewtonSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override;
      void generate_initial_iterate(const Subproblem& subproblem, Iterate& initial_iterate, Evaluations& evaluations) override;
      void compute_least_squares_multipliers(const Subproblem& subproblem, Iterate& iterate, Evaluations& evaluations) override;

      [[nodiscard]] Direction& solve(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
         double trust_region_radius, const Vector<double>& initial_point, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] bool has_second_order_corrections() const override;
      const Direction& compute_second_order_correction(const Subproblem& subproblem, const Iterate& current_iterate,
         const Vector<double>& constraints_SOC) override;

      [[nodiscard]] const SolverWorkspace& get_workspace() const override;

   protected:
      Direction direction;
      InverseLBFGSHessian& hessian_model;
      NewtonWorkspace workspace;

      Vector<double> rhs;
   };
} // namespace

#endif // UNO_INVERSENEWTONSOLVER_H