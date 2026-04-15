// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NEWTONSOLVER_H
#define UNO_NEWTONSOLVER_H

#include "SubproblemSolver.hpp"
#include "SolverWorkspace.hpp"
#include "linear_algebra/Vector.hpp"

namespace uno {
   class NewtonWorkspace: public SolverWorkspace {
   public:
      NewtonWorkspace() = default;

      double compute_hessian_quadratic_form(const Subproblem& /*subproblem*/, const Vector<double>& /*vector*/) const override {
         // no explicit Hessian. Since the (inverse) Hessian model is positive definite, the predicted reduction can
         // be kept first order
         return 0.;
      }
   };

   // forward declaration
   class InverseLBFGSHessian;

   class NewtonSolver: public SubproblemSolver {
   public:
      explicit NewtonSolver(InverseLBFGSHessian& hessian_model);
      ~NewtonSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& initial_point, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   protected:
      InverseLBFGSHessian& hessian_model;
      NewtonWorkspace workspace;

      Vector<double> rhs;
   };
} // namespace

#endif // UNO_NEWTONSOLVER_H