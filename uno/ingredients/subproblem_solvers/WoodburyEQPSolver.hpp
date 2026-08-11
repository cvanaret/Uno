// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_WOODBURYEQPSOLVER_H
#define UNO_WOODBURYEQPSOLVER_H

#include <memory>
#include "DirectSymmetricIndefiniteLinearSolver.hpp"
#include "SubproblemSolver.hpp"
#include "optimization/Direction.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class DenseMatrix;
   class Options;
   class DirectQuasiNewtonHessian;

   // The WoodburyEQPSolver is a special case of EQPSolver with a quasi-Newton Hessian model
   // When the Hessian approximation is given by a low-rank correction to a diagonal part:
   // H = δ I + E P Eᵀ,
   // the linear system
   // (δ I + E P Eᵀ     Jᵀ) (d_x)   = - (r_x)
   // (J                0 ) (-d_y)      (r_y)
   // must be decomposed as:
   // (δ I + E P Eᵀ     Jᵀ) = (δ I   Jᵀ) + (E) P (Eᵀ  0)    (*)
   // (J                0 )   (J     0 )   (0)
   // solved the following way:
   // 1. factorize the matrix with the diagonal part only:
   // (δ I   Jᵀ)
   // (J     0 )
   // 2. use the Woodbury formula to invert (*) symbolically
   class WoodburyEQPSolver: public SubproblemSolver {
   public:
      WoodburyEQPSolver(const DirectQuasiNewtonHessian& hessian_model, const Options& options);
      ~WoodburyEQPSolver() override = default;

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
      const DirectQuasiNewtonHessian& hessian_model;
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> linear_solver;
      bool analysis_performed{false};

      void compute_low_rank_correction(const Subproblem& subproblem, Vector<double>& b) const;
      [[nodiscard]] static bool solve_dense_indefinite_system(DenseMatrix<double>& T, const Vector<double>& c, Vector<double>& d);
   };
} // namespace

#endif // UNO_WOODBURYEQPSOLVER_H