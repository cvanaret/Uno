// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_WOODBURYEQPSOLVER_H
#define UNO_WOODBURYEQPSOLVER_H

#include <memory>
#include "DirectSymmetricIndefiniteLinearSolver.hpp"
#include "SubproblemSolver.hpp"

namespace uno {
   // forward declarations
   template <typename ElementType>
   class DenseMatrix;
   class Options;
   class QuasiNewtonHessian;

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
      WoodburyEQPSolver(const QuasiNewtonHessian& hessian_model, const Options& options);
      ~WoodburyEQPSolver() override = default;

      void initialize_memory(const Subproblem& subproblem) override;

      void solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius, const Vector<double>& initial_point,
         Direction& direction, Evaluations& current_evaluations, const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   protected:
      const QuasiNewtonHessian& hessian_model;
      std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> linear_solver;
      bool analysis_performed{false};

      void compute_low_rank_correction(const Subproblem& subproblem, LinearSystem& linear_system, Vector<double>& b) const;
      [[nodiscard]] static bool solve_dense_indefinite_system(DenseMatrix<double>& T, const Vector<double>& c, Vector<double>& d);
   };
} // namespace

#endif // UNO_WOODBURYEQPSOLVER_H