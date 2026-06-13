// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MMASOLVER_H
#define UNO_MMASOLVER_H

#include <vector>
#include "SCPSolver.hpp"

namespace uno {
   class Options;

   class MMASolver : public SCPSolver {
   public:
      MMASolver(size_t number_variables, size_t number_constraints, const Options& options);
      ~MMASolver() override = default;

   protected:
      void compute_diagonal_hessian(const Vector<double>& initial_point, Evaluations& current_evaluations, std::vector<double>& Dx) override;

   private:
      std::vector<double> x_old1;
      std::vector<double> x_old2;
      std::vector<double> lower_asymptotes;
      std::vector<double> upper_asymptotes;

      double asyinit;
      double asyincr;
      double asydecr;
      double external_move_limit;
      double internal_limit;
      size_t max_inner_iterations;

      void update_asymptotes(const Vector<double>& current_x, const std::vector<double>& lower_bounds, const std::vector<double>& upper_bounds);
   };
} // namespace

#endif // UNO_MMASOLVER_H
