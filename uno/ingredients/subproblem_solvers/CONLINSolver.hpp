// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CONLINSOLVER_H
#define UNO_CONLINSOLVER_H

#include <vector>
#include "SCPSolver.hpp"

namespace uno {
   class Options;

   class CONLINSolver : public SCPSolver {
   public:
      CONLINSolver(size_t number_variables, size_t number_constraints, const Options& options);
      ~CONLINSolver() override = default;

   protected:
      void compute_diagonal_hessian(const Vector<double>& initial_point, Evaluations& current_evaluations, std::vector<double>& Dx) override;
   };
} // namespace

#endif // UNO_CONLINSOLVER_H
