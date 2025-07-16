// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNOSOLVERWRAPPER_H
#define UNO_UNOSOLVERWRAPPER_H

#include <functional>
#include "Uno.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/PointerWrapper.hpp"

namespace uno {
   // forward declarations
   class Options;

   using objective_function_type = std::function<double(PointerWrapper<Vector<double>>)>;
   using constraint_functions_type = std::function<void(PointerWrapper<Vector<double>>, PointerWrapper<Vector<double>>)>;
   using objective_gradient_type = std::function<void(PointerWrapper<Vector<double>>, PointerWrapper<SparseVector<double>>)>;
   using jacobian_type = std::function<void(PointerWrapper<Vector<double>> /*x*/,
      PointerWrapper<SymmetricMatrix<size_t, double>> /*jacobian*/)>;
   using hessian_type = std::function<void(PointerWrapper<Vector<double>> /*x*/, double objective_multiplier,
      PointerWrapper<Vector<double>> /*y*/, PointerWrapper<SymmetricMatrix<size_t, double>> /*hessian*/)>;

   class UnoSolverWrapper {
   public:
      UnoSolverWrapper(bool constrained_model, const Options& options);

      void solve(size_t number_variables, size_t number_constraints, const objective_function_type& evaluate_objective,
         const constraint_functions_type& evaluate_constraints, const objective_gradient_type& evaluate_objective_gradient,
         const jacobian_type& evaluate_jacobian, const hessian_type& evaluate_hessian, const Options& options);

   protected:

   };
} // namespace

#endif // UNO_UNOSOLVERWRAPPER_H