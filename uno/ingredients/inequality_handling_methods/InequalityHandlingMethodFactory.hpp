// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INEQUALITYHANDLINGMETHODFACTORY_H
#define UNO_INEQUALITYHANDLINGMETHODFACTORY_H

#include <memory>
#include <vector>

namespace uno {
   // forward declaration
   class Options;
   class InequalityHandlingMethod;

   class InequalityHandlingMethodFactory {
      public:
         static std::unique_ptr<InequalityHandlingMethod> create(size_t number_variables, size_t number_constraints, size_t number_objective_gradient_nonzeros,
               size_t number_jacobian_nonzeros, size_t number_hessian_nonzeros, const Options& options);

         static std::vector<std::string> available_strategies();
   };
} // namespace

#endif // UNO_INEQUALITYHANDLINGMETHODFACTORY_H
