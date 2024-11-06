// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HESSIANMODELFACTORY_H
#define UNO_HESSIANMODELFACTORY_H

#include <memory>
#include <string>

namespace uno {
   // forward declarations
   class HessianModel;
   class Options;

   class HessianModelFactory {
   public:
      static std::unique_ptr<HessianModel> create(const std::string& hessian_model, size_t dimension, size_t maximum_number_nonzeros,
            bool convexify, const Options& options);
   };
} // namespace

#endif // UNO_HESSIANMODELFACTORY_H,