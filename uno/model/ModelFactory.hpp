// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MODELFACTORY_H
#define UNO_MODELFACTORY_H

#include <memory>
#include "Model.hpp"

namespace uno {
   // forward declarations
   class Iterate;
   class Options;

   class ModelFactory {
   public:
      static std::unique_ptr<Model> reformulate(std::unique_ptr<Model> model, Iterate& initial_iterate, const Options& options);
   };
} // namespace

#endif // UNO_MODELFACTORY_H
