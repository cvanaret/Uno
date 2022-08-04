// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MODELFACTORY_H
#define UNO_MODELFACTORY_H

#include <memory>
#include "Model.hpp"
#include "Iterate.hpp"
#include "interfaces/AMPL/AMPLModel.hpp"
#include "tools/Options.hpp"

class ModelFactory {
public:
   static std::unique_ptr<Model> reformulate(const AMPLModel& ampl_model, Iterate& first_iterate, const Options& options);
};

#endif // UNO_MODELFACTORY_H
