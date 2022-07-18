// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#ifndef UNO_MODELFACTORY_H
#define UNO_MODELFACTORY_H

#include <memory>
#include "Model.hpp"
#include "tools/Options.hpp"

class ModelFactory {
public:
   static std::unique_ptr<Model> create(const std::string& problem_name);
};

#endif // UNO_MODELFACTORY_H
