// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "Uno_C_API.h"

// TODO
class CModel {
public:
   CModel() = default;
   ~CModel() = default;
};

void* create_uno_model() {
   return new CModel();
}

void destroy_uno_model(void* uno_model) {
   assert(uno_model != nullptr);
   delete static_cast<CModel*>(uno_model);
}