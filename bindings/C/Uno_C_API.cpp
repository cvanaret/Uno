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

// current version is 2.0.1
void uno_get_version(int* major, int* minor, int* patch) {
   *major = 2;
   *minor = 0;
   *patch = 1;
}

void* uno_create_model() {
   return new CModel();
}

void uno_destroy_model(void* uno_model) {
   assert(uno_model != nullptr);
   delete static_cast<CModel*>(uno_model);
}