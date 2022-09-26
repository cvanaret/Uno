// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include "SparseVectorWrapper.hpp"

SparseVector<double>* SparseVector_new(size_t capacity) {
   return new SparseVector<double>(capacity);
}

void SparseVector_delete(SparseVector<double>* vector) {
   delete vector;
}

void SparseVector_insert(SparseVector<double>* vector, size_t key, double value) {
   vector->insert(key, value);
}

void SparseVector_clear(SparseVector<double>* vector) {
   vector->clear();
}

void SparseVector_display(SparseVector<double>* vector) {
   vector->for_each([](size_t key, double value) {
      std::cout << "(" << key << ", " << value << ")\n";
   });
}
