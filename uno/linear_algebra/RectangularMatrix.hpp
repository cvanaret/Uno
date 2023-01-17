// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RECTANGULARMATRIX_H
#define UNO_RECTANGULARMATRIX_H

#include "SparseVector.hpp"

template <typename T>
using RectangularMatrix = std::vector<SparseVector<T>>;

#endif // UNO_RECTANGULARMATRIX_H