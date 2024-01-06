// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RECTANGULARMATRIX_H
#define UNO_RECTANGULARMATRIX_H

#include "SparseVector.hpp"

// TODO use more appropriate sparse representation
template <typename ElementType>
using RectangularMatrix = std::vector<SparseVector<ElementType>>;

#endif // UNO_RECTANGULARMATRIX_H