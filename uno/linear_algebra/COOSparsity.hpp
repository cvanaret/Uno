// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COOSPARSITY_H
#define UNO_COOSPARSITY_H

#include <vector>
#include "../interfaces/C/uno_int.h"

namespace uno {
   class COOSparsity {
   public:
      std::vector<uno_int> row_indices;
      std::vector<uno_int> column_indices;

      COOSparsity(size_t number_nonzeros): row_indices(number_nonzeros), column_indices(number_nonzeros) {
      }
   };
} // namespace

#endif // UNO_COOSPARSITY_H