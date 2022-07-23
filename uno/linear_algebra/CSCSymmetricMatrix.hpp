// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CSCSYMMETRICMATRIX_H
#define UNO_CSCSYMMETRICMATRIX_H

#include "SymmetricMatrix.hpp"

class CSCSymmetricMatrix : public SymmetricMatrix {
   // Compressed Sparse Column
public:
   CSCSymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization);

   void reset() override;
   void for_each(const std::function<void (size_t, size_t, double)>& f) const override;
   void for_each(size_t column_index, const std::function<void (size_t, double)>& f) const;
   void insert(double term, size_t row_index, size_t column_index) override;
   void finalize_column(size_t column_index) override;
   [[nodiscard]] double smallest_diagonal_entry() const override;
   void set_regularization(const std::function<double(size_t index)>& regularization_function) override;

   void print(std::ostream& stream) const override;

protected:
   std::vector<size_t> column_starts{};
   std::vector<size_t> row_indices{};
   size_t current_column{0};
};

#endif // UNO_CSCSYMMETRICMATRIX_H