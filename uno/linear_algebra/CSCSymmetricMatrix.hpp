// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project root for details.

#ifndef UNO_CSCSYMMETRICMATRIX_H
#define UNO_CSCSYMMETRICMATRIX_H

#include "SymmetricMatrix.hpp"

class CSCSymmetricMatrix : public SymmetricMatrix {
   // Compressed Sparse Column
public:
   CSCSymmetricMatrix(size_t dimension, size_t maximum_number_nonzeros, size_t padding_size = 0);

   void reset() override;
   void for_each(const std::function<void (size_t, size_t, double)>& f) const override;
   void for_each(size_t column_index, const std::function<void (size_t, double)>& f) const;
   void insert(double term, size_t row_index, size_t column_index) override;
   void pop() override;
   void finalize(size_t column_index) override;
   void add_identity_multiple(double multiple, size_t number_variables) override;
   [[nodiscard]] double smallest_diagonal_entry() const override;

   void force_explicit_diagonal_elements();
   void remove_variables(const std::vector<int>& variable_indices);

   static CSCSymmetricMatrix identity(size_t dimension);

   void print(std::ostream& stream) const override;

protected:
   std::vector<size_t> column_starts{};
   std::vector<size_t> row_indices{};
   size_t current_column{0};
   // when elements are inserted one by one, keep track of the current column
   std::vector<size_t> remaining_column_padding;
};

#endif // UNO_CSCSYMMETRICMATRIX_H