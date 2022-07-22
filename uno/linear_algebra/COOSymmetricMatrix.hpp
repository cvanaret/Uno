// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COOSYMMETRICMATRIX_H
#define UNO_COOSYMMETRICMATRIX_H

#include <functional>
#include "SymmetricMatrix.hpp"

class COOSymmetricMatrix: public SymmetricMatrix {
   // Coordinate list
public:
   COOSymmetricMatrix(size_t dimension, size_t original_capacity, bool use_regularization);

   void reset() override;
   void for_each(const std::function<void (size_t, size_t, double)>& f) const override;
   void insert(double term, size_t row_index, size_t column_index) override;
   void finalize(size_t column_index) override;
   [[nodiscard]] double smallest_diagonal_entry() const override;
   void set_regularization(const std::function<double(size_t index)>& f) override;

   void print(std::ostream& stream) const override;

protected:
   std::vector<size_t> row_indices;
   std::vector<size_t> column_indices;

   [[nodiscard]] size_t find(size_t row_index, size_t column_index) const;
   void initialize_regularization();
};

#endif // UNO_COOSYMMETRICMATRIX_H