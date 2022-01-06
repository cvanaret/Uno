#ifndef UNO_COOSYMMETRICMATRIX_H
#define UNO_COOSYMMETRICMATRIX_H

#include "SymmetricMatrix.hpp"

class COOSymmetricMatrix: public SymmetricMatrix {
   // Coordinate list
public:
   COOSymmetricMatrix(size_t dimension, size_t capacity);

   void reset() override;
   void for_each(const std::function<void (size_t, size_t, double)>& f) const override;
   void insert(double term, size_t row_index, size_t column_index) override;
   void pop() override;
   void add_identity_multiple(double multiple) override;
   [[nodiscard]] double smallest_diagonal_entry() const override;

   void print(std::ostream& stream) const override;

protected:
   std::vector<size_t> row_indices;
   std::vector<size_t> column_indices;

   [[nodiscard]] size_t find(size_t row_index, size_t column_index) const;
};

#endif // UNO_COOSYMMETRICMATRIX_H