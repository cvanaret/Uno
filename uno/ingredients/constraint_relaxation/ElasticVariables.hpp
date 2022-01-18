#ifndef UNO_ELASTICVARIABLES_H
#define UNO_ELASTICVARIABLES_H

#include "linear_algebra/SparseVector.hpp"

struct ElasticVariables {
   SparseVector<size_t> positive;
   SparseVector<size_t> negative;
   explicit ElasticVariables(size_t capacity): positive(capacity), negative(capacity) {}
   [[nodiscard]] size_t size() const { return this->positive.size() + this->negative.size(); }
};

#endif // UNO_ELASTICVARIABLES_H
