// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_ELASTICVARIABLES_H
#define UNO_ELASTICVARIABLES_H

#include "linear_algebra/SparseVector.hpp"

namespace uno {
   // forward declaration
   class Model;

   struct ElasticVariablesSizes {
      size_t positive{0};
      size_t negative{0};
   };

   class ElasticVariables {
   public:
      SparseVector<size_t> positive{};
      SparseVector<size_t> negative{};

      ElasticVariables(size_t number_positive_variables, size_t number_negative_variables);
      [[nodiscard]] size_t size() const { return this->positive.size() + this->negative.size(); }

      static ElasticVariables generate(const Model& model);

   protected:
      static ElasticVariablesSizes count(const Model& model);
   };
} // namespace

#endif // UNO_ELASTICVARIABLES_H
