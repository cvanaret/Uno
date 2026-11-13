// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PROBLEMTYPE_H
#define UNO_PROBLEMTYPE_H

namespace uno {
   enum class ProblemType{LINEAR, QUADRATIC, NONLINEAR};

   inline std::string to_string(ProblemType problem_type) {
      if (problem_type == ProblemType::LINEAR) {
         return "LP";
      }
      else if (problem_type == ProblemType::QUADRATIC) {
         return "QP";
      }
      else {
         return "NLP";
      }
   }
} // namespace

#endif // UNO_PROBLEMTYPE_H