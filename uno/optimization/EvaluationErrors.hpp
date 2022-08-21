// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_EVALUATIONERRORS_H
#define UNO_EVALUATIONERRORS_H



struct EvaluationError : public std::exception {
   [[nodiscard]] const char* what() const noexcept override = 0;
};

struct GradientEvaluationError : EvaluationError {
   [[nodiscard]] const char* what() const noexcept override {
      return "A numerical error was encountered while evaluating a gradient";
   }
};

struct FunctionEvaluationError : EvaluationError {
   [[nodiscard]] const char* what() const noexcept override {
      return "A numerical error was encountered while evaluating a function";
   }
};

#endif // UNO_EVALUATIONERRORS_H
