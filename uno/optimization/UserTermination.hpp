// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_USERTERMINATION_H
#define UNO_USERTERMINATION_H

namespace uno {
   // thrown by a user callback (a model evaluation or a solver callback) to stop Uno immediately: Uno returns the last
   // accepted iterate with the USER_TERMINATION status.
   // It is a request, not an error: unlike an EvaluationError, it is not handled by the globalization mechanisms (no
   // backtracking), and it deliberately does not derive from std::exception, so that the handlers of errors
   // (catch (const std::exception&)) let it through, wherever the callback is called
   struct UserTermination {
      [[nodiscard]] const char* what() const noexcept {
         return "The user requested the termination";
      }
   };
} // namespace

#endif // UNO_USERTERMINATION_H
