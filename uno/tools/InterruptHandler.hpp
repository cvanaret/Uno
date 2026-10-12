// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INTERRUPTHANDLER_H
#define UNO_INTERRUPTHANDLER_H

namespace uno {
   class InterruptGuard {
   public:
      InterruptGuard();
      ~InterruptGuard();
      InterruptGuard(const InterruptGuard&) = delete;
      InterruptGuard& operator=(const InterruptGuard&) = delete;

      [[nodiscard]] static bool interrupt_requested();

   private:
      using Handler = void (*)(int);
      Handler previous_handler;
   };
} // namespace

#endif //UNO_INTERRUPTHANDLER_H