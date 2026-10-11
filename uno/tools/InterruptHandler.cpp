// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "InterruptHandler.hpp"
#include <atomic>
#include <csignal>

namespace {
   std::atomic<bool> interrupted{false};
   static_assert(std::atomic<bool>::is_always_lock_free, "signal-safe flag requires lock-free atomic");
}

// C linkage: strictly required for a signal handler
extern "C" void uno_sigint_handler(int /*signal*/) {
   interrupted.store(true, std::memory_order_relaxed);
   // restore default: a second Ctrl+C kills the process (Windows does this implicitly anyway)
   std::signal(SIGINT, SIG_DFL);
}

namespace uno {
   InterruptGuard::InterruptGuard() {
      interrupted.store(false, std::memory_order_relaxed);
      this->previous_handler = std::signal(SIGINT, uno_sigint_handler);
   }

   InterruptGuard::~InterruptGuard() {
      if (this->previous_handler != SIG_ERR) {
         std::signal(SIGINT, this->previous_handler);
      }
   }

   bool InterruptGuard::interrupt_requested() {
      return interrupted.load(std::memory_order_relaxed);
   }
}