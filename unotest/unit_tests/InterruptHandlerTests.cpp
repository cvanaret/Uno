#include <gtest/gtest.h>
#include <csignal>
#include "tools/InterruptHandler.hpp"

namespace {
   volatile std::sig_atomic_t sentinel_called = 0;
   extern "C" void sentinel_handler(int) { sentinel_called = 1; }
}

TEST(InterruptGuard, RaiseSetsFlag) {
   uno::InterruptGuard guard;
   EXPECT_FALSE(uno::InterruptGuard::interrupt_requested());
   std::raise(SIGINT);
   EXPECT_TRUE(uno::InterruptGuard::interrupt_requested());
   // do NOT raise again here: the handler has reset SIGINT to SIG_DFL
}

TEST(InterruptGuard, NewGuardResetsFlag) {
   { uno::InterruptGuard guard; std::raise(SIGINT); }
   uno::InterruptGuard guard;
   EXPECT_FALSE(uno::InterruptGuard::interrupt_requested());
}

TEST(InterruptGuard, RestoresPreviousHandler) {
   auto old = std::signal(SIGINT, sentinel_handler);
   { uno::InterruptGuard guard; }
   std::raise(SIGINT);
   EXPECT_EQ(sentinel_called, 1);
   std::signal(SIGINT, old);
}