// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Runtime behavior of HSL_RUNTIME_LOADING, without the real (licensed) libhsl:
// - no library: the HSL solvers are not advertised and requesting one fails with a clear exception (no crash);
// - a partial library (the stub, which exports MA27 only): only MA27 is advertised and usable;
// - the library name is resolved as documented: libhsl_path option > UNO_HSL_LIBRARY environment variable.

#include <algorithm>
#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
#include <gtest/gtest.h>
#ifndef _WIN32
#include <dlfcn.h>
#endif
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/HSL/HSLLoader.hpp"
#include "ingredients/subproblem_solvers/MA27/MA27Solver.hpp"
#include "ingredients/subproblem_solvers/MA57/MA57Solver.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "options/Options.hpp"

extern "C" bool LIBHSL_isfunctional();

using namespace uno;

namespace {
   // stub library built alongside the tests (path injected by CMake)
   const std::string stub_library = UNO_TEST_HSL_STUB_PATH;
   const std::string missing_library = "uno-nonexistent-directory/libhsl-missing";

   // set (or unset) an environment variable for the lifetime of the object, then restore it
   class ScopedEnvironmentVariable {
   public:
      ScopedEnvironmentVariable(const char* name, const std::optional<std::string>& value): name(name) {
         if (const char* previous = std::getenv(name); previous != nullptr) {
            this->previous_value = previous;
         }
         set(value);
      }
      ~ScopedEnvironmentVariable() {
         set(this->previous_value);
      }
      ScopedEnvironmentVariable(const ScopedEnvironmentVariable&) = delete;
      ScopedEnvironmentVariable& operator=(const ScopedEnvironmentVariable&) = delete;

   private:
      const char* name;
      std::optional<std::string> previous_value{};

      void set(const std::optional<std::string>& value) const {
#ifdef _WIN32
         _putenv_s(this->name, value.has_value() ? value->c_str() : "");
#else
         if (value.has_value()) {
            setenv(this->name, value->c_str(), 1);
         }
         else {
            unsetenv(this->name);
         }
#endif
      }
   };

   Options make_options(const std::string& linear_solver, const std::optional<std::string>& libhsl_path = std::nullopt) {
      Options options;
      options.set_string("linear_solver", linear_solver);
      options.set_bool("MA57_use_scaling", true);
      if (libhsl_path.has_value()) {
         options.set_string("libhsl_path", *libhsl_path);
      }
      return options;
   }

   bool is_available(const std::string& solver) {
      const auto solvers = SymmetricIndefiniteLinearSolverFactory::available_solvers();
      return std::find(solvers.begin(), solvers.end(), solver) != solvers.end();
   }

   // the factory rejects an HSL solver whose symbols cannot be resolved, with an explicit message
   void expect_not_functional(const Options& options) {
      try {
         SymmetricIndefiniteLinearSolverFactory::create(options);
         FAIL() << "the factory created " << options.get_string("linear_solver") << " without the HSL routines";
      }
      catch (const std::invalid_argument& exception) {
         EXPECT_NE(std::string(exception.what()).find("not functional"), std::string::npos) << exception.what();
      }
   }

   // number of calls to the stub's ma27id (checks that Uno called into the dlopen'd library)
   std::optional<int> stub_ma27id_calls() {
#ifdef _WIN32
      return std::nullopt;
#else
      // RTLD_NOLOAD: only succeeds if Uno has loaded the stub
      void* handle = dlopen(stub_library.c_str(), RTLD_NOW | RTLD_NOLOAD);
      if (handle == nullptr) {
         return std::nullopt;
      }
      // get the symbol from the stub lib
      const auto* calls = static_cast<const int*>(dlsym(handle, "uno_hsl_stub_ma27id_calls"));
      const std::optional<int> result = (calls != nullptr) ? std::optional<int>{*calls} : std::nullopt;
      dlclose(handle);
      return result;
#endif
   }
} // namespace

// each test starts and ends with no HSL library loaded (the loader caches the first successful load)
class HSLRuntimeLoading: public testing::Test {
protected:
   void SetUp() override {
      unload_hsl_library();
   }

   void TearDown() override {
      unload_hsl_library();
   }
};

TEST_F(HSLRuntimeLoading, MissingLibraryIsReportedCleanly) {
   const ScopedEnvironmentVariable environment("UNO_HSL_LIBRARY", missing_library);
   EXPECT_FALSE(load_hsl_library());
   EXPECT_FALSE(LIBHSL_isfunctional());
   EXPECT_FALSE(ma27_symbols_available());
   EXPECT_FALSE(ma57_symbols_available());
   EXPECT_FALSE(ma86_symbols_available());

   // the HSL solvers are not advertised
   EXPECT_FALSE(is_available("MA27"));
   EXPECT_FALSE(is_available("MA57"));
   EXPECT_FALSE(is_available("MA86"));

   // requesting one fails cleanly, also with an explicit (missing) libhsl_path
   for (const std::string solver: {"MA27", "MA57", "MA86"}) {
      SCOPED_TRACE(solver);
      expect_not_functional(make_options(solver));
      expect_not_functional(make_options(solver, missing_library));
   }

   // constructing a solver directly throws instead of calling a null function pointer
   EXPECT_THROW(MA27Solver{}, std::runtime_error);
   EXPECT_THROW(MA57Solver{true}, std::runtime_error);
}

TEST_F(HSLRuntimeLoading, PartialLibraryOnlyProvidesItsSolvers) {
   const ScopedEnvironmentVariable environment("UNO_HSL_LIBRARY", stub_library);
   ASSERT_TRUE(load_hsl_library()) << "could not load the stub " << stub_library;
   EXPECT_TRUE(ma27_symbols_available());
   EXPECT_FALSE(ma57_symbols_available());
   EXPECT_FALSE(ma86_symbols_available());

   // only MA27 is advertised...
   EXPECT_TRUE(is_available("MA27"));
   EXPECT_FALSE(is_available("MA57"));
   EXPECT_FALSE(is_available("MA86"));

   // ... and only MA27 can be created
   const auto solver = SymmetricIndefiniteLinearSolverFactory::create(make_options("MA27"));
   EXPECT_NE(solver, nullptr);
   expect_not_functional(make_options("MA57"));
   expect_not_functional(make_options("MA86"));
}

TEST_F(HSLRuntimeLoading, SolverCallsIntoTheLoadedLibrary) {
   const ScopedEnvironmentVariable environment("UNO_HSL_LIBRARY", stub_library);
   // the stub stays mapped once loaded (RTLD_NODELETE): its counter may include the calls of previous tests
   const int previous_calls = stub_ma27id_calls().value_or(0);
   const auto solver = SymmetricIndefiniteLinearSolverFactory::create(make_options("MA27"));
   ASSERT_NE(solver, nullptr);
#ifdef _WIN32
   GTEST_SKIP() << "the stub call counter is only read on POSIX systems";
#else
   // the MA27 constructor calls ma27id once, through the function pointer resolved in the stub
   EXPECT_EQ(stub_ma27id_calls(), std::optional<int>{previous_calls + 1});
#endif
}

TEST_F(HSLRuntimeLoading, EnvironmentVariableIsHonored) {
   // no libhsl_path option: the library is taken from UNO_HSL_LIBRARY
   {
      const ScopedEnvironmentVariable environment("UNO_HSL_LIBRARY", stub_library);
      EXPECT_TRUE(is_available("MA27"));
      EXPECT_NE(SymmetricIndefiniteLinearSolverFactory::create(make_options("MA27")), nullptr);
   }
   unload_hsl_library();
   {
      const ScopedEnvironmentVariable environment("UNO_HSL_LIBRARY", missing_library);
      EXPECT_FALSE(is_available("MA27"));
   }
}

TEST_F(HSLRuntimeLoading, LibhslPathOptionTakesPrecedence) {
   // the libhsl_path option wins over a (missing) library in UNO_HSL_LIBRARY
   const ScopedEnvironmentVariable environment("UNO_HSL_LIBRARY", missing_library);
   EXPECT_NE(SymmetricIndefiniteLinearSolverFactory::create(make_options("MA27", stub_library)), nullptr);
   EXPECT_TRUE(ma27_symbols_available());
}

TEST_F(HSLRuntimeLoading, FirstSuccessfulLoadWins) {
   const ScopedEnvironmentVariable environment("UNO_HSL_LIBRARY", missing_library);
   ASSERT_TRUE(load_hsl_library(stub_library));
   // a later request for another library keeps the loaded one (and warns)
   EXPECT_TRUE(load_hsl_library(missing_library));
   EXPECT_TRUE(ma27_symbols_available());
}

TEST_F(HSLRuntimeLoading, LibraryCanBeReloaded) {
   const ScopedEnvironmentVariable environment("UNO_HSL_LIBRARY", missing_library);
   ASSERT_TRUE(load_hsl_library(stub_library));
   unload_hsl_library();
   EXPECT_FALSE(ma27_symbols_available());
   ASSERT_TRUE(load_hsl_library(stub_library));
   EXPECT_TRUE(ma27_symbols_available());
}
