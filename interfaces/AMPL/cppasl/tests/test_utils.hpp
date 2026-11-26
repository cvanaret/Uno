#pragma once
// Minimal self-contained test harness and helpers (no external dependency).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include "cppasl/nl_model.hpp"
#include "nl_writer.hpp"

namespace testing {

   struct TestCase {
      const char* name;
      std::function<void()> body;
   };
   std::vector<TestCase>& registry();
   extern int failed_checks;
   extern int total_checks;

   struct Registrar {
      Registrar(const char* name, std::function<void()> body) { registry().push_back({name, std::move(body)}); }
   };

   inline bool report(bool condition, const char* expression, const char* file, int line, const std::string& details = {}) {
      ++total_checks;
      if (!condition) {
         ++failed_checks;
         std::fprintf(stderr, "  FAILED %s:%d: %s %s\n", file, line, expression, details.c_str());
      }
      return condition;
   }

   inline bool close(double a, double b, double tolerance) {
      return std::fabs(a - b) <= tolerance * std::max(1., std::max(std::fabs(a), std::fabs(b)));
   }

   /// writes the problem (optionally to $CPPASL_TEST_OUTPUT/<name>.nl for external cross-checks) and reads it back
   cppasl::NlModel write_and_read(const nlwriter::NlProblem& problem, nlwriter::Format format, const std::string& name,
      const cppasl::ReaderOptions& options = {});

   /// finite-difference checks of all derivatives, and consistency of the Hessian, HVP and linear operators;
   /// returns the largest relative error
   double check_derivatives(const cppasl::NlModel& model, const std::vector<double>& x, double objective_multiplier,
      const std::vector<double>& multipliers, double tolerance, const char* file, int line);

   /// dense lower-triangular Lagrangian Hessian from the sparse values
   std::vector<double> dense_hessian(const cppasl::NlModel& model, const std::vector<double>& values);

   inline std::vector<double> random_vector(std::mt19937_64& random, std::size_t size, double low = -1., double high = 1.) {
      std::uniform_real_distribution<double> distribution(low, high);
      std::vector<double> v(size);
      for (double& value: v) value = distribution(random);
      return v;
   }
} // namespace testing

#define TEST_CASE(name) \
   static void name(); \
   static testing::Registrar name##_registrar(#name, name); \
   static void name()
#define CHECK(condition) testing::report((condition), #condition, __FILE__, __LINE__)
#define CHECK_CLOSE(a, b, tolerance) \
   testing::report(testing::close((a), (b), (tolerance)), #a " ~ " #b, __FILE__, __LINE__, \
      (std::ostringstream() << "(" << (a) << " vs " << (b) << ")").str())
#define CHECK_DERIVATIVES(model, x, sigma, y, tolerance) \
   testing::check_derivatives((model), (x), (sigma), (y), (tolerance), __FILE__, __LINE__)
