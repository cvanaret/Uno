#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include "cppasl/nl_header.hpp"
#include "cppasl/problem_classification.hpp"

namespace cppasl {

   namespace detail {
      struct ModelData;
   }

   struct ReaderOptions {
      /// Objective whose Hessian enters the Lagrangian Hessian (structure and values); -1 for none.
      int hessian_objective_index{0};
      /// If true (default), every polynomial part of degree <= 2 is stored as a lower-triangular COO Hessian and
      /// evaluated by sparse matrix-vector products. If false, degree-2 terms are differentiated by the AD tapes
      /// (useful to cross-validate the two code paths).
      bool detect_quadratic_structure{true};
   };

   enum class SuffixTarget { Variables = 0, Constraints = 1, Objectives = 2, Problem = 3 };

   struct Suffix {
      std::string name;
      SuffixTarget target{SuffixTarget::Variables};
      bool is_real{false};
      std::vector<double> values; ///< dense (one entry per target item), integer suffixes are stored exactly
   };

   /// Algebraic structure of one objective or constraint
   enum class FunctionStructure { Linear, Quadratic, Nonlinear };

   /// Read-only view of the Hessian H of the quadratic part (1/2 x^T H x) of one function, lower triangle, COO.
   /// Diagonal entries come first, then the strictly lower entries sorted by (column, row).
   struct QuadraticFormView {
      std::size_t number_nonzeros{0};
      std::size_t number_diagonal_nonzeros{0};
      const int* row_indices{nullptr};
      const int* column_indices{nullptr};
      const double* values{nullptr};
   };

   class EvaluationWorkspace;
   class EvaluationKernels;

   /// An optimization model read from an .nl file.
   /// The model is immutable once read: all evaluation methods are const and thread-safe provided that each thread
   /// uses its own EvaluationWorkspace (the ASL2 "EvalWorkspace" design). Copies share the underlying data.
   ///
   /// Each function is decomposed at read time into
   ///    f(x) = constant + linear^T x + 1/2 x^T H x + sum_e scale_e * element_e(x)
   /// where H is stored as a lower COO matrix and the elements are the non-quadratic parts of the top-level sum,
   /// compiled into flat AD tapes (partially separable structure).
   ///
   /// Conventions: indices are 0-based; the Jacobian is stored row-wise (CSR); the Lagrangian Hessian
   /// sigma * H_f + sum_i y_i * H_ci is the lower triangle stored column-wise (CSC, row >= column).
   class NlModel {
   public:
      static NlModel read_file(const std::string& file_name, const ReaderOptions& options = {});
      static NlModel read_from_memory(const char* contents, std::size_t size, const ReaderOptions& options = {});

      [[nodiscard]] const NlHeader& header() const;
      [[nodiscard]] int number_variables() const;
      [[nodiscard]] int number_constraints() const;
      [[nodiscard]] int number_objectives() const;

      // bounds and initial point (infinite bounds are +/- std::numeric_limits<double>::infinity())
      [[nodiscard]] const std::vector<double>& variable_lower_bounds() const;
      [[nodiscard]] const std::vector<double>& variable_upper_bounds() const;
      [[nodiscard]] const std::vector<double>& constraint_lower_bounds() const;
      [[nodiscard]] const std::vector<double>& constraint_upper_bounds() const;
      [[nodiscard]] const std::vector<double>& initial_primal_point() const;
      [[nodiscard]] const std::vector<double>& initial_dual_point() const;
      /// variable complementary to each constraint (-1 if none)
      [[nodiscard]] const std::vector<int>& complementary_variables() const;
      [[nodiscard]] bool is_integer_variable(int variable_index) const;
      [[nodiscard]] bool is_maximization(int objective_index = 0) const;
      [[nodiscard]] const std::vector<Suffix>& suffixes() const;

      // algebraic structure
      [[nodiscard]] FunctionStructure objective_structure(int objective_index = 0) const;
      [[nodiscard]] FunctionStructure constraint_structure(int constraint_index) const;
      [[nodiscard]] QuadraticFormView objective_quadratic_form(int objective_index = 0) const;
      [[nodiscard]] QuadraticFormView constraint_quadratic_form(int constraint_index) const;
      /// Problem class (LP/QP/QCQP/NLP) and convexity (LAPACK Cholesky on the quadratic forms, see options)
      [[nodiscard]] ProblemClassification classify(const ClassificationOptions& options = {}) const;

      // Jacobian sparsity (CSR)
      [[nodiscard]] std::size_t number_jacobian_nonzeros() const;
      [[nodiscard]] const std::vector<std::size_t>& jacobian_row_starts() const;
      [[nodiscard]] const std::vector<int>& jacobian_row_indices() const; ///< COO rows, convenience
      [[nodiscard]] const std::vector<int>& jacobian_column_indices() const;
      // Lagrangian Hessian sparsity (lower triangle, CSC)
      [[nodiscard]] int hessian_objective_index() const;
      [[nodiscard]] std::size_t number_hessian_nonzeros() const;
      [[nodiscard]] const std::vector<std::size_t>& hessian_column_starts() const;
      [[nodiscard]] const std::vector<int>& hessian_row_indices() const;
      [[nodiscard]] const std::vector<int>& hessian_column_indices() const; ///< COO columns, convenience

      // evaluations
      [[nodiscard]] double evaluate_objective(EvaluationWorkspace& workspace, const double* x, int objective_index = 0) const;
      void evaluate_objective_gradient(EvaluationWorkspace& workspace, const double* x, double* gradient,
         int objective_index = 0) const;
      void evaluate_constraints(EvaluationWorkspace& workspace, const double* x, double* constraints) const;
      void evaluate_jacobian(EvaluationWorkspace& workspace, const double* x, double* jacobian_values) const;
      /// values of sigma * Hessian(objective) + sum_i y_i * Hessian(constraint i), in the order of the structure
      void evaluate_lagrangian_hessian(EvaluationWorkspace& workspace, const double* x, double objective_multiplier,
         const double* constraint_multipliers, double* hessian_values) const;
      /// result = (sigma * Hessian(objective) + sum_i y_i * Hessian(constraint i)) * vector (matrix-free)
      void evaluate_lagrangian_hessian_vector_product(EvaluationWorkspace& workspace, const double* x,
         double objective_multiplier, const double* constraint_multipliers, const double* vector, double* result) const;

      // linear operators on evaluated values
      /// result = J * vector (size m)
      void multiply_jacobian(const double* jacobian_values, const double* vector, double* result) const;
      /// result = J^T * vector (size n)
      void multiply_jacobian_transpose(const double* jacobian_values, const double* vector, double* result) const;
      /// result = W * vector where W is given by its lower triangle (hessian_values)
      void multiply_hessian(const double* hessian_values, const double* vector, double* result) const;

      /// Writes an AMPL solution file (text .sol, as ASL's write_sol): message, options, duals y (size m, may be
      /// null), primal x (size n, may be null), the solve result code (AMPL's solve_result_num) and output suffixes
      /// (dense values, only the nonzero entries are written).
      void write_solution(const std::string& file_name, const std::string& message, const double* x, const double* y,
         int solve_result_code = 0, const std::vector<Suffix>& output_suffixes = {}) const;

   private:
      explicit NlModel(std::shared_ptr<const detail::ModelData> data);
      std::shared_ptr<const detail::ModelData> data;
      friend class EvaluationWorkspace;
   };

   /// Per-thread scratch memory for the evaluations (values of the AD tapes are cached per function and point).
   class EvaluationWorkspace {
   public:
      explicit EvaluationWorkspace(const NlModel& model);

   private:
      friend class NlModel;
      friend class EvaluationKernels;
      std::shared_ptr<const detail::ModelData> data;
      std::vector<double> tape_values;          ///< forward values of all element tapes
      std::vector<double> defined_gradients;    ///< gradients of the shared defined variables
      std::vector<double> defined_values, defined_partials, defined_adjoints; ///< sweeps of their tapes
      std::vector<double> defined_tangents, defined_direction_tangents, defined_second_order_adjoints, defined_second_partials;
      std::vector<double> defined_weights; ///< Hessian: sum of the multiplier-weighted adjoints of each defined variable
      std::uint64_t defined_point_version{0}, defined_gradient_version{0};
      std::vector<double> tape_partials;        ///< first partial derivatives stored by the forward sweeps (2 per instruction)
      std::vector<std::uint64_t> function_point_version; ///< version of x at which a function's tapes were evaluated
      std::vector<double> current_point;
      std::uint64_t point_version{1};
      // single-element scratch
      std::vector<double> partials;             ///< 5 per instruction (see tape.hpp)
      std::vector<double> adjoints;
      std::vector<double> tangents;
      std::vector<double> second_order_adjoints;
      std::vector<double> local_direction;
      std::vector<double> element_hessian;      ///< packed lower triangle of one element Hessian
      std::vector<std::int32_t> direction_of_local, local_of_direction; ///< group Hessians: operand variables
   };

} // namespace cppasl
