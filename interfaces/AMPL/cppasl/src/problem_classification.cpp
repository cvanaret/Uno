// Classification LP / QP / QCQP / NLP and convexity detection of the quadratic forms.
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "cppasl/nl_model.hpp"
#include "cppasl/problem_classification.hpp"
#include "model_data.hpp"

#ifdef CPPASL_USE_LAPACK
extern "C" void dpotrf_(const char* uplo, const int* n, double* a, const int* lda, int* info, std::size_t uplo_length);
#endif

namespace cppasl {

   std::string ProblemClassification::to_string() const {
      std::string type_name;
      switch (this->type) {
         case ProblemType::LinearProgram: type_name = "LP"; break;
         case ProblemType::QuadraticProgram: type_name = "QP"; break;
         case ProblemType::QuadraticallyConstrainedQuadraticProgram: type_name = "QCQP"; break;
         case ProblemType::NonlinearProgram: type_name = "NLP"; break;
      }
      if (this->has_integer_variables) type_name = "MI" + type_name;
      switch (this->convexity) {
         case Convexity::Convex: return "convex " + type_name;
         case Convexity::Nonconvex: return "nonconvex " + type_name;
         default: return type_name;
      }
   }

   namespace {
      /// In-place Cholesky factorization of a dense column-major lower triangle. Returns false if not positive definite.
      bool cholesky_factorize(std::vector<double>& matrix, int dimension) {
#ifdef CPPASL_USE_LAPACK
         int info = 0;
         const char lower = 'L';
         dpotrf_(&lower, &dimension, matrix.data(), &dimension, &info, 1);
         return info == 0;
#else
         const auto n = static_cast<std::size_t>(dimension);
         for (std::size_t j = 0; j < n; ++j) {
            double* column_j = matrix.data() + j * n;
            for (std::size_t k = 0; k < j; ++k) { // left-looking update: column_j -= L(j,k) L(:,k)
               const double ljk = matrix[j + k * n];
               if (ljk == 0.) continue;
               const double* column_k = matrix.data() + k * n;
               for (std::size_t i = j; i < n; ++i) column_j[i] -= ljk * column_k[i];
            }
            if (!(column_j[j] > 0.)) return false;
            const double pivot = std::sqrt(column_j[j]);
            column_j[j] = pivot;
            for (std::size_t i = j + 1; i < n; ++i) column_j[i] /= pivot;
         }
         return true;
#endif
      }

      struct DisjointSets {
         std::vector<int> parent;
         explicit DisjointSets(std::size_t size): parent(size) { std::iota(this->parent.begin(), this->parent.end(), 0); }
         int find(int i) {
            while (this->parent[static_cast<std::size_t>(i)] != i) {
               int& p = this->parent[static_cast<std::size_t>(i)];
               p = this->parent[static_cast<std::size_t>(p)]; // path halving
               i = p;
            }
            return i;
         }
         void unite(int i, int j) {
            i = this->find(i);
            j = this->find(j);
            if (i != j) this->parent[static_cast<std::size_t>(std::max(i, j))] = std::min(i, j);
         }
      };

      /// Is sign * H positive semidefinite? H is a lower COO matrix. The support is split into its connected
      /// components; each one is tested by cheap certificates (negative diagonal, 2x2 minors, diagonal dominance),
      /// then by a Cholesky factorization of the shifted dense block.
      Convexity positive_semidefiniteness(const QuadraticFormView& form, double sign, const ClassificationOptions& options) {
         if (form.number_nonzeros == 0) return Convexity::Convex;
         // compress the support
         std::vector<int> support;
         support.reserve(2 * form.number_nonzeros);
         for (std::size_t q = 0; q < form.number_nonzeros; ++q) {
            support.push_back(form.row_indices[q]);
            support.push_back(form.column_indices[q]);
         }
         std::sort(support.begin(), support.end());
         support.erase(std::unique(support.begin(), support.end()), support.end());
         const std::size_t k = support.size();
         auto local = [&](int variable) {
            return static_cast<int>(std::lower_bound(support.begin(), support.end(), variable) - support.begin());
         };
         std::vector<int> local_rows(form.number_nonzeros), local_columns(form.number_nonzeros);
         std::vector<double> diagonal(k, 0.), off_diagonal_sums(k, 0.);
         double maximum_diagonal = 0.;
         for (std::size_t q = 0; q < form.number_nonzeros; ++q) {
            local_rows[q] = local(form.row_indices[q]);
            local_columns[q] = local(form.column_indices[q]);
            const double value = sign * form.values[q];
            if (local_rows[q] == local_columns[q]) {
               diagonal[static_cast<std::size_t>(local_rows[q])] += value;
            }
         }
         for (double d: diagonal) maximum_diagonal = std::max(maximum_diagonal, std::fabs(d));
         const double shift = options.semidefiniteness_tolerance * std::max(1., maximum_diagonal);

         // cheap certificates of indefiniteness: negative diagonal, or a 2x2 principal minor < 0
         for (double d: diagonal) if (d < -shift) return Convexity::Nonconvex;
         DisjointSets components(k);
         for (std::size_t q = 0; q < form.number_nonzeros; ++q) {
            const int i = local_rows[q], j = local_columns[q];
            if (i == j) continue;
            const double h = sign * form.values[q];
            const double dii = diagonal[static_cast<std::size_t>(i)] + shift, djj = diagonal[static_cast<std::size_t>(j)] + shift;
            if (h * h > dii * djj) return Convexity::Nonconvex;
            off_diagonal_sums[static_cast<std::size_t>(i)] += std::fabs(h);
            off_diagonal_sums[static_cast<std::size_t>(j)] += std::fabs(h);
            components.unite(i, j);
         }
         // certificate of semidefiniteness: diagonal dominance
         bool is_dominant = true;
         for (std::size_t i = 0; i < k && is_dominant; ++i) is_dominant = diagonal[i] + shift >= off_diagonal_sums[i];
         if (is_dominant) return Convexity::Convex;

         // group the local indices by component
         std::vector<int> component_of(k), position_in_component(k);
         std::vector<std::vector<int>> members;
         {
            std::vector<int> component_id(k, -1);
            for (std::size_t i = 0; i < k; ++i) {
               const auto root = static_cast<std::size_t>(components.find(static_cast<int>(i)));
               if (component_id[root] < 0) {
                  component_id[root] = static_cast<int>(members.size());
                  members.emplace_back();
               }
               component_of[i] = component_id[root];
               position_in_component[i] = static_cast<int>(members[static_cast<std::size_t>(component_id[root])].size());
               members[static_cast<std::size_t>(component_id[root])].push_back(static_cast<int>(i));
            }
         }
         // entries grouped by component
         std::vector<std::vector<std::size_t>> entries(members.size());
         for (std::size_t q = 0; q < form.number_nonzeros; ++q) {
            entries[static_cast<std::size_t>(component_of[static_cast<std::size_t>(local_rows[q])])].push_back(q);
         }
         Convexity result = Convexity::Convex;
         std::vector<double> dense;
         for (std::size_t c = 0; c < members.size(); ++c) {
            const std::size_t size = members[c].size();
            if (size == 1) continue; // 1x1: the diagonal was checked
            bool component_dominant = true;
            for (int i: members[c]) {
               component_dominant = component_dominant && diagonal[static_cast<std::size_t>(i)] + shift >= off_diagonal_sums[static_cast<std::size_t>(i)];
            }
            if (component_dominant) continue;
            if (size > static_cast<std::size_t>(options.maximum_dense_dimension)) {
               result = Convexity::Unknown;
               continue;
            }
            dense.assign(size * size, 0.);
            for (std::size_t q: entries[c]) {
               const auto i = static_cast<std::size_t>(position_in_component[static_cast<std::size_t>(local_rows[q])]);
               const auto j = static_cast<std::size_t>(position_in_component[static_cast<std::size_t>(local_columns[q])]);
               // lower triangle, column-major (the COO is lower: row >= column, and the order is preserved)
               dense[std::max(i, j) + std::min(i, j) * size] += sign * form.values[q];
            }
            for (std::size_t i = 0; i < size; ++i) dense[i + i * size] += shift;
            if (!cholesky_factorize(dense, static_cast<int>(size))) return Convexity::Nonconvex;
         }
         return result;
      }

      Convexity combine(Convexity a, Convexity b) {
         if (a == Convexity::Nonconvex || b == Convexity::Nonconvex) return Convexity::Nonconvex;
         if (a == Convexity::Unknown || b == Convexity::Unknown) return Convexity::Unknown;
         return Convexity::Convex;
      }
   } // namespace

   ProblemClassification NlModel::classify(const ClassificationOptions& options) const {
      const detail::ModelData& model = *this->data;
      ProblemClassification classification;
      classification.has_integer_variables =
         std::any_of(model.variable_is_integer.begin(), model.variable_is_integer.end(), [](char c) { return c != 0; });

      // the classified objective: the one of the Lagrangian Hessian, otherwise the first one
      const int objective = (model.hessian_objective_index >= 0) ? model.hessian_objective_index : ((model.number_objectives > 0) ? 0 : -1);
      bool has_elements = false, has_quadratic_constraints = false, has_quadratic_objective = false;
      if (objective >= 0) {
         has_elements = model.objective_block(objective).has_elements();
         has_quadratic_objective = model.objective_block(objective).has_quadratic_part();
      }
      has_elements = has_elements || !model.nonlinear_constraint_indices.empty();
      has_quadratic_constraints = !model.quadratic_constraint_indices.empty();

      if (has_elements) {
         classification.type = ProblemType::NonlinearProgram;
         classification.convexity = Convexity::Unknown;
         return classification;
      }
      classification.type = has_quadratic_constraints ? ProblemType::QuadraticallyConstrainedQuadraticProgram
         : (has_quadratic_objective ? ProblemType::QuadraticProgram : ProblemType::LinearProgram);

      // minimize f: f convex; maximize f: f concave. g(x) <= u: g convex; g(x) >= l: g concave; both: g affine
      Convexity convexity = Convexity::Convex;
      if (has_quadratic_objective) {
         const double sign = this->is_maximization(objective) ? -1. : 1.;
         convexity = positive_semidefiniteness(this->objective_quadratic_form(objective), sign, options);
      }
      for (int i: model.quadratic_constraint_indices) {
         if (convexity == Convexity::Nonconvex) break;
         const bool has_lower = std::isfinite(model.constraint_lower_bounds[static_cast<std::size_t>(i)]);
         const bool has_upper = std::isfinite(model.constraint_upper_bounds[static_cast<std::size_t>(i)]);
         if (has_lower && has_upper) convexity = Convexity::Nonconvex;
         else if (has_upper) convexity = combine(convexity, positive_semidefiniteness(this->constraint_quadratic_form(i), 1., options));
         else if (has_lower) convexity = combine(convexity, positive_semidefiniteness(this->constraint_quadratic_form(i), -1., options));
      }
      classification.convexity = convexity;
      return classification;
   }

} // namespace cppasl
