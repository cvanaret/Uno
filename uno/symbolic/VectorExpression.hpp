// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTOREXPRESSION_H
#define UNO_VECTOREXPRESSION_H

#include "symbolic_traits.hpp"

namespace uno {
   // A symbolic sparse vector: it has size() structural entries; entry k has
   // global index indices[k] and value component_function(indices[k]).
   // `indices` is typically a std::vector<size_t> or a uno::Range.
   //
   // Storage policy:
   //   - indices: reference for lvalues (no copy of a large vector),
   //              owned by value for rvalues (e.g. a temporary Range) — no dangling.
   //   - component_function: always owned by value (a lambda is cheap to move),
   //              so it never dangles even if built from a local.
   template <typename Indices, typename Callable>
   class VectorExpression {
      // type of indices.begin(); const because we only ever read through const members
      using inner_iterator = decltype(std::declval<const std::remove_reference_t<Indices>&>().begin());

   public:
      class iterator {
      public:
         using value_type        = std::pair<size_t, double>;
         using reference         = std::pair<size_t, double>;
         using pointer           = void;
         using difference_type   = std::ptrdiff_t;
         using iterator_category = std::input_iterator_tag;

         iterator(const VectorExpression& expression, inner_iterator inner):
            expression(&expression), inner(std::move(inner)) { }

         [[nodiscard]] std::pair<size_t, double> operator*() const {
            const size_t global_index = *this->inner;              // element of `indices`
            return {global_index, this->expression->component_function(global_index)};
         }

         iterator& operator++() {
            ++this->inner;
            return *this;
         }

         [[nodiscard]] friend bool operator==(const iterator& a, const iterator& b) {
            return a.inner == b.inner;
         }
         [[nodiscard]] friend bool operator!=(const iterator& a, const iterator& b) {
            return a.inner != b.inner;
         }

      private:
         const VectorExpression* expression;
         inner_iterator inner;
      };

      // for algorithms that query the element type
      using value_type = double;

      // CTAD picks Indices/Callable via the deduction guide below:
      //   lvalue arg -> Indices = const T&  (Indices&& collapses to const T&, stored as reference)
      //   rvalue arg -> Indices = T         (Indices&& is T&&, moved into an owned member)
      VectorExpression(Indices&& indices, Callable&& component_function):
         indices(std::forward<Indices>(indices)),
         component_function(std::forward<Callable>(component_function)) { }

      [[nodiscard]] size_t size() const noexcept {
         return this->indices.size();
      }

      // Evaluate the component function at a GLOBAL index (the index space the
      // callable expects), not at a structural position. Provided for parity with
      // the original; norms consume this expression through begin()/end().
      [[nodiscard]] double operator[](size_t global_index) const {
         return this->component_function(global_index);
      }

      [[nodiscard]] iterator begin() const noexcept { return iterator(*this, this->indices.begin()); }
      [[nodiscard]] iterator end()   const noexcept { return iterator(*this, this->indices.end()); }

   private:
      Indices indices;                       // const T& for lvalues, T for rvalues
      std::decay_t<Callable> component_function;   // owned by value
   };

   // The crucial piece: first parameter is a forwarding reference, so value
   // category is preserved into `Indices` (reference for lvalues, value for rvalues).
   template <typename Indices, typename Callable>
   VectorExpression(Indices&&, Callable&&) -> VectorExpression<Indices, Callable>;
} // namespace

#endif // UNO_VECTOREXPRESSION_H