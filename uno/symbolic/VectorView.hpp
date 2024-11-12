// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VECTORVIEW_H
#define UNO_VECTORVIEW_H

namespace uno {
   // span of an arbitrary container: allocation-free view of a certain length
   template <typename Expression>
   class VectorView {
   public:
      /*
      class iterator {
      public:
         using value_type = std::pair<size_t, typename std::remove_reference_t<Expression>::value_type>;

         iterator(const VectorView<Expression>& view, size_t index) : view(view), index(index) { }

         value_type operator*() {
            // copy the element in the pair. Cheap only for trivial types
            return {this->index, this->view[this->index]};
         }
         
         iterator& operator++() { this->index++; return *this; }

         friend bool operator!=(const iterator& a, const iterator& b) { return &a.view != &b.view || a.index != b.index; };

      private:
         const VectorView<Expression>& view;
         size_t index;
      };
   */
      using value_type = typename std::remove_reference_t<Expression>::value_type;

      VectorView(Expression&& expression, size_t start, size_t end):
            expression(std::forward<Expression>(expression)), start(start), end(std::min(end, expression.size())) {
         if (end < start) {
            throw std::runtime_error("The view ends before its starting point");
         }
      }

      // preconditions: expression != nullptr, i < length
      [[nodiscard]] value_type& operator[](size_t index) noexcept { return this->expression[this->start + index]; }
      [[nodiscard]] value_type operator[](size_t index) const noexcept { return this->expression[this->start + index]; }
      [[nodiscard]] size_t size() const noexcept { return this->end - this->start; }

      // [[nodiscard]] iterator begin() const noexcept { return iterator(*this, 0); }
      // [[nodiscard]] iterator end() const noexcept { return iterator(*this, this->length); }

   protected:
      Expression expression;
      const size_t start;
      const size_t end;
   };

   // free function
   template <typename Expression>
   VectorView<Expression> view(Expression&& expression, size_t start, size_t end) {
      return {std::forward<Expression>(expression), start, end};
   }

   template <typename Expression>
   std::ostream& operator<<(std::ostream& stream, const VectorView<Expression>& view) {
      for (size_t index: Range(view.size())) {
         stream << view[index] << " ";
      }
      stream << '\n';
      return stream;
   }
} // namespace

#endif //UNO_VECTORVIEW_H
