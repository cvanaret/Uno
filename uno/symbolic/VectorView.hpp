// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_VIEW_H
#define UNO_VIEW_H

// span of an arbitrary container: allocation-free view of a certain length
template <typename Expression>
class VectorView {
public:
   using value_type = typename std::remove_reference_t<Expression>::value_type;

   class Iterator {
   public:
      Iterator(const VectorView<Expression>& view, size_t index) : view(view), index(index) {
      }
      std::pair<size_t, value_type> operator*() {
         // copy the element in the pair. Cheap only for trivial types
         return {this->index, this->view[this->index]};
      }
      // prefix increment
      Iterator& operator++() { this->index++; return *this; }

      friend bool operator== (const Iterator& a, const Iterator& b) { return &(a.view) == &(b.view) && a.index == b.index; };
      friend bool operator!= (const Iterator& a, const Iterator& b) { return &(a.view) != &(b.view) || a.index != b.index; };

   private:
      const VectorView<Expression>& view;
      size_t index;
   };

   VectorView(Expression&& expression, size_t length) noexcept;

   // preconditions: expression != nullptr, i < length
   [[nodiscard]] const value_type& operator[](size_t index) const noexcept { return this->expression[index]; }
   [[nodiscard]] size_t size() const noexcept { return this->length; }

   [[nodiscard]] Iterator begin() const noexcept { return Iterator(*this, 0); }
   [[nodiscard]] Iterator end() const noexcept { return Iterator(*this, this->length); }

protected:
   Expression expression;
   const size_t length;
};

template <typename Expression>
VectorView<Expression>::VectorView(Expression&& expression, size_t length) noexcept:
      expression(std::forward<Expression>(expression)), length(std::min(length, expression.size())) {
}

// free function
template <typename Expression>
VectorView<Expression> view(Expression&& expression, size_t length) {
   return {std::forward<Expression>(expression), length};
}

#endif //UNO_VIEW_H