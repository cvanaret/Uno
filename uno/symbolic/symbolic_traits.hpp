// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMBOLIC_TRAITS_H
#define UNO_SYMBOLIC_TRAITS_H

#include <type_traits>
#include <utility>

namespace uno {
   // empty marker base for all symbolic expression templates
   struct SymbolicExpression {
   };

   template<typename T>
   inline constexpr bool is_symbolic_v = std::is_base_of_v<SymbolicExpression, std::remove_cv_t<std::remove_reference_t<T>>>;

   template<typename T, typename = void>
   struct has_value_type : std::false_type {
   };

   template<typename T>
   struct has_value_type<T, std::void_t<typename T::value_type> > : std::true_type {
   };

   template<class T>
   using storage_t = std::conditional_t<std::is_lvalue_reference_v<T>, T, std::decay_t<T> >;

   namespace detail {
      // forwards a member out of an rvalue expression: an owned member (Declared is not a reference) is moved,
      // a reference member refers to an object the expression does not own, so it stays an lvalue
      template <typename Declared, typename Member>
      constexpr decltype(auto) forward_member(Member& member) noexcept {
         if constexpr (std::is_lvalue_reference_v<Declared>) {
            return static_cast<Declared>(member);
         }
         else {
            return std::move(member);
         }
      }
   } // namespace detail

   #define UNO_FORWARD_ACCESSOR(name, member)                                 \
   constexpr decltype(auto) name() & noexcept {                               \
      return (member);                                                        \
   }                                                                          \
                                                                              \
   constexpr decltype(auto) name() const& noexcept {                          \
      return (member);                                                        \
   }                                                                          \
                                                                              \
   constexpr decltype(auto) name() && noexcept {                              \
      return ::uno::detail::forward_member<decltype(member)>(member);         \
   }                                                                          \
                                                                              \
   constexpr decltype(auto) name() const&& noexcept {                         \
      return ::uno::detail::forward_member<decltype(member)>(member);         \
   }
} // namespace

#endif // UNO_SYMBOLIC_TRAITS_H