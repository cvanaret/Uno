// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMBOLIC_TRAITS_H
#define UNO_SYMBOLIC_TRAITS_H

#include <type_traits>
#include <utility>

namespace uno {
    template<typename T, typename = void>
    struct has_value_type: std::false_type {};

    template<typename T>
    struct has_value_type<T, std::void_t<typename T::value_type>>: std::true_type {};

    template<class T>
    using storage_t = std::conditional_t<std::is_lvalue_reference_v<T>, T, std::decay_t<T>>;

    #define UNO_FORWARD_ACCESSOR(name, member)                     \
    constexpr decltype(auto) name() & noexcept {                   \
      return (member);                                            \
    }                                                              \
    \
    constexpr decltype(auto) name() const& noexcept {              \
      return (member);                                            \
    }                                                              \
    \
    constexpr decltype(auto) name() && noexcept {                  \
      return std::move(member);                                   \
    }                                                              \
    \
    constexpr decltype(auto) name() const&& noexcept {             \
      return std::move(member);                                   \
    }
} // namespace

#endif // UNO_SYMBOLIC_TRAITS_H