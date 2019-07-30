#pragma once
#ifndef T3TYPETRAITS_H
#define T3TYPETRAITS_H

#include <cstdint>
#include <type_traits>

namespace t3 {

template <typename...> using TryToInstantiate = void;
using DisregardThis = void;

template <template <typename...> class Expression, typename Attempt,
          typename... Ts>
struct IsDetectedImpl : std::false_type {};

template <template <typename...> class Expression, typename... Ts>
struct IsDetectedImpl<Expression, TryToInstantiate<Expression<Ts...>>, Ts...>
    : std::true_type {};

template <template <typename...> class Expression, typename... Ts>
constexpr auto IsDetected =
    IsDetectedImpl<Expression, DisregardThis, Ts...>::value;

// is_assignable yields true for stron typedef
template <typename ReturnType, typename Accessor>
using AccessExpression =
    decltype(std::declval<ReturnType &>() = std::declval<decltype(
                 std::declval<Accessor const>()[std::declval<uint64_t>()])>());

template <typename ValueType, typename Mutator>
using MutateExpression =
    decltype(std::declval<decltype(
                 std::declval<Mutator &>()[std::declval<uint64_t>()]) &>() =
                 std::declval<ValueType>());

template <typename ReturnType, typename Accessor, typename... Accessors>
constexpr auto HasAccessor =
    HasAccessor<ReturnType, Accessor> &&HasAccessor<ReturnType, Accessors...>;

template <typename ReturnType, typename Accessor>
constexpr auto HasAccessor<ReturnType, Accessor> =
    IsDetected<AccessExpression, ReturnType, Accessor>;

template <typename ValueType, typename Mutator, typename... Mutators>
constexpr auto HasMutator =
    HasMutator<ValueType, Mutator> &&HasMutator<ValueType, Mutators...>;

template <typename ValueType, typename Mutator>
constexpr auto HasMutator<ValueType, Mutator> =
    IsDetected<MutateExpression, ValueType, Mutator>;

} // namespace t3

#endif // T3TYPETRAITS_H
