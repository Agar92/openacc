#pragma once
#ifndef T3ARRAYVIEW_H
#define T3ARRAYVIEW_H

#include "T3Defs.h"
#include "T3TypeTraits.h"

namespace t3 {

template <typename T, typename MemberGetter> class ArrayStructureMemberView;
// for type deduction before C++17
template <typename T, typename MemberGetter>
auto makeAOSMemberView(T *data, MemberGetter getter);

template <typename T> class ConstantView;
template <typename T> auto makeConstantView(T constant);

template <typename T> class ConstantView {
  friend auto makeConstantView<T>(T constant);

public:
  auto operator[](size_t) const { return fConstant; }

private:
  ConstantView(T constant) : fConstant(constant) {}
  T const fConstant;
};

template <typename T, typename MemberGetter> class ArrayStructureMemberView {
  friend auto makeAOSMemberView<T, MemberGetter>(T *data, MemberGetter getter);

public:
  auto &operator[](size_t ind) { return fGetter(fData[ind]); }
  auto operator[](size_t ind) const { return fGetter(fData[ind]); }

private:
  ArrayStructureMemberView(T *data, MemberGetter getter)
      : fData(data), fGetter(getter) {}
  T *const fData;
  MemberGetter fGetter;
};

template <typename T> auto makeConstantView(T constant) {
  return ConstantView<T>(constant);
}

template <typename T, typename MemberGetter>
auto makeAOSMemberView(T *data, MemberGetter getter) {
  return ArrayStructureMemberView<T, MemberGetter>(data, getter);
}

} // namespace t3

#endif // T3ARRAYVIEW_H
