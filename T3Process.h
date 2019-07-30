#pragma once
#ifndef T3PROCESS_H
#define T3PROCESS_H

#include <array>
#include <cstdint>
#include <memory>
#include <tuple>
#include <utility>

#include "T3Defs.h"
#include "T3LorentzVector.h"
#include "T3TypeTraits.h"

namespace t3 {

template <typename ProcessImpl> class Process {
public:
  using Base_t = ProcessImpl;
  Process(Process &&) = default;
  Process(Process const &) = default;
  Process(ProcessImpl &&aProcessImpl)
      : fProcessImpl(std::make_unique<ProcessImpl>(std::move(aProcessImpl))) {}
  template <typename... Args>
  Process(Args &&... args)
      : fProcessImpl(
            std::make_unique<ProcessImpl>(std::forward<Args>(args)...)) {}
  Process &operator=(Process &&) = default;
  template <typename Floating>
  auto GetCS(Floating e, PDG_t incPDG, MatID_t matID) const;
  template <typename ArrayViewE, typename ArrayViewPDG, typename ArrayViewMat,
            typename ArrayViewOut>
  auto GetCS(ArrayViewE const e, ArrayViewPDG const incPDG,
             ArrayViewMat const matID, uint64_t N, ArrayViewOut output) const
      -> bool;
  template <typename ArrayViewP, typename ArrayViewPDG, typename ArrayViewMat,
            typename ArrayViewRng, typename TupleOfArrayViewsPDG,
            typename TupleOfArrayViewsP>
  auto GetFS(ArrayViewP const p, ArrayViewPDG const incPDG,
             ArrayViewMat const matID, ArrayViewRng randomEngine, uint64_t N,
             TupleOfArrayViewsPDG &outPDG, TupleOfArrayViewsP &outP) const
      -> bool;

private:
  std::unique_ptr<ProcessImpl> fProcessImpl;

private:
  template <typename ArrayViewP, typename ArrayViewPDG, typename ArrayViewMat,
            typename ArrayViewRng, typename TupleOfArrayViewsPDG,
            typename TupleOfArrayViewsP, size_t... Is>
  auto GetFSImpl(ArrayViewP const p, ArrayViewPDG const incPDG,
                 ArrayViewMat const matID, ArrayViewRng randomEngine,
                 uint64_t N, TupleOfArrayViewsPDG &outPDG,
                 TupleOfArrayViewsP &outP, std::index_sequence<Is...>) const
      -> bool;
};

//______________________________________________________________________________

template <class ProcessImpl>
template <typename Floating>
auto Process<ProcessImpl>::GetCS(Floating e, PDG_t incPDG,
                                 MatID_t matID) const {
  return fProcessImpl->GetCS(e, incPDG, matID);
}

template <class ProcessImpl>
template <typename ArrayViewE, typename ArrayViewPDG, typename ArrayViewMat,
          typename ArrayViewOut>
auto Process<ProcessImpl>::GetCS(ArrayViewE const e, ArrayViewPDG const incPDG,
                                 ArrayViewMat const matID, uint64_t N,
                                 ArrayViewOut output) const -> bool {
  static_assert(t3::HasAccessor<double, ArrayViewE> ||
                    std::is_floating_point<ArrayViewE>::value,
                "e must be a floating point value or have operator[](size_t) "
                "returning a value convertible to double");
  static_assert(t3::HasAccessor<t3::PDG_t, ArrayViewPDG> ||
                    std::is_same<t3::PDG_t, ArrayViewPDG>::value,
                "incPDG must be a PDG_t value or have operator[](size_t) "
                "returning PDG_t");
  static_assert(t3::HasAccessor<uint64_t, ArrayViewMat> ||
                    std::is_same<t3::MatID_t, ArrayViewMat>::value,
                "matID must have operator[](size_t) returning MatID_t");
  static_assert(t3::HasMutator<double, ArrayViewOut>,
                "output must have operator[](size_t) returning a reference to "
                "a type assignable from double");
//
#ifdef OPENACC
#pragma acc parallel loop
#else
#pragma omp parallel for simd
#endif
//
  for (auto ind = 0u; ind < N; ++ind) {
    auto const ei = [e, ind]() {
      if constexpr (std::is_floating_point<ArrayViewE>::value) {
        return e;
      } else {
        return e[ind];
      }
    }();
    auto const incPDGi = [incPDG, ind]() {
      if constexpr (std::is_same<t3::PDG_t, ArrayViewPDG>::value) {
        return incPDG;
      } else {
        return incPDG[ind];
      }
    }();
    auto const matIDi = [matID, ind]() {
      if constexpr (std::is_same<t3::MatID_t, ArrayViewMat>::value) {
        return matID;
      } else {
        return matID[ind];
      }
    }();
    output[ind] = fProcessImpl->GetCS(ei, incPDGi, matIDi);
  }
  return true;
}

template <class ProcessImpl>
template <typename ArrayViewP, typename ArrayViewPDG, typename ArrayViewMat,
          typename ArrayViewRng, typename TupleOfArrayViewsPDG,
          typename TupleOfArrayViewsP>
auto Process<ProcessImpl>::GetFS(ArrayViewP const p, ArrayViewPDG const incPDG,
                                 ArrayViewMat const matID,
                                 ArrayViewRng randomEngine, uint64_t N,
                                 TupleOfArrayViewsPDG &outPDG,
                                 TupleOfArrayViewsP &outP) const -> bool {
  static_assert(t3::HasAccessor<t3::PDG_t, ArrayViewPDG> ||
                    std::is_same<t3::PDG_t, ArrayViewPDG>::value,
                "incPDG must be a PDG_t value or have operator[](size_t) "
                "returning PDG_t");
  static_assert(t3::HasAccessor<uint64_t, ArrayViewMat> ||
                    std::is_same<t3::MatID_t, ArrayViewMat>::value,
                "matID must have operator[](size_t) returning MatID_t");

  return GetFSImpl(
      p, incPDG, matID, randomEngine, N, outPDG, outP,
      std::make_index_sequence<std::tuple_size<TupleOfArrayViewsPDG>::value>{});
}

template <class ProcessImpl>
template <typename ArrayViewP, typename ArrayViewPDG, typename ArrayViewMat,
          typename ArrayViewRng, typename TupleOfArrayViewsPDG,
          typename TupleOfArrayViewsP, size_t... Is>
auto Process<ProcessImpl>::GetFSImpl(
    ArrayViewP const p, ArrayViewPDG const incPDG, ArrayViewMat const matID,
    ArrayViewRng randomEngine, uint64_t N, TupleOfArrayViewsPDG &outPDG,
    TupleOfArrayViewsP &outP, std::index_sequence<Is...>) const -> bool {
  static_assert(t3::HasMutator<t3::PDG_t, decltype(std::get<Is>(outPDG))...>,
                "output must have operator[](size_t) returning a reference to "
                "a type assignable from double");

//
#ifdef OPENACC
#pragma acc parallel loop
#else 
#pragma omp parallel for simd
#endif  
//
  for (auto ind = 0u; ind < N; ++ind) {
    auto const &pi = p[ind];
    auto const incPDGi = [incPDG, ind]() {
      if constexpr (std::is_same<t3::PDG_t, ArrayViewPDG>::value) {
        return incPDG;
      } else {
        return incPDG[ind];
      }
    }();
    auto const matIDi = [matID, ind]() {
      if constexpr (std::is_same<t3::MatID_t, ArrayViewMat>::value) {
        return matID;
      } else {
        return matID[ind];
      }
    }();
    auto &randomEnginei = randomEngine[ind];
    std::tie(std::get<Is>(outPDG)[ind]..., std::get<Is>(outP)[ind]...) =
        fProcessImpl->GetFS(pi, incPDGi, matIDi, randomEnginei);
  }
  return true;
}

} // namespace t3

#endif // T3PROCESS_H
