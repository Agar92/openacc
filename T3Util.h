#pragma once
#ifndef T3UTIL_H
#define T3UTIL_H

#include "T3Defs.h"
#include "T3TypeTraits.h"

#include <iostream>
#include <string>
#include <thread>

#define SNITCHDEFALTCONSTRUCTORS(CLASSNAME)                                    \
  CLASSNAME(CLASSNAME &&) { ConstructorSnitch::reportMove(#CLASSNAME); }       \
  CLASSNAME(CLASSNAME const &) { ConstructorSnitch::reportCopy(#CLASSNAME); }  \
  CLASSNAME() { ConstructorSnitch::reportParameterless(#CLASSNAME); }          \
  CLASSNAME &operator=(CLASSNAME &&) = default;                                \
  CLASSNAME &operator=(CLASSNAME const &) = default;

namespace t3 {

auto ReadColumns(std::string const &file, uint64_t nColumns = 2)
    -> t3::vector<t3::vector<double>>;

template <typename RandomEngine>
auto generateSubCanonical(RandomEngine &engine);

template <bool reverse = false, uint8_t debugLevel = 0, typename KeyGetter,
          typename OutputMutator>
auto BucketSort(KeyGetter aKeyGetter, OutputMutator aMutator, uint64_t N,
                size_t maxKey) {
  static_assert(
      t3::HasAccessor<size_t, KeyGetter>,
      "matID must have operator[](size_t) returning type assignable to size_t");
  static_assert(t3::HasMutator<size_t, OutputMutator>,
                "output must have operator[](size_t) returning a reference to "
                "a type assignable from size_t");
  auto const numKeys = maxKey + 1;
  auto getKey = [aKeyGetter, maxKey](size_t ind) {
    auto const val = aKeyGetter[ind];
    if constexpr (reverse)
      return maxKey - val;
    else
      return val;
  };
  static auto const numThreads = std::thread::hardware_concurrency();
  uint64_t constexpr bigNumber = 16384ul;
  uint64_t const numChunks = [N]() {
    if (N < bigNumber)
      return 1u;
    else
      return numThreads;
  }();
  if constexpr (debugLevel > 0) {
    std::cout << "t3::BucketSort: N = " << N << ", numThreads = " << numThreads
              << ", numChunks = " << numChunks << ", numKeys = " << numKeys
              << std::endl;
  }
  t3::vector<size_t> shifts(numChunks * numKeys, 0u);
  auto const divider = N / numChunks;
  auto const R = N % numChunks;
  auto const firstRChunksSize = divider + 1;
  auto const otherChunksSize = divider;
  auto const chunkStart = [firstRChunksSize, otherChunksSize,
                           R](size_t chunkNum) {
    if (chunkNum <= R) {
      return firstRChunksSize * chunkNum;
    } else {
      return firstRChunksSize * R + otherChunksSize * (chunkNum - R);
    }
  };
  auto const chunkEnd = [chunkStart, firstRChunksSize, otherChunksSize,
                         R](size_t chunkNum) {
    return chunkStart(chunkNum) +
           (chunkNum < R ? firstRChunksSize : otherChunksSize);
  };
  auto indexInShifts = [numKeys](auto keyVal, auto chunkNum) {
    return chunkNum * numKeys + keyVal;
  };

#pragma omp parallel for
  for (auto chunkNum = 0u; chunkNum < numChunks; ++chunkNum) {
    auto const start = chunkStart(chunkNum);
    auto const end = chunkEnd(chunkNum);
    // no point in simd - more memory access than computation
    for (auto ind = start; ind < end; ++ind) {
      auto const keyVal = getKey(ind);
      ++shifts[indexInShifts(keyVal, chunkNum)];
    }
  }

  if (debugLevel > 1) {
    std::cout << "t3::BucketSort: counts:";
    for (auto keyVal = 0u; keyVal < numKeys; ++keyVal) {
      std::cout << "\n";
      for (auto chunkNum = 0u; chunkNum < numChunks; ++chunkNum) {
        std::cout << shifts[indexInShifts(keyVal, chunkNum)] << " ";
      }
    }
    std::cout << std::endl;
  }

  // FIXME loop dependency
  auto accumulate = 0u;
  for (auto keyVal = 0u; keyVal < numKeys; ++keyVal) {
    for (auto chunkNum = 0u; chunkNum < numChunks; ++chunkNum) {
      auto const value = shifts[indexInShifts(keyVal, chunkNum)];
      shifts[indexInShifts(keyVal, chunkNum)] = accumulate;
      accumulate += value;
    }
  }

  if (debugLevel > 1) {
    std::cout << "t3::BucketSort: shifts:";
    for (auto keyVal = 0u; keyVal < numKeys; ++keyVal) {
      std::cout << "\n";
      for (auto chunkNum = 0u; chunkNum < numChunks; ++chunkNum) {
        std::cout << shifts[indexInShifts(keyVal, chunkNum)] << " ";
      }
    }
    std::cout << std::endl;
  }

#pragma omp parallel for
  for (auto chunkNum = 0u; chunkNum < numChunks; ++chunkNum) {
    auto const start = chunkStart(chunkNum);
    auto const end = chunkEnd(chunkNum);
    // no point in simd - more memory access than computation
    for (auto ind = start; ind < end; ++ind) {
      auto const keyVal = getKey(ind);
      auto &shift = shifts[indexInShifts(keyVal, chunkNum)];
      aMutator[shift] = ind;
      ++shift;
    }
  }

  if (debugLevel > 2) {
    std::cout << "t3::BucketSort: outputs:";
    for (auto ind = 0u; ind < N; ++ind) {
      std::cout << " " << aMutator[ind];
    }
    std::cout << std::endl;
    std::cout << "t3::BucketSort: keys:";
    for (auto ind = 0u; ind < N; ++ind) {
      std::cout << " " << aKeyGetter[aMutator[ind]];
    }
    std::cout << std::endl;
  }
}

namespace ConstructorSnitch {
auto reportCopy(std::string name) -> void;
auto reportMove(std::string name) -> void;
auto reportParameterless(std::string name) -> void;
} // namespace ConstructorSnitch

template <typename Floating, typename RandomEngine>
auto GenerateSubCanonical(RandomEngine &engine) {
  auto const randomNumber = engine();
  auto maxRandomNumber = engine.max();
  auto maxRandomNumberAsFloating =
      static_cast<Floating>(maxRandomNumber);
  auto const result =
      static_cast<Floating>(randomNumber) / maxRandomNumberAsFloating;
  return result;
}

} // namespace t3

#endif // T3UTIL_H
