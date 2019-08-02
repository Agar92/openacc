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
#include "T3Particle.h"

namespace t3 {

template <typename ProcessImpl> class Process {
public:
  using Base_t = ProcessImpl;
  //Process(Process &&) = default;
  //Process(Process const &) = default;
  //Process(ProcessImpl &&aProcessImpl)
  //    : fProcessImpl(ProcessImpl(aProcessImpl)) {}
  template <typename... Args>
  Process(Args... args)
      : fProcessImpl(
            ProcessImpl(args...)) {}
  //Process &operator=(Process &) = default;
  template <typename Floating, typename ParticleTable>
  auto GetCS(Particle<Floating> * particles, MatID_t matID, uint64_t N,
             Floating * outputCSArray, ParticleTable * aParticleTable) const
      -> bool;
  //1 secondary particle:
  template <typename Floating>
  auto GetFS(Particle<Floating> * particles, MatID_t matID, uint64_t N,
             PDG_t * outPDG1, LorentzVector<Floating> * outP1,
             Floating * csBorderData) const
      -> bool;
  //2 secondary particles:
  template <typename Floating>
  auto GetFS(Particle<Floating> * particles, MatID_t matID, uint64_t N,
             PDG_t * outPDG1, LorentzVector<Floating> * outP1, PDG_t * outPDG2,
             LorentzVector<Floating> * outP2, Floating * csBorderData) const
      -> bool;

private:
  ProcessImpl fProcessImpl;

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
template <typename Floating, typename ParticleTable>
auto Process<ProcessImpl>::GetCS(Particle<Floating> * particles,
                                 MatID_t matID, uint64_t N,
                                 Floating * outputCSArray,
                                 ParticleTable * aParticleTable) const
  -> bool {         
//
#ifdef OPENACC
#pragma acc parallel loop
#else
#pragma omp parallel for simd
#endif
//
  for (auto ind = 0u; ind < N; ++ind) {
    outputCSArray[ind] = fProcessImpl.GetCS(particles[ind].p.E() - aParticleTable->GetMass(particles[ind].pdg), particles[ind].pdg, matID);
  }
  return true;
}

template <class ProcessImpl>
template <typename Floating>
auto Process<ProcessImpl>::GetFS(
    Particle<Floating> * particles, MatID_t matID, uint64_t N, PDG_t * outPDG1,
    LorentzVector<Floating> * outP1, Floating * csBorderData) const
    -> bool {
//
#ifdef OPENACC
#pragma acc parallel loop
#else 
#pragma omp parallel for simd
#endif  
//
  for (auto ind = 0u; ind < N; ++ind) {
    auto const fProcessImplTemporary=fProcessImpl.GetFS( particles[ind].p, particles[ind].pdg, matID, particles[ind].rs, csBorderData, ind);
    outPDG1[ind]=fProcessImplTemporary.pdg1;
    outP1[ind]=fProcessImplTemporary.P1;
  }
  return true;
}
  
template <class ProcessImpl>
template <typename Floating>
auto Process<ProcessImpl>::GetFS(
    Particle<Floating> * particles, MatID_t matID, uint64_t N, PDG_t * outPDG1,
    LorentzVector<Floating> * outP1, PDG_t * outPDG2,
    LorentzVector<Floating> * outP2, Floating * csBorderData) const
    -> bool {
//
#ifdef OPENACC
#pragma acc parallel loop
#else 
#pragma omp parallel for simd
#endif  
//
  for (auto ind = 0u; ind < N; ++ind) {
    auto const fProcessImplTemporary=fProcessImpl.GetFS(particles[ind].p, particles[ind].pdg, matID, particles[ind].rs, csBorderData, ind);
    outPDG1[ind]=fProcessImplTemporary.pdg1;
    outPDG2[ind]=fProcessImplTemporary.pdg2;
    outP1[ind]=fProcessImplTemporary.P1;
    outP2[ind]=fProcessImplTemporary.P2;
  }
  return true;
}
} // namespace t3

#endif // T3PROCESS_H
