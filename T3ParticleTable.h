#pragma once
#ifndef T3PARTICLETABLE_H
#define T3PARTICLETABLE_H

#include <atomic>
#include <map>
#include <vector>
#include <iostream>

#include "T3Defs.h"

namespace t3 {
using namespace units;
class ParticleTable {
public:
  ParticleTable();
  auto IsNucleus(PDG_t aPDG) const;
  auto GetMass(PDG_t aPDG) const;
  auto GetZ(PDG_t aPDG) const;
  auto GetA(PDG_t aPDG) const;
  std::pair<uint64_t, uint64_t> GetZA(PDG_t aPDG) const;
  static constexpr auto makePDGfromZandA(uint64_t, uint64_t);
private:
  static std::atomic<bool> wasInitialized;
  static std::map<PDG_t, double> fMasses;
  static auto constexpr photonMass = 0.0;
  static auto constexpr electronMass = 0.511 * MeV;
  static auto constexpr protonMass = 938.272 * MeV;
  static auto constexpr neutronMass = 939.566 * MeV;
  static auto constexpr aeMass = 931.494 * MeV;
  static auto constexpr nucleusPDGStartAfter = t3::PDG_t(1000000000);
};

inline std::pair<uint64_t, uint64_t> ParticleTable::GetZA(PDG_t aPDG) const {
  if (aPDG < nucleusPDGStartAfter) {
    return std::make_pair(0u, 0u);
  }
  auto const baryonNumber = (aPDG / 10) % 1000;
  auto const protonNumber = (aPDG / 10000) % 1000;
  return std::make_pair(std::move(protonNumber), std::move(baryonNumber));
}

inline auto ParticleTable::IsNucleus(PDG_t aPDG) const { return aPDG > nucleusPDGStartAfter; }

inline auto ParticleTable::GetZ(PDG_t aPDG) const { return GetZA(aPDG).first; }

inline auto ParticleTable::GetA(PDG_t aPDG) const { return GetZA(aPDG).second; }

inline constexpr auto ParticleTable::makePDGfromZandA(uint64_t protonNumber,
                                                      uint64_t baryonNumber) {
  return static_cast<PDG_t>(static_cast<int64_t>(
      1000000000u + 10000u * protonNumber + baryonNumber * 10u));
}

inline auto ParticleTable::GetMass(PDG_t aPDG) const
{
  if(aPDG!=makePDGfromZandA(1, 2) || aPDG!=makePDGfromZandA(22, 48)){}
  if(aPDG > nucleusPDGStartAfter)
  {
    if(aPDG == makePDGfromZandA(1,2)) return 1875.6 * MeV;
    else if(aPDG == makePDGfromZandA(22,48))
    {
       int Z=GetZ(aPDG);
       double A=GetA(aPDG);
       return A * aeMass - Z * electronMass - 48.4917 * MeV;
    }
    else
    {
      auto const pair = GetZA(aPDG);
      auto const protonNumber=pair.first;
      auto const baryonNumber=pair.second;
      auto const neutronNumber = baryonNumber - protonNumber;
      return protonNumber * protonMass + neutronNumber * neutronMass;
    }
  }
  return fMasses.at(aPDG);
}

} // namespace t3

#endif // T3PARTICLETABLE_H
