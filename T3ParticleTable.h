#pragma once
#ifndef T3PARTICLETABLE_H
#define T3PARTICLETABLE_H

#include <map>
#include <iostream>
#include "T3Globals.h"
#include "T3Defs.h"

namespace t3 {
using namespace units;
class ParticleTable {
public:
  ParticleTable();
  /*ParticleTable(const ParticleTable & aParticleTable);*/
  auto IsNucleus(PDG_t aPDG) const;
  auto GetMass(PDG_t aPDG) const;
  auto GetZ(PDG_t aPDG) const;
  auto GetA(PDG_t aPDG) const;
  std::pair<uint64_t, uint64_t> GetZA(PDG_t aPDG) const;
  PDG_t makePDGfromZandA(uint64_t, uint64_t) const;
private:
  std::map<PDG_t, double> fMasses;
  double photonMass = 0.0;
  double electronMass = 0.511 * MeV;
  double protonMass = 938.272 * MeV;
  double neutronMass = 939.566 * MeV;
  double aeMass = 931.494 * MeV;
  double nucleusPDGStartAfter = t3::PDG_t(1000000000);
};

inline std::pair<uint64_t, uint64_t> ParticleTable::GetZA(PDG_t aPDG) const {
  if (aPDG < nucleusPDGStartAfter) {
    return std::make_pair(0u, 0u);
  }
  auto const baryonNumber = (aPDG / 10) % 1000;
  auto const protonNumber = (aPDG / 10000) % 1000;
  return std::make_pair(protonNumber, baryonNumber);
}

inline auto ParticleTable::IsNucleus(PDG_t aPDG) const { return aPDG > nucleusPDGStartAfter; }

inline auto ParticleTable::GetZ(PDG_t aPDG) const { return GetZA(aPDG).first; }

inline auto ParticleTable::GetA(PDG_t aPDG) const { return GetZA(aPDG).second; }

inline PDG_t ParticleTable::makePDGfromZandA(uint64_t protonNumber,
                                                      uint64_t baryonNumber) const {
  return static_cast<PDG_t>(static_cast<int64_t>(
      1000000000u + 10000u * protonNumber + baryonNumber * 10u));
}

inline auto ParticleTable::GetMass(PDG_t aPDG) const
{
  //\\//std::cout<<"aPDG="<<aPDG<<std::endl;
//#1
//THIS GIVES AN ERROR:
  /*
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
  */
//#2
//THIS WORKS CORRECTLY:
  /*
  auto m=0.0;
  switch(aPDG)
  {
  case PDG_t(22):
    m=0.0;
    break;
  case PDG_t(11):
    m=electronMass;
    break;
  case PDG_t(-11):
    m=electronMass;
    break;
  case PDG_t(2212):
    m=protonMass;
    break;
  case PDG_t(2112):
    m=neutronMass;
    break;
  default:
    m=0.0;
    break;
  }
  return m;
  */
//#3
//THIS WORKS CORRECTLY:  
  //*
  auto m=0.0;
  if(aPDG==PDG_t(22)) m=0.0;
  else if(aPDG==PDG_t(11)) m=electronMass;
  else if(aPDG==PDG_t(-11)) m=electronMass;
  else if(aPDG==PDG_t(2212)) m=protonMass;
  else if(aPDG==PDG_t(2112)) m=neutronMass;
  return m;
  //*/
}
  
} // namespace t3
#endif // T3PARTICLETABLE_H
