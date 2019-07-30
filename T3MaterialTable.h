#pragma once
#ifndef T3MATERIALTABLE_H
#define T3MATERIALTABLE_H

#include <atomic>
#include <utility>
#include <vector>

#include "T3Defs.h"

namespace t3 {

class MaterialTable {

public:
  MaterialTable();
  //make here return by reference.
  auto const & GetIsotopes(MatID_t matID) const  { return fMaterials.at(matID); }
  auto GetFractions(MatID_t matID) const { return fFractions.at(matID); }
  auto GetConcentrations(MatID_t matID) const { return fConcentrations.at(matID); }
  auto GetDensity(MatID_t matID) const { return fDensities.at(matID); }
  auto GetNumberOfIsotopes(MatID_t matID) const { return fMaterials.at(matID).size(); }
  auto GetNMat() const { return fMaterials.size(); }

private:
  static std::atomic<bool> wasInitialized;
  static std::vector<std::vector<PDG_t>> fMaterials;
  static std::vector<std::vector<double>> fFractions;
  static std::vector<double> fConcentrations;
  static std::vector<double> fDensities;
};
} // namespace t3

#endif // T3MATERIALTABLE_H
