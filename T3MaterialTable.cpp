#include "T3Defs.h"
#include "T3MaterialTable.h"
#include <mutex>
#include "T3ParticleTable.h"
#include <iostream>

namespace t3 {

decltype(MaterialTable::fMaterials) MaterialTable::fMaterials{};
decltype(MaterialTable::fFractions) MaterialTable::fFractions{};
decltype(MaterialTable::fConcentrations) MaterialTable::fConcentrations{};
decltype(MaterialTable::fDensities) MaterialTable::fDensities{};
decltype(MaterialTable::wasInitialized) MaterialTable::wasInitialized{false};

MaterialTable::MaterialTable() {
  if (wasInitialized)
    return;
  else {
    static std::mutex initializationMutex;
    std::lock_guard<decltype(initializationMutex)> guard(initializationMutex);
    if (wasInitialized)
      return;
    fMaterials = {
        std::vector<PDG_t>{PDG_t{1000010020}},                    // 0 deuterium
        std::vector<PDG_t>{PDG_t{1000220480}},                    // 1 22Ti48
        std::vector<PDG_t>{PDG_t{1000010020}, PDG_t{1000220480}}, // 2 TiD_2
        std::vector<PDG_t>{PDG_t{1000040090}},                    // 3 4Be9
        std::vector<PDG_t>{PDG_t{1000922350},
                          PDG_t{1000922380}}}; // 4//mix of U235 and U238.
    fFractions = {std::vector<double>{1.}, std::vector<double>{1.},
                  std::vector<double>{2. / 3., 1 - 2. / 3.},
                  std::vector<double>{1.},
                  std::vector<double>{0.0072, 1 - 0.0072}};
    fDensities = {0.000169 * gcm3, 4.11 * gcm3, 3.91 * gcm3, 1.85 * gcm3, 19.1 * gcm3};
    fConcentrations.resize(GetNMat());
    ParticleTable aParticleTable = ParticleTable();
    for(size_t i=0; i<fMaterials.size(); ++i)
    {
      auto Atot=0.0;
      for(size_t j=0; j<fMaterials.at(i).size(); ++j)
      {
        Atot += fFractions.at(i).at(j) * aParticleTable.GetMass(fMaterials.at(i).at(j));
      }
      fConcentrations.at(i) = fDensities.at(i) / Atot;
    }
    wasInitialized=true;
  }
}
} // namespace t3
