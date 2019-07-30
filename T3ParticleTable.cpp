#include "T3ParticleTable.h"
#include <mutex>

namespace t3 {

decltype(ParticleTable::fMasses) ParticleTable::fMasses{};
decltype(ParticleTable::wasInitialized) ParticleTable::wasInitialized{false};

ParticleTable::ParticleTable() {
  if (wasInitialized) {
    return;
  }
  static std::mutex initializationMutex;
  std::lock_guard<decltype(initializationMutex)> lock(initializationMutex);
  // FIXME dummy
  //Our units: energy is in Mev, cross section - in mm^2.
  //PDG codes:gamma=22;electron=11;positron=-11;proton=2212;neutron=2112;deuteron=1000010020;
  //tritium=1000010030; He3=1000020030;He42 (alpha-particle)=1000020040.
  fMasses[PDG_t(22)]   = 0.0;
  fMasses[PDG_t(11)]   = electronMass;
  fMasses[PDG_t(2212)] = protonMass;
  fMasses[PDG_t(2112)] = neutronMass;
  wasInitialized = true;
}

} // namespace t3
