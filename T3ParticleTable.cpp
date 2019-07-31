#include "T3ParticleTable.h"
#include <mutex>

namespace t3 {
ParticleTable::ParticleTable() {
  fMasses[PDG_t(22)]   = 0.0;
  fMasses[PDG_t(11)]   = electronMass;
  fMasses[PDG_t(2212)] = protonMass;
  fMasses[PDG_t(2112)] = neutronMass;
}
} // namespace t3
