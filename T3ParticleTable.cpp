#include "T3ParticleTable.h"

namespace t3 {
ParticleTable::ParticleTable() {
  fMasses[PDG_t(22)]   = 0.0;         //photon
  fMasses[PDG_t(11)]   = electronMass;//electron
  fMasses[PDG_t(-11)]  = electronMass;//positron
  fMasses[PDG_t(2212)] = protonMass;  //proton
  fMasses[PDG_t(2112)] = neutronMass; //neutron
}
/*
ParticleTable::ParticleTable(const ParticleTable & a)
{
  //this->fMasses[PDG_t(22)]   = a.fMasses[PDG_t(22)];
  //this->fMasses[PDG_t(11)]   = a.fMasses[PDG_t(11)];
  //this->fMasses[PDG_t(2212)] = a.fMasses[PDG_t(2212)];
  //this->fMasses[PDG_t(2112)] = a.fMasses[PDG_t(2112)];
  this->fMasses[PDG_t(22)]   = 0.0;         //photon
  this->fMasses[PDG_t(11)]   = electronMass;//electron
  this->fMasses[PDG_t(-11)]  = electronMass;//positron
  this->fMasses[PDG_t(2212)] = protonMass;  //proton
  this->fMasses[PDG_t(2112)] = neutronMass; //neutron
  this->photonMass = a.photonMass;
  this->electronMass = a.electronMass;
  this->protonMass = a.protonMass;
  this->neutronMass = a.neutronMass;
  this->aeMass = a.aeMass;
  this->nucleusPDGStartAfter = a.nucleusPDGStartAfter;
}
*/
} // namespace t3
