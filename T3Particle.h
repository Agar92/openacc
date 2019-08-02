#pragma once
#ifndef T3PARTICLE_H
#define T3PARTICLE_H

#include <cmath>
#include "T3Defs.h"
#include "T3RNG.h"
#include "T3LorentzVector.h"
#include "T3ParticleTable.h"

namespace t3 {

using RandGen = RNDGenerator;

template <typename Floating> struct Particle {
  Particle()
      : Particle<Floating>(t3::LorentzVector<Floating>(0.0,0.0,0.0,0.0),
                           t3::LorentzVector<Floating>(1.0/sqrt(3.0),1.0/sqrt(2.0),1.0/sqrt(6.0),units::G),
                           0.0, t3::PDG_t(2112), 1.0, 0, 1, 1) {}

  Particle(t3::LorentzVector<Floating> _r, t3::LorentzVector<Floating> _p, Floating _de, t3::PDG_t _pdg,
           Floating _wt, RandGen::result_type _rs, int _ir, int _id);
  t3::LorentzVector<Floating> r;// r = {rx, ry, rz, t}
  t3::LorentzVector<Floating> p;// p = {px, py, pz, en}
  Floating de;
  Floating wt;
  RandGen rs;
  t3::PDG_t pdg;
  int ir;
  void SetEtot(Floating Etot, ParticleTable * aParticleTable)
  {
    const auto m = aParticleTable->GetMass(pdg);
    const auto pnew = sqrt(Etot * Etot - m*m);
    p.SetPxPyPzE(pnew*vx(), pnew*vy(), pnew*vz(), Etot);
  }
  Floating GetdE() const { return de; }
  Floating GetEtot() const { return p.E(); }
  Floating vx() const { return p.x() / p.R(); }
  Floating vy() const { return p.y() / p.R(); }
  Floating vz() const { return p.z() / p.R(); }
  int ix() const { return static_cast<int>(std::floor(r.x() / ag)); }
  int jy() const { return static_cast<int>(std::floor(r.y() / ag)); }
  int kz() const { return static_cast<int>(std::floor(r.z() / ag)); }
  int id;
  auto GenerateCanonical() {
     return t3::GenerateSubCanonical<Floating>(rs);
  }
};

template <typename Floating>
Particle<Floating>::Particle(t3::LorentzVector<Floating> _r, t3::LorentzVector<Floating> _p, Floating _de,
                             t3::PDG_t _pdg, Floating _wt, RandGen::result_type _rs, int _ir, int _id)
    : r{_r}, p{_p}, de(_de), pdg(_pdg), wt(_wt), rs(_rs), ir(_ir), id(_id){}

} // namespace t3

#endif // T3PARTICLE_H
