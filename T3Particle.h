#pragma once
#ifndef T3PARTICLE_H
#define T3PARTICLE_H

#include <cmath>
#include "T3Globals.h"
#include "T3RNG.h"
#include "T3ParticleTable.h"
#include "T3LorentzVector.h"

namespace t3 {

template <typename Floating> struct Particle {
Particle()
  : Particle<Floating>(LorentzVector<Floating>(0.0,0.0,0.0,0.0),
                       LorentzVector<Floating>(1.0/sqrt(3.0),1.0/sqrt(2.0),1.0/sqrt(6.0),G),
                       0.0, t3::PDG_t(2112), 1.0, 0, 1, 1) {}

  Particle(LorentzVector<Floating> _r, LorentzVector<Floating> _p, Floating _de, t3::PDG_t _pdg,
         Floating _wt, RandGen::result_type _rs, int _ir, int _id);
///\\\///float4 r;// r = {rx, ry, rz, t}
LorentzVector<Floating> r;// r = {rx, ry, rz, t}
///\\\///float4 p;// p = {px, py, pz, en}
LorentzVector<Floating> p;// p = {px, py, pz, en}
Floating de;
Floating wt;
RandGen rs;
PDG_t pdg;
int ir;
//var1, var2 ARE BECAUSE WITHOUT THEM HAVE COMPILER ERROR:
//PGCC-S-1000-Call in OpenACC region to procedure 'memcpy' which has no acc routine information
//the structure must be 512 bit or so!!!
int var1;
int var2;
//int var3;
//int var4;
/*
void SetEtot(Floating Elskin, ParticleTable & aParticleTable, PDG_t * extrapdg, Floating * extram, int i)
{
  extrapdg[i]=pdg;
  const auto m =aParticleTable.GetMass(pdg);
  extram[i]=m;
  const auto plsnew = std::sqrt(Elskin*(2*m+Elskin));
  SetPxPyPzE(plsnew*vx(), plsnew*vy(), plsnew*vz(), Elskin+m);
}
*/
void SetEtot(Floating Etot, ParticleTable & aParticleTable, PDG_t * extrapdg, Floating * extram, int i)
{
  extrapdg[i]=pdg;
  const Floating m =aParticleTable.GetMass(pdg);
  extram[i]=m;
  if(Etot<m) return;
  const auto plsnew = std::sqrt(Etot*Etot-m*m);
  SetPxPyPzE(plsnew*vx(), plsnew*vy(), plsnew*vz(), Etot);
}
void SetPxPyPzE(Floating px, Floating py, Floating pz, Floating E){p.SetPxPyPzE(px, py, pz, E);}
Floating m() const {return p.m();}
Floating GetdE() const { return de; }
Floating GetEtot() const { return p.E(); }
Floating R() const { return std::sqrt(r.x()*r.x()+r.y()*r.y()+r.z()*r.z());}
Floating P() const { return std::sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());}
Floating vx() const { return p.x() / P(); }
Floating vy() const { return p.y() / P(); }
Floating vz() const { return p.z() / P(); }
int ix() const { return static_cast<int>(std::floor(r.x() / ag)); }
int jy() const { return static_cast<int>(std::floor(r.y() / ag)); }
int kz() const { return static_cast<int>(std::floor(r.z() / ag)); }
int id;
auto GenerateCanonical() {
   return t3::GenerateSubCanonical<Floating>(rs);
}
};

template <typename Floating>
Particle<Floating>::Particle(LorentzVector<Floating> _r, LorentzVector<Floating> _p, Floating _de,
                             t3::PDG_t _pdg, Floating _wt, RandGen::result_type _rs, int _ir, int _id)
    : r{_r}, p{_p}, de(_de), pdg(_pdg), wt(_wt), rs(_rs), ir(_ir), id(_id){}


} // namespace t3

#endif // T3PARTICLE_H
