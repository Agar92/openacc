#pragma once
#ifndef T3MULTIPLESCATTERINGCSIMPL_H
#define T3MULTIPLESCATTERINGCSIMPL_H

#include "T3ParticleTable.h"
#include "T3MaterialTable.h"
#include "T3Util.h"
#include "T3Defs.h"
#include <cmath>

namespace t3 {
using namespace units;
template <typename Floating> class MultipleScatteringCS {
public:
  MultipleScatteringCS() = default;
  inline Floating
    GetCS(Floating e, PDG_t,
          PDG_t targetPDG) const;
  inline Floating
  GetCS(Floating e, PDG_t,
        MatID_t matID) const;
  ParticleTable aParticleTable;
  MaterialTable aMaterialTable;
};
template <typename Floating>
Floating MultipleScatteringCS<Floating>::GetCS(Floating e,
                                               PDG_t incPDG, PDG_t targetPDG) const {
  auto const Edisplace_deuteron = 10 * eV;
  auto const Edisplace_titan    = 25 * eV;
  if (incPDG == PDG_t(22) || incPDG == PDG_t(2112)) return 0.;
  auto const deuteronPDG = aParticleTable.makePDGfromZandA(1, 2);
  auto const titanPDG    = aParticleTable.makePDGfromZandA(22, 48);
  auto const c=1.0;
  auto const alpha=1.0/137;
  auto const ec=std::sqrt(alpha);
  auto const me=aParticleTable.GetMass(PDG_t(11));
  auto const h=1.0;
  auto const a0=h*h/me/ec/ec;
  auto const CTF=1.0/2.0*std::pow(3*M_PI/4,2.0/3.0);
  auto const m=aParticleTable.GetMass(incPDG);
  auto const M=aParticleTable.GetMass(targetPDG);
  auto const z=aParticleTable.GetZ(incPDG);
  auto const Z=aParticleTable.GetZ(targetPDG);
  auto const Elsinc=m+e;
  auto const s=m*m+2*Elsinc*M+M*M;
  auto const mr=m*M/std::sqrt(s);
  auto const pls2=e*(2*m+e);
  auto const rat=m/M;
  auto const rat2=rat*rat;
  auto const pcm2=pls2/(1.0+rat2+2*std::sqrt(rat2+pls2/M/M));
  auto const pcm=std::sqrt(pcm2);
  auto const aI=CTF*a0/std::sqrt(std::pow(z, 2.0/3.0) + std::pow(Z, 2.0/3.0));
  auto const beta_r=1.0;
  auto const c12=alpha*z*Z/beta_r;
  auto const c122=c12*c12;
  auto const const1=h/2/aI/pcm;
  auto const const2=const1*const1;
  auto const AS=const2*(1.13+3.76*c122);
  auto tmin=0.0;
  if(targetPDG == deuteronPDG)   tmin=2*Edisplace_deuteron*M;
  else if(targetPDG == titanPDG) tmin=2*Edisplace_titan*M;
  auto const tmax_real = 4 * pcm2;
  auto tmax = 0.0;
  if(tmax_real<=tmin) tmax = tmax_real;
  else                tmax = tmin;
  auto const ca0=mr*alpha*z*Z;
  auto const pcm4=pcm2*pcm2;
  auto const ca=ca0*ca0*M_PI/4/pcm4/pcm2;
  auto const ASC1=tmax/4/pcm2+AS;
  auto const ASC=tmax/AS/ASC1;
  auto const hc=200.0 * MeV * fm;
  auto cs=ca*ASC*hc*hc;
  if(incPDG == targetPDG) cs *= 1.0;
  return cs;
}

template <typename Floating>
Floating MultipleScatteringCS<Floating>::GetCS(Floating e,
                                               PDG_t incPDG, MatID_t matID) const {
  if (incPDG == PDG_t(22) || incPDG == PDG_t(2112)) return 0.;
  size_t const numi = aMaterialTable.GetNumberOfIsotopes(matID);
  Floating csIsotope[numi];
  for(size_t i=0; i<numi; ++i)
  {
    csIsotope[i] = GetCS(e,incPDG,aMaterialTable.GetIsotopes(matID)[i]);
  }
  auto cs=0.0;
  for(size_t i=0; i<numi; ++i)
  {
    cs += aMaterialTable.GetFractions(matID)[i]*csIsotope[i];
  }
  return cs;
}
  
} // namespace t3

#endif // T3MULTIPLESCATTERINGCSIMPL_H
