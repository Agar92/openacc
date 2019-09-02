#pragma once
#ifndef T3MULTIPLESCATTERINGCSIMPL_H
#define T3MULTIPLESCATTERINGCSIMPL_H

#include <cmath>
#include "T3Defs.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"

namespace t3 {
using namespace units;
template <typename Floating> class MultipleScatteringCS/*RutherfordCS*/ {
public:
  MultipleScatteringCS() = default;
  inline Floating
  GetCS(Floating e/*LS kinetic energy of the inc particle in MeV*/, PDG_t incPDG,
        PDG_t targetPDG, const ParticleTable & aParticleTable) const; // PDG_t-PDG code of incident particle; targetPDG - PDG code of the target isotope.
  inline Floating
  GetCS(Floating e, PDG_t,
        MatID_t matID, Floating * csIsotope,
        const ParticleTable & aParticleTable, const MaterialTable & aMaterialTable,
        int particle_index) const; // PDG_t - PDG code of incident particle; MatID_t - target material number.
  Floating const Edisplace_deuteron = 10 * eV;//displacement energy for deuteron is 10 eV.
  Floating const Edisplace_titan    = 25 * eV;//displacement energy for Ti48 is 25 eV.
};
//returns the microscopic  cross section in mm^2
template <typename Floating>
Floating MultipleScatteringCS<Floating>::GetCS(Floating e/*LS kinetic energy of the inc particle in MeV*/,
                                               PDG_t incPDG, PDG_t targetPDG/*PDG of target isotope of the material*/,
                                               const ParticleTable & aParticleTable) const {
  //e - LS kinetic energy of the incident particle
  //if the particle is a gamma or a neutron, the Rutherford cross section is 0.0.
  //if (incPDG == PDG_t(22) || incPDG == PDG_t(2112)) return 0.;
  const Floating m=aParticleTable.GetMass(incPDG);
  const Floating Tls=e-m;
  const auto deuteronPDG = aParticleTable.makePDGfromZandA(1, 2);
  const auto titanPDG    = aParticleTable.makePDGfromZandA(22, 48);

  auto const c=1.0;
  auto const alpha=1.0/137;
  auto const ec=sqrt(alpha);
  auto const me=aParticleTable.GetMass(PDG_t(11));
  auto const h=1.0;
  auto const a0=h*h/me/ec/ec;
  auto const CTF=1.0/2.0*pow(3*M_PI/4,2.0/3.0);
  auto const M=aParticleTable.GetMass(targetPDG);
  auto const z=aParticleTable.GetZ(incPDG);
  auto const Z=aParticleTable.GetZ(targetPDG);
  auto const Elsinc=m+e;
  auto const s=m*m+2*Elsinc*M+M*M;
  auto const mr=m*M/sqrt(s);
  auto const pls2=e*(2*m+e);
  auto const rat=m/M;
  auto const rat2=rat*rat;
  auto const pcm2=pls2/(1.0+rat2+2*sqrt(rat2+pls2/M/M));
  auto const pcm=sqrt(pcm2);
  auto const aI=CTF*a0/sqrt(pow(z, 2.0/3.0) + pow(Z, 2.0/3.0));

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
  if(incPDG == targetPDG) cs *= 1.0/*2.0*/;

  //returns isotope-isotope cross section in mm^2.
// \\ //  
  cs=1.0;
// \\ //  
  return cs;
}


template <typename Floating>
Floating MultipleScatteringCS<Floating>::GetCS(Floating e,
                                               PDG_t incPDG, MatID_t matID,
                                               Floating * csBorderDataCS,
                                               const ParticleTable & aParticleTable,
                                               const MaterialTable & aMaterialTable,
                                               int particle_index) const {
  //\\//if (incPDG == PDG_t(22) || incPDG == PDG_t(2112)) return 0.;
  const int numi = aMaterialTable.GetNumberOfIsotopes(matID);
  const int MaxNumberOfIsotopes=aMaterialTable.GetMaxNumberOfIsotopes();
  
  //this does not work:   Floating csIsotope[numi];
  //C array or std::vector does not work, they give an error:
  //PGCC-S-0155-Accelerator region ignored; ...
  for(size_t i=0; i<numi; ++i)
  {
    csBorderDataCS[particle_index*MaxNumberOfIsotopes+i*numi] = GetCS(e,incPDG,aMaterialTable.GetIsotopes(matID,i), aParticleTable);
  }
  auto cs=0.0;
  for(size_t i=0; i<numi; ++i)
  {
    cs += aMaterialTable.GetFractions(matID,i)*csBorderDataCS[i];
  }
// \\ // 
  cs=1.0;
// \\ //   
  return cs;
}

} // namespace t3
#endif // T3MULTIPLESCATTERINGCSIMPL_H
