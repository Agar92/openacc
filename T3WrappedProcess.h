#pragma once
#ifndef T3WRAPPEDPROCESS_H
#define T3WRAPPEDPROCESS_H

//#include "T3DoNothingProcessImpl.h"
//#include "T3DuplicateProcessImpl.h"
//#include "T3IsotropicAbsoluteElasticFSImpl.h"
//#include "T3IsotropicElasticFSImpl.h"
#include "T3Process.h"
//#include "T3ProcessWithTabulatedCSImpl.h"
//added
//#include "T3RutherfordCSImpl.h"
//#include "T3RutherfordFSImpl.h"
//#include "T3InelasticddCSImpl.h"
//#include "T3InelasticddFSImpl.h"
#include "T3MultipleScatteringCSImpl.h"
#include "T3MultipleScatteringFSImpl.h"
//end of added

namespace t3 {
//CSImpl - cross section with GetCS() function implementation (for example, T3RutherfordCSImpl.h)
//FSImpl - final state with GetFS() function implementation (for example, T3IsotropicAbsoluteElasticFSImpl.h)
//ProcessImplFromCSandFS<CSImpl, FSImpl> forms a class with GetCS() and GetFS() for a single particle.
//Process<ProcessImplFromCSandFS<CSImpl, FSImpl>> forms a class with GetCS() and GetFS() for a vector of particles.
//=>ProcessFromCSandFS is a class with GetCS() and GetFS() for a vector of particles.
template <typename CSImpl, typename FSImpl>
using ProcessFromCSandFS = Process<ProcessImplFromCSandFS<CSImpl, FSImpl>>;

//using DoNothingProcess = Process<DoNothingProcessImpl<double>>;
//using DuplicateProcess = Process<DuplicateProcessImpl<double>>;
//using TabulatedProcess = Process<ProcessWithTabulatedCS<double>>;
using MultipleScatteringProcess = Process<MultipleScatteringCS<double>>;

//using ConstantCSIsotropicAbsoluteElasticProcess =
//  ProcessFromCSandFS<ConstantCS<double>, IsotropicAbsoluteElasticFS<double>>;
////////added
//  using CSRutherfordProcess =
//      ProcessFromCSandFS<RutherfordCS<double>, RutherfordFS<true, double>>;
  
    using CSMultipleScatteringProcess =
        ProcessFromCSandFS<MultipleScatteringCS<double>, MultipleScatteringFS<true, double>>;
  
//using CSInelasticddProcess =
//    ProcessFromCSandFS<InelasticddCS<double>, InelasticddFS<double>>;
//end of added
//using ConstantCSIsotropicElasticProcessWithoutRecoil =
//    ProcessFromCSandFS<ConstantCS<double>, IsotropicElasticFS<false, double>>;

//using ConstantCSIsotropicElasticProcessWithRecoil =
//    ProcessFromCSandFS<ConstantCS<double>, IsotropicElasticFS<true, double>>;

} // namespace t3

#endif // T3WRAPPEDPROCESS_H
