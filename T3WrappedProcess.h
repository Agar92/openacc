#pragma once
#ifndef T3WRAPPEDPROCESS_H
#define T3WRAPPEDPROCESS_H

#include "T3Process.h"
#include "T3MultipleScatteringCSImpl.h"
#include "T3MultipleScatteringFSImpl.h"

namespace t3 {
template <typename CSImpl, typename FSImpl>
using ProcessFromCSandFS = Process<ProcessImplFromCSandFS<CSImpl, FSImpl>>;

using MultipleScatteringProcess = Process<MultipleScatteringCS<double>>;

using CSMultipleScatteringProcess =
        ProcessFromCSandFS<MultipleScatteringCS<double>, MultipleScatteringFS<true, double>>;
} // namespace t3
#endif // T3WRAPPEDPROCESS_H
