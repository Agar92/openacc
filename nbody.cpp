#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <cmath>
#include <omp.h>
#include <unistd.h>
#include <vector>
#include <cfloat>
#include <algorithm>

#include "T3DataHolder.h"

using namespace dataholder;

//////////////////////////////////////////////////////////////////////////////
// Program main
//////////////////////////////////////////////////////////////////////////////
int main(int, char **) {  
  auto begin = std::chrono::steady_clock::now();

  DataHolder<FloatingType> * data;
#ifdef OPENACC  
  cudaMalloc((void**) & data, sizeof(DataHolder<FloatingType>));
#else  
  data = new DataHolder<FloatingType>;
#endif  
  if (report)
    std::cout << "size of DataHolder is ~"
              << sizeof(DataHolder<FloatingType>) / float(1ul << 30ul)
              << "GB, size of d is ~" << sizeof(data) / float(1ul << 30ul) << "GB"
              << std::endl;

  data->InitParticle();
  for (unsigned int step = 1; GetNumOfAliveParticles() > 0; ++step) {
    data->Propagate();
    data->Compress();
    data->React();
    data->Inject();
    if (report) {
      std::cout << step << "   " << GetNumOfAliveParticles() << "    "
                << GetNoNew() <<  "   " << GetNumOfInjectedParticles() <<std::endl;
    }
  }
  data->Histogram_theta();

  auto end = std::chrono::steady_clock::now();
  auto elapsed_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

#ifdef OPENACC
  cudaFree(data);
#endif
  
  delete data;
}


