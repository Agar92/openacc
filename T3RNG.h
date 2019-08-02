#pragma once
#ifndef T3RNG_H
#define T3RNG_H

#include <cmath>
#include <cstdint>
//------------------------//
//Random number generator.//
//------------------------//

namespace t3 {

class RNDGenerator
{
private:
  unsigned int fseed;
public:
  void seed(unsigned int seed){fseed=seed;}
  using result_type=unsigned int;
  RNDGenerator():fseed(1){}
  RNDGenerator(unsigned int seed):fseed(seed){}
  result_type min() const {return 0u;}
  result_type max() const {return 0xFFFFFFFF;}
  result_type Rand32(unsigned int xn)
  {
    u_quad_t a=0x5DEECE66D;
    u_quad_t c=0xB;
    return (unsigned int)((a*xn+c) & 0xFFFFFFFF);
  }
  result_type operator()()
  {
    fseed=Rand32(fseed);
    return fseed;
  }
};

} // namespace t3

#endif // T3RNG_H
