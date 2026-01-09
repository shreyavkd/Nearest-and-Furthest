#ifndef RNG_H
#define RNG_H

#include <random>

namespace utils
{
namespace rand
{
    void seedRand(unsigned long long seed);
    double randDouble(double min, double max);
}
}

#endif