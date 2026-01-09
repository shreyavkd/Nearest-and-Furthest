#include "rng.h"

namespace utils
{
namespace rand
{
    static std::mt19937_64 engine;

    void seedRand(unsigned long long seed)
    {
        engine.seed(seed);
    }

    double randDouble(double min, double max)
    {
        std::uniform_real_distribution<double> d(min, max);
        return d(engine);
    }

}
}