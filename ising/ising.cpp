//
// Created by mitja on 13.4.2019.
//

#include <chrono>
#include "Ising2D.h"

static const size_t steps = N * 10;

static const double J = 1.;

static const double iBeta = 0;
static const double dBeta = 0.1;
static const uint32_t nBeta = 1000;

int main()
{
    auto start = std::chrono::high_resolution_clock::now();

    Ising2D lattice(randomState());

    double m = 0.;
    for (size_t n = 0; n < 1e7; ++n)
    {
        m += lattice.magnetization();
    }
    std::cout << m << std::endl;

    auto finish = std::chrono::high_resolution_clock::now();
    std::cout   << "Finished in "
                << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count()
                << "ms" << std::endl;
    return 0;
}
