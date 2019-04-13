//
// Created by mitja on 13.4.2019.
//

#ifndef VAJA_IV_1_ISING2D_H
#define VAJA_IV_1_ISING2D_H

#include <bitset>
#include <random>

static const uint32_t L = 1 << 7;
static const uint32_t N = L * L;

class Ising2D
{
public:
    Ising2D();
    Ising2D(unsigned long init);
    inline double magnetization(){return (double)(2*m_spins.count() - N) / N;}
    void step(double J, double h, double beta);
private:
    std::bitset<N> m_spins;
    std::uniform_int_distribution<uint32_t> m_spot;
    std::uniform_real_distribution<double> m_chi;
    std::mt19937 m_twister;
};


#endif //VAJA_IV_1_ISING2D_H
