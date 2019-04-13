//
// Created by mitja on 13.4.2019.
//

#ifndef VAJA_IV_1_ISING2D_H
#define VAJA_IV_1_ISING2D_H

#include <bitset>
#include <random>
#include <functional>

static const uint32_t L = 1 << 7;
static const uint32_t N = L * L;

class Ising2D
{
public:
    Ising2D(double j, double h);
    Ising2D(double j, double h, const std::string& init);
    inline double magnetization(){return (2*(double)m_spins.count() - N) / N;}
    inline double energy(){return m_energy;}
    inline const std::bitset<N>& lattice() const {return m_spins;}
    void step(double beta);
private:
    std::bitset<N> m_spins;

    std::uniform_int_distribution<uint32_t> m_spot;
    std::uniform_real_distribution<double> m_chi;
    std::mt19937 m_twister;

    double m_j;
    double m_h;
    double m_energy;
};

inline std::string randomState()
{
    std::mt19937 twister(std::random_device{}());
    std::uniform_int_distribution<char> bits('0', '1');
    auto bit = std::bind(bits, twister);

    std::string state;
    for (uint32_t n = 0; n < N; ++n) state += bit();

    return state;
}

#endif //VAJA_IV_1_ISING2D_H
