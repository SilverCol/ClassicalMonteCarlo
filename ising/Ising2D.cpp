//
// Created by mitja on 13.4.2019.
//

#include "Ising2D.h"

Ising2D::Ising2D(double j, double h) :
m_spins(),
m_spot(0, N - 1),
m_chi(0., 1.),
m_twister(std::random_device{}()),
m_j(j),
m_h(h),
m_energy(-j*2*N - std::abs(h)*N)
{
    if (h >= 0) m_spins = ~m_spins;
    if (j < 0)
    {
        if (L % 2 == 0) std::cerr << "Should be using odd dimension for AFM." << std::endl;
        for (uint32_t i = 0; i < N; i += 2)
        {
            m_spins[i].flip();
            if (i % L == L - 2) ++i;
            else if (i % L == L - 1) --i;
        }
    }
    else if (L % 2 != 0) std::cerr << "Using odd dimension for FM." << std::endl;
}

Ising2D::Ising2D(double j, double h, const std::string& init) :
m_spins(init),
m_spot(0, N - 1),
m_chi(0., 1.),
m_twister(std::random_device{}()),
m_j(j),
m_h(h),
m_energy(-h*(2*m_spins.count() - N))
{
    // add bond contribution to energy
    for (uint32_t r = 0; r < N; ++r)
    {
        // right neighbour
        uint32_t rn = r + 1;
        if (rn % L == 0) rn -= L;
        if (m_spins[r] ^ m_spins[rn]) m_energy += j;
        else m_energy -= j;

        // lower neighbour
        rn = r + L;
        if (rn >= N) rn = r % L;
        if (m_spins[r] ^ m_spins[rn]) m_energy += j;
        else m_energy -= j;
    }
}

void Ising2D::step(double beta)
{
    // random spot
    const uint32_t r = m_spot(m_twister);

    // energy difference
    double change = m_h;

    // right neighbour
    uint32_t rn = r + 1;
    if (rn % L == 0) rn -= L;
    if (m_spins.test(rn)) change += m_j;
    else change -= m_j;

    // left neighbour
    if (r % L == 0) rn = r + L - 1;
    else rn = r - 1;
    if (m_spins.test(rn)) change += m_j;
    else change -= m_j;

    // lower neighbour
    rn = r + L;
    if (rn >= N) rn %= L;
    if (m_spins.test(rn)) change += m_j;
    else change -= m_j;

    // upper neighbour
    if (r / L == 0) rn = N - L + r;
    else rn = r - L;
    if (m_spins.test(rn)) change += m_j;
    else change -= m_j;

    // on spot
    if (m_spins.test(r)) change *= 2;
    else change *= -2;

    // to flip or not to flip
    if (change < 0)
    {
        m_spins.flip(r);
        m_energy += change;
    }
    else if (m_chi(m_twister) < std::exp(-beta * change))
    {
        m_spins.flip(r);
        m_energy += change;
    }
}
