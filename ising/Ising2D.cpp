//
// Created by mitja on 13.4.2019.
//

#include "Ising2D.h"

Ising2D::Ising2D():
m_spins(),
m_spot(0, N - 1),
m_chi(0., 1.),
m_twister(std::random_device{}())
{}

Ising2D::Ising2D(const std::string& init):
m_spins(init),
m_spot(0, N - 1),
m_chi(0., 1.),
m_twister(std::random_device{}())
{}

void Ising2D::step(double J, double h, double beta)
{
    // random spot
    const uint32_t r = m_spot(m_twister);

    // energy difference
    double change = h;

    // right neighbour
    uint32_t rn = r + 1;
    if (rn % L == 0) rn -= L;
    if (m_spins.test(rn)) change += J;
    else change -= J;

    // left neighbour
    if (r % L == 0) rn = r + L - 1;
    else rn = r - 1;
    if (m_spins.test(rn)) change += J;
    else change -= J;

    // lower neighbour
    rn = r + L;
    if (rn >= N) rn %= L;
    if (m_spins.test(rn)) change += J;
    else change -= J;

    // upper neighbour
    if (r / L == 0) rn = N - L + r;
    else rn = r - L;
    if (m_spins.test(rn)) change += J;
    else change -= J;

    // on spot
    if (m_spins.test(r)) change *= 2;
    else change *= -2;

    // to flip or not to flip
    if (change < 0) m_spins.flip(r);
    else if (m_chi(m_twister) < std::exp(-beta * change)) m_spins.flip(r);
}
