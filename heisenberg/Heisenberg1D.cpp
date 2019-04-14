//
// Created by mitja on 14.4.2019.
//

#include "Heisenberg1D.h"

Heisenberg1D::Heisenberg1D(double j, double h):
m_spins(N, {0.0, 0.0, 1.0}),
m_spot(0, N - 1),
m_chi(0.0, 1.0),
m_phi(0.0, 6.283185307179586),
m_z(-1.0, 1.0),
m_twister(std::random_device{}()),
m_j(j),
m_h(h),
m_energy(-(j + h)*N),
m_magnet({0.0, 0.0, 1.0})
{
    if (h < 0) // downwards magnetic
    {
        for (auto& spin : m_spins) spin.z = -1.0;
        m_magnet.z = -1.0;
    }
    if (j < 0) // AFM
    {
        // init to Neel
        for (uint32_t j = 0; j < N; j += 2) m_spins[j].z *= -1.0;
        m_energy = j*N;
        m_magnet.z = 0.0;
    }
}

Heisenberg1D::Heisenberg1D(double j): // for constant 0 magnetization
m_spins(N, {0.0, 0.0, 1.0}),
m_spot(0, N - 1),
m_chi(0.0, 1.0),
m_phi(0.0, 6.283185307179586),
m_z(-1.0, 1.0),
m_twister(std::random_device{}()),
m_j(j),
m_h(),
m_energy(j*N),
m_magnet({0.0, 0.0, 0.0})
{
    for (uint32_t j = 0; j < N; j += 2) m_spins[j].z = -1.0;
}

Spin Heisenberg1D::randomSpin()
{
    double z = m_z(m_twister);
    double phi = m_phi(m_twister);

    double v = std::sqrt(1.0 - z*z);

    return {v*std::cos(phi), v*std::cos(phi), z};
}

void Heisenberg1D::step(double beta)
{
    // random spot, and new spin -> change in spin
    const uint32_t j = m_spot(m_twister);
    const Spin s = randomSpin();
    const Spin ds = m_spins[j] - s;

    // energy difference
    double change = 0;

    // right neighbour
    uint32_t jn = (j + 1) % N;
    change += m_spins[jn] * ds;

    // left neighbour
    if (j != 0) jn = j - 1;
    else jn = N;
    change += m_spins[jn] * ds;

    change *= m_j;
    change += m_h * ds.z;

    // to flip or not to flip
    if (change < 0)
    {
        m_spins[j] = s;
        m_energy += change;
        m_magnet -= ds / N;
    }
    else if (m_chi(m_twister) < std::exp(-beta * change))
    {
        m_spins[j] = s;
        m_energy += change;
        m_magnet -= ds / N;
    }
}

void Heisenberg1D::stepCarefully(double beta) // for constant 0 magnetization
{
    // random spot, and new spin -> change in spin
    const uint32_t j = m_spot(m_twister);
    const Spin s = randomSpin();
    const Spin ds = m_spins[j] - s;

    // new neighbour to compensate
    const uint32_t jp = (j + 1) % N;
    const Spin sp = m_spins[jp] + ds;

    // energy difference
    double change = m_spins[j] * m_spins[jp] - s * sp;
    // right neighbour
    uint32_t jpn = (jp + 1) % N;
    // left neighbour
    uint32_t jn;
    if (j != 0) jn = j - 1;
    else jn = N;
    change += (m_spins[jn] - m_spins[jpn]) * ds;
    change *= m_j;

    // to flip or not to flip
    if (change < 0)
    {
        m_spins[j] = s;
        m_spins[jp] = sp;
        m_energy += change;
    }
    else if (m_chi(m_twister) < std::exp(-beta * change))
    {
        m_spins[j] = s;
        m_spins[jp] = sp;
        m_energy += change;
    }
}
