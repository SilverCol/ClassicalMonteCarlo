//
// Created by mitja on 14.4.2019.
//

#ifndef VAJA_IV_1_HEISENBERG1D_H
#define VAJA_IV_1_HEISENBERG1D_H

#include <vector>
#include <random>

static const uint32_t N = 1 << 10;

struct Spin
{
    double x;
    double y;
    double z;
};

Spin operator+(const Spin& s1, const Spin& s2) {return {s1.x + s2.x, s1.y + s2.y, s1.z + s2.z};}
Spin operator-(const Spin& s1, const Spin& s2) {return {s1.x - s2.x, s1.y - s2.y, s1.z - s2.z};}
double operator*(const Spin& s1, const Spin& s2) {return s1.x * s2.x + s1.y * s2.y + s1.z * s2.z;}

class Heisenberg1D
{
public:
    Heisenberg1D(double j, double h);
    Heisenberg1D(double j);
    void step(double beta);
    void stepCarefully(double beta);
private:
    Spin randomSpin();

    std::vector<Spin> m_spins;

    std::uniform_int_distribution<uint32_t> m_spot;
    std::uniform_real_distribution<double> m_chi;
    std::uniform_real_distribution<double> m_phi;
    std::uniform_real_distribution<double> m_z;
    std::mt19937 m_twister;

    double m_j;
    double m_h;
    double m_energy;
};


#endif //VAJA_IV_1_HEISENBERG1D_H
