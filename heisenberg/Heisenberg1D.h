//
// Created by mitja on 14.4.2019.
//

#ifndef VAJA_IV_1_HEISENBERG1D_H
#define VAJA_IV_1_HEISENBERG1D_H

#include <vector>
#include <random>
#include "pcg_random.hpp"

static const uint32_t N = 1 << 10;

struct Spin
{
    double x;
    double y;
    double z;
};

inline Spin operator+(const Spin& s1, const Spin& s2) {return {s1.x + s2.x, s1.y + s2.y, s1.z + s2.z};}
inline Spin operator-(const Spin& s1, const Spin& s2) {return {s1.x - s2.x, s1.y - s2.y, s1.z - s2.z};}
inline Spin operator-(const Spin& s) {return {-s.x, -s.y, -s.z};}
inline void operator-=(Spin& s1, const Spin& s2) {s1 = s1 - s2;}
inline double operator*(const Spin& s1, const Spin& s2) {return s1.x * s2.x + s1.y * s2.y + s1.z * s2.z;}
inline Spin operator*(double a, const Spin& s) {return {a * s.x, a * s.y, a * s.z};}
inline Spin vector_product(const Spin& s1, const Spin& s2)
{return {s1.y*s2.z - s1.z*s2.y, s1.z*s2.x - s1.x*s2.z, s1.x*s2.y - s1.y*s2.x};}
inline Spin operator/(const Spin& s1, double d) {return {s1.x / d, s1.y / d, s1.z / d};}
inline double abs(const Spin& s){return std::sqrt(s*s);}

class Heisenberg1D
{
public:
    Heisenberg1D(double j, double h);
    explicit Heisenberg1D(double j);

    inline double magnetization(){return abs(m_magnet);}
    inline double energy(){return m_energy;}
    inline double correlation(uint32_t r){return m_spins[0].z * m_spins[r].z;}

    void step(double beta);
    void stepCarefully(double beta);

private:
    Spin randomSpin();
    std::tuple<Spin, Spin, Spin> carefulSpin(const Spin& s1, const Spin& s2);

    std::vector<Spin> m_spins;

    std::uniform_int_distribution<uint32_t> m_spot;
    std::uniform_real_distribution<double> m_chi;
    std::uniform_real_distribution<double> m_phi;
    std::uniform_real_distribution<double> m_z;
    pcg32_fast m_twister;

    double m_j;
    double m_h;
    double m_energy;
    Spin m_magnet;
};


#endif //VAJA_IV_1_HEISENBERG1D_H
