//
// Created by mitja on 13.4.2019.
//

#include<chrono>
#include <iostream>
#include <fstream>
#include "Ising2D.h"

static const size_t steps = N * 100;
static const uint8_t modulo = 1;

static const double J = 1.;
static const double H = 0.;

static const double iBeta = .6;
static const double dBeta = -.0004;
static const uint32_t nBeta = 1000;

void writeBinary(std::vector<double>& data, const std::string& file)
{
    std::ofstream output(file, std::ios::binary);
    for (double& x : data)
    {
        output.write(reinterpret_cast<char*>(&x), sizeof(x));
    }
    output.close();
}

int main()
{
    auto start = std::chrono::high_resolution_clock::now();

    // Ising2D lattice(J, H, randomState());
    Ising2D lattice(J, H);

    std::vector<double> output;

    std::cout << "Iterating over betas." << std::endl;
    double beta = iBeta;
    for (uint32_t n = 0; n < nBeta; ++n)
    {
        std::cout << "beta " << beta << std::endl;
        std::vector<double> magnetVector;
        std::vector<double> energyVector;

        for (size_t i = 0; i < steps; ++i)
        {
            lattice.step(beta);
            if (i % modulo == 0)
            {
                magnetVector.push_back(lattice.magnetization());
                energyVector.push_back(lattice.energy());
            }
        }

        output.push_back(beta); // beta

        output.push_back(
                modulo * std::accumulate(magnetVector.begin(), magnetVector.end(), 0.0) / steps
                ); //mean mag
        output.push_back(
                modulo * std::inner_product(magnetVector.begin(), magnetVector.end(), magnetVector.begin(), 0.0) / steps
                ); //mean mag sqr

        output.push_back(
                modulo * std::accumulate(energyVector.begin(), energyVector.end(), 0.0) / steps
        ); //mean enery
        output.push_back(
                modulo * std::inner_product(energyVector.begin(), energyVector.end(), energyVector.begin(), 0.0) / steps
        ); //mean energy sqr

        beta += dBeta;
    }

    std::string file("../../ising/data/");
    if (dBeta < 0) file.append("heating");
    else file.append("cooling");
    file.append(std::to_string(N));
    file.append("_");
    file.append(std::to_string(J));
    file.append("_");
    file.append(std::to_string(H));
    file.append(".bin");

    std::cout << "Writing to file: " << file << std::endl;
    writeBinary(output, file);

    auto finish = std::chrono::high_resolution_clock::now();
    std::cout   << "Finished in "
                << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count()
                << "ms" << std::endl;
    return 0;
}
