//
// Created by mitja on 13.4.2019.
//

#include<chrono>
#include <iostream>
#include <fstream>
#include <map>
#include "Ising2D.h"

static const size_t steps = N * 100;
static const uint8_t modulo = 10;

static const double J = 1.;
static const double H = 0.;

static const double iBeta = 1;
static const double dBeta = -.01;
static const uint32_t nBeta = 100;

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

    std::map<double, std::pair<std::vector<double>, std::vector<double> > > observables;

    std::cout << "Iterating over betas." << std::endl;
    double beta = iBeta;
    for (uint32_t n = 0; n < nBeta; ++n)
    {
        std::cout << "beta " << beta << std::endl;
        observables.emplace(beta, std::make_pair(std::vector<double> (), std::vector<double> ()));
        std::vector<double>& magnetVector = observables.at(beta).first;
        std::vector<double>& energyVector = observables.at(beta).second;

        for (size_t i = 0; i < steps; ++i)
        {
            lattice.step(beta);
            if (i % modulo == 0)
            {
                magnetVector.push_back(lattice.magnetization());
                energyVector.push_back(lattice.energy());
            }
        }

        beta += dBeta;
    }

    std::cout << "Calculating means." << std::endl;
    std::vector<double> output;
    for (auto& observable : observables)
    {
        output.push_back(observable.first); // beta

        std::vector<double>& magnetVector = observable.second.first;
        output.push_back(
                modulo * std::accumulate(magnetVector.begin(), magnetVector.end(), 0.0) / steps
                ); //mean mag
        output.push_back(
                modulo * std::inner_product(magnetVector.begin(), magnetVector.end(), magnetVector.begin(), 0.0) / steps
                ); //mean mag sqr

        std::vector<double>& energyVector = observable.second.second;
        output.push_back(
                modulo * std::accumulate(energyVector.begin(), energyVector.end(), 0.0) / steps
        ); //mean enery
        output.push_back(
                modulo * std::inner_product(energyVector.begin(), energyVector.end(), energyVector.begin(), 0.0) / steps
        ); //mean energy sqr
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
