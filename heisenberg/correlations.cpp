//
// Created by mitja on 13.4.2019.
//

#include<chrono>
#include <iostream>
#include <fstream>
#include <map>
#include "Heisenberg1D.h"

static const size_t steps = 1 << 20;
static const size_t modulo = 1 << 10;

static const double J = 1.;
static const double H = 0.;

static const double iBeta = 100;
static const double dBeta = -1;
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

    Heisenberg1D chain(J);

    std::map<double, std::vector<std::vector<double> > > correlations;

    std::cout << "Iterating over betas." << std::endl;
    double beta = iBeta;
    for (uint32_t n = 0; n < nBeta; ++n)
    {
        std::cout << "beta " << beta << std::endl;
        correlations.emplace(beta, std::vector<std::vector<double> >(N));
        std::vector<std::vector<double> >& correlationTable = correlations.at(beta);

        for (size_t i = 0; i < steps; ++i)
        {
            chain.stepCarefully(beta);
            if (i % modulo == 0)
            {
                for (uint32_t j = 0; j < N; ++j)
                {
                    correlationTable[j].push_back(chain.correlation(j));
                }
            }
        }

        beta += dBeta;
    }

    std::cout << "Calculating means." << std::endl;
    std::vector<double> output;
    for (auto& entry : correlations)
    {
        output.push_back(entry.first); // beta

        std::vector<std::vector<double> >& correlationTable = entry.second;
        for (const auto& spot : correlationTable)
        {
            output.push_back( modulo * std::accumulate(spot.begin(), spot.end(), 0.0) / steps ); // mean
        }
    }

    std::string file("../../heisenberg/data/");
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
