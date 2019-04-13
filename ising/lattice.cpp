//
// Created by mitja on 13.4.2019.
//

#include<chrono>
#include <iostream>
#include <fstream>
#include "Ising2D.h"

static const size_t steps = N * 100;

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

    std::vector<double> betas;
    std::string latticeProgression;

    std::cout << "Iterating over betas." << std::endl;
    double beta = iBeta;
    for (uint32_t n = 0; n < nBeta; ++n)
    {
        std::cout << "beta " << beta << std::endl;
        betas.push_back(beta);

        for (size_t i = 0; i < steps; ++i)
        {
            lattice.step(beta);
            if (i % 16384 == 0)
            {
                for (uint32_t r = 0; r < N; ++r)
                {
                    if (lattice.lattice().test(r)) latticeProgression += '1';
                    else latticeProgression += '0';
                }
            }
        }

        beta += dBeta;
    }

    std::string bFile("../data/betas");
    bFile.append(std::to_string(N));
    bFile.append("_");
    bFile.append(std::to_string(J));
    bFile.append("_");
    bFile.append(std::to_string(H));
    bFile.append(".bin");
    std::cout << "Writing betas to file: " << bFile << std::endl;
    writeBinary(betas, bFile);

    std::string lFile("../data/lattice");
    lFile.append(std::to_string(N));
    lFile.append("_");
    lFile.append(std::to_string(J));
    lFile.append("_");
    lFile.append(std::to_string(H));
    lFile.append(".bin");
    std::cout << "Writing lattice progression to file: " << lFile << std::endl;

    std::ofstream output(lFile, std::ios::binary);
    output.write(latticeProgression.data(), sizeof(latticeProgression.data()));
    output.close();

    auto finish = std::chrono::high_resolution_clock::now();
    std::cout   << "Finished in "
                << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count()
                << "ms" << std::endl;
    return 0;
}
