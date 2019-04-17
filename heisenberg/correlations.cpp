//
// Created by mitja on 13.4.2019.
//

#include<chrono>
#include <iostream>
#include <fstream>
#include "Heisenberg1D.h"

static const size_t steps = 1 << 15;
static const size_t modulo = 1 << 10;

static const double J = 1.;
static const double H = 0.;

static const double iBeta = 20;
static const double dBeta = .1;
static const uint32_t nBeta = 200;

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

    std::vector<double> output1;

    std::cout << "Iterating over betas." << std::endl;
    double beta = iBeta;
    for (uint32_t n = 0; n < nBeta; ++n)
    {
        std::cout << "beta " << beta << std::endl;
        std::vector<std::vector<double> > correlationTableX(N);
        std::vector<std::vector<double> > correlationTableY(N);
        std::vector<std::vector<double> > correlationTableZ(N);

        for (size_t i = 0; i < steps; ++i)
        {
            chain.stepCarefully(beta);
            if (i % modulo == 0)
            {
                for (uint32_t j = 0; j < N; ++j)
                {
                    correlationTableX[j].push_back(chain.correlationX(j));
                    correlationTableY[j].push_back(chain.correlationY(j));
                    correlationTableZ[j].push_back(chain.correlationZ(j));
                }
            }
        }

        output1.push_back(beta); // beta

        for (size_t j = 0; j < correlationTableX.size(); ++j)
        {
            output1.push_back( modulo * std::accumulate(correlationTableX[j].begin(), correlationTableX[j].end(), 0.0) / steps ); // mean
            output1.push_back( modulo * std::accumulate(correlationTableY[j].begin(), correlationTableY[j].end(), 0.0) / steps ); // mean
            output1.push_back( modulo * std::accumulate(correlationTableZ[j].begin(), correlationTableZ[j].end(), 0.0) / steps ); // mean
        }

        beta -= dBeta;
    }
    
    std::vector<double> output2;

    std::cout << "Iterating over betas." << std::endl;
    for (uint32_t n = 0; n < nBeta; ++n)
    {
        std::cout << "beta " << beta << std::endl;
        std::vector<std::vector<double> > correlationTableX(N);
        std::vector<std::vector<double> > correlationTableY(N);
        std::vector<std::vector<double> > correlationTableZ(N);

        for (size_t i = 0; i < steps; ++i)
        {
            chain.stepCarefully(beta);
            if (i % modulo == 0)
            {
                for (uint32_t j = 0; j < N; ++j)
                {
                    correlationTableX[j].push_back(chain.correlationX(j));
                    correlationTableY[j].push_back(chain.correlationY(j));
                    correlationTableZ[j].push_back(chain.correlationZ(j));
                }
            }
        }

        output2.push_back(beta); // beta

        for (size_t j = 0; j < correlationTableX.size(); ++j)
        {
            output2.push_back( modulo * std::accumulate(correlationTableX[j].begin(), correlationTableX[j].end(), 0.0) / steps ); // mean
            output2.push_back( modulo * std::accumulate(correlationTableY[j].begin(), correlationTableY[j].end(), 0.0) / steps ); // mean
            output2.push_back( modulo * std::accumulate(correlationTableZ[j].begin(), correlationTableZ[j].end(), 0.0) / steps ); // mean
        }

        beta += dBeta;
    }
    
    std::string file1("../../heisenberg/data/heating");
    file1.append(std::to_string(N));
    file1.append("_");
    file1.append(std::to_string(J));
    file1.append("_");
    file1.append(std::to_string(H));
    file1.append(".bin");
    std::cout << "Writing to file: " << file1 << std::endl;
    writeBinary(output1, file1);

    std::string file2("../../heisenberg/data/cooling");
    file2.append(std::to_string(N));
    file2.append("_");
    file2.append(std::to_string(J));
    file2.append("_");
    file2.append(std::to_string(H));
    file2.append(".bin");
    std::cout << "Writing to file: " << file2 << std::endl;
    writeBinary(output2, file2);
   
    auto finish = std::chrono::high_resolution_clock::now();
    std::cout   << "Finished in "
                << std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count()
                << "ms" << std::endl;
    return 0;
}
