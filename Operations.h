#pragma once

/*
* @file Operations.h
* @brief Declares the operation subclasses, including their act() and measure() functions
* @author Dominik Kuczynski
*/

#include<vector>
#include<complex>

#include"QuantumCircuit.h"

class Hadamard : public Operation{
    public:
        Hadamard(int q);

        /*
        * @brief Applies the operation matrix to the given statevector
        *
        */
        void act(std::vector<std::complex<double> >& statevector) override;

    private:
        const static std::vector<std::vector<std::complex<double> > > operationMatrix;
};

class Measure : public Operation{
    public:
        Measure(int q, int* c);

        void act(std::vector<std::complex<double> >& statevector) override;
};

class CNot : public Operation{
    public:
        CNot(int qControl, int qTarget);

        void act(std::vector<std::complex<double> >& statevector) override;

    private:
        const static std::vector<std::vector<std::complex<double> > > operationMatrix;
};