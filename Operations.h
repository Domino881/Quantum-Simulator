#pragma once

/*
* @file Operations.h
* @brief Declares the operation subclasses, including their act() and measure() functions
* @author Dominik Kuczynski
*/

#include<vector>
#include<complex>

#include"QuantumCircuit.h"
#include"CMatrix.h"

class singleQubitGate : public Operation{
    public:
        /*
        * @brief Constructs a single qubit (no cbit) gate with a given operation matrix
        * @param label label of the gate, e.g. "h" for Hadamard
        * @param q qubit on which the gate is to act
        */
        singleQubitGate(const std::string& label, const int& q,
                        const CMatrix& operationMatrix);

        ~singleQubitGate() {};

        /*
        * @brief Applies the operation matrix to the given statevector
        * @returns The total operation matrix (to be cached)
        */
        CMatrix act(std::vector<std::complex<double> >& statevector) override;

    private:
        const CMatrix operationMatrix;

        // Cache total ops ?? 
        static std::map<std::string, CMatrix> cachedTotalOps;
};

class ControlledGate : public Operation{
    public:
        ControlledGate(const std::string& label, const int& qControl, const int& qTarget,
                       const CMatrix& operationMatrix);
        ~ControlledGate() {};

        /*
        * @brief Acts with the controlled-NOT operation on the given statevector
        * @returns The total operation matrix (to be cached)
        */
        CMatrix act(std::vector<std::complex<double> >& statevector) override;

    private:
        const CMatrix operationMatrix;
};

class SwapGate : public Operation{
    public:
        SwapGate(const std::string& label, const int& q1, const int& q2);
        ~SwapGate() {};

        /*
        * @brief Swaps the states of two qubits
        * @returns The total operation matrix (to be cached)
        */
        CMatrix act(std::vector<std::complex<double> >& statevector) override;
    private:
        const CMatrix operationMatrix;
};

class Measure : public Operation{
    public:
        Measure(int q, int* c);
        ~Measure() {};

        /*
        * @brief Measures the given qubits to the classical register
        * @returns A CMatrix with dimension 0
        */
        CMatrix act(std::vector<std::complex<double> >& statevector) override;
};
