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
        */
        void act(std::vector<std::complex<double> >& statevector) override;

    private:
        const CMatrix operationMatrix;
};

class Measure : public Operation{
    public:
        Measure(int q, int* c);
        ~Measure() {};

        /*
        * @brief Measures the given qubits to the classical register
        */
        void act(std::vector<std::complex<double> >& statevector) override;
};

class ControlledGate : public Operation{
    public:
        ControlledGate(const std::string& label, const int& qControl, const int& qTarget,
                       const CMatrix& operationMatrix);
        ~ControlledGate() {};

        /*
        * @brief Acts with the controlled-NOT operation on the given statevector
        */
        void act(std::vector<std::complex<double> >& statevector) override;

    private:
        const CMatrix operationMatrix;
};