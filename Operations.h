#ifndef OPERATIONS_H
#define OPERATIONS_H

#include"CMatrix.h"
#include"QuantumCircuit.h"
#include<vector>
#include<complex>

class Hadamard : public Operation{
    public:
        Hadamard(int q);

        /*
        * @brief Applies the operation matrix to the given statevector
        *
        */
        void act(std::vector<std::complex<double> >& statevector) override;

    private:
        const static std::vector<std::vector<std::complex<double> > > op_matrix;
};

class Measure : public Operation{
    public:
        Measure(int q, int c);

        void act(std::vector<std::complex<double> >& statevector) override {};

        /*
        * @brief Measures the qubit given by this->qubits[0] and stores the result to cbit.
        * @param quantumRegister a list of pointers to all qubits in the circuit
        */
        void measure(std::vector<std::complex<double> >& sv, int& cbit) const override;
};

class CNot : public Operation{
    public:
        CNot(int q_control, int q_target);

        void act(std::vector<std::complex<double> >& statevector) override;

    private:
        const static std::vector<std::vector<std::complex<double> > > op_matrix;
};

#endif