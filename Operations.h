#ifndef OPERATIONS_H
#define OPERATIONS_H

#include"CMatrix.h"
#include"Qubit.h"
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
        const std::vector<std::vector<std::complex<double> > > op_matrix;
};

class Measure : public Operation{
    public:
        Measure(int q, int c);

        void act(std::vector<std::complex<double> >& statevector) override {};

        void measure(std::vector<std::shared_ptr<Qubit> >& qubits, int& cbit) const override;
};

#endif