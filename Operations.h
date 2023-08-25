#ifndef OPERATIONS_H
#define OPERATIONS_H

#include"CMatrix.h"
#include"Qubit.h"
#include"QuantumCircuit.h"
#include<vector>
#include<complex>
using namespace std;
using namespace complex_literals;

class Hadamard : public Operation{
    public:
        Hadamard(int q);

        /*
        * @brief Applies the operation matrix to the given statevector
        *
        */
        void act(vector<complex<double> >& statevector);

    private:
        const vector<vector<complex<double> > > op_matrix={{sqrt(0.5), sqrt(0.5)},
                                                           {sqrt(0.5),-sqrt(0.5)}};
};

#endif