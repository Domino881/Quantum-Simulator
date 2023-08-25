#ifndef OPERATIONS_H
#define OPERATIONS_H

#include"CMatrix.h"
#include"QuantumCircuit.h"
#include<vector>
#include<complex>
using namespace std;
using namespace complex_literals;

class Hadamard : public Operation{
    public:
        Hadamard(int q);

        void act();

    private:
        const complex<double> op_matrix[2][2] =     {{sqrt(0.5), sqrt(0.5)},
                                                    {sqrt(0.5),-sqrt(0.5)}};
};

#endif