#ifndef QUBIT_H
#define QUBIT_H

#include"CMatrix.h"
#include<complex>
#include<queue>
#include<ctime>
#include<cstdlib>

class Qubit{
    public:
        Qubit(): statevector({1.f,0.f}) { std::srand(std::time(nullptr)); };

        /*
        * @brief Uses a pseudo-random number to determine the state of the qubit after measuring.
        *  @returns Either 1 or 0, with probabilities dependent on the qubit's statevector
        * 
        * */
        bool measure();

        std::vector<std::complex<double> > statevector;
};

#endif