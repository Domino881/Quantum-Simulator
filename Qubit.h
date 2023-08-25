#ifndef QUBIT_H
#define QUBIT_H

#include"CMatrix.h"
#include<complex>
#include<queue>
using namespace std;
using namespace complex_literals;

class Qubit{
    public:
        Qubit(): statevector({1,0}) {};

        /*
        * @brief Uses a pseudo-random number to determine the state of the qubit after measuring.
        *  @returns Either 1 or 0, with probabilities dependent on the qubit's statevector
        * 
        * */
        const bool measure();

        vector<complex<double> > statevector;
};

#endif