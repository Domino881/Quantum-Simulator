#include"CMatrix.h"
#include"Qubit.h"
#include<complex>
#include<queue>
#include<cassert>
#include<ctime>
#include<iostream>
#include<cstdlib>
#include<cmath>

bool Qubit::measure(){
    double probabilities[2] = {std::abs(this->statevector[0])*std::abs(this->statevector[0]), 
                       std::abs(this->statevector[1])*std::abs(this->statevector[1])};

    const double tol = 1e-5;
    assert(std::abs(1-probabilities[0]-probabilities[1]) < tol);

    double r = (double)std::rand()/RAND_MAX;

    printf("probabilities: %.2f   %.2f, r=%.3f, result=%d\n", probabilities[0], probabilities[1], r, r>probabilities[0]);
    return r > probabilities[0];
}