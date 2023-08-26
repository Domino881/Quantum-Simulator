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
    double probs[2] = {std::abs(this->statevector[0])*std::abs(this->statevector[0]), 
                       std::abs(this->statevector[1])*std::abs(this->statevector[1])};

    const double tol = 1e-5;
    assert(std::abs(1-probs[0]-probs[1]) < tol);

    double r = (double)std::rand()/RAND_MAX;

    printf("probabilities: %.2f   %.2f, r=%.3f, result=%d\n", probs[0], probs[1], r, r>probs[0]);
    return r > probs[0];
}