#include"CMatrix.h"
#include"Qubit.h"
#include<complex>
#include<queue>
#include<cassert>
#include<ctime>
#include<iostream>
#include<cmath>
using namespace std;
using namespace complex_literals;

const bool Qubit::measure(){
    double probs[2] = {abs(this->statevector[0])*abs(this->statevector[0]), abs(this->statevector[1])*abs(this->statevector[1])};
    assert(abs(1-probs[0]-probs[1]) < 1e-5);

    srand(time(0));
    double r = (double)rand()/RAND_MAX;

    return r > probs[0];
}