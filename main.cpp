#include<cstdio>
#include<iostream>
#include<complex>
#include<bitset>
#include<cmath>

#include"QuantumCircuit.h"
using namespace std::complex_literals;

int main(){

    QuantumCircuit qc(5, 5);
    qc.h(0);
    qc.h(1);
    qc.h(2);
    qc.h(3);
    qc.x(4);

    double alpha = 2.0*M_PI/3.0;

    for(int i=0; i<4; i++){
        for(int j=0; j<1<<i; j++){
            qc.cPhase(alpha, i, 4);
        }
    }
    qc.barrier();

    qc.swap(0,3);
    qc.swap(1,2);

    qc.h(0);

    qc.cPhase(-M_PI/2.0, 0,1);
    qc.h(1);

    qc.cPhase(-M_PI/4.0, 0,2);
    qc.cPhase(-M_PI/2.0, 1,2);
    qc.h(2);

    qc.cPhase(-M_PI/8.0, 0,3);
    qc.cPhase(-M_PI/4.0, 1,3);
    qc.cPhase(-M_PI/2.0, 2,3);
    qc.h(3);

    qc.measure(0,0);
    qc.measure(1,1);
    qc.measure(2,2);
    qc.measure(3,3);

    int shots = 1000;
    qc.run(shots);
    qc.draw();
    // qc.debug_print();

    auto results = qc.getCounts();
    for(auto x : results){
        printf("%s:", x.first.c_str());
        for(int i=0;i<(50*x.second) /shots; i++)printf("#");
        printf(" %d\n", x.second);
    }

}