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



    qc.measure(0,0);
    qc.measure(1,1);

    qc.run(10);
    qc.debug_print();
    qc.draw();

    auto results = qc.getCounts();
    for(auto x : results){
        printf("%s: %d\n", x.first.c_str(), x.second);
    }

}