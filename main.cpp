#include<cstdio>
#include<iostream>
#include<complex>
#include<bitset>

#include"QuantumCircuit.h"
using namespace std::complex_literals;

int main(){

    QuantumCircuit qc(3, 3);
    qc.measure(1,1);
    qc.h(1);
    qc.h(0);
    qc.cx(0,2);
    qc.cPhase(0.5, 1,0);

    qc.run();
    qc.draw();

    auto results = qc.getCounts();
    for(auto x : results){
        printf("%s: %d\n", x.first.c_str(), x.second);
    }

}