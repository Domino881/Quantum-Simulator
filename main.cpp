#include<cstdio>
#include<iostream>
#include<complex>
#include<bitset>

#include"QuantumCircuit.h"
using namespace std::complex_literals;

int main(){

    QuantumCircuit qc(2, 2);
    qc.x(0);
    qc.h(1);
    qc.cPhase(3, 0,1);
    qc.h(1);
    qc.measure(0,0);
    qc.measure(1,1);

    qc.run(5000);
    qc.draw();

    auto results = qc.getCounts();
    for(auto x : results){
        printf("%s: %d\n", x.first.c_str(), x.second);
    }

}