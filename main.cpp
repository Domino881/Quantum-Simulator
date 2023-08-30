#include<cstdio>
#include<iostream>
#include<complex>
#include<bitset>

#include"QuantumCircuit.h"
using namespace std::complex_literals;

int main(){

    QuantumCircuit qc(3, 3);
    qc.x(0);
    qc.cx(0,2);
    qc.measure(0,0);
    qc.measure(1,1);
    qc.measure(2,2);

    qc.run();
    qc.draw();
    qc.debug_print();
    auto results = qc.getCounts();
    for(auto x : results){
        printf("%s: %d\n", x.first.c_str(), x.second);
    }

}