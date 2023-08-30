#include<cstdio>
#include<iostream>
#include<complex>
#include<bitset>

#include"QuantumCircuit.h"
using namespace std::complex_literals;
#define numCbits 3
#define numQubits 3

int main(){

    QuantumCircuit qc(numQubits, numCbits);
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
        const char* bitmask = std::bitset<10>(x.first).to_string().c_str();
        printf("%s: %d\n", bitmask+10-numCbits, x.second);
    }

}