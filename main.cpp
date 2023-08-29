#include<cstdio>
#include<iostream>
#include<complex>
#include<bitset>

#include"QuantumCircuit.h"
using namespace std::complex_literals;

int main(){

    const int numQubits = 2;
    const int numCbits = 2;

    QuantumCircuit qc(numQubits, numCbits);
    qc.h(0);
    qc.cx(0,1);
    qc.measure(0,0);
    qc.measure(1,1);

    qc.run(1000);
    qc.draw();
    auto results = qc.getCounts();
    for(auto x : results){
        const char* bitmask = std::bitset<10>(x.first).to_string().c_str();
        printf("%-.2s: %d\n", bitmask+10-numCbits, x.second);
    }

}