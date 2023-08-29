#include<cstdio>
#include<iostream>
#include<complex>
#include<bitset>

#include"QuantumCircuit.h"
using namespace std::complex_literals;

int main(){

    QuantumCircuit qc(2,5);
    qc.h(0);
    qc.measure(0,0);

    qc.debug_print();
    qc.run();
    qc.draw();
    auto results = qc.getCounts();
    for(auto x : results){
        printf("%-.2s: %d\n", std::bitset<10>(x.first).to_string().c_str(), x.second);
    }
    // vector<complex<double> > diag = {1,2,3};
    // CMatrix a = diagMatrix(diag);
    // vector<complex<double> > b = {1,2,3};

    // vector<complex<double> > r = matmul(a,b);
    // for(auto x: r)printf("%.2f ", real(x));
    

}