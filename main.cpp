#include<cstdio>
#include<iostream>
#include<complex>

#include"QuantumCircuit.h"
using namespace std::complex_literals;

int main(){

    QuantumCircuit qc(2,5);
    qc.h(0);
    qc.measure(0, 0);
    qc.h(0);
    qc.measure(0, 1);
    qc.h(0);
    qc.measure(0, 2);
    qc.h(0);
    qc.measure(0, 3);
    qc.h(0);
    qc.measure(0, 4);

    qc.run();
    qc.debug_print();
    qc.draw();
    // vector<complex<double> > diag = {1,2,3};
    // CMatrix a = diagMatrix(diag);
    // vector<complex<double> > b = {1,2,3};

    // vector<complex<double> > r = matmul(a,b);
    // for(auto x: r)printf("%.2f ", real(x));
    

}