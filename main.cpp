#include<cstdio>
#include<string>
#include<vector>
#include<iostream>
#include<complex>
#include"QuantumCircuit.h"
using namespace std;
using namespace complex_literals;

int main(){

    QuantumCircuit qc(2,5);
    qc.h(0);
    qc.h(1);
    qc.measure(0,0);
    qc.h(0);
    qc.measure(0,1);
    qc.h(0);
    qc.measure(0,2);
    qc.h(0);
    qc.measure(0,3);
    qc.measure(0,4);


    qc.constructDag();
    qc.debug_print();
    qc.run();
    qc.debug_print();
    // vector<complex<double> > diag = {1,2,3};
    // CMatrix a = diagMatrix(diag);
    // vector<complex<double> > b = {1,2,3};

    // vector<complex<double> > r = matmul(a,b);
    // for(auto x: r)printf("%.2f ", real(x));
    

}