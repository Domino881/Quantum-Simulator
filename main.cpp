#include<cstdio>
#include<string>
#include<vector>
#include<iostream>
#include<complex>
#include"QuantumCircuit.h"
using namespace std;
using namespace complex_literals;

int main(){

    QuantumCircuit qc(2,1);
    qc.h(0);
    qc.h(1);

    qc.debug_print();

}