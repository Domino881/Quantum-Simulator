#include"Operations.h"
#include"CMatrix.h"
#include"QuantumCircuit.h"
#include<vector>
#include<complex>
using namespace std;
using namespace complex_literals;

Hadamard::Hadamard(int q){
    this->qubits = {q};
    this->name='h';
}