#include"Operations.h"
#include"CMatrix.h"
#include"Qubit.h"
#include"QuantumCircuit.h"
#include<vector>
#include<complex>
using namespace std;
using namespace complex_literals;

Hadamard::Hadamard(int q){
    this->qubits = {q};
    this->name='h';
    this->cbits = {};
}

void Hadamard::act(vector<complex<double> >& statevector){
    CMatrix op(this->op_matrix);
    statevector = matmul(op, statevector);
}