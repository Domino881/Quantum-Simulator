#include"Operations.h"
#include"CMatrix.h"
#include"Qubit.h"
#include"QuantumCircuit.h"
#include<vector>
#include<complex>

Hadamard::Hadamard(int q): op_matrix({{std::sqrt(0.5), -std::sqrt(0.5)},
                                    {std::sqrt(0.5), -std::sqrt(0.5)}}){
    this->qubits = {q};
    this->name='h';
    this->cbits = {};
}

void Hadamard::act(std::vector<std::complex<double> >& statevector){
    CMatrix op(this->op_matrix);
    statevector = matmul(op, statevector);
}

Measure::Measure(int q, int c){
    this->qubits = {q};
    this->name = 'm';
    this->cbits = {c};
}

void Measure::measure(std::vector<std::shared_ptr<Qubit> >& qubits, int& cbit) const{
    cbit = qubits[this->qubits[0]]->measure();
    qubits[this->qubits[0]]->statevector = {(double)cbit==0, (double)cbit==1};
}