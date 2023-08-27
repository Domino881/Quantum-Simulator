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
// TODO:::: fix the ordering cause its backwards!!!
void Hadamard::act(std::vector<std::complex<double> >& statevector){

    int num_qubits = std::log2(statevector.size());

    CMatrix op(this->op_matrix);
    CMatrix id = diagMatrix(2);

    std::vector<CMatrix*> ops;
    for(int i=0;i<this->qubits[0];i++){
        ops.push_back(&id);
    }
    ops.push_back(&op);

    for(int i=0;i<num_qubits - this->qubits[0] - 1; i++){
        ops.push_back(&id);
    }
    CMatrix total_op = kroneckerProduct(ops);

    statevector = matmul(total_op, statevector);
}



Measure::Measure(int q, int c){
    this->qubits = {q};
    this->name = 'm';
    this->cbits = {c};
}
void Measure::measure(std::vector<std::complex<double> > sv, int& cbit) const{
}