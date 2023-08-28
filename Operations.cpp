#include"Operations.h"
#include"CMatrix.h"
#include"Qubit.h"
#include"QuantumCircuit.h"
#include<vector>
#include<complex>
#include<cassert>
#include<cstdlib>

const std::vector<std::vector<std::complex<double> > > Hadamard::op_matrix = {
                                                                            {std::sqrt(0.5), std::sqrt(0.5)},
                                                                            {std::sqrt(0.5), -std::sqrt(0.5)}
                                                                            };

Hadamard::Hadamard(int q){
    this->qubits = {q};
    this->name="h";
    this->cbits = {};
}

void Hadamard::act(std::vector<std::complex<double> >& statevector){

    int num_qubits = std::log2(statevector.size());

    CMatrix op(this->op_matrix);
    CMatrix id = diagMatrix(2);

    std::vector<CMatrix*> ops;
    for(int i=0;i<num_qubits - this->qubits[0] - 1; i++){
        ops.push_back(&id);
    }
    ops.push_back(&op);
    for(int i=0;i<this->qubits[0];i++){
        ops.push_back(&id);
    }

    CMatrix total_op = kroneckerProduct(ops);

    statevector = matmul(total_op, statevector);
}



Measure::Measure(int q, int c){
    std::srand(std::time(nullptr));
    this->qubits = {q};
    this->name = "m";
    this->cbits = {c};
}
void Measure::measure(std::vector<std::complex<double> >& sv, int& cbit) const{
    int q = this->qubits[0];
    double sqrt_prob[2] = {0.f,0.f};

    long long int bs = 1<<q;
    for(unsigned i=0; i<sv.size(); i++){
        sqrt_prob[bs&i > 0] += std::abs(sv[i]);
    }
    double r = std::rand() / RAND_MAX;
    bool result = r < std::pow(sqrt_prob[1], 2);

    for(unsigned i=0; i<sv.size(); i++){
        if( (bs&i) != result ){sv[i]=0; continue;}

        sv[i] /= sqrt_prob[result];
    }
    cbit = (int)result;
}



const std::vector<std::vector<std::complex<double> > > CNot::op_matrix = {
                                                                         {1.f, 0.f, 0.f, 0.f},
                                                                         {0.f, 0.f, 0.f, 1.f},
                                                                         {0.f, 0.f, 1.f, 0.f},
                                                                         {0.f, 1.f, 0.f, 0.f}
                                                                         };
// TODO: fix for non-consecutive qubits
CNot::CNot(int q_control, int q_target){
    this->qubits = {q_control, q_target};
    assert(q_target-q_control == 1);

    this->name="cx";
    this->cbits = {};
}
void CNot::act(std::vector<std::complex<double> >& statevector){
    int num_qubits = std::log2(statevector.size());

    CMatrix op(this->op_matrix);
    CMatrix id = diagMatrix(2);

    std::vector<CMatrix*> ops;
    for(int i=0;i<num_qubits - this->qubits[1] - 1; i++){
        ops.push_back(&id);
    }
    ops.push_back(&op);
    for(int i=0;i<this->qubits[0];i++){
        ops.push_back(&id);
    }

    CMatrix total_op = kroneckerProduct(ops);
    statevector = matmul(total_op, statevector);
}