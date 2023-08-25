#include"QuantumCircuit.h"
#include"Operations.h"
#include"CMatrix.h"
#include<complex>
#include<queue>
#include<cassert>
#include<ctime>
#include<iostream>
#include<cmath>
using namespace std;
using namespace complex_literals;

const bool Qubit::measure(){
    double probs[2] = {abs(this->statevector[0])*abs(this->statevector[0]), abs(this->statevector[1])*abs(this->statevector[1])};
    assert(abs(1-probs[0]-probs[1]) < 1e-5);

    srand(time(0));
    double r = (double)rand()/RAND_MAX;

    return r > probs[0];
}

Operation::Operation(){
    this->next = {};
    this->dependencies = 0;
    this->id = -1;
}

bool DagCompare::operator()(const Operation* a, const Operation* b) const{
    assert(a->id != -1 && b->id != -1);
    if(a->dependencies == b->dependencies)
        return a->id < b->id;

    return a->dependencies < b->dependencies;
}

QuantumCircuit::QuantumCircuit(int num_qbits, int num_cbits){
    this->id_counter=0;

    vector<double> cr;
    cr.resize(num_cbits);
    this->ClassicalRegister = cr;

    for(int i=0; i<num_qbits; i++){
        Qubit* q = new Qubit;
        this->QuantumRegister.push_back(q);
    }
}

void QuantumCircuit::h(int q){
    assert(q<this->QuantumRegister.size());
    Hadamard h(q);

    h.id = this->id_counter++;
    this->operations.push_back(h);
}

void QuantumCircuit::debug_print() const{
    printf("q register: ");
    for(auto x : this->QuantumRegister){
        printf("Qubit([%.1f, %.1f]) ", pow(abs(x->statevector[0]),2),pow(abs(x->statevector[1]),2));
    }
    printf("\nc register: ");
    for(auto x : this->ClassicalRegister)printf("%.2f ", x);

    printf("\noperations: ");
    for(auto x : operations){
        printf("(%c, [", x.name);
        for(auto y: x.qubits)printf("%p ", y);
        printf("],[");
        for(auto y: x.cbits)printf("%d ", y);
        printf("]) ");
    }
    cout<<"\n";
}

void QuantumCircuit::constructDag(){
    vector<Operation* > nextOp(this->QuantumRegister.size(), nullptr);

    for(int i=this->operations.size()-1; i>=0; i--){
        Operation& op = this->operations[i];
        for(auto q : op.qubits){
            if(nextOp[q] == nullptr)continue;

            op.next.push_back(nextOp[q]);
            nextOp[q]->dependencies++;

            nextOp[q] = &op;
        }

    }
}