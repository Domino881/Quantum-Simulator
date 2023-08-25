#include"QuantumCircuit.h"
#include"Operations.h"
#include"CMatrix.h"
#include"Qubit.h"
#include<complex>
#include<queue>
#include<cassert>
#include<ctime>
#include<iostream>
#include<cmath>
using namespace std;
using namespace complex_literals;

Operation::Operation(){
    this->next = {};
    this->dependencies = 0;
    this->id = -1;
}

bool DagCompare::operator()(const Operation* a, const Operation* b) const{
    assert(a->id != -1 && b->id != -1);
    if(a->dependencies == b->dependencies)
        return a->id > b->id;

    return a->dependencies > b->dependencies;
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
    Hadamard* h = new Hadamard(q);

    h->id = this->id_counter++;
    this->operations.push_back(h);
}

void QuantumCircuit::debug_print() const{
    printf("\n%-13s", "q register: ");
    for(auto x : this->QuantumRegister){
        printf("Qubit([%.1f, %.1f])  ", pow(abs(x->statevector[0]),2),pow(abs(x->statevector[1]),2));
    }
    printf("\n%-13s", "c register: ");
    for(auto x : this->ClassicalRegister)printf("%.2f  ", x);

    printf("\n%-13s", "operations: ");
    for(auto x : operations){
        printf("(%c%d, [", x->name, x->id);
        for(auto y: x->qubits)printf("%d ", y);
        printf("],[");
        // for(auto y: x->cbits)printf("%d ", y);
        printf("])  ");
    }

    printf("\n%-13s", "dag: ");
    priority_queue<Operation*, vector<Operation*>, DagCompare> dag_copy;
    dag_copy = this->dag;
    while(!dag_copy.empty()){
        Operation* op = dag_copy.top();
        dag_copy.pop();
        printf("(%c%d, [", op->name, op->id);
        for(auto y: op->qubits)printf("%d ", y);
        printf("],[");
        for(auto y: op->cbits)printf("%d ", y);
        printf("], %d)  ", op->dependencies);
    }
    printf("\n");
}

void QuantumCircuit::constructDag(){
    vector<Operation* > nextOp(this->QuantumRegister.size(), nullptr);

    for(int i=this->operations.size()-1; i>=0; i--){
        Operation* op = this->operations[i];
        for(auto q : op->qubits){
            if(nextOp[q] == nullptr){
                nextOp[q] = op;
                continue;
            }

            op->next.push_back(nextOp[q]);
            nextOp[q]->dependencies++;

            nextOp[q] = op;
        }
    }
    for(int i=0;i<this->operations.size();i++){
        this->dag.push(this->operations[i]);
    }
}

void QuantumCircuit::run(){
    assert(!this->dag.empty());

    while(!this->dag.empty()){
        Operation* op = this->dag.top();
        vector<complex<double> >* sv;
        //TODO *TODO* _TODO_ fix for multiple qubit gates!
        sv = &(this->QuantumRegister[op->qubits[0]]->statevector);
        op->act(*sv);
        this->dag.pop();
    }
}