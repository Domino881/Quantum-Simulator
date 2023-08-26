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
#include<memory>

Operation::Operation(){
    this->next = {};
    this->dependencies = 0;
    this->id = -1;
}

bool DagCompare::operator()(const std::shared_ptr<Operation> a, const std::shared_ptr<Operation> b) const{
    assert(a->id != -1 && b->id != -1);
    if(a->dependencies == b->dependencies)
        return a->id > b->id;

    return a->dependencies > b->dependencies;
}




QuantumCircuit::QuantumCircuit(int num_qbits, int num_cbits){
    this->id_counter=0;

    std::vector<int> cr;
    cr.resize(num_cbits, 0);
    this->ClassicalRegister = cr;

    for(int i=0; i<num_qbits; i++){
        auto q = std::make_shared<Qubit>();
        this->QuantumRegister.push_back(q);
    }
}

void QuantumCircuit::h(int q){
    assert(q<this->QuantumRegister.size());
    std::shared_ptr<Operation> h = std::make_shared<Hadamard>(q);

    h->id = this->id_counter++;
    this->operations.push_back(h);
}

void QuantumCircuit::measure(int q, int c){
    assert(q<this->QuantumRegister.size());
    assert(c<this->ClassicalRegister.size());

    std::shared_ptr<Operation> measure = std::make_shared<Measure>(q,c);
    measure->id = this->id_counter++;
    this->operations.push_back(measure);
}

void QuantumCircuit::debug_print() const{
    printf("\n%-13s", "q register: ");
    for(auto x : this->QuantumRegister){
        printf("Qubit([%.1f, %.1f])  ", pow(abs(x->statevector[0]),2),pow(abs(x->statevector[1]),2));
    }
    printf("\n%-13s", "c register: ");
    for(auto x : this->ClassicalRegister)printf("%d  ", x);

    printf("\n%-13s", "operations: ");

    const auto ops = &(this->operations);
    for(int i=0;i<ops->size();i++){
        printf("(%c%d, [", ops->at(i)->name, ops->at(i)->id);
        for(auto y: ops->at(i)->qubits)printf("%d ", y);
        printf("],[");
        // for(auto y: x->cbits)pruniqueintf("%d ", y);
        printf("])  ");
    }

    printf("\n%-13s", "dag: ");
    auto dag_copy = this->dag;
    while(!dag_copy.empty()){
        auto op = dag_copy.top();
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
    std::vector<std::shared_ptr<Operation> > nextOp(this->QuantumRegister.size(), nullptr);

    for(int i=this->operations.size()-1; i>=0; i--){
        for(auto q : this->operations[i]->qubits){
            if(nextOp[q] == nullptr){
                nextOp[q] = this->operations[i];
                continue;
            }

            this->operations[i]->next.push_back(nextOp[q]);
            nextOp[q]->dependencies++;

            nextOp[q] = this->operations[i];
        }
    }
    for(int i=0;i<this->operations.size();i++){
        this->dag.push(this->operations[i]);
    }
}

void QuantumCircuit::run(){
    assert(!this->dag.empty());

    while(!this->dag.empty()){
        auto op = this->dag.top();

        // Measurements treated separately
        if(op->name == 'm'){

            //TODO make this the entire statevector!!
            int* cbit = &(this->ClassicalRegister[op->cbits[0]]);
            op->measure(this->QuantumRegister, *cbit);
        }

        //TODO *TODO* _TODO_ fix for multiple qubit gates!
        op->act(this->QuantumRegister[op->qubits[0]]->statevector);
        this->dag.pop();
    }
}