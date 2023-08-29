#include"QuantumCircuit.h"

#include"Operations.h"

#include<cstdio>
#include<cassert>
#include<cmath>
#include<string>
#include<bitset>

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




QuantumCircuit::QuantumCircuit(int numQubits, int num_cbits): numQubits(numQubits){
    this->idCounter=0;

    std::vector<int> cr;
    cr.resize(num_cbits, 0);
    this->classicalRegister = cr;

    // for(int i=0; i<numQubits; i++){
    //     auto q = std::make_shared<Qubit>();
    //     this->QuantumRegister.push_back(q);
    // }

    std::vector<std::complex<double> > init_sv(1<<numQubits, 0.f);
    init_sv[0] = 1.f;
    this->totalStatevector = init_sv;
}

void QuantumCircuit::h(int q){
    assert(q<this->numQubits);
    std::shared_ptr<Operation> h = std::make_shared<Hadamard>(q);

    h->id = this->idCounter++;
    this->operations.push_back(h);
}

void QuantumCircuit::measure(int q, int c){
    assert(q<this->numQubits);
    assert(c<(int)this->classicalRegister.size());

    std::shared_ptr<Operation> measure = std::make_shared<Measure>(q,c);
    measure->id = this->idCounter++;
    this->operations.push_back(measure);
}

void QuantumCircuit::cx(int q_control, int q_target){
    assert(q_control<this->numQubits && q_target<this->numQubits);

    std::shared_ptr<Operation> cx = std::make_shared<CNot>(q_control, q_target);
    cx->id = this->idCounter++;
    this->operations.push_back(cx);
}

void QuantumCircuit::debug_print() const{
    const static int bs_size = 2;
    printf("\n%-13s", "q register: ");
    for(int b=0;b<(1<<this->numQubits); b++){
        // assuming max 4 qubits!!
        std::bitset<bs_size> bs(b);
        printf("[%s]: %.2f  ", bs.to_string().c_str(), pow(abs(this->totalStatevector[b]),2));
    }
    printf("\n%-13s", "c register: ");
    for(auto x : this->classicalRegister)printf("%d  ", x);

    printf("\n%-13s", "operations: ");

    const auto ops = &(this->operations);
    for(unsigned i=0;i<ops->size();i++){
        printf("(%s%d, [", ops->at(i)->name.c_str(), ops->at(i)->id);
        for(auto y: ops->at(i)->qubits)printf("%d ", y);
        printf("],[");
        for(auto y: ops->at(i)->cbits)printf("%d ", y);
        printf("])  ");
    }

    printf("\n%-13s", "dag: ");
    auto dag_copy = this->dag;
    while(!dag_copy.empty()){
        auto op = dag_copy.top();
        dag_copy.pop();
        printf("(%s%d, [", op->name.c_str(), op->id);
        for(auto y: op->qubits)printf("%d ", y);
        printf("],[");
        for(auto y: op->cbits)printf("%d ", y);
        printf("], %d)  ", op->dependencies);
    }
    printf("\n");
}
void QuantumCircuit::draw() const{
    std::vector<std::vector<char> > blocks{};
    for(int i=0;i<this->numQubits;i++){
        std::vector<char> v(this->operations.size() * 5 + 4, '-');
        blocks.push_back(v);
        std::vector<char> w(this->operations.size() * 5 + 4, ' ');
        blocks.push_back(w);
    }

    for(unsigned col=0;col<this->operations.size();col++){
        if(this->operations[col]->name == "h"){
            blocks[2*this->operations[col]->qubits[0]][2 + 5*col - 1] = '[';
            blocks[2*this->operations[col]->qubits[0]][2 + 5*col] = 'H';
            blocks[2*this->operations[col]->qubits[0]][2 + 5*col + 1] = ']';
        }
        if(this->operations[col]->name == "cx"){
            blocks[2*this->operations[col]->qubits[0]][2 + 5*col] = 'o';

            int mid = (this->operations[col]->qubits[0] + this->operations[col]->qubits[1]);
            blocks[mid][2 + 5*col] = '|';

            blocks[2*this->operations[col]->qubits[1]][2 + 5*col-1] = '[';
            blocks[2*this->operations[col]->qubits[1]][2 + 5*col] = 'X';
            blocks[2*this->operations[col]->qubits[1]][2 + 5*col+1] = ']';
        }
        if(this->operations[col]->name == "m"){
            blocks[2*this->operations[col]->qubits[0]][2 + 5*col - 1] = '(';
            blocks[2*this->operations[col]->qubits[0]][2 + 5*col] = 'M';
            blocks[2*this->operations[col]->qubits[0]][2 + 5*col + 1] = ')';
        }
    }
    printf("\n");
    for(int row=0; row<2*this->numQubits; row++){
        if(row%2==0)printf("|0> ");
        else printf("    ");

        for(unsigned col=0; col<blocks[row].size(); col++){
            printf("%c", blocks[row][col]);
        }
        printf("\n");
    }
}

void QuantumCircuit::constructDag(){
    std::vector<std::shared_ptr<Operation> > nextOp(this->numQubits, nullptr);

    for(int i=(int)this->operations.size()-1; i>=0; i--){
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
    for(unsigned i=0;i<this->operations.size();i++){
        this->dag.push(this->operations[i]);
    }
}

void QuantumCircuit::run(){
    assert(!this->dag.empty());

    while(!this->dag.empty()){
        auto op = this->dag.top();

        // Measurements treated separately
        if(op->name == "m"){

            // The whole quantum register is added because of possible entanglement
            op->measure(this->totalStatevector, this->classicalRegister[op->cbits[0]]);
        }

        op->act(this->totalStatevector);
        this->dag.pop();
    }
    // controversial:
    this->operations.clear();
}