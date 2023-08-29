#include"QuantumCircuit.h"

#include"Operations.h"

#include<cstdio>
#include<cassert>
#include<cmath>
#include<string>
#include<bitset>
#include<ctime>
#include<random>

Operation::Operation(){
    this->next = {};
    this->dependencies = 0;
    this->id = -1;
}

QuantumCircuit::QuantumCircuit(int numQubits, int num_cbits): numQubits(numQubits){
    this->idCounter=0;
    std::srand(std::time(nullptr));

    std::vector<int> cr;
    cr.resize(num_cbits, 0);
    this->classicalRegister = cr;

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

    std::shared_ptr<Operation> measure = std::make_shared<Measure>(q,&(this->classicalRegister[c]));
    measure->id = this->idCounter++;
    this->operations.push_back(measure);
}

void QuantumCircuit::cx(int qControl, int qTarget){
    assert(qControl<this->numQubits && qTarget<this->numQubits);

    std::shared_ptr<Operation> cx = std::make_shared<CNot>(qControl, qTarget);
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

        for(auto y: ops->at(i)->cbits)printf("%lld ", (y - &this->classicalRegister[0]));

        printf("])");
        printf(i%5!=4?"  ":"\n             ");
    }

    printf("\n%-13s", "sortedDag: ");

    for(unsigned i=0;i<this->sortedDag.size();i++){
        auto op = this->sortedDag[i];
        printf("(%s%d, [", op->name.c_str(), op->id);

        for(auto y: op->qubits)printf("%d ", y);
        printf("],[");

        for(auto y: op->cbits)printf("%lld ", (y - &this->classicalRegister[0]));
        printf("], %d)", op->dependencies);

        printf(i%5!=4?" > ":" >\n             ");
    }
    printf("\n");
}
void QuantumCircuit::draw() const{
    std::vector<std::vector<char> > blocks{};
    for(int i=0;i<this->numQubits;i++){
        std::vector<char> v(2, '-');
        blocks.push_back(v);
        std::vector<char> w(this->operations.size() * 5 + 4, ' ');
        blocks.push_back(w);
    }

    auto copysortedDag = this->sortedDag;
    assert(this->sortedDag.size() == this->operations.size());

    for(unsigned i=0;i<this->sortedDag.size();i++){
        auto op = this->sortedDag[i];
        auto copyQubits = op->qubits;

        if(op->name == "h"){
            blocks[2*copyQubits[0]].push_back('[');
            blocks[2*copyQubits[0]].push_back('H');
            blocks[2*copyQubits[0]].push_back(']');
            blocks[2*copyQubits[0]].push_back('-');
            blocks[2*copyQubits[0]].push_back('-');
        }
        if(op->name == "cx"){
            int maxSize = std::max(blocks[2*copyQubits[0]].size(), blocks[2*copyQubits[1]].size());

            while((int)blocks[2*copyQubits[0]].size() < maxSize)
                blocks[2*copyQubits[0]].push_back('-');
            while((int)blocks[2*copyQubits[1]].size() < maxSize)
                blocks[2*copyQubits[1]].push_back('-');

            blocks[2*copyQubits[0]].push_back('-');
            blocks[2*copyQubits[0]].push_back('*');
            blocks[2*copyQubits[0]].push_back('-');
            blocks[2*copyQubits[0]].push_back('-');
            blocks[2*copyQubits[0]].push_back('-');

            blocks[2*copyQubits[1]].push_back('[');
            blocks[2*copyQubits[1]].push_back('X');
            blocks[2*copyQubits[1]].push_back(']');
            blocks[2*copyQubits[1]].push_back('-');
            blocks[2*copyQubits[1]].push_back('-');

            blocks[copyQubits[0]+copyQubits[1]][1+maxSize] = '|';
        }
        if(op->name == "m"){
            blocks[2*copyQubits[0]].push_back('[');
            blocks[2*copyQubits[0]].push_back('M');
            blocks[2*copyQubits[0]].push_back(']');
            blocks[2*copyQubits[0]].push_back('-');
            blocks[2*copyQubits[0]].push_back('-');
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
    // Clears the sortedDag
    this->sortedDag.clear();
    
    std::queue<std::shared_ptr<Operation> > queue{};

    for(auto op : this->operations){
        if(op->dependencies == 0){
            queue.push(op);
        }
    }
    while(!queue.empty()){
        auto op = queue.front();
        queue.pop();

        for(auto nextOp : op->next){
            if(--nextOp->dependencies == 0){
                queue.push(nextOp);
            }
        }
        this->sortedDag.push_back(op);
    }
}

void QuantumCircuit::run(){
    if(this->sortedDag.size() < this->operations.size())
        this->constructDag();

    auto sortedDagCopy = this->sortedDag;
    for(auto op : sortedDag){

        op->act(this->totalStatevector);

    }
}