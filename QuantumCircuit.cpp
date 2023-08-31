#include"QuantumCircuit.h"

#include"Operations.h"
#include"CMatrix.h"

#include<cstdio>
#include<cassert>
#include<cmath>
#include<string>
#include<bitset>
#include<ctime>
#include<random>

using VCD = std::vector<std::complex<double> >; 

std::map<std::string, const CMatrix> QuantumCircuit::opMatrices = {
    {"hadamard", CMatrix({{std::sqrt(0.5), std::sqrt(0.5)}, {std::sqrt(0.5), -std::sqrt(0.5)}})},

    {"NOT", CMatrix({{std::sqrt(0.f), std::sqrt(1.f)}, {std::sqrt(1.f), -std::sqrt(0.f)}})}
};

Operation::Operation(const std::string& label, const std::vector<int>& qubits, const std::vector<int*>& cbits): 
    label(label), qubits(qubits), cbits(cbits), next({}), dependencies(0), id(-1) {}

QuantumCircuit::QuantumCircuit(int numQubits, int num_cbits): numQubits(numQubits){
    std::srand(std::time(nullptr));
    this->idCounter=0;

    std::vector<int> cr(num_cbits, 0);
    this->classicalRegister = cr;

    VCD ket0(1<<numQubits, 0.f);
    ket0[0] = 1.f;
    this->totalStatevector = ket0;
}

void QuantumCircuit::h(int q){
    assert(q<this->numQubits);

    std::shared_ptr<Operation> h = std::make_shared<singleQubitGate>("h", q, this->opMatrices["hadamard"]);

    h->id = this->idCounter++;
    this->operations.push_back(h);
}
void QuantumCircuit::x(int q){
    assert(q<this->numQubits);

    std::shared_ptr<Operation> x = std::make_shared<singleQubitGate>("x", q, this->opMatrices["NOT"]);
    x->id = this->idCounter++;
    this->operations.push_back(x);
}

void QuantumCircuit::cx(int qControl, int qTarget){
    assert(qControl<this->numQubits && qTarget<this->numQubits);

    std::shared_ptr<Operation> cx = std::make_shared<ControlledGate>("cx", qControl, qTarget, this->opMatrices["NOT"]);
    cx->id = this->idCounter++;
    this->operations.push_back(cx);
}

void QuantumCircuit::phase(double lambda, int q){
    assert(q<this->numQubits);

    using namespace std::complex_literals;
    const CMatrix pMatrix({{1.f, 0.f},
                           {0.f, std::exp(1i * lambda)}});
    std::shared_ptr<Operation> p = std::make_shared<singleQubitGate>("p", q, pMatrix);
    p->id = this->idCounter++;
    this->operations.push_back(p);
}

void QuantumCircuit::cPhase(double lambda, int qControl, int qTarget){
    assert(qControl<this->numQubits && qTarget<this->numQubits);

    using namespace std::complex_literals;
    const CMatrix cpMatrix({{1.f, 0.f},
                            {0.f, std::exp(1i * lambda)}});
    std::shared_ptr<Operation> cp = std::make_shared<ControlledGate>("cp", qControl, qTarget, cpMatrix);
    cp->id = this->idCounter++;
    this->operations.push_back(cp);
}

void QuantumCircuit::swap(int q1, int q2){

}

void QuantumCircuit::measure(int q, int c){
    assert(q<this->numQubits);
    assert(c<(int)this->classicalRegister.size());

    std::shared_ptr<Operation> measure = std::make_shared<Measure>(q,&(this->classicalRegister[c]));
    measure->id = this->idCounter++;
    this->operations.push_back(measure);
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
        printf("(%s%d, [", ops->at(i)->label.c_str(), ops->at(i)->id);

        for(auto y: ops->at(i)->qubits)printf("%d ", y);
        printf("],[");

        for(auto y: ops->at(i)->cbits)printf("%lld ", (y - &this->classicalRegister[0]));

        printf("])");
        printf(i%5!=4?"  ":"\n             ");
    }

    printf("\n%-13s", "sortedDag: ");

    for(unsigned i=0;i<this->sortedDag.size();i++){
        auto op = this->sortedDag[i];
        printf("(%s%d, [", op->label.c_str(), op->id);

        for(auto y: op->qubits)printf("%d ", y);
        printf("],[");

        for(auto y: op->cbits)printf("%lld ", (y - &this->classicalRegister[0]));
        printf("], %d)", op->dependencies);

        printf(i%5!=4?" > ":" >\n             ");
    }
    printf("\n");
}

void QuantumCircuit::draw(){
    if (this->sortedDag.size() < this->operations.size())
        this->constructDag();

    std::vector<std::vector<char> > blocks;
    std::vector<unsigned> blockSizes(2 * this->numQubits, 1);

    for (int i = 0;i < this->numQubits;i++) {
        std::vector<char> v(this->operations.size() * 5 + 4, ' ');
        v[0] = '-'; v[1] = '-';
        blocks.push_back(v);
        std::vector<char> w(this->operations.size() * 5 + 4, ' ');
        blocks.push_back(w);
    }
    for (unsigned i = 0;i < this->sortedDag.size();i++) {
        auto& op = this->sortedDag[i];
        auto& copyQubits = op->qubits;

        //single-qubit gates
        if (copyQubits.size() == 1) {
            for (char c : {'[', (char)(op->label[0] - 32), ']', '-', '-'}) {
                blocks[2 * copyQubits[0]][++blockSizes[2 * copyQubits[0]]] = c;
            }
        }
        //controlled gates
        else if (op->label[0] == 'c') {
            int minQ = 2 * std::min(copyQubits[0], copyQubits[1]);
            int diffQ = 2 * std::abs(copyQubits[0] - copyQubits[1]);
            unsigned maxSize = 0;
            for (int j = 0; j < diffQ; j++) {
                maxSize = std::max(maxSize, blockSizes[j]);
            }
            for (int j = 0; j <= diffQ; j++) {
                if (j % 2 == 0) {
                    while (blockSizes[minQ + j] < maxSize) {
                        if (blocks[minQ + j][++blockSizes[minQ + j]] == ' ')
                            blocks[minQ + j][blockSizes[minQ + j]] = '-';
                    }
                }
                blocks[minQ + j][2 + maxSize] = '|';
            }

            for (char c : {'-', '*', '-', '-', '-'}) {
                blocks[2 * copyQubits[0]][++blockSizes[2 * copyQubits[0]]] = c;
            }

            if (op->label == "cx") {
                for (char c : { '(', '+', ')', '-', '-' }) {
                    blocks[2 * copyQubits[1]][++blockSizes[2 * copyQubits[1]]] = c;
                }
            }
            if (op->label == "cp") {
                for (char c : {'[', 'P', ']', '-', '-'}) {
                    blocks[2 * copyQubits[1]][++blockSizes[2 * copyQubits[1]]] = c;
                }
            }
        }
    }
    unsigned maxSize = 0;
    for (int j = 0; j < this->numQubits; j++) {
        maxSize = std::max(blockSizes[2*j], maxSize);
    }

    for (int j=0;j<this->numQubits; j++){
        while(blockSizes[2*j] < maxSize){
            if(blocks[2*j][++blockSizes[2*j]] == ' ')
                blocks[2 * j][blockSizes[2 * j]] = '-';
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

void QuantumCircuit::reset(){
    VCD ket0(1<<(this->numQubits), 0.f);
    ket0[0] = 1.f;

    this->totalStatevector = ket0;
}

void QuantumCircuit::constructDag(){
    // Clears the sortedDag
    this->sortedDag.clear();
    
    std::queue<std::shared_ptr<Operation> > queue{};

    // Pushes operations with no dependencies into the queue
    for(auto op : this->operations){
        if(op->dependencies == 0){
            queue.push(op);
        }
    }
    
    // Executes the operation from front of the queue and updates dependencies
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

void QuantumCircuit::run(int shots){
    if(this->sortedDag.size() < this->operations.size())
        this->constructDag();

    for(int shot=0;shot<shots;shot++){
        // Zero out the statevector
        this->reset();

        for(auto op : this->sortedDag){
            op->act(this->totalStatevector);
        }

        //add counts of current c-register state ('bitmask') to counts
        unsigned numCbits = this->classicalRegister.size();
        std::string bitmask(numCbits, '0');
        for(unsigned j = 0;j < numCbits;j++){
            bitmask[numCbits-1 - j] = (char)(this->classicalRegister[j] + '0');
        }

        if(!counts.count(bitmask))
            counts[bitmask] = 1;
        else   
            counts[bitmask]++;
    }
}