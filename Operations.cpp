#include"Operations.h"

#include<cassert>
#include<cstdlib>
#include<complex>

#include"CMatrix.h"
#include"CMatrix.h"

const std::vector<std::vector<std::complex<double> > > Hadamard::operationMatrix = {
                                                                            {std::sqrt(0.5), std::sqrt(0.5)},
                                                                            {std::sqrt(0.5), -std::sqrt(0.5)}
                                                                            };

Hadamard::Hadamard(int q){
    this->qubits = {q};
    this->name="h";
    this->cbits = {};
}

void Hadamard::act(std::vector<std::complex<double> >& statevector){

    int numQubits = std::log2(statevector.size());

    CMatrix op(this->operationMatrix);
    CMatrix id = identityMatrix(2);

    // Constructs an operator matrix that will act on the total statevector.
    // This matrix is of the form I x I x ... x H x ... x I where H is on the position determined
    // by this->qubits[0] (counting from the back) and x denotes the kronecker product.
    std::vector<CMatrix*> ops;
    for(int i=0;i<numQubits - this->qubits[0] - 1; i++){
        ops.push_back(&id);
    }
    ops.push_back(&op);
    for(int i=0;i<this->qubits[0];i++){
        ops.push_back(&id);
    }

    CMatrix totalOp = kroneckerProduct(ops);

    statevector = matmul(totalOp, statevector);
}



Measure::Measure(int q, int* c){
    this->qubits = {q};
    this->name = "m";
    this->cbits = {c};
}
void Measure::act(std::vector<std::complex<double> >& statevector){
    int q = this->qubits[0];
    double probability[2] = {0.f,0.f};

    // An int of the form 00...1...00 in binary, where 1 is on the qth position.
    // Then, if binary-ANDed with an int corresponding to a basis state of the s.vector,
    // it gives 0 if the basis state included the qth qubit being 0 and >0 otherwise.
    int bs = 1<<q;

    for(unsigned i=0; i<statevector.size(); i++){
        probability[(bs&i) > 0] += std::pow(std::abs(statevector[i]),2);
    }
    double r = (double)std::rand() / RAND_MAX;
    // printf("probability: %.2f, r: %.2f\n", probability[1], r);
    bool result = r <= probability[1];

    for(unsigned i=0; i<statevector.size(); i++){
        // States corresponding to the qubit's state opposite to measured have 0 probability.
        if( (bs&i) != result ){statevector[i]=0; continue;}

        // Remaining states become normalized.
        statevector[i] /= std::sqrt(probability[result]);
    }
    *this->cbits[0] = (int)result;
}



const std::vector<std::vector<std::complex<double> > > CNot::operationMatrix = {
                                                                         {0.f, 1.f},
                                                                         {1.f, 0.f}
                                                                         };
// TODO: fix for non-consecutive qubits
CNot::CNot(int qControl, int qTarget){
    this->qubits = {qControl, qTarget};
    // assert(qTarget-qControl == 1);

    this->name="cx";
    this->cbits = {};
}
void CNot::act(std::vector<std::complex<double> >& statevector){
    int numQubits = std::log2(statevector.size());

    // op = X, and the total operation matrix is written as I x |0><0| + X x |1><1| 
    // (in the case of consecutive qubits)

    // In the general case, totalOp becomes:
    // I x ... x I x |0><0| x ... x I x ... x I    +    I x ... x I x |1><1| x ... x X x ... x I 
    //                 ^            ^                                   ^            ^
    //              qControl     qTarget                             qControl     qTarget
    CMatrix op(this->operationMatrix);
    CMatrix id = identityMatrix(2);

    std::vector<CMatrix*> ops0(numQubits, &id);
    CMatrix ket0bra0({{1.f,0.f},{0.f,0.f}});
    ops0[numQubits-1 - this->qubits[0]] = &ket0bra0;

    std::vector<CMatrix*> ops1(numQubits, &id);
    CMatrix ket1bra1({{0.f,0.f},{0.f,1.f}});
    ops1[numQubits-1 - this->qubits[0]] = &ket1bra1;

    ops1[numQubits-1 - this->qubits[1]] = &op;

    auto opCombined0 = kroneckerProduct(ops0);
    auto opCombined1 = kroneckerProduct(ops1);

    CMatrix totalOp = opCombined0 + opCombined1;
    statevector = matmul(totalOp, statevector);
}


const std::vector<std::vector<std::complex<double> > > Not::operationMatrix = {
                                                                            {0, 1},
                                                                            {1, 0}
                                                                            };
Not::Not(int q){
    this->qubits = {q};
    this->name="x";
    this->cbits = {};
}

void Not::act(std::vector<std::complex<double> >& statevector){
    // For documentation, see Hadamard::act()

    int numQubits = std::log2(statevector.size());

    CMatrix op(this->operationMatrix);
    CMatrix id = identityMatrix(2);

    std::vector<CMatrix*> ops;
    for(int i=0;i<numQubits - this->qubits[0] - 1; i++){
        ops.push_back(&id);
    }
    ops.push_back(&op);
    for(int i=0;i<this->qubits[0];i++){
        ops.push_back(&id);
    }

    CMatrix totalOp = kroneckerProduct(ops);

    statevector = matmul(totalOp, statevector);
}