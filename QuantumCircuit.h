#ifndef QUANTUMCIRCUIT_H
#define QUANTUMCIRCUIT_H

#include"CMatrix.h"
#include<complex>
#include<queue>
using namespace std;
using namespace complex_literals;

class Qubit{
    public:
        Qubit(): statevector({1,0}) {};

        /*
        * @brief Uses a pseudo-random number to determine the state of the qubit after measuring.
        *  @returns Either 1 or 0, with probabilities dependent on the qubit's statevector
        * 
        * */
        const bool measure();

        vector<complex<double> > statevector;
};

class Operation{
    public:
        Operation();

        vector<int> qubits;
        vector<int> cbits;

        vector<Operation*> next;
        int dependencies;
        int id;
        char name;
};

struct DagCompare{
    /*
    * @brief Sorts operations topologically
    */
    bool operator()(const Operation* a, const Operation* b) const;
};

class QuantumCircuit{
    /*
    *
    * 
    * */
    public:
        /*
        * @brief Creates a new quantum circuit.
        * @param num_qubits number of qubits in the circuit
        * @param num_cbits number of classical bits in the circuit
        */
        QuantumCircuit(int num_qubits, int num_cbits);

        /*
        * @brief Adds a hadamard gate to the circuit
        * @param q target qubit for the gate
        */
        void h(int q);

        void swap(int q1, int q2);

        void cx(int q_control, int q_target);

        /*
        * brief Prints out all internal variables
        */
        void debug_print() const;

        void constructDag();

    private:
        vector<Qubit* > QuantumRegister;
        vector<double> ClassicalRegister;

        vector<Operation> operations;
        priority_queue<Operation*, vector<Operation*>, DagCompare> dag;

        int id_counter;

};

#endif