#ifndef QUANTUMCIRCUIT_H
#define QUANTUMCIRCUIT_H

#include"CMatrix.h"
#include<complex>
#include<queue>
#include"Qubit.h"
using namespace std;
using namespace complex_literals;

class Operation{
    public:
        Operation();

        /*
        * Virtual - to be overriden in every derived object 
        * (like Hadamard::act)
        */
        virtual void act(vector<complex<double> >& statevector) = 0;

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
        vector<double> ClassicalRegister;

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

        /*
        * @brief Constructs a directed acyclic graph of operations, based on their order
        *
        */
        void constructDag();

        /*
        * @brief Executes the gates in the circuit
        *
        */
        void run();

    private:
        vector<Qubit* > QuantumRegister;

        vector<Operation*> operations;
        priority_queue<Operation*, vector<Operation*>, DagCompare> dag;

        int id_counter;

};

#endif