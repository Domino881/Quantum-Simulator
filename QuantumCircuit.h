#ifndef QUANTUMCIRCUIT_H
#define QUANTUMCIRCUIT_H

#include"CMatrix.h"
#include"Qubit.h"
#include<complex>
#include<queue>
#include<vector>
#include<memory>

class Operation{
    public:
        Operation();

        /*
        * Virtual - to be overriden in every derived object 
        * (like Hadamard::act)
        */
        virtual void act(std::vector<std::complex<double> >& statevector) = 0;
        virtual void measure(std::vector<std::shared_ptr<Qubit> >& qubits, int& cbit) const {};

        std::vector<int> qubits;
        std::vector<int> cbits;

        std::vector<std::shared_ptr<Operation> > next;
        int dependencies;
        int id;
        char name;
};

struct DagCompare{
    /*
    * @brief Sorts operations topologically
    */
    bool operator()(const std::shared_ptr<Operation> a, const std::shared_ptr<Operation> b) const;
};

class QuantumCircuit{
    /*
    *
    * 
    * */
    public:
        std::vector<int> ClassicalRegister;

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

        void measure(int q, int c);

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
        std::vector<std::shared_ptr<Qubit> > QuantumRegister;

        std::vector<std::shared_ptr<Operation> > operations;
        std::priority_queue<std::shared_ptr<Operation>, std::vector<std::shared_ptr<Operation>>, DagCompare> dag;

        int id_counter;

};

#endif