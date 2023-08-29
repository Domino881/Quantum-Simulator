#pragma once

/*
* @file QuantumCircuit.h
* @brief Declares the base Operation class, as well as the QuantumCircuit class with its core functionality
* @author Dominik Kuczynski
*/

#include<ctime>
#include<memory>
#include<vector>
#include<queue>
#include<complex>

class Operation{
    public:
        Operation();

        /*
        * Virtual - to be overriden in every derived object 
        * (like Hadamard::act)
        */
        virtual void act(std::vector<std::complex<double> >& statevector) = 0;

        // The qubits affected by / needed for the operation
        std::vector<int> qubits;

        // The classical bits affected by / needed for the operation
        std::vector<int*> cbits;

        // The operations following this one on this->qubits
        std::vector<std::shared_ptr<Operation> > next;

        // The number of operations that need do be done before this one
        int dependencies;
        int id;
        std::string name;
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
        std::vector<int> classicalRegister;

        /*
        * @brief Creates a new quantum circuit.
        * @param num_qubits number of qubits in the circuit
        * @param num_cbits number of classical bits in the circuit
        */
        QuantumCircuit(int numQubits, int numCbits);

        /*
        * @brief Adds a hadamard gate to the circuit
        * @param q target qubit for the gate
        */
        void h(int q);

        /*
        * @brief Adds a measurement to the circuit
        * @param q qubit to be measured
        * @param c classical bit to store the result
        */
        void measure(int q, int c);

        void cx(int qControl, int qTarget);

        void swap(int q1, int q2);

        /*
        * brief Prints out all internal variables
        */
        void debug_print() const;
        void draw() const;

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
        const int numQubits;
        // The total statevector of the system - represented by the Kronecker (tensor) product of the 
        // states of the qubits.
        std::vector<std::complex<double> > totalStatevector;

        // Chronological list of operations
        std::vector<std::shared_ptr<Operation> > operations;

        // Directed acyclic graph to represent operation dependencies
        std::priority_queue<std::shared_ptr<Operation>, std::vector<std::shared_ptr<Operation>>, DagCompare> dag;

        int idCounter;

};