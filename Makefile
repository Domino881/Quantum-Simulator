# the compiler: gcc for C program, define as g++ for C++
CC = g++

# windows remove command
RM = del

# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS  = -g -Wall

all: main

main: CMatrix.o QuantumCircuit.o Operations.o Qubit.o
	$(CC) $(CFLAGS) -o main main.cpp CMatrix.o QuantumCircuit.o Operations.o Qubit.o

QuantumCircuit.o: Operations.h QuantumCircuit.h Qubit.h

Operations.o: Operations.h Qubit.h

Qubit.o: CMatrix.h Qubit.h

CMatrix.o: CMatrix.h

clean:
	$(RM) *.o
	$(RM) *.exe