# the compiler: gcc for C program, define as g++ for C++
CC = g++

# windows remove command
RM = del

# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS  = -g -Wall

all: main

main: main.o CMatrix.o QuantumCircuit.o Operations.o
	$(CC) $(CFLAGS) -o main main.o CMatrix.o QuantumCircuit.o Operations.o

CMatrix.o: CMatrix.h

QuantumCircuit.o: Operations.h QuantumCircuit.h

Operations.o: Operations.h
