# the compiler: gcc for C program, define as g++ for C++
CC = g++

# windows remove command
RM = del

# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
CPPFLAGS  = -g -Wall

all: main

main: CMatrix.o QuantumCircuit.o Operations.o
	$(CC) $(CPPFLAGS) -o main main.cpp CMatrix.o QuantumCircuit.o Operations.o 

QuantumCircuit.o: Operations.h QuantumCircuit.h 

Operations.o: Operations.h 

CMatrix.o: CMatrix.h

clean:
	$(RM) *.o
	$(RM) *.exe