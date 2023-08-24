# the compiler: gcc for C program, define as g++ for C++
CC = g++

# windows remove command
RM = del

# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS  = -g -Wall

all: main

main: main.o mymath.o
	$(CC) $(CFLAGS) -o main main.o mymath.o

main.o: main.cpp CMatrix.h
	$(CC) $(CFLAGS) -c main.cpp

mymath.o: CMatrix.cpp CMatrix.h
	$(CC) $(CFLAGS) -c CMatrix.cpp

clean:
	$(RM) main