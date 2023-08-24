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

main.o: main.cpp mymath.h
	$(CC) $(CFLAGS) -c main.cpp

mymath.o: mymath.cpp mymath.h
	$(CC) $(CFLAGS) -c mymath.cpp

clean:
	$(RM) main