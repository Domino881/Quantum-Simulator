#include<cstdio>
#include<string>
#include<vector>
#include<iostream>
#include<complex>
#include"CMatrix.h"
using namespace std;
using namespace complex_literals;

int main(){

    CMatrix d(3);
    d(0,0) = 3; d(0,1) = 1; d(0,2) = 2;
    d.print(); printf("\n");

    CMatrix e(2);
    e(0,1)=-1; e(1,0)=1;

    printf("d is unit? %d\n", d.isUnitary());
    printf("e is unit? %d\n", e.isUnitary());
    printf("id is unit? %d\n", diagMatrix(4).isUnitary());
}