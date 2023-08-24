#include<cstdio>
#include<string>
#include<vector>
#include<iostream>
#include<complex>
#include"mymath.h"
using namespace std;
using namespace complex_literals;

int main(){

    CMatrix d(3);

    d(0,1) = 3.0 + 1i;
    printf("%.2f\n", real(d(0,1)));
    d.print();
    d.t().compconj().print();

}