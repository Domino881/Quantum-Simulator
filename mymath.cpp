#include"mymath.h"
#include<complex>
#include<vector>
using namespace std;
using namespace complex_literals;

CMatrix::CMatrix(int dim){
    vector<complex<double> > row;
    for (int i=0; i<dim; i++) row.push_back(0.0);

    for (int i=0; i<dim; i++) this->matrix.push_back(row);
}

void CMatrix::print(){
    int dim = this->matrix.size();
    for(int i=0; i<dim; i++){
        if(i==0)printf("/ ");
        else if(i==dim-1)printf("\\ ");
        else printf("| ");
        
        for(int j=0; j<dim; j++){
            printf("%.2f%+.2fi ", real(this->matrix[i][j]), imag(this->matrix[i][j]));
        }

        if(i==0)printf("\\\n");
        else if(i==dim-1)printf("/\n");
        else printf("|\n");
    }
}

const int CMatrix::dim(){
    return this->matrix.size();
}

complex<double>& CMatrix::operator()(int a, int b){
    return this->matrix[a][b];
}

CMatrix CMatrix::t(){
    int dim = this->dim();
    CMatrix transposed(dim);

    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            transposed(i,j) = this->operator()(j,i);
        }
    }
    return transposed;
}

CMatrix CMatrix::compconj(){
    CMatrix conjug = *this;
    for(int i=0; i<conjug.dim(); i++){
        for(int j=0; j<conjug.dim(); j++){
            conjug(i,j) = conj(conjug(i,j));
        }
    }
    return conjug;
}