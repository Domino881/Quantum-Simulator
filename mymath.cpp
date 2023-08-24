#include"mymath.h"
#include<complex>
#include<vector>
#include<cassert>
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

const int CMatrix::dim() const{
    return this->matrix.size();
}

complex<double>& CMatrix::operator()(int a, int b){
    return this->matrix[a][b];
}
const complex<double>& CMatrix::operator()(int a, int b) const{
    return this->matrix[a][b];
}

CMatrix CMatrix::t() const{
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

CMatrix diagMatrix(vector< complex<double> >& diagonal){
    CMatrix r(diagonal.size());
    for(unsigned i=0;i<diagonal.size();i++){
        r(i,i)=diagonal[i];
    }
    return r;
}
CMatrix diagMatrix(int dim){
    vector<complex<double> > v;
    for(int i=0;i<dim;i++)v.push_back(1.0);
    return diagMatrix(v);
}

CMatrix matmul(const CMatrix& a, const CMatrix& b){
    assert(a.dim() == b.dim());

    CMatrix res(a.dim());
    for(int i=0; i<a.dim(); i++){
        for(int j=0; j<a.dim(); j++){
            for(int k=0; k<a.dim(); k++){
                res(i,j) += a(i,k) * b(k,j);
            }
        }
    }
    return res;
}

CMatrix CMatrix::operator+(const CMatrix& a){
    CMatrix r = *this;
    for(int i=0; i<r.dim(); i++){
        for(int j=0; j<r.dim(); j++){
            r(i,j) += a(i,j);
        }
    }
    return r;
}

CMatrix CMatrix::operator-() const{
    CMatrix r = *this;
    for(int i=0; i<r.dim(); i++){
        for(int j=0; j<r.dim(); j++){
            r(i,j) = -r(i,j);
        }
    }
    return r;
}

const bool CMatrix::isUnitary(double epsilon) const{
    CMatrix dag = this->dagger();
    CMatrix t = matmul(*this, dag);

    CMatrix res = diagMatrix(t.dim()) - t;
    for(int i=0;i<t.dim();i++){
        for(int j=0;j<t.dim();j++){
            if(abs(res(i,j)) >= epsilon) return false;
        }
    }
    return true;
}