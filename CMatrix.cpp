#include"CMatrix.h"
#include<complex>
#include<vector>
#include<cassert>
#include<memory>

CMatrix::CMatrix(const int dim){
    std::vector<std::complex<double> > row;
    for (int i=0; i<dim; i++) row.push_back(0.0);

    for (int i=0; i<dim; i++) this->matrix.push_back(row);
}
CMatrix::CMatrix(const std::vector<std::vector<std::complex<double> > >& m){
    assert(m.size() == m[0].size());
    auto matrix = std::make_unique<CMatrix>(m.size());
    for(int i=0; i<m.size(); i++){
        for(int j=0; j<m.size(); j++){
            matrix->matrix[i][j] = m[i][j];
        }
    }
    this->matrix = matrix->matrix;
}

void CMatrix::print() const{
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

int CMatrix::dim() const{
    return this->matrix.size();
}

std::complex<double> CMatrix::operator()(const int& a, const int& b){
    return this->matrix[a][b];
}
const std::complex<double> CMatrix::operator()(const int& a, const int& b) const{
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

CMatrix diagMatrix(std::vector< std::complex<double> >& diagonal){
    CMatrix r(diagonal.size());
    for(unsigned i=0;i<diagonal.size();i++){
        r(i,i)=diagonal[i];
    }
    return r;
}
CMatrix diagMatrix(const int& dim){
    std::vector<std::complex<double> > v;
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
std::vector<std::complex<double> > matmul(const CMatrix& a, const std::vector<std::complex<double> >& b){
    assert(a.dim() == b.size());
    std::vector<std::complex<double> > res(b.size(), 0);
    for(int j=0;j<b.size();j++){
        for(int k=0; k<a.dim(); k++){
            res[j] += a(j,k) * b[k];
        }
    }
    return res;
}

CMatrix CMatrix::operator+(const CMatrix& a) const{
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

bool CMatrix::isUnitary(const double epsilon) const{
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