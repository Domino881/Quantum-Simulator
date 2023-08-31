#include"CMatrix.h"

#include<cassert>
#include<memory>

CMatrix::CMatrix(const unsigned dim){
    std::vector<std::complex<double> > row;
    for (unsigned i=0; i<dim; i++) row.push_back(0.0);

    for (unsigned i=0; i<dim; i++) this->matrix.push_back(row);
}
CMatrix::CMatrix(const std::vector<std::vector<std::complex<double> > >& m){
    assert(m.size() == m[0].size());
    auto matrix = std::make_unique<CMatrix>(m.size());
    for(unsigned i=0; i<m.size(); i++){
        for(unsigned j=0; j<m.size(); j++){
            matrix->matrix[i][j] = m[i][j];
        }
    }
    this->matrix = matrix->matrix;
}

void CMatrix::print() const{
    unsigned dim = this->dim();

    for(unsigned i=0; i<dim; i++){
        if(i==0)printf("/ ");
        else if(i==dim-1)printf("\\ ");
        else printf("| ");
        
        for(unsigned j=0; j<dim; j++){
            printf("%.2f%+.2fi ", real(this->matrix[i][j]), imag(this->matrix[i][j]));
        }

        if(i==0)printf("\\\n");
        else if(i==dim-1)printf("/\n");
        else printf("|\n");
    }
}

unsigned CMatrix::dim() const{
    return this->matrix.size();
}

std::complex<double>& CMatrix::operator()(const int& i, const int& j){
    return this->matrix[i][j];
}
const std::complex<double> CMatrix::operator()(const int& i, const int& j) const{
    return this->matrix[i][j];
}

CMatrix CMatrix::t() const{
    unsigned dim = this->dim();
    CMatrix transposed(dim);

    for(unsigned i=0; i<dim; i++){
        for(unsigned j=0; j<dim; j++){
            transposed(i,j) = this->operator()(j,i);
        }
    }
    return transposed;
}

CMatrix CMatrix::compconj(){
    CMatrix conjug = *this;
    for(unsigned i=0; i<conjug.dim(); i++){
        for(unsigned j=0; j<conjug.dim(); j++){
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
CMatrix identityMatrix(const int& dim){
    CMatrix id(dim);
    for(int i=0;i<dim;i++)id(i,i)=1.f;
    return id;
}

CMatrix matmul(const CMatrix& a, const CMatrix& b){
    assert(a.dim() == b.dim());

    CMatrix res(a.dim());
    for(unsigned i=0; i<a.dim(); i++){
        for(unsigned j=0; j<a.dim(); j++){
            for(unsigned k=0; k<a.dim(); k++){
                res(i,j) += a(i,k) * b(k,j);
            }
        }
    }
    return res;
}
std::vector<std::complex<double> > matmul(const CMatrix& a, const std::vector<std::complex<double> >& b){
    assert(a.dim() == b.size());

    std::vector<std::complex<double> > res(b.size(), 0);

    for(unsigned j=0;j<b.size();j++){
        for(unsigned k=0; k<a.dim(); k++){
            res[j] += a(j,k) * b[k];
        }
    }
    return res;
}

CMatrix CMatrix::operator+(const CMatrix& m) const{
    CMatrix r = *this;

    for(unsigned i=0; i<r.dim(); i++){
        for(unsigned j=0; j<r.dim(); j++){
            r(i,j) += m(i,j);
        }
    }
    return r;
}

CMatrix CMatrix::operator-() const{
    CMatrix r = *this;

    for(unsigned i=0; i<r.dim(); i++){
        for(unsigned j=0; j<r.dim(); j++){
            r(i,j) = -r(i,j);
        }
    }
    return r;
}

bool CMatrix::isUnitary(const double epsilon) const{
    CMatrix dag = this->dagger();
    CMatrix t = matmul(*this, dag);

    CMatrix res = identityMatrix(t.dim()) - t;

    for(unsigned i=0;i<t.dim();i++){
        for(unsigned j=0;j<t.dim();j++){
            if(abs(res(i,j)) >= epsilon) return false;
        }
    }
    return true;
}

CMatrix CMatrix::pow(int p) const{
    assert(p>0);

    CMatrix res = *this;
    p--;
    while(p){
        if(p%2 == 0){
            res = matmul(res, res);
            p /= 2;
        }
        else{
            res = matmul(res, *this);
            p--;
        }
    }
    return res;
}

CMatrix kroneckerProduct(const CMatrix& a, const CMatrix& b){
    CMatrix res(a.dim() * b.dim());

    for(unsigned i=0;i<res.dim();i++){
        for(unsigned j=0; j<res.dim(); j++){
            res(i,j) = a(i/b.dim(), j/b.dim()) * b(i%b.dim(), j%b.dim());
        }
    }
    return res;
}
CMatrix kroneckerProduct(const std::vector<CMatrix*>& v){
    CMatrix res = *v[0];
    for(unsigned i=1; i<v.size(); i++){
        res = kroneckerProduct(res, *v[i]);
    }
    return res;
}

std::vector<std::complex<double> > kroneckerProduct(const std::vector<std::complex<double> >& a, const std::vector<std::complex<double> >& b){
    std::vector<std::complex<double> > res(a.size() * b.size());
    for(unsigned i=0; i<res.size(); i++){
        res[i] = a[i/b.size()] * b[i%b.size()];
    }
    return res;
}
std::vector<std::complex<double> > kroneckerProduct(const std::vector<std::vector<std::complex<double> > >& v){
    std::vector<std::complex<double> > res=v[0];
    for(unsigned i=1; i<v.size(); i++){
        res = kroneckerProduct(res, v[i]);
    }
    return res;
}