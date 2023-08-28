#ifndef CMATRIX_H
#define CMATRIX_H

#include<complex>
#include<iostream>
#include<vector>
#include<cassert>

class CMatrix{
    private:
        std::vector< std::vector<std::complex<double> > > matrix;

    public:
        /*
        * @brief Constructs a complex square matrix, filled with zeros
        * @param dim dimension of the matrix
        * @returns A zero matrix of dimension dim
        */
        CMatrix(const unsigned dim);

        /*
        * @brief Constructs a CMatrix object from a std::vector matrix 
        *
        */
        CMatrix(const std::vector<std::vector<std::complex<double> > >& m);

        /*
        * @brief Prints the complex matrix in a nice format
        */
        void print() const;

        /*
        * @returns The dimension of the matrix
        */
        unsigned dim() const;

    //--------------OPERATORS--------------
        std::complex<double>& operator()(const int& a, const int& b);
        const std::complex<double> operator()(const int& a, const int& b) const;

        CMatrix operator+(const CMatrix& a) const;
        CMatrix operator-() const;
        CMatrix operator-(const CMatrix& a) const{
            return (*this)+(-a);
        }

        /*
        * @returns The transpose of the matrix
        */
        CMatrix t() const;

        /*
        * @returns A matrix with entries conjugate to the original
        */
        CMatrix compconj();

        /*
        * @returns The dagger of the matrix
        */
        CMatrix dagger() const{
            return (this->t()).compconj();
        };

        /*
        * @returns true if matrix is unitary
        * @param epsilon maximum allowed error from the identity
        */
        bool isUnitary(const double epsilon=1e-5) const;

};

/*
* @brief Returns a diagonal CMatrix object
* @param diagonal std::vector of entries along the diagonal
*/
CMatrix diagMatrix(const std::vector< std::complex<double> >& diagonal);

/*
* @brief Returns the identity matrix of a given size
* @param dim dimension of the identity matrix
*/
CMatrix diagMatrix(const int& dim);

/*
* @brief Performs matrix multiplication
*  @param a,b matrices to be multiplied
* @returns The product of a and b
*/
CMatrix matmul(const CMatrix& a, const CMatrix& b);

/*
* @brief Multiplies the vector b on the left with the matrix a
* @returns The product ab
*/
std::vector<std::complex<double> > matmul(const CMatrix& a, const std::vector<std::complex<double> >& b);

CMatrix kroneckerProduct(const CMatrix& a, const CMatrix& b);
CMatrix kroneckerProduct(const std::vector<CMatrix*>& v);

std::vector<std::complex<double> > kroneckerProduct(const std::vector<std::complex<double> >& a, const std::vector<std::complex<double> >& b);
std::vector<std::complex<double> > kroneckerProduct(const std::vector<std::vector<std::complex<double> > >& v);
#endif