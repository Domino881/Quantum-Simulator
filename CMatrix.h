#ifndef CMATRIX_H
#define CMATRIX_H

#include<complex>
#include<iostream>
#include<vector>
#include<cassert>
using namespace std;
using namespace complex_literals;

class CMatrix{
    private:
        vector< vector<complex<double> > > matrix;

    public:
        /*
        * @brief Constructs a complex square matrix, filled with zeros
        * @param dim dimension of the matrix
        * @returns A zero matrix of dimension dim
        */
        CMatrix(int dim);

        /*
        * @brief Constructs a CMatrix object from a std::vector matrix 
        *
        */
        CMatrix(const vector<vector<complex<double> > >& m);

        /*
        * @brief Prints the complex matrix in a nice format
        */
        void print();

        /*
        * @returns The dimension of the matrix
        */
        const int dim() const;

    //--------------OPERATORS--------------
        complex<double>& operator()(int a, int b);
        const complex<double>& operator()(int a, int b) const;

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
        const bool isUnitary(double epsilon=1e-5) const;

};

/*
* @brief Returns a diagonal CMatrix object
* @param diagonal vector of entries along the diagonal
*/
CMatrix diagMatrix(vector< complex<double> >& diagonal);

/*
* @brief Returns the identity matrix of a given size
* @param dim dimension of the identity matrix
*/
CMatrix diagMatrix(int dim);

/*
* @brief Performs matrix multiplication
*  @param a,b matrices to be multiplicated
* @returns The product of a and b
*/
CMatrix matmul(const CMatrix& a, const CMatrix& b);

vector<complex<double> > matmul(const CMatrix& a, const vector<complex<double> >& b);

#endif