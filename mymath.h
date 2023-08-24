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
        */
        CMatrix(int dim);

        /*
        * @brief Prints the complex matrix in a nice format
        */
        void print();

        /*
        * @brief Returns the dimension of the matrix
        */
        const int dim() const;

        complex<double>& operator()(int a, int b);
        const complex<double>& operator()(int a, int b) const;

        /*
        * @brief Returns true if matrix is unitary
        * @param epsilon maximum allowed error from the identity
        */
        const bool isUnitary(double epsilon=1e-5) const;

        /*
        * @brief Returns the transpose of the matrix
        */
        CMatrix t() const;

        /*
        * @brief Returns a matrix with entries conjugate to the original
        */
        CMatrix compconj();

        /*
        * @brief Returns the dagger of the matrix
        */
        CMatrix dagger() const{
            return (this->t()).compconj();
        };

        CMatrix operator+(const CMatrix& a);
        CMatrix operator-() const;
        CMatrix operator-(const CMatrix& a){
            return (*this)+(-a);
        }
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

CMatrix matmul(CMatrix& a, CMatrix& b);