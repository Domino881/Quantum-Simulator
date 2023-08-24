#include<complex>
#include<iostream>
#include<vector>
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
        const int dim();

        complex<double>& operator()(int a, int b);

        /*
        * @brief Returns true if matrix is unitary
        */
        const bool isUnitary();

        /*
        * @brief Returns the transpose of the matrix
        */
        CMatrix t();

        /*
        * @brief Returns a matrix with entries conjugate to the original
        */
        CMatrix compconj();

        /*
        * @brief Returns the dagger of the matrix
        */
        CMatrix dagger();
};