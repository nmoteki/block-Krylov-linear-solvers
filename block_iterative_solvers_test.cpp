#include <iostream>
#include <cmath>
#include <string>
#include <random>
#include <chrono>
#include <fstream>
#include <tuple>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "block_iterative_solvers.hpp"
#include "linear_algebra_addon.hpp"
using namespace std;
using namespace Eigen;
using namespace std::chrono;

int main(int argc, char *argv[]){

    double tol= 1e-8;
    int itermax= 1000;

    std::srand((unsigned int) time(0));
    int L=atoi(argv[1]); // dimension of matrix A

    MatrixXcd A;


    switch (argc)
    {
    case 2:
        {
        // load A from a file in matrix market format .mtx
        cout << "please enter filename for A (.mtx): " ;
        string matfname;
        cin >> matfname;
        SparseMatrix<complex<double>> mat;
        //Matrix collections is available in   http://www.cise.ufl.edu/research/sparse/matrices/
        if(!loadMarket(mat,matfname)){
            cerr << "can't open file" << matfname << endl;
            exit(EXIT_FAILURE);
        }
        cout << "A was loaded from " <<  matfname << endl;
        A= MatrixXcd(mat);

        break;
        }
    case 3:
        {
        // Set A as a random matrix
        int n=atoi(argv[2]); // number of RHS vector
        A= MatrixXcd::Random(n,n);
        A= A.transpose()*A; // complex-symmetric
        cout << "A is a " << n << " x " << n << " square random matrix" << endl;
        break;
        }
    default:
        {

        break;
        }
    }

    MatrixXcd B= MatrixXcd::Random(A.rows(),L);


    {
        time_point<steady_clock> start= steady_clock::now();
        cout << "solving AX=B using bl_cocg_rq ..." << endl;
            MatrixXcd X= bl_cocg_rq(A, B, tol, itermax); // only applicable to complex-symmetric matrix A
        cout << " bl_cocg_rq relative error: " << (A*X-B).norm()/B.norm() << endl << endl;
        auto end= steady_clock::now();
        cout << "bl_cocg_rq computation time= " << duration_cast<seconds>((end - start)).count()  << " sec" << endl << endl<< endl;
    }


    {
        time_point<steady_clock> start= steady_clock::now();
        cout << "solving AX=B using bl_bicg_rq ..." << endl;
            MatrixXcd X= bl_bicg_rq(A, B, tol, itermax); // applicable to general matrix A but needs A.adjoint()
        cout << " bl_bicg_rq relative error: " << (A*X-B).norm()/B.norm() << endl << endl;
        auto end= steady_clock::now();
        cout << "bl_bicg_rq computation time= " << duration_cast<seconds>((end - start)).count() << " sec" << endl << endl<< endl;
    }

    {
        time_point<steady_clock> start= steady_clock::now();
        cout << "solving AX=B using bl_bicr_rq ..." << endl;
            MatrixXcd X= bl_bicr_rq(A, B, tol, itermax); // applicable to general matrix A but needs A.adjoint()
        cout << " bl_bicr_rq relative error: " << (A*X-B).norm()/B.norm() << endl << endl;
        auto end= steady_clock::now();
        cout << "bl_bicr_rq computation time= " << duration_cast<seconds>((end - start)).count() << " sec" << endl << endl<< endl;
    }

    {
        time_point<steady_clock> start= steady_clock::now();
        cout << "solving AX=B using bl_bicgstab_rq ..." << endl;
            MatrixXcd X= bl_bicgstab_rq(A, B, tol, itermax); // only applicable to complex-symmetric matrix A
        cout << " bl_bicgstab_rq relative error: " << (A*X-B).norm()/B.norm() << endl << endl;
        auto end= steady_clock::now();
        cout << "bl_bicgstab_rq computation time= " << duration_cast<seconds>((end - start)).count() << " sec" << endl << endl<< endl;
    }

    {
        time_point<steady_clock> start= steady_clock::now();
        cout << "solving AX=B using bl_bicgstab ..." << endl;
            MatrixXcd X= bl_bicgstab(A, B, tol, itermax); // only applicable to complex-symmetric matrix A
        cout << " bl_bicgstab relative error: " << (A*X-B).norm()/B.norm() << endl << endl;
        auto end= steady_clock::now();
        cout << "bl_bicgstab computation time= " << duration_cast<seconds>((end - start)).count() << " sec" << endl << endl<< endl;
    }

    {
        time_point<steady_clock> start= steady_clock::now();
        cout << "solving AX=B using bl_bicggr ..." << endl;
            MatrixXcd X= bl_bicggr(A, B, tol, itermax); // only applicable to complex-symmetric matrix A
        cout << " bl_bicggr relative error: " << (A*X-B).norm()/B.norm() << endl << endl;
        auto end= steady_clock::now();
        cout << "bl_bicggr computation time= " << duration_cast<seconds>((end - start)).count() << " sec" << endl << endl<< endl;
    }








}
