#include <iostream>
#include <tuple>
#include <Eigen/Dense>
#include "linear_algebra_addon.hpp"
using namespace Eigen;
using namespace std;

MatrixXcd bl_bicgstab_rq(const MatrixXcd& A, const MatrixXcd& B, const double& tol, const int& itermax)
{
  //Originally derived based on
  // Rashedi et al. 2016, On short recurrence Krylov type methods for linear systems with many right-hand sides, J. Comput. Appl. Math.
  // Guennouni et al 2003, A block version of BICGSTAB for linear systems with multiple right-hand sides, Electronic Transaction on Numerical Anaysis

  int N= B.rows();
  int L= B.cols();
  MatrixXcd thinQ(MatrixXcd::Identity(N,L));
  double Bnorm= B.norm();
  MatrixXcd X= MatrixXcd::Zero(N,L); // Initial guess of X (zeros)
  MatrixXcd R0= B-A*X;
  MatrixXcd R;
  MatrixXcd C;
  tie(R,C)= qr_reduced(R0);
  MatrixXcd Y0= R0;
  double norm0= R0.norm();
  double normr= norm0;
  MatrixXcd P= R;


  for(int k= 0; k < itermax; ++k){
    MatrixXcd AP= A*P;
    FullPivLU<MatrixXcd> lu(Y0.adjoint()*AP);
    MatrixXcd alpha= lu.solve(Y0.adjoint()*R);
    MatrixXcd S= R-AP*alpha;
    MatrixXcd T= A*S;
    MatrixXcd SC= S*C;
    MatrixXcd TC= T*C;
    complex<double> w= (SC.adjoint()*TC).trace()/(TC.adjoint()*TC).trace();
    X= X+P*alpha*C+w*SC;
    R= S-w*T;
    MatrixXcd U;
    tie(R,U)= qr_reduced(R);
    C= U*C;


    double err= C.norm()/Bnorm;
    cout << "bl_bicgstab_rq: " << "iter= " << k << " relative err= " << err << endl;
    if(err < tol) break;


    MatrixXcd beta= -lu.solve(Y0.adjoint()*T*U.inverse());
    P= R+(P-AP*w)*beta;


  }
  if((A*X-B).norm()/Bnorm > 10*tol){
      cerr << "bl_bicgstab_rq did not converge to solution within error tolerance !" << endl;
     // exit(EXIT_FAILURE);
  }

  return X;
}
