#include <iostream>
#include <tuple>
#include <Eigen/Dense>
#include "linear_algebra_addon.hpp"
using namespace Eigen;
using namespace std;

MatrixXcd bl_cocg_rq(const MatrixXcd& A, const MatrixXcd& B, const double& tol, const int& itermax)
{
// only for complex-symmetric matrix
// Gu et al 2016, arXiv, Block variants of COCG and COCR methods
// for solving complex symmetric linear systems with multiple right-hand sides

  double Bnorm= B.norm();
  MatrixXcd X= MatrixXcd::Zero(B.rows(),B.cols()); // Initial guess of X (zeros)

  MatrixXcd Q;
  MatrixXcd xi;
  tie(Q,xi)= qr_reduced(B-A*X);
  MatrixXcd S= Q;

  for(int k= 0; k < itermax; ++k){
      MatrixXcd AS= A*S;
      MatrixXcd alpha= (S.transpose()*AS).fullPivLu().solve(Q.transpose()*Q);
      X= X+S*alpha*xi;
      MatrixXcd Qnew;
      MatrixXcd tau;
      tie(Qnew,tau)= qr_reduced(Q-AS*alpha);
      xi= tau*xi;
      double err= xi.norm()/Bnorm;
      cout << "bl_cocg_rq: " << "iter= " << k << " relative err= " << err << endl;
      if(err < tol) break;
      MatrixXcd beta= (Q.transpose()*Q).fullPivLu().solve(tau.transpose()*Qnew.transpose()*Qnew);
      Q= Qnew;
      S= Q+S*beta;
  }

  if((A*X-B).norm()/Bnorm > 10*tol){
      cerr << "bl_cocg_rq did not converge to solution within error tolerance !" << endl;
     // exit(EXIT_FAILURE);
  }

  return X;
}
