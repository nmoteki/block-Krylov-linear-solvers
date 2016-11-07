#include <iostream>
#include <tuple>
#include <Eigen/Dense>
#include "linear_algebra_addon.hpp"
using namespace Eigen;
using namespace std;

MatrixXcd bl_bicg_rq(const MatrixXcd& A, const MatrixXcd& B, const double& tol, const int& itermax)
{
    // Rashedi et al 2016, On short recurrence Krylov type methods for linear systems with many right-hand sides
  double Bnorm= B.norm();
  MatrixXcd X= MatrixXcd::Zero(B.rows(),B.cols()); // Initial guess of X (zeros)
  MatrixXcd R= B-A*X;
  MatrixXcd R_hat= R; // or Rhat= R.conjugate();
  MatrixXcd Q;
  MatrixXcd C;
  tie(Q,C)= qr_reduced(R);
  MatrixXcd Q_hat;
  MatrixXcd C_hat;
  tie(Q_hat,C_hat)= qr_reduced(R_hat);
  MatrixXcd V= Q;
  MatrixXcd V_hat= Q_hat;

  for(int k= 0; k < itermax; ++k){
      MatrixXcd W= A*V;
      MatrixXcd W_hat= A.adjoint()*V_hat;

      MatrixXcd alpha= (V_hat.adjoint()*W).fullPivLu().solve(Q_hat.adjoint()*Q);
      MatrixXcd alpha_hat= (V.adjoint()*W_hat).fullPivLu().solve(Q.adjoint()*Q_hat);

      X= X+V*alpha*C;

      MatrixXcd Qnew;
      MatrixXcd S;
      tie(Qnew,S)= qr_reduced(Q-W*alpha);
      C= S*C;

      MatrixXcd Qnew_hat;
      MatrixXcd S_hat;
      tie(Qnew_hat,S_hat)= qr_reduced(Q_hat-W_hat*alpha_hat);
      C_hat= S_hat*C_hat;

      double err= C.norm()/Bnorm;
      cout << "bl_bicg_rq: " << "iter= " << k << " relative err= " << err << endl;
      if(err < tol) break;

      MatrixXcd beta= (Q_hat.adjoint()*Q).fullPivLu().solve(S_hat.adjoint()*Qnew_hat.adjoint()*Qnew);
      MatrixXcd beta_hat= (Q.adjoint()*Q_hat).fullPivLu().solve(S.adjoint()*Qnew.adjoint()*Qnew_hat);

      Q= Qnew;
      Q_hat= Qnew_hat;

      V= Q+V*beta;
      V_hat= Q_hat+V_hat*beta_hat;
  }

  if((A*X-B).norm()/Bnorm > 10*tol){
      cerr << "bl_bicg_rq did not converge to solution within error tolerance !" << endl;
     // exit(EXIT_FAILURE);
  }

  return X;
}
