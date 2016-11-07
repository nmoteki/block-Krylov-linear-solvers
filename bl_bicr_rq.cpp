#include <iostream>
#include <tuple>
#include <Eigen/Dense>
#include "linear_algebra_addon.hpp"
using namespace Eigen;
using namespace std;

MatrixXcd bl_bicr_rq(const MatrixXcd& A, const MatrixXcd& B, const double& tol, const int& itermax)
{
    // Tadano et al 2014, Improvement of the accuracy of the approximate solution of the block BiCR method
  double Bnorm= B.norm();
  MatrixXcd X= MatrixXcd::Zero(B.rows(),B.cols()); // Initial guess of X (zeros)
  MatrixXcd R= B-A*X;
  MatrixXcd Q;
  MatrixXcd xi;
  tie(Q,xi)= qr_reduced(R);

  MatrixXcd R_til= R; // or R_til= R.conjugate();
  MatrixXcd Q_til;
  MatrixXcd xi_til;
  tie(Q_til,xi_til)= qr_reduced(R_til);

  MatrixXcd S= Q;
  MatrixXcd S_til= Q_til;
  MatrixXcd U= A*Q;
  MatrixXcd U_til= A.adjoint()*Q_til;
  MatrixXcd V_til= U_til;

  for(int k= 0; k < itermax; ++k){

      MatrixXcd alpha= (U_til.adjoint()*U).fullPivLu().solve(V_til.adjoint()*Q);
      MatrixXcd alpha_til= (U.adjoint()*U_til).fullPivLu().solve(Q.adjoint()*V_til);

      X= X+S*alpha*xi;

      MatrixXcd Qnew;
      MatrixXcd tau;
      tie(Qnew,tau)= qr_reduced(Q-U*alpha);
      MatrixXcd Qnew_til;
      MatrixXcd tau_til;
      tie(Qnew_til,tau_til)= qr_reduced(Q_til-U_til*alpha_til);

      xi= tau*xi;
      MatrixXcd Vnew_til= A.adjoint()*Qnew_til;

      double err= xi.norm()/Bnorm;
      cout << "bl_bicr_rq: " << "iter= " << k << " relative err= " << err << endl;
      if(err < tol) break;

      MatrixXcd beta= (V_til.adjoint()*Q).fullPivLu().solve(tau_til.adjoint()*Vnew_til.adjoint()*Qnew);
      MatrixXcd beta_til= (Q.adjoint()*V_til).fullPivLu().solve(tau.adjoint()*Qnew.adjoint()*Vnew_til);

      Q= Qnew;
      Q_til= Qnew_til;
      V_til= Vnew_til;
      S= Q+S*beta;
      S_til= Q_til+S_til*beta_til;
      U_til= V_til+U_til*beta_til;
      U= A*S;
  }

  if((A*X-B).norm()/Bnorm > 10*tol){
      cerr << "bl_bicr_rq did not converge to solution within error tolerance !" << endl;
     // exit(EXIT_FAILURE);
  }

  return X;
}
