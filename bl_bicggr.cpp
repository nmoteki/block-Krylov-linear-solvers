#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

MatrixXcd bl_bicggr(const MatrixXcd& A, const MatrixXcd& B, const double& tol, const int& itermax)
{
// Block BiCGGR [Tadano etal 2009 JSIAM letters]
  double Bnorm= B.norm();
  MatrixXcd X= MatrixXcd::Zero(B.rows(),B.cols()); // Initial guess of X (zeros)
  MatrixXcd R= B-A*X;
  MatrixXcd P= R;
  MatrixXcd V= A*R;
  MatrixXcd W= V;
  MatrixXcd R0til= R; //MatrixXcd::Random(n,L);
  MatrixXcd R0til_H= R0til.adjoint();
  for(int k= 0; k < itermax; ++k){
    MatrixXcd alfa= (R0til_H*V).fullPivLu().solve(R0til_H*R);
    complex<double> qsi= (W.adjoint()*R).trace()/(W.adjoint()*W).trace();
    MatrixXcd S= P-qsi*V;
    MatrixXcd U= S*alfa;
    MatrixXcd Y= A*U;
    X= X+qsi*R+U;
    MatrixXcd Rnew= R-qsi*W-Y;
    double err= Rnew.norm()/Bnorm;
    cout << "bl_bicggr: " << "iter= " << k << " relative err= " << err << endl;
    if(err < tol) break;
    W= A*Rnew;
    MatrixXcd gamma = (R0til_H*R).fullPivLu().solve(R0til_H*Rnew/qsi);
    R=Rnew;
    P= R+U*gamma;
    V= W+Y*gamma;

  }
  if((A*X-B).norm()/Bnorm > 10*tol){
      cerr << "bl_bicggr did not converge to solution within error tolerance !" << endl;
     // exit(EXIT_FAILURE);
  }

  return X;
}
