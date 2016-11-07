#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

MatrixXcd bl_bicgstab(const MatrixXcd& A, const MatrixXcd& B, const double& tol, const int& itermax)
{
// Block BiCGSTAB [Tadano etal 2009 JSIAM letters]
  double Bnorm= B.norm();
  MatrixXcd X= MatrixXcd::Zero(B.rows(),B.cols()); // Initial guess of X (zeros)
  MatrixXcd R= B-A*X;
  MatrixXcd P= R;
  MatrixXcd R0til= R; //MatrixXcd::Random(B.rows(),B.cols());
  MatrixXcd R0til_H= R0til.adjoint();
  for(int k= 0; k < itermax; ++k){
      MatrixXcd V= A*P;
      FullPivLU<MatrixXcd> lu(R0til_H*V);
      MatrixXcd alfa= lu.solve(R0til_H*R);
      MatrixXcd T= R-V*alfa;
      MatrixXcd Z= A*T;
      complex<double> qsi= (Z.adjoint()*T).trace()/(Z.adjoint()*Z).trace();
      X= X+P*alfa+qsi*T;
      R= T-qsi*Z;
      double err= R.norm()/Bnorm;
      cout << "bl_bicgstab: " << "iter= " << k << " relative err= " << err << endl;
      if(err < tol) break;
      MatrixXcd beta= lu.solve(-R0til_H*Z);
      P= R+(P-qsi*V)*beta;
  }
  if((A*X-B).norm()/Bnorm > 10*tol){
      cerr << "bl_bicgstab did not converge to solution within error tolerance !" << endl;
     // exit(EXIT_FAILURE);
  }

  return X;
}
