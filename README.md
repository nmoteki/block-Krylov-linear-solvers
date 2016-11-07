##  "block-Krylov linear solvers" -- C++ library for solving a systems of linear equations AX=B using the block-Krylov type methods
---

### General Descriptions:
  A collection of Block-Krylov type solvers for a system of linear equations with multiple right-hand sides AX=B. Matrix elements are assumed to be complex number of double precision. Mathematical algorithms are implemented by C++ language using Eigen. The list of solvers is:
  - **bl_cocg_rq** : Block Condugate Orthogonal Conjugate Gradient (COCG) method with residual orthonormalization [Gu et al. 2016]. Only for  complex-symmetric matricies (i.e., A == A^T). Only one Matrix-Vector product per iteration.
  - **bl_bicg_rq** : Block Bi-Conjugate Grandient (BiCG) method with residual orthonormalization [Rashedi et al. 2016]. Applicable to any non-singular matricies. The Hermitian of A is also used in the solver.
  - **bl_bicr_rq** : Block Bi-Conjugate Residual (BiCR) method with residual orthonormalization [Tadano et al. 2014]. Applicable to any non-singular matricies. The Hermitian of A is also used in the solver.
  - **bl_bicgstab_rq** : Block Bi-Conjugate Gradient stabilzied (BiCGSTAB) method with residual orthonormalization [Rashedi et al. 2016;  Personal communication with Dr. Andreas Frommer]. Applicable to any non-singular matricies.
  - **bl_bicgstab** : Block Bi-Conjugate Gradient stabilzied (BiCGSTAB) method [Guennouni et al. 2003; Tadano et al. 2009]. Applicable to any non-singular matricies.
  - **bl_bicggr** : Block Bi-Conjugate Gradient Gap-Reducing (BiCGGR) method [Tadano et al. 2009]. Applicable to any non-singular matricies.


### Prerequisites:
  1. A C++ compiler supporting C++11. The author tested the code using the gcc version 6.2.0.
  2. A C++ library for linear algebra "Eigen" with version 3.2 or newer. The author tested the code using the Eigen version 3.2.10.

### Test run:
  1. Modify the compilation parameters in "Makefile" according to your environment.
  In particular,please change the "CXX" (= C++ complier) and "INCLUDES" (= absolute path of Eigen folder).

  2. Compilation and linking are performed by

            make

  3. A test run of the solvers for a system AX=B is performed by

          ./block_iterative_solvers_test _number_of_RHS_columns_ [_number_of_rows_]

    If you skip the *number_of_rows*, you will be prompted to enter a filename of matrix A. If *number_of_rows* is entered, a random square matrix of size (*number_of_rows*, *number_of_rows*) is assumed for A. The matrix B is set to be a random rectangular matrix of size (*number_of_rows*, *number_of_RHS_columns*).  


### Interface of solver:
Function interface is common to all the solvers.

    Eigen::MatrixXcd _solver_name_(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B, const double& tol, const int& itermax);

The "Eigen::MatrixXcd" is a Eigen class of a dynamically-allocated matrix with type std::complex<double>. Each solver returns a solution matrix X of size (*number_of_rows*, *number_of_RHS_columns*) in "Eigen::MatrixXcd" format.

### Tolerance:
Every solver internally assures ||AX-B|| / ||B|| < 10*tol* , before returning the numerical solution X. Here, the "|| ||" meant Frobenius norm.

### Maximum number of iterations:
The *itermax* sets the maximum number of iterations. If relative error doesn't decrease below *tol* until the maximum number of iteration, solver is aborted with an error message "*solver_name* did not converge to solution within error tolerance !".

### Theory
  For mathematical backgrounds, please see the references sited in the comment lines of each .cpp file.

### Contact
nobuhiro.moteki@gmail.com

### Licence
Please see LICENSE.txt
