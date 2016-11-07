Eigen::MatrixXcd bl_cocg_rq(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B, const double& tol, const int& itermax);
Eigen::MatrixXcd bl_bicg_rq(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B, const double& tol, const int& itermax);
Eigen::MatrixXcd bl_bicr_rq(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B, const double& tol, const int& itermax);
Eigen::MatrixXcd bl_bicgstab_rq(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& b, const double& tol, const int& itermax);
Eigen::MatrixXcd bl_bicgstab(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B, const double& tol, const int& itermax);
Eigen::MatrixXcd bl_bicggr(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B, const double& tol, const int& itermax);
