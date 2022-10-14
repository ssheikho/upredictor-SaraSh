#ifndef FIT_FUNCTIONS_H
#define FIT_FUNCTIONS_H

#include <Eigen/Core>
//#include <Eigen>

//Points are assumed to be 1 pt per column in Cartesian format 2D
Eigen::MatrixXd buildConicConstraintMat(Eigen::MatrixXd inPts);
Eigen::MatrixXd buildConicConstraintMatHom(Eigen::MatrixXd inPtsHom);
Eigen::MatrixXd arrangeConicVectIntoImplicit(Eigen::MatrixXd inVect);
Eigen::MatrixXd aThirtyThreeFromConicSolution(Eigen::MatrixXd inVect);
Eigen::MatrixXd twoDHomToThreeDHom(Eigen::MatrixXd inMat, double zOff);
Eigen::MatrixXd cart3DRotMatToHomRT(Eigen::MatrixXd inR);


std::pair<Eigen::MatrixXd, Eigen::MatrixXd> KmeansCluster
	(Eigen::MatrixXd Xmat, int k);

#endif
