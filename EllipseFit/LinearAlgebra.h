#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include <stdio.h>
#include <sys/time.h>
#include <math.h> 

#include <Eigen/Dense>

#include <map>

Eigen::MatrixXd centerPts 
	(Eigen::MatrixXd inPts3dRowWise);
void centerRetPts 
	(Eigen::MatrixXd &inPts3dRowWise);

Eigen::MatrixXd buildTensorProduct
	(Eigen::MatrixXd inA, Eigen::MatrixXd inB);

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> 
	solveEigensystem(Eigen::MatrixXd a);

Eigen::MatrixXd buildCovMatPCA
	(Eigen::MatrixXd inPts3dRowWise);

// diagEigMat: the diagonal matrix  with the 3  eigenvalues of covMat
Eigen::MatrixXd getDiagEigMat
	(std::pair<Eigen::MatrixXd, Eigen::MatrixXd> eigenSys);

// OrthEigVMat:orthogonal matrix row-stacking the corresponding eigenvectors
Eigen::MatrixXd getOrthEigVectsMatRowStacked
	(std::pair<Eigen::MatrixXd, Eigen::MatrixXd> eigenSys);

Eigen::MatrixXd allignPtsToPCA
	(Eigen::MatrixXd inPts3dRowWise
		, Eigen::MatrixXd EigVectsMatRowStacked);

/** Eqn 23 - Quinn **/
Eigen::MatrixXd buildPCAGausNoiseVarEigVects		
	(int nRows, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
			 eigenSys);

double getInvFdist (int nRows, double alpha);

/** Eqn 18 - Quinn **/
Eigen::MatrixXd buildEigVectsErrs
	(int nRows, Eigen::MatrixXd noiseVar);

//alglib::ae_int_t df1, alglib::ae_int_t df2,
Eigen::MatrixXd getHypErrShellAxes
	(Eigen::MatrixXd eigVectsErrs
		, Eigen::MatrixXd diagEigValsMat);


Eigen::MatrixXd getFixedLEllipsoidAxes
	(Eigen::MatrixXd HypErrShellAxes);

Eigen::MatrixXd getHypErrBounds
	(Eigen::MatrixXd ErrHyperboloidAxes
		, Eigen::MatrixXd inPts3DPCAaligned);

/** Eqn 35 - Quinn **/
std::pair<double, double> getMaxMinAngularErrs
	(Eigen::MatrixXd hErrAxes);


std::pair<Eigen::MatrixXd, Eigen::MatrixXd> solveEigensystem(Eigen::MatrixXd a);
Eigen::MatrixXd buildCovMatPCA(Eigen::MatrixXd inPts3D);

#endif
