#include "LinearAlgebra.h"

#include <iostream>

#include "linalg.h"
#include "specialfunctions.h"

using namespace Eigen;
using namespace std;

//n.3
Eigen::MatrixXd centerPts (Eigen::MatrixXd inPts3dRowWise){
	size_t n = inPts3dRowWise.rows();

	MatrixXd retPts3dCentered = MatrixXd::Zero(n,3);
	retPts3dCentered = inPts3dRowWise.rowwise() 
		- inPts3dRowWise.colwise().mean();
	return retPts3dCentered;
}


void centerRetPts (Eigen::MatrixXd &inPts3dRowWise){
	inPts3dRowWise = inPts3dRowWise.rowwise() 
		- inPts3dRowWise.colwise().mean();
}


Eigen::MatrixXd buildTensorProduct
	(Eigen::MatrixXd inA, Eigen::MatrixXd inB) {
	//inA (m.n); inB (p.q)
	int m, n, p, q;
	m = inA.rows();
	n = inA.cols();
	p = inB.rows();
	q = inB.cols();
	
	Eigen::MatrixXd retMat(m*p,n*q);
	for (int rowA = 0; rowA < m; rowA++) {
		for (int colA = 0; colA < n; colA++) {
			retMat.block(rowA*p,colA*q,p,q) 
				= inA(rowA,colA)*inB;
		}
	}
	return retMat;
}
//Eigen: V is a matrix with the 
//eigenvectors as its columns
std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
	solveEigensystem(Eigen::MatrixXd a) {

	EigenSolver<MatrixXd> eigensolver(a);

	//cout << eigensolver.eigenvalues() << endl;
	//cout << eigensolver.eigenvectors() << endl;

	//SelfAdjointEigenSolver<MatrixXd> eigensolver(a);

	std::map<double, MatrixXd> retVal;

	for(size_t i = 0; i < a.rows(); i++)
		retVal.insert(
			pair<double, MatrixXd>(
				eigensolver.eigenvalues()(i).real()
				, eigensolver.eigenvectors().real().block(0, i, a.rows(), 1)));

	MatrixXd eigenValues(a.rows(), 1);
	MatrixXd eigenVectors(a.rows(), a.rows());
	size_t i = 0;

	for(std::map<double,
		Eigen::MatrixXd>::reverse_iterator iter =
			retVal.rbegin(); iter != retVal.rend(); iter++) {
		
		eigenValues(i,0) = iter->first;

		//col-stacking Eig. vectors 
		eigenVectors.block(0,i,a.rows(),1) = iter->second;
		i++;
	}

	return std::pair<MatrixXd, MatrixXd>(eigenValues, eigenVectors);
}

//X is a n by 3 matrix in which  all the 3D coordinates of the points are row-stacked.
// X_mean is the mean coordinate, i.e., the centroid of the points. 
MatrixXd buildCovMatPCA (MatrixXd inPts3dRowWiseCentered){
	
	size_t n = inPts3dRowWiseCentered.rows();
	
	MatrixXd covMat = MatrixXd::Zero(3,3);	
	covMat = inPts3dRowWiseCentered.transpose()
		* inPts3dRowWiseCentered
		/ (n-1.0);

	return covMat;
}

// diagEigMat: the diagonal matrix  with the 3  eigenvalues of covMat
Eigen::MatrixXd getDiagEigMat
	(std::pair<Eigen::MatrixXd, Eigen::MatrixXd> eigenSys) {

	MatrixXd retDiagEigMat = MatrixXd::Identity(3,3);
	for (int i =0; i<3; i++)
		retDiagEigMat(i,i) = eigenSys.first(i,0);
	
	return retDiagEigMat;
}

// OrthEigVMat:orthogonal matrix row-stacking the corresponding eigenvectors
Eigen::MatrixXd getOrthEigVectsMatRowStacked
	(std::pair<Eigen::MatrixXd, Eigen::MatrixXd> eigenSys) {

	return eigenSys.second.transpose();
}

// nx3
Eigen::MatrixXd allignPtsToPCA
	(Eigen::MatrixXd inPts3dRowWiseCentered
		, Eigen::MatrixXd EigVectsMatRowStacked) {
/*
	Eigen::MatrixXd EigVectsMatRowStacked = 
		getOrthEigVectsMatRowStacked(inPts3dRowWiseCentered);
*/
	// nx3
	MatrixXd ptsPCAaligned = 
		MatrixXd::Zero(inPts3dRowWiseCentered.rows(),3);
	ptsPCAaligned = inPts3dRowWiseCentered
		* EigVectsMatRowStacked.transpose();

	return ptsPCAaligned;
}

/** Eqn 23 - Quinn **/
//inPts3dRowWise (nx3) 
Eigen::MatrixXd buildPCAGausNoiseVarEigVects		
	(int nRows
		, std::pair<Eigen::MatrixXd, Eigen::MatrixXd> eigenSys)
{

	MatrixXd retGaussE = MatrixXd::Zero(3,1);

	for (int i = 0; i<3; i++)	{
		retGaussE (i,0) = sqrt(2.0 
			* eigenSys.first(i,0) / double(nRows-2)) 
			* eigenSys.first(2,0) ;
	}	
	return retGaussE;
}

//alglib::ae_int_t df1, alglib::ae_int_t df2,
//Eigen::MatrixXd buildErrorBoundsPCAaligned(Eigen::MatrixXd inPts3dRowWise, double alpha)	{
//inPts3dRowWise (nx3) 
double getInvFdist (int nRows, double alpha){

	alglib::ae_int_t df1 = (alglib::ae_int_t) 2;
	alglib::ae_int_t df2 = 
		(alglib::ae_int_t) nRows-df1;

	double FisherTestDensity = 
		alglib::invfdistribution(df1, df2, alpha);
	
	return FisherTestDensity;

}

//inPts3dRowWise (nx3)
Eigen::MatrixXd buildEigVectsErrs
	(int nRows, Eigen::MatrixXd noiseVar){

	//MatrixXd noiseVar = MatrixXd::Zero(3,1);
	
	double FisherTestDensit = 
		getInvFdist(nRows, 0.05);

	Eigen::MatrixXd retEigVectsErrs = 
		MatrixXd::Zero(3,1);
	retEigVectsErrs = noiseVar * FisherTestDensit;

	return retEigVectsErrs;
}

Eigen::MatrixXd getHypErrShellAxes
	(Eigen::MatrixXd eigVectsErrs
		, Eigen::MatrixXd diagEigValsMat){

	//h
	Eigen::MatrixXd hyperboloidAxes = 
		MatrixXd:: Zero(3,1);
/*	
	//e_λ
	Eigen::MatrixXd eigVectsErrs = 
		MatrixXd::Zero(3,1);
	eigVectsErrs = buildEigVectsErrs(noiseVar);
	//λ
	Eigen::MatrixXd diagEigValsMat = 
		MatrixXd::Zero(3,3);
	diagEigValsMat = getDiagEigMat(inPts3dRowWise);
*/

	hyperboloidAxes(0,0) = 
		diagEigValsMat(0,0) - eigVectsErrs(0,0);
	hyperboloidAxes(1,0) =
		diagEigValsMat(1,1) - eigVectsErrs (1,0);
	hyperboloidAxes(2,0) = 
		diagEigValsMat(2,2) + eigVectsErrs (2,0);
	
	return hyperboloidAxes;
}

Eigen::MatrixXd getFixedLEllipsoidAxes
	(Eigen::MatrixXd HypErrShellAxes){
	
	/***Fixed-Length Normal Vector A.2
		Error space - ellipsoid major axes with a center
    offset √2h3 from the origin along the 3 axis. ***/    
	Eigen::MatrixXd fixedLEllipsoidAxes  = 
		MatrixXd::Zero(3,1);
	
	fixedLEllipsoidAxes(0,0) = HypErrShellAxes(2,0) 
		* HypErrShellAxes(2,0) /HypErrShellAxes(0,0);
	fixedLEllipsoidAxes(1,0) = HypErrShellAxes(2,0) 
		* HypErrShellAxes(2,0) /HypErrShellAxes(1,0);
	fixedLEllipsoidAxes(2,0) = HypErrShellAxes(2,0);
	
	return fixedLEllipsoidAxes;
}

/*inPts3DPCAaligned (nx3)*/
Eigen::MatrixXd getHypErrBounds
	(Eigen::MatrixXd ErrHyperboloidAxes
		, Eigen::MatrixXd inPts3DPCAaligned) {
	
	int nRows = inPts3DPCAaligned.rows();
	MatrixXd retHypErrBounds = 
		MatrixXd::Zero(nRows,2);
/*
	//Hyperblois Axes    
	Eigen::MatrixXd ErrHyperboloidAxes =
		MatrixXd::Zero(3,1);
	ErrHyperboloidAxes = 
		getHypErrShellAxes(inPts3dRowWise);


	//inPts3DPCAaligned (nx3)
	Eigen::MatrixXd inPts3DPCAaligned = 
		allignPtsToPCA(inPts3dRowWise
			, EigVectsMatRowStacked);
*/
	for (int i=0; i<nRows; i++){
		retHypErrBounds(i,0) = ErrHyperboloidAxes(2,0) 
			* sqrt(pow((inPts3DPCAaligned(i,0)
					/ErrHyperboloidAxes(0,0)),2.0) + 1);

		retHypErrBounds(i,1) = ErrHyperboloidAxes(2,0) 
			* sqrt(pow((inPts3DPCAaligned(i,1)
					/ErrHyperboloidAxes(1,0)),2.0) + 1);
	}
	return retHypErrBounds;
}


/** Eqn 35 - Quinn **/
std::pair<double, double> getMaxMinAngularErrs
	(Eigen::MatrixXd hErrAxes){
	
	double minTheta = 2.0 
		* atan(hErrAxes(2,0)/hErrAxes(1,0));
	double maxTheta = 2.0 
		* atan(hErrAxes(2,0)/hErrAxes(0,0));

	return std::pair<double, double>
		(maxTheta, minTheta);

}

