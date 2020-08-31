#include "Plane.h"
#include "UBCUtil.h"
#include "LinearAlgebra.h"

/* The homogeneous representation of the plane is the
 4-vector π = (π 1 , π 2 , π 3 , π 4 )^T
π 1 X 1 + π 2 X 2 + π 3 X 3 + π 4 X 4 = 0
The first 3 components of π correspond to the plane normal of Euclidean geometry

– using inhomogeneous notation (3.2) becomes the familiar plane equation written in 
3-vector notation as n. X + d = 0, and d = π 4, where 
n = (π 1 , π 2 , π 3 )^T, X = ( X , Y , Z )^T, X 4 = 1

In this form d/||n|| is the distance of the plane from the origin. */

// (nX3) inPtsAlongRows
Plane::Plane(MatrixXd inPts) :
	_plane(svdSolve((inPts).transpose()))
	, _planeNormalVect(_plane.block(0,0,3,1))
	, _norm(_planeNormalVect.norm())
	, _planeN(_plane / _norm)
	//, _planeNormalVectN(_planeN.block(0,0,3,1))
	, _planeOrigin(inPts.rowwise().mean().block(0,0,3,1)) 
{
	// Eigen System Solver //
	Eigen::MatrixXd inPtsAlongRowsC = 
		centerPts((inPts).transpose().leftCols(3));

	MatrixXd covMat = buildCovMatPCA
		(inPtsAlongRowsC);
	_eigenSys = solveEigensystem(covMat); 			

	Eigen::MatrixXd EigVectsMatRowStacked =
		getOrthEigVectsMatRowStacked (_eigenSys);	
	Eigen::MatrixXd diagEigValsMat = getDiagEigMat(_eigenSys);
	
	/*	 The nominal Plane in PCA 	*/	
	_planeZoff = EigVectsMatRowStacked.row(2).dot(_planeOrigin);
	//flip the coordinate system 
	Vector3d zAxis (0.0, 0.0, 1.0) ;  
	if (EigVectsMatRowStacked.row(2).dot(zAxis) < 0.0) {
		EigVectsMatRowStacked.row(0) = 
			- EigVectsMatRowStacked.row(0);
		EigVectsMatRowStacked.row(2) = 
			- EigVectsMatRowStacked.row(2);
	}

	_planeNormalVectN = 
		EigVectsMatRowStacked.row(2).transpose();
	_planeV1N =	EigVectsMatRowStacked.row(0).transpose();
	_planeV2N =	EigVectsMatRowStacked.row(1).transpose();

	// plane angles wrt the global coordinate reference 
	Vector3d crossProduct = _planeNormalVectN.cross(zAxis);
	double dotProduct = _planeNormalVectN.dot(zAxis);
	
	_planeNtoZ =
			atan2(crossProduct.norm(), dotProduct);
	_planeYawZ = atan2(_planeNormalVectN(0,0)
		, _planeNormalVectN(1,0)) - M_PI/2.0;
	//only for phase 1 - quick hack
	while (_planeYawZ < -M_PI) _planeYawZ+= M_PI;

	_planeDipX = acos(_planeNormalVectN(2,0)
		/_planeNormalVectN.norm());
	// rake - only for the purpose of reporting error 
	_planeTiltY = acos(crossProduct.dot(_planeV2N));

}
Plane::~Plane() {}

Vector4d Plane::getPlane() {
	return _plane;
}

Vector3d Plane::getPlaneNormalVect() {
	return _planeNormalVect;
}

double Plane::getNorm() {
	return _norm;
}


Vector4d Plane::getPlaneN() {
	return _planeN;
}

Vector3d Plane::getPlaneNormalVectN() {
	return _planeNormalVectN;
}

Vector3d Plane::getPlaneV1N() {
	return _planeV1N;
}

Vector3d Plane::getPlaneV2N() {
	return _planeV2N;
}

Vector3d Plane::getPlaneOrigin() {
	return _planeOrigin;
}

double Plane::getPlaneZoff() {
	return 	_planeZoff;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
 Plane::getPlaneEigenSys() {
	return _eigenSys;
}

double Plane::getPlaneNtoZ() {
	return _planeNtoZ;
}

double Plane::getPlaneYawZ() {
	return _planeYawZ;
}

double Plane::getPlaneDipX() {
	return _planeDipX;
}

double Plane::getPlaneTiltY() {
	return _planeTiltY;
}




