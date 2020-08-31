#ifndef BASIC_FORMULAS_H
#define BASIC_FORMULAS_H

#include <Eigen/SVD>
#include "Eigen/Core"
#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include <string>
#include <fstream>
#include <cmath>        // std::abs
#include <math.h>       

using namespace Eigen;

/*
Eigen Block of size (p,q), starting at (i,j)	
matrix.block(i,j,p,q);
matrix.block<p,q>(i,j);
*/

using namespace Eigen;

template<typename T>
T simpleThirdOrder(T constTerm, T aConst, T bConst, T cConst, T parametricTerm) {
	return constTerm +
		parametricTerm * aConst +
		parametricTerm * parametricTerm * bConst +
		parametricTerm * parametricTerm * parametricTerm * cConst;
}

template<typename T>
T simpleForthOrder(T constTerm, T aConst, T bConst, T cConst
	, T dConst, T parametricTerm) {
	return constTerm +
		parametricTerm * aConst +
		parametricTerm * parametricTerm * bConst +
		parametricTerm * parametricTerm * parametricTerm * cConst +
		parametricTerm * parametricTerm * parametricTerm 
			* parametricTerm * dConst;
}

template<typename T>
Eigen::Matrix<T, 4, 4> InverseTmat(
	Eigen::Matrix<T, 4, 4> inTmat) {
	Eigen::Matrix<T, 4, 4> TmatInverse 
		= Eigen::Matrix<T, 4, 4>::Identity(4,4);
	Eigen::Matrix<T, 3, 3> rotMat = inTmat.block(0,0,3,3);
	Eigen::Matrix<T, 3, 1> transMat = inTmat.block(0,3,3,1);
	TmatInverse.block(0,0,3,3) = rotMat.transpose();
	TmatInverse.block(0,3,3,1) = -rotMat.transpose() * transMat;
	return TmatInverse;
}

//ensure theta is within [-pi,pi], 
//	&& theta(t) is a monotonic function 
template<typename T>
Eigen::Matrix<T, Dynamic, 1> fixThetas(
	Eigen::Matrix<T, Dynamic, 1> inMat) {
	Eigen::Matrix<T, Dynamic, 1> retVal
		= Eigen::Matrix<T,Dynamic, 1>::
			Zero(inMat.rows(), 1);

	T curTheta = inMat(0,0);
	while(curTheta >= T(-M_PI)) curTheta -= T(M_PI);
	while(curTheta <= T(M_PI)) curTheta += T(M_PI);
	retVal(0,0) = curTheta;

	for(size_t i = 1; i < inMat.rows(); i++) {
		T curTheta = inMat(i,0);
		while(curTheta >= T(-M_PI)) curTheta -= T(2.0 * M_PI);
		while(curTheta <= T(M_PI)) curTheta += T(2.0 * M_PI);

		T testA = curTheta + T(2.0 * M_PI);
		T testB = curTheta - T(2.0 * M_PI);
		
		T distNull = curTheta - retVal(i - 1, 0);
		T distA = testA - retVal(i - 1, 0);
		T distB = testB - retVal(i - 1, 0);

		retVal(i,0) = distNull < distA ? curTheta : testA;
		retVal(i,0) = distNull < distB ? curTheta : testB;
	}

	return retVal;
}


/*
Extracting a rotation can also be done in different ways: MotionBuilder does it by creating a rotated referential from the three markers. First, a vector is defined from the first marker to the second, and a second vector from the first marker to the third. The cross-product of these vectors gives the Z axis of the rotated referential. The sum of these two vectors is perpendicular to the Z axis (by property of cross-product) and is used as the X axis of the rotated referential. The remaining Y axis can easily be computed as the cross-product of Z axis and X axis
*/
template<typename T>
/** Wam base reference sys via computing the 
				bfp of the chest **/
Eigen::Matrix<T, 3, 3> buildRefFramefrom3Pts
	(Eigen::Matrix<T, 3, 1> pt1, Eigen::Matrix<T, 3, 1> pt2
		, Eigen::Matrix<T, 3, 1> pt3) {

	Eigen::Matrix<T, 3, 3> retT;

	/* First, a vector is defined from the first marker to the
		 second, and a second vector from the first marker to the 
			third.*/
	Eigen::Vector3d V1 = (pt1-pt2),
		V2 = (pt3-pt2);
	/* The cross-product of these vectors gives the 
			Y axis of the rotated wam - pointing to back. */
	retT.row(1) = -(((V1).cross(V2)).normalized().transpose());
	
	/* The sum of these two vectors is perpendicular to the
		 Z axis (by property of cross-product) and is used as
		 the X axis of the rotated referential.*/
	retT.row(0) = ((V1+V2).normalized().transpose());

 	/* The remaining Z axis can easily be computed as the
		 cross-product of X axis and Y axis */
	retT.row(2) = ((retT.row(0).transpose()).cross(
		retT.row(1).transpose())).normalized().transpose();

	return retT.transpose();
	}


//(3).	get Euler Parameters from RotMat

// the Euler-angles of the rotation matrix *this using the convention defined by the respective rotation axis as an integer in {2,1,2}. 
/*
template <typename U>
inline void EulerZYZfromRmat	(U* EulAngles 
	,	Eigen::Matrix<U, 3, 3>	Rmat)	{
	// If we choose the value for θ given by Equation (2.29),
	// 	then s θ > 0, and
	U thetaA = ceres::atan2(Rmat(2,2), U(ceres::sqrt(
		1.0 - Rmat(2,2) * Rmat(2,2))));
	U phiA = ceres::atan2(Rmat(0,2), Rmat(1,2));
	U gammaA = ceres::atan2(-Rmat(2,0), Rmat(2,1));
	// If we choose the value for θ given by Equation (2.30), 
	//	then s θ < 0, and
	U thetaB = ceres::atan2(Rmat(2,2), U(-ceres::sqrt(
		1.0 - Rmat(2,2) * Rmat(2,2))));



	EulAngles[1] = ceres::atan2(ceres::sqrt(Rmat(0,2) *
		Rmat(0,2) + Rmat(1,2) * Rmat(1,2)), Rmat(2,2));

	EulAngles[2] = ceres::atan2(Rmat(2,1), -Rmat(2,0));
}
*/

template<typename T>
Eigen::Matrix<T, Dynamic, Dynamic> linearInterp(
	Eigen::Matrix<T, Dynamic, 1> from
	, Eigen::Matrix<T, Dynamic, 1> to
	, size_t nSlices) {
	Eigen::Matrix<T, Dynamic, Dynamic> along = to - from;
	T alongNorm = along.norm();
	Eigen::Matrix<T, Dynamic, Dynamic> alongN = along / alongNorm;
	T sliceSize = alongNorm / T(nSlices);
	Eigen::Matrix<T, Dynamic, Dynamic> slice = sliceSize * alongN;

	Eigen::Matrix<T, Dynamic, Dynamic> outPts(from.rows(), nSlices + 1);
	for(size_t i = 0; i <= nSlices; i++) {
		outPts.block(0, i, from.rows(), 1) = slice * T(i) + from; 
	}

	return outPts;
}


template<typename T>
Eigen::Matrix<T, Dynamic, Dynamic> EllipticalInterp(
	Eigen::Matrix<T, Dynamic, 1> from
	, Eigen::Matrix<T, Dynamic, 1> to
	, size_t nSlices) {

	Eigen::Matrix<T, Dynamic, Dynamic> vect = to - from;
	T vectNorm = vect.norm();
	Eigen::Matrix<T, Dynamic, Dynamic> vectN = vect / vectNorm;
	T sliceSize = M_PI / T(nSlices);
	
	T a = vectNorm / T(2.0);
	T b = vectNorm / T(5.0);
	Eigen::Matrix<T, Dynamic, Dynamic> 
		thetaVect(1, nSlices + 1);
	
	// from collected data people travel from 1/2pi to 3/4pi
	// with an overall deltaTheta of pi
	for(size_t i = 0; i <= nSlices; i++)
		thetaVect(0, i) = M_PI/T(2.0) + sliceSize; 		
		
	for(size_t i = 0; i <= nSlices; i++)
		thetaVect(1, i); 		
		
	Eigen::Matrix<double,Dynamic,Dynamic> 
		pointOnEllipsePseudo2D(3,nSlices);
	for(size_t i = 0; i <= nSlices; i++) {
		pointOnEllipsePseudo2D(0,0) = a * cos(thetaVect(1, i));
		pointOnEllipsePseudo2D(1,0) = b * sin(thetaVect(1, i));
		pointOnEllipsePseudo2D(2,0) = 1.0;
	}
		
	Eigen::Matrix<T, Dynamic, Dynamic> slice = sliceSize * vectN;

	Eigen::Matrix<T, Dynamic, Dynamic> outPts(from.rows(), nSlices + 1);
	for(size_t i = 0; i <= nSlices; i++) {
		outPts.block(0, i, from.rows(), 1) = slice * T(i) + from; 
	}

	return outPts;
}



#endif
