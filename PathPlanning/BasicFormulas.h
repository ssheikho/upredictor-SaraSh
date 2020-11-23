#ifndef BASIC_FORMULAS_H
#define BASIC_FORMULAS_H


#include "UBCUtil.h"

//#include "LinearAlgebra.h"

#include "ParseMathematica.h"
#include "ParseCSV.h"
#include <string>
#include <iostream>
//#include "FitFunctions.h"
//#include "Ellipse3D.h"
//#include "Plane.h"
//#include "LinearAlgebra.cpp"
//#include "RigidTrans2D.h"

//#include "BaysFitFunctions.h"
//#include "EllipseConicConstraints.h"

#include "DifferentialIKTypes.h"

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
/** chest reference sys via computing the 
				bfp of the chest and according to
				maya local frame of spine **/
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
			Z axis - pointing to back. */
	retT.row(1) = (((V1).cross(V2)).normalized().transpose());
	
	/* The sum of these two vectors is perpendicular to the
		 Z axis (by property of cross-product) and is used as
		 the X axis of the rotated referential - pointing upward.*/
	retT.row(0) = ((V1+V2).normalized().transpose());

 	/* The remaining Y axis can easily be computed as the
		 cross-product Z x X */
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


//ensure theta is within [-pi,pi], 
//	&& theta(t) is a monotonic function 
template<typename T>
Eigen::Matrix<T, Dynamic, Dynamic> fixThetas(Eigen::Matrix<T, Dynamic, Dynamic> a) {
	
	Eigen::Matrix<T, Dynamic, Dynamic> 
		retVal(a.rows(), a.cols());

	T curTheta = a(0,0);
	while(curTheta >= T(-M_PI)) curTheta -= T(M_PI);
	while(curTheta <= T(M_PI)) curTheta += T(M_PI);
	retVal(0,0) = curTheta;

	for(size_t i = 1; i < a.cols(); i++) {
		T curTheta = a(0,i);
		while(curTheta >= T(-M_PI)) curTheta -= T(2.0 * M_PI);
		while(curTheta <= T(M_PI)) curTheta += T(2.0 * M_PI);

		T testA = curTheta + T(2.0 * M_PI);
		T testB = curTheta - T(2.0 * M_PI);
		
		T distNull = curTheta - retVal(0, i - 1);
		T distA = testA - retVal(0, i - 1);
		T distB = testB - retVal(0, i - 1);

		retVal(0, i) = distNull < distA ? curTheta : testA;
		retVal(0, i) = distNull < distB ? curTheta : testB;
	}

	return retVal;
}


template<typename T>
Eigen::Matrix<T, Dynamic, Dynamic> EllipseFit(
	Eigen::Matrix<T, Dynamic, Dynamic> inPts) {
	Eigen::Matrix<T, Dynamic, Dynamic> outPts(
		inPts.rows(), inPts.cols());

	int nPts = inPts.rows();
	/********************************************************
	Eigen::Matrix<T, 3, Dynamic> 
		inPtsRPi = inPts.block(Markers::RPi,0,3,nPts),
		inPtsRTh = inPts.block(Markers::RTh,0,3,nPts),
		inPtsRWr = inPts.block(Markers::RWr,0,3,nPts),
		inPtsRLA = inPts.block(Markers::RLA,0,3,nPts),
		inPtsREl = inPts.block(Markers::REl,0,3,nPts),
		inPtsRUA = inPts.block(Markers::RUA,0,3,nPts),
		inPtsRSh = inPts.block(Markers::RSh,0,3,nPts);

	Eigen::Matrix<T, Dynamic, Dynamic> thetasEf(1,nPts);
	//(0) Shoulder, not modeling as ellipse, just copied 
	//			from human data
	outPts.block(Markers::RSh,0,3,nPts) = inPtsRSh;
	//(a) ShoToUa
	Ellipse3D *e3dUA = new Ellipse3D(cartToHom(inPtsRUA));
	thetasEf = e3dUA->findThetas();
	outPts.block(Markers::RUA,0,3,nPts) = homToCart(
		e3dUA->getPointAtThetasH(fixThetas(thetasEf)));
	//(b) shoToEl
	Ellipse3D e3dEl(cartToHom(inPtsREl));
	thetasEf = e3dEl.findThetas();
	outPts.block(Markers::REl,0,3,nPts) = homToCart(
		e3dEl.getPointAtThetasH(fixThetas(thetasEf)));
	//(c) shoToLa
	Ellipse3D e3dLA(cartToHom(inPtsRLA));
	thetasEf = e3dLA.findThetas();
	outPts.block(Markers::RLA,0,3,nPts) = homToCart(
		e3dLA.getPointAtThetasH(fixThetas(thetasEf)));
	//(d) shoToWr
	Ellipse3D e3dWr(cartToHom(inPtsRWr));
	thetasEf = e3dWr.findThetas();
	outPts.block(Markers::RWr,0,3,nPts) = homToCart(
		e3dWr.getPointAtThetasH(fixThetas(thetasEf)));
	//(e) shoToTh
	Ellipse3D e3dTh(cartToHom(inPtsRTh));
	thetasEf = e3dTh.findThetas();
	outPts.block(Markers::RTh,0,3,nPts) = homToCart(
		e3dTh.getPointAtThetasH(fixThetas(thetasEf)));
	//(e) shoToPi
	Ellipse3D e3dPi(cartToHom(inPtsRPi));
	thetasEf = e3dPi.findThetas();
	outPts.block(Markers::RPi,0,3,nPts) = homToCart(
		e3dPi.getPointAtThetasH(fixThetas(thetasEf)));
*/
	return outPts;
}



#endif
