#ifndef NAIVE_IK_H
#define NAIVE_IK_H

#include "BasicFormulas.h"
#include "ForwardKin.h"
#include <UBCUtil.h>

#include "ceres/ceres.h"

#include <Eigen/SVD>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>       

using namespace ceres;
using namespace Eigen;
using namespace std;

template <typename T>
struct NaiveIK {
	NaiveIK( 
		ForwardKin<T> &fk
		/* 5 input markers*/
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsRUAinBase
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsRELinBase
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsRLAinBase
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsRWinBase
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsRTHinBase
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		) : _fk(fk)
		, _inPtsRUAinBase(inPtsRUAinBase)
		, _inPtsRELinBase(inPtsRELinBase)
		, _inPtsRLAinBase(inPtsRLAinBase)
		, _inPtsRWinBase(inPtsRWinBase)
		, _inPtsRTHinBase(inPtsRTHinBase)
		, _origin(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(4,1))
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		 {
		_origin << T(0.0) /*T(-0.045)*/, T(0.0), T(0.0), T(1.0);


		//marker UA pos on the wam -> J1,J2, move in -x3 direction
		T nJoints = _fk.getNJoints();
		T polynomialOrder = T(3);

		size_t nParams =  nJoints * (polynomialOrder + 1);
		double startConfig[nParams]; 
		//initialize 
		int ctr = 0;
		for(size_t joint = 0; joint < 3; joint++) {
			startConfig[ctr] = 0.0; 
			ctr++;
			for(size_t term = 0; term < polynomialOrder; term++) {
				startConfig[ctr] = 0.0;
				ctr++;
			}
		}

		Problem problem;
		//IKfromUAtoSh::IKfromUAtoShCostFunction *cfUa;// = 
		//	IKfromUAtoSh<double>::CreateThirdOrderSevenParam (
		//	_fk, _inPtsRUAinBase, _parametricTerm
		//	, _jointMinAngles, _jointMaxAngles);
		//problem.AddResidualBlock(cfUa, NULL, startConfig);


	} 
	
	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_inPtsRUAinBase, _inPtsRELinBase, _inPtsRLAinBase
		, _inPtsRWinBase, _inPtsRTHinBase, _origin;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
};

/* 
IK CHAIN 1 
	INPUT1: (observedUaPt)		-> IK1 -> OUTPUT1: J1,J2,J3'	


template <typename T>
struct IKfromUAtoSh {

	IKfromUAtoSh( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedUaPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		) :
		_fk(fk)
		, _observedUaPt(observedUaPt)
		, _origin(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(4,1))
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		 {

		//marker Elb pos on the wam, move in -x3 direction
		_origin << T(0.0), T(0.0), T(0.0), T(1.0);
	}

	template <typename U>
	bool operator()(const U* const candidateParam, U* residuals) const {
		U *thetas = new U[_fk.getNJoints()];
	
		//calculate thetas 1 to 7
		for(size_t i = 0; i < 3; i++) {
			size_t iTimesFour = i * 4;
			thetas[i] = 
				simpleThirdOrder(
					candidateParam[iTimesFour]
					, candidateParam[iTimesFour + 1]
					, candidateParam[iTimesFour + 2]
					, candidateParam[iTimesFour + 3]
					, _parametricTerm);
		}
		thetas[3] = U(0.0);
		for(size_t i = 0; i < 7; i++)
			_fk.getDH(i).setTheta(thetas[i]);

		//predicting each marker position wrt the base  
		//of the robot
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> predElPt =
				homToCart(_fk.getMat(4) * _origin);

		
			residuals[0] = U(predElPt(0,0) - _observedUaPt(0,0));
			residuals[1] = U(predElPt(1,0) - _observedUaPt(1,0));
			residuals[2] = U(predElPt(2,0) - _observedUaPt(2,0));



		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* CreateThirdOrderSevenParam(
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedUaPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		 ) {
		
		return new ceres::NumericDiffCostFunction
			<IKfromUAtoSh, ceres::RIDDERS,3,12>(
			new IKfromUAtoSh (fk, observedUaPt
				, parametricTerm, jointMinAngles, jointMaxAngles));
	}

	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_observedUaPt, _origin;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
};
*/
/*
IK CHAIN 2 IKfitFromElbToLa

	INPUT2: (OUTPUT1(J1,J2,J3')	& observedLaPtInBase & observed)		-> IK2 -> OUTPUT2: J3,J4,J5'
*/
#endif
