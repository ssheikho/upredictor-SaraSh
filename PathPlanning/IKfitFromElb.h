#ifndef IK_FIT_FROM_ELB_H
#define IK_FIT_FROM_ELB_H

#include "BasicFormulas.h"
#include "ForwardKin.h"
#include <UBCUtil.h>

#include "ceres/ceres.h"

#include <iostream>
#include <vector>

template <typename T>
struct IKfitFromElb {

	IKfitFromElb( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedElPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		) :
		_fk(fk)
		, _observedElPt(observedElPt)
		, _origin(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(4,1))
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		 {

		//marker Elb pos on the wam, move in -x3 direction
		_origin << T(0.0) /*T(0.045)*/, T(0.0), T(0.0), T(1.0);
	}

	template <typename U>
	bool operator()(const U* const candidateParam, U* residuals) const {
		U *thetas = new U[3];
	
		//calculate thetas 1 to 3
		for(size_t i = 0; i < 3; i++) {
			size_t iTimesFour = i * 4;
			thetas[i] = U(
				simpleThirdOrder(
					candidateParam[iTimesFour]
					, candidateParam[iTimesFour + 1]
					, candidateParam[iTimesFour + 2]
					, candidateParam[iTimesFour + 3]
					, _parametricTerm));
		}

		//if joints within limit, set to wam		
		bool retVal = true;
		for(size_t i = 0; i < 3; i++) {
			retVal =  retVal && (thetas[i] >= _jointMinAngles[i]);
			retVal =  retVal && (thetas[i] <= _jointMaxAngles[i]);
			if (retVal)
				_fk.getDH(i).setTheta(thetas[i]);
		}

		//predicting each marker position wrt the base  
		//of the robot
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> predElPt =
				homToCart(_fk.getMat(3) * _origin);

		residuals[0] = U(predElPt(0,0)) - U(_observedElPt(0,0));
		residuals[1] = U(predElPt(1,0) - _observedElPt(1,0));
		residuals[2] = U(predElPt(2,0) - _observedElPt(2,0));

		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* CreateThirdOrderSevenParam(
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedElPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		 ) {
		
		return new ceres::NumericDiffCostFunction
			<IKfitFromElb, ceres::RIDDERS,3,12>(
			new IKfitFromElb (fk, observedElPt
				, parametricTerm, jointMinAngles, jointMaxAngles));
	}

	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_observedElPt, _origin;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
};

#endif 
