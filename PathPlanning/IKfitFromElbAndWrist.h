#ifndef IK_FIT_FROM_ELB_AND_WRIST_H
#define IK_FIT_FROM_ELB_AND_WRIST_H

#include "BasicFormulas.h"
#include "ForwardKin.h"
#include <UBCUtil.h>

#include "ceres/ceres.h"

#include <iostream>
#include <vector>

template <typename T>
struct IKfitFromElbAndWrist {

	IKfitFromElbAndWrist( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedElPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedWrPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		) :
		_fk(fk)
		, _observedElPt(observedElPt)
		, _observedWrPt(observedWrPt)
		, _origin(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(4,1))
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		 {
		//marker Elb pos on the wam, move in -x3 direction
		_origin << -0.045, 0.0, 0.0, 1.0;
	}

	template <typename U>
	bool operator()(const U* const candidateParam, U* residuals) const {
		U *thetas = new U[_fk.getNJoints()];
	
		//calculate thetas 1 to 7
		for(size_t i = 0; i < _fk.getNJoints(); i++) {
			size_t iTimesFour = i * 4;
			thetas[i] = 
				simpleThirdOrder(
					candidateParam[iTimesFour]
					, candidateParam[iTimesFour + 1]
					, candidateParam[iTimesFour + 2]
					, candidateParam[iTimesFour + 3]
					, _parametricTerm);
		}
		for(size_t i = 0; i < _fk.getNJoints(); i++)
			_fk.getDH(i).setTheta(thetas[i]);

		//predicting each marker position wrt the base  
		//of the robot
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> predElPt =
				homToCart(_fk.getMat(_fk.getNJoints()-4) * _origin);
		//Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> predWrPt =
		//		homToCart(_fk.getMat(_fk.getNJoints()-3) * _origin);

		
			residuals[0] = U(100.0) * (predElPt(0,0) - _observedElPt(0,0));
			residuals[1] = U(100.0) * (predElPt(1,0) - _observedElPt(1,0));
			residuals[2] = U(100.0) * (predElPt(2,0) - _observedElPt(2,0));
/*
			residuals[3] = U(120) * (predWrPt(0,0) - _observedWrPt(0,0));
			residuals[4] = U(120) * (predWrPt(1,0) - _observedWrPt(1,0));
			residuals[5] = U(120) * (predWrPt(2,0) - _observedWrPt(2,0));
*/


		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* CreateThirdOrderSevenParam(
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedElPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedWrPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		 ) {
		
		return new ceres::NumericDiffCostFunction
			<IKfitFromElbAndWrist, ceres::RIDDERS,3,28>(
			new IKfitFromElbAndWrist (fk, observedElPt, observedWrPt
				, parametricTerm, jointMinAngles, jointMaxAngles));
	}

	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_observedElPt, _observedWrPt, _origin;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
};

#endif 
