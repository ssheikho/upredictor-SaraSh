#ifndef IK_FIT_FROM_LA_TO_W_H
#define IK_FIT_FROM_LA_TO_W_H

#include "BasicFormulas.h"
#include "ForwardKin.h"
#include <UBCUtil.h>

#include "ceres/ceres.h"

#include <iostream>
#include <vector>
/**
solve Joint 5 from Wr position,
**/

template <typename T>
struct IKfitFromLaToW {

	IKfitFromLaToW( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedWrPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedThPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> startConfigToLa
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		) :
		_fk(fk)
		, _observedWrPt(observedWrPt)
		, _observedThPt(observedThPt)
		, _origin(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(4,1))
		, _WrtoThtransWam(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(4,1))
		, _startConfigToLa(startConfigToLa)
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		 {

		//marker Elb pos on the wam, move in +z6 direction
		_origin << T(0.0), T(0.0), T(0.0), T(1.0);
		_WrtoThtransWam << T(0.0),  T(0.0),  T(0.06), T(1.0);
	}

	template <typename U>
	bool operator()(const U* const candidateParam, U* residuals) const {
		U *thetas = new U[_fk.getNJoints()];
	
		//calculate thetas 1 to 4 from elb & La
		for(size_t i = 0; i < 4; i++) 
			thetas[i] =  U(_startConfigToLa(i,0));

		//calculate wrist thetas i.e. 5,6,(7?)
		for(size_t i = 0; i < 3; i++) {
			size_t iTimesFour = i * 4;
			thetas[i+4] = 
				simpleThirdOrder(
					candidateParam[iTimesFour]
					, candidateParam[iTimesFour + 1]
					, candidateParam[iTimesFour + 2]
					, candidateParam[iTimesFour + 3]
					, _parametricTerm);
		}
		for(size_t i = 0; i < _fk.getNJoints(); i++)
			_fk.getDH(i).setTheta(thetas[i]);
/*
		//predicting each WRIST marker position wrt the base  
		//of the robot
		// Wr position is in the origin of theta5,
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> predWrPt =
				homToCart(_fk.getMat(_fk.getNJoints()-2) * _origin);

		residuals[0] = U(predWrPt(0,0) - _observedWrPt(0,0));
		residuals[1] = U(predWrPt(1,0) - _observedWrPt(1,0));
		residuals[2] = U(predWrPt(2,0) - _observedWrPt(2,0));

*/
		//Lenght  is always const in +z
		//move up to J6, and the move by 0.06 in Z6 direction
		//Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> predThPt =
		//		homToCart(_fk.getMat(_fk.getNJoints()-1) * _WrtoThtransWam);

		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> predThPt =
				homToCart(_fk.getMat(_fk.getNJoints()) * _origin);

		
			residuals[0] = U(predThPt(0,0) - _observedThPt(0,0));
			residuals[1] = U(predThPt(1,0) - _observedThPt(1,0));
			residuals[2] = U(predThPt(2,0) - _observedThPt(2,0));



		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* CreateThirdOrderSevenParam(
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedWrPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedThPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> startConfigToLa
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		 ) {
		
		return new ceres::NumericDiffCostFunction
			<IKfitFromLaToW, ceres::RIDDERS,3,12>(
			new IKfitFromLaToW (fk, observedWrPt, observedThPt, startConfigToLa
				, parametricTerm, jointMinAngles, jointMaxAngles));
	}

	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_observedWrPt, _observedThPt, _origin
		, _WrtoThtransWam, _startConfigToLa;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
};

#endif 
