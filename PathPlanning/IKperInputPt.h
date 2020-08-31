#ifndef IK_PER_INPUT_PT_H
#define IK_PER_INPUT_PT_H

#include "BasicFormulas.h"
#include "ForwardKin.h"
#include <UBCUtil.h>

#include "ceres/ceres.h"

#include <iostream>
#include <vector>

using ceres::AutoDiffCostFunction;
using ceres::DynamicAutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;
using ceres::CauchyLoss;
using namespace Eigen;

template <typename T>
struct IKperInputPt {

	IKperInputPt( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedElPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedWrPt
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		) :
		_fk(fk)
		, _observedElPt(observedElPt)
		, _observedWrPt(observedWrPt)
		, _origin(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(4,1))
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		, _a(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _alpha(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _d(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		 {

		_origin << T(-0.045), T(0.0), T(0.0), T(1.0);
		
		_a << T(0.0), T(0.0), T(0.045), T(-0.045), T(0.0)
				, T(0.0), T(0.0);
		_alpha << T(-M_PI/2.0), T(M_PI/2.0), T(-M_PI/2.0)
					 , T(M_PI/2.0), T(-M_PI/2.0), T(M_PI/2.0), T(0.0);
		_d << T(0.0), T(0.0), T(0.55), T(0.0), T(0.3)
			 , T(0.0), T(0.06);
	}

	template <typename U>
	bool operator()(const U* const thetas, U* residuals) const {
/*
		U *thetas = new U[7];
		//calculate thetas 1 to 7
		for(size_t i = 0; i < 7; i++) {
			size_t iTimesFour = i * 4;
			thetas[i] = candidateParam[iTimesFour] +
				U(_parametricTerm) * candidateParam[iTimesFour + 1] +
				U(_parametricTerm * _parametricTerm) * 
					 candidateParam[iTimesFour + 2] 
				+	U(_parametricTerm * _parametricTerm *
					_parametricTerm) * candidateParam[iTimesFour + 3];
		}
	
*/
		//ELBOW
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> fkMatEl =
			Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
				::Identity(4,4);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			dhTransformEl(4, 4);
		for(size_t i = 0; i < 2; i++) {
			//compute dh matrix (i-1)T(i)
			dhTransformEl(0,0) = U(cos(thetas[i]));
			dhTransformEl(0,1) = U(-cos(_alpha(i,0)) * sin(thetas[i]));
			dhTransformEl(0,2) = U(sin(_alpha(i,0)) * sin(thetas[i]));
			dhTransformEl(0,3) = U(_a(i,0) * cos(thetas[i]));

			dhTransformEl(1,0) = U(sin(thetas[i]));
			dhTransformEl(1,1) = U(cos(_alpha(i,0)) * cos(thetas[i]));
			dhTransformEl(1,2) = U(-cos(thetas[i]) * sin(_alpha(i,0)));
			dhTransformEl(1,3) = U(_a(i,0) * sin(thetas[i]));

			dhTransformEl(2,0) = U(0.0);
			dhTransformEl(2,1) = U(sin(_alpha(i,0)));
			dhTransformEl(2,2) = U(cos(_alpha(i,0)));
			dhTransformEl(2,3) = U(_d(i,0));

			dhTransformEl(3,0) = U(0.0);
			dhTransformEl(3,1) = U(0.0);
			dhTransformEl(3,2) = U(0.0);
			dhTransformEl(3,3) = U(1.0);

			fkMatEl = fkMatEl * dhTransformEl;
		}

		//predicting each ELBOW marker position wrt the base  
		//of the robot
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			 predElPt(3,1);

		predElPt(0,0) = U(fkMatEl(0,0) * _origin(0,0)) + fkMatEl(0,3); 
		predElPt(1,0) = U(fkMatEl(1,0) * _origin(0,0)) + fkMatEl(1,3); 
		predElPt(2,0) = U(fkMatEl(2,0) * _origin(0,0)) + fkMatEl(2,3); 
 
		/** WRIST
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> fkMatW =
			Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
				::Identity(4,4);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			dhTransformW(4, 4);

		for(size_t i = 0; i < 4; i++) {
			//compute dh matrix (i-1)T(i)
			dhTransformW(0,0) =  U(cos(thetas[i]));
			dhTransformW(0,1) =  U(-cos(_alpha(i,0)) * sin(thetas[i]));
			dhTransformW(0,2) =  U(sin(_alpha(i,0)) * sin(thetas[i]));
			dhTransformW(0,3) =  U(_a(i,0) * cos(thetas[i]));

			dhTransformW(1,0) = U(sin(thetas[i]));
			dhTransformW(1,1) = U(cos(_alpha(i,0)) * cos(thetas[i]));
			dhTransformW(1,2) =	U(-cos(thetas[i]) * sin(_alpha(i,0)));
			dhTransformW(1,3) = U(_a(i,0) * sin(thetas[i]));

			dhTransformW(2,0) = U(0.0);
			dhTransformW(2,1) = U(sin(_alpha(i,0)));
			dhTransformW(2,2) = U(cos(_alpha(i,0)));
			dhTransformW(2,3) = U(_d(i,0));

			dhTransformW(3,0) = U(0.0);
			dhTransformW(3,1) = U(0.0);
			dhTransformW(3,2) = U(0.0);
			dhTransformW(3,3) = U(1.0);
		
			fkMatW = fkMatW * dhTransformW;
		}

		//predicting WRIST marker position wrt the base  
		//of the robot
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> predWrPt(3,1);
		predWrPt(0,0) = fkMatW(0,3); 
		predWrPt(1,0) = fkMatW(1,3); 
		predWrPt(2,0) = fkMatW(2,3); 
**/		
		residuals[0] = predElPt(0,0) - _observedElPt(0,0);
		residuals[1] = predElPt(1,0) - _observedElPt(1,0);
		residuals[2] = predElPt(2,0) - _observedElPt(2,0);
/*
		residuals[3] = predWrPt(0,0) - _observedWrPt(0,0);
		residuals[4] = predWrPt(1,0) - _observedWrPt(1,0);
		residuals[5] = predWrPt(2,0) - _observedWrPt(2,0);
*/
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* CreateThirdOrderSevenParam(
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedElPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedWrPt
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		 ) {
 	return (new ceres::AutoDiffCostFunction<IKperInputPt,3, 2>
		( new IKperInputPt (fk, observedElPt, observedWrPt
				, jointMinAngles, jointMaxAngles)));
	}

	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_observedElPt, _observedWrPt, _origin;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_a, _alpha, _d;

};

#endif 
