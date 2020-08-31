#ifndef IK_FIT_FROM_ELB_TO_LA_H
#define IK_FIT_FROM_ELB_TO_LA_H

#include "BasicFormulas.h"
#include "ForwardKin.h"
#include <UBCUtil.h>

#include "IKfitFromElb.h"
#include "ceres/ceres.h"

#include <iostream>
#include <vector>

template <typename T>
struct IKfit {

	IKfit(ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsVec
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		) : 
		_fk(fk)
		, _inPtsVec(inPtsVec)
		 {
	}
	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_inPtsVec, _origin, _observedElPt;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;

};


template <typename T>
struct IKfitFromElbToLa {

	IKfitFromElbToLa( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedLaPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> startConfigEl
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		) :
		_fk(fk)
		, _observedLaPt(observedLaPt)
		, _origin(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(4,1))
		, _startConfigEl(startConfigEl)
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		 {

		//marker Elb pos on the wam, move in -x3 direction
		_origin << T(0.0), T(0.0),  T(0.0), T(1.0);
		//wrist position in J4(i.e. elbow joint)
		//_origin << T(0.0), T(0.0),  T(0.3), T(1.0);
	}

	template <typename U>
	bool operator()(const U* const candidateParam, U* residuals) const {
		U *thetas = new U[5];
	
 /* 
	ceres::Problem problemEl;
	T* startConfigEl = new T[12]; //12

	ceres::CostFunction* cfShtoEl = IKfitFromElb<T>::
		CreateThirdOrderSevenParam (_fk, _observedElPt
		, _parametricTerm, _jointMinAngles, _jointMaxAngles);
	problemEl.AddResidualBlock(cfShtoEl, NULL, startConfigEl);

	size_t ctr = 0;
	size_t polynomialOrder = 3;
	for(size_t joint = 0; joint < 3; joint++) {
		problemEl.SetParameterLowerBound(startConfigEl, ctr
			, _jointMinAngles[joint]);
		problemEl.SetParameterUpperBound(startConfigEl, ctr
			, _jointMaxAngles[joint]);
		ctr++;
		for(size_t term = 0; term < polynomialOrder; term++) {
			problemEl.SetParameterLowerBound(
				startConfigEl, ctr, _jointMinAngles[joint]);
			problemEl.SetParameterUpperBound(
				startConfigEl, ctr, _jointMaxAngles[joint]);
			ctr++;
		}
	}	

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;


	//options.use_nonmonotonic_steps = true;
	//options.check_gradients = true;
	options.max_num_iterations = 2000;
	options.function_tolerance = 1.0 * pow(10.0,-50);
	options.gradient_tolerance = 1.0 * pow(10.0,-50);
	options.parameter_tolerance =  1.0 * pow(10.0,-50);
	ceres::Solve(options, &problemEl, &summary);	

	// Fit nth order polynomial to thetas[0], thetas[1] & thetas[2] 
	//from ELBOW in joint space
	for (size_t j = 0; j < 2; j++){
		size_t jTimesFour = j * 4;
		thetas[j] = simpleThirdOrder(
					startConfigEl[jTimesFour]
					, startConfigEl[jTimesFour + 1]
					, startConfigEl[jTimesFour + 2]
					, startConfigEl[jTimesFour + 3]
					, _parametricTerm);
	}
*/
	thetas[0] =  _startConfigEl(0,0);
	thetas[1] =  _startConfigEl(1,0);
	//calculate thetas 3 AND 4 5 for LA
	for(size_t i = 0; i < 3; i++) {

		size_t iTimesFour = i * 4;
		thetas[i+2] = 
			simpleThirdOrder(
				candidateParam[iTimesFour]
				, candidateParam[iTimesFour + 1]
				, candidateParam[iTimesFour + 2]
				, candidateParam[iTimesFour + 3]
				, _parametricTerm);

		}

		
		//IF J4 GOES OUT OF BOUND, REVERSE J3 at time 0
		// wich would be a0 theta = a0 +a1t+...];
		bool retVal = true;
			for(size_t i = 0; i < 5; i++) {
				retVal =  retVal && (thetas[i] >= _jointMinAngles[i]);
				retVal =  retVal && (thetas[i] <= _jointMaxAngles[i]);
			}
			//if joints within limit, set to wam			
			if (retVal) {
				for(size_t i = 0; i < 5; i++) 
					_fk.getDH(i).setTheta(thetas[i]);	
			}
	
			//predicting each marker position wrt the base  
			//of the robot
			Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> predLaPt =
				homToCart(_fk.getMat(5) * _origin);
	
			residuals[0] = U(predLaPt(0,0) - _observedLaPt(0,0));
			residuals[1] = U(predLaPt(1,0) - _observedLaPt(1,0));
			residuals[2] = U(predLaPt(2,0) - _observedLaPt(2,0));
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* CreateThirdOrderSevenParam(
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedLaPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedElPtb
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		 ) {
		
		return new ceres::NumericDiffCostFunction
			<IKfitFromElbToLa, ceres::RIDDERS,3,12>(
			new IKfitFromElbToLa (fk, observedLaPt, observedElPtb
				, parametricTerm, jointMinAngles, jointMaxAngles));
	}

	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_observedLaPt, _origin, _startConfigEl;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
};

#endif 
