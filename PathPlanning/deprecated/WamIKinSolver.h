									
	/*** Wam IK SOLVER ---
		(2.A.1) Calculate thetas_i from startConfigs and t_i:  
				Ttot(s) = 1(s)/100(f) *nRows(f) 	
				t_slizeSize = nRows(f)
		  	thetas_i = f(startConfig_i) 
					= a0 + a1*t + a2*t^2 + a3*t^3
	**/  
											

#ifndef WAM_IKIN_SOLVER_H
#define WAM_IKIN_SOLVER_H

#include "BasicFormulas.h"
#include "ForwardKin.h"
#include <UBCUtil.h>
#include "BasePoseSolver.h"

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include <Eigen/SVD>
#include <string>
#include <fstream>
#include <math.h>       

#include <iostream>
#include <vector>

using namespace Eigen;
using namespace ceres;
using namespace std;

using ceres::AutoDiffCostFunction;
using ceres::DynamicAutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;
using ceres::CauchyLoss;

template <typename T>
struct WamIKinSolver {

	WamIKinSolver( 
		// input raw human data from Vicon
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsAlongRowsMB
		, ForwardKin<T> &fk
		/* 5 input markers*/
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsRUAinBase
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsRELinBase
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsRLAinBase
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsRWinBase
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsRTHinBase
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &JsWam
		) : _inPtsAlongRowsMB(inPtsAlongRowsMB)
		, _fk(fk)
		, _inPtsRUAinBase(inPtsRUAinBase)
		, _inPtsRELinBase(inPtsRELinBase)
		, _inPtsRLAinBase(inPtsRLAinBase)
		, _inPtsRWinBase(inPtsRWinBase)
		, _inPtsRTHinBase(inPtsRTHinBase)
		, _origin(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(4,1))
		,_JsWam(JsWam) 
		{

		_origin << T(0.0), T(0.0), T(0.0), T(1.0);
		int nRows = inPtsRUAinBase.cols();

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  
			inPtsRW = 0.01 
				* _inPtsAlongRowsMB.leftCols(3).transpose()
			, inPtsREl = 0.01 
				* _inPtsAlongRowsMB.block(0,9,nRows,3).transpose()
			, inPtsRSh = 0.01 
				* _inPtsAlongRowsMB.block(0,15,nRows,3).transpose();

		// sholder-to-elbow & -elbow-to-wrist vectors
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>  
			shoToElVs = (inPtsREl - inPtsRSh).colwise().normalized()
			, elToWrVs = (inPtsRW - inPtsREl).colwise().normalized();


		//marker UA pos on the wam -> J1,J2, move in -x3 direction
		double deltaSize = 1.0/double(nRows);
		T nJoints = _fk.getNJoints();
		T polynomialOrder = T(3);
		size_t nParams =  nJoints * (polynomialOrder + 1);
		double startConfig[nParams] = {T(0.0)}; 
	
		_jointMinAngles.push_back(-2.6);
		_jointMinAngles.push_back(-2.0);
		_jointMinAngles.push_back(-2.8);
		_jointMinAngles.push_back(-0.9);
		_jointMinAngles.push_back(-4.76);
		_jointMinAngles.push_back(-1.6);
		_jointMinAngles.push_back(-3.0);

		_jointMaxAngles.push_back(2.6);
		_jointMaxAngles.push_back(2.0);
		_jointMaxAngles.push_back(2.8);
		_jointMaxAngles.push_back(3.1);
		_jointMaxAngles.push_back(1.24);
		_jointMaxAngles.push_back(1.6);
		_jointMaxAngles.push_back(3.0);
		
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> shoToUa(4,1);
		shoToUa << T(0.0), T(0.0), T(0.55), T(1.0);

		int toDofUa = 2, toDofEl = 3;
		size_t nParamsUa =  toDofUa*(polynomialOrder + 1),
			nParamsEl =  polynomialOrder + 1;
		double startConfigUa[nParamsUa] = {T(0.0)}
			, startConfigEl[nParamsEl] = {T(0.0)}
			, startConfigLa[nParamsEl] = {T(0.0)};

		Problem problemUa, problemEl, problemLa;
		Solver::Summary summaryUa, summaryEl, summaryLa;
		Solver::Options options;
		options.use_nonmonotonic_steps = true;
		//options.check_gradients = true;
		options.max_num_iterations = 2000;
		options.function_tolerance = 1.0 * pow(10.0,-100);
		options.gradient_tolerance = 1.0 * pow(10.0,-100);
		options.parameter_tolerance =  1.0 * pow(10.0,-100);
		options.max_num_consecutive_invalid_steps = 100;
			
		/* IK CHAIN 1 
			INPUT1: (observedUaPt)	-> IK1 
				-> OUTPUT1: J1,J2		*/
		for(size_t i = 0; i < nRows; i++) {
			double multiplier = double(i) * deltaSize; //0.01;
		/*
			//cf - EL	without offset; i.e. Theta_1, Theta_2
			CostFunction *cfUa = IKfromUAtoSh::
				CreateThirdOrderSevenParam (_fk,  _inPtsRUAinBase.col(i)
				, shoToUa, multiplier, _jointMinAngles, _jointMaxAngles);
			problemUa.AddResidualBlock(cfUa, new ceres::CauchyLoss(0.5) 
				, startConfigUa);
	*/
			CostFunction *cfEl = IKfromEltoSh::
				CreateThirdOrderSevenParam (_fk,  _inPtsRUAinBase.col(i)
				, _inPtsRELinBase.col(i), multiplier
				, _jointMinAngles, _jointMaxAngles);
			problemEl.AddResidualBlock(cfEl, new ceres::CauchyLoss(0.5) 
				, startConfigUa, startConfigEl);
		}		

		//set the joint limits into ceres for Theta_1, Theta_2
		int ctr = 0;
		for(size_t joint = 0; joint < toDofUa; joint++) {
			problemEl.SetParameterLowerBound(startConfigUa, ctr
				, _jointMinAngles[joint]);
			problemEl.SetParameterUpperBound(startConfigUa, ctr
				, _jointMaxAngles[joint]);
			ctr++;
			for(size_t term = 0; term < polynomialOrder; term++) {
				problemEl.SetParameterLowerBound(
					startConfigUa, ctr, _jointMinAngles[joint]);
				problemEl.SetParameterUpperBound(
					startConfigUa, ctr, _jointMaxAngles[joint]);
				ctr++;
			}
		}	

		//set the joint limits into ceres for Theta_3 
		ctr = 0;
		problemEl.SetParameterLowerBound(startConfigEl, ctr
			, _jointMinAngles[2]);
		problemEl.SetParameterUpperBound(startConfigEl, ctr
			, _jointMaxAngles[2]);
		ctr++;
		for(size_t term = 0; term < polynomialOrder; term++) {
			problemEl.SetParameterLowerBound(
				startConfigEl, ctr, _jointMinAngles[2]);
			problemEl.SetParameterUpperBound(
				startConfigEl, ctr, _jointMaxAngles[2]);
		}
		ceres::Solve(options, &problemEl, &summaryEl);
		cout << summaryEl.FullReport() << endl;

		// Fit nth order polynomial to THETA1, THETA2 & THETA3
		for(size_t i = 0; i < nRows; i++) {
			double multiplier = double(i) * deltaSize; 
			for(size_t j = 0; j < toDofUa; j++) {
				size_t jTimesFour = j * 4;
				_JsWam(i,j) = simpleThirdOrder(
					startConfigUa[jTimesFour]
					, startConfigUa[jTimesFour + 1]
					, startConfigUa[jTimesFour + 2]
					, startConfigUa[jTimesFour + 3]
					, multiplier);
			
			}	
			_JsWam(i,toDofUa) = simpleThirdOrder(
				startConfigEl[0]
				, startConfigEl[1]
				, startConfigEl[2]
				, startConfigEl[3]
				, multiplier);
		}
		
	/*
	IK CHAIN 1.2 IKfromUaToEl
	
		INPUT2: ((J1,J2)	& observedElPtInBase)		
			-> IK2 -> OUTPUT2: J3 
	

		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> JsUaI
			=	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
				::Zero(1,7);
		for(size_t i = 0; i < nRows; i++) {
			double multiplier = double(i) * deltaSize; //0.01;
		
			//Elbow offset Vector
			Vector3d inElOffsetVinBase 
				= _inPtsRELinBase.col(i) - _inPtsRUAinBase.col(i);
			//inEl = outUA + inElOffset
				_fk.setThetaVect(_JsWam.row(i).transpose());
			Vector3d inPtEl = homToCart(_fk.getMat(toDofUa) * shoToUa)
				+ inElOffsetVinBase;

			JsUaI(0,0) = _JsWam(i,0);
			JsUaI(0,1) = _JsWam(i,1);

			//cf - EL	with offset; i.e. Theta_3
			CostFunction *cfUaToEl = IKfromUaToEl::
				CreateThirdOrderSevenParam (_fk,  JsUaI, inPtsRELinBase.col(i)
				, multiplier, _jointMinAngles, _jointMaxAngles);
			problemEl.AddResidualBlock(cfUaToEl, new ceres::CauchyLoss(0.5) 
				, startConfigEl);
		}		

		//set the joint limits into ceres for Theta_3
		 ctr = 0;
		problemEl.SetParameterLowerBound(startConfigEl, ctr
			, _jointMinAngles[2]);
		problemEl.SetParameterUpperBound(startConfigEl, ctr
			, _jointMaxAngles[2]);
		ctr++;
		for(size_t term = 0; term < polynomialOrder; term++) {
			problemEl.SetParameterLowerBound(
				startConfigEl, ctr, _jointMinAngles[2]);
			problemEl.SetParameterUpperBound(
				startConfigEl, ctr, _jointMaxAngles[2]);
			ctr++;
		}	
		ceres::Solve(options, &problemEl, &summaryEl);
		cout << summaryEl.FullReport() << endl;

		// Fit nth order polynomial to THETA3
		for(size_t i = 0; i < nRows; i++) {
			double multiplier = double(i) * deltaSize; 
			_JsWam(i,toDofUa) = simpleThirdOrder(
				startConfigEl[0]
				, startConfigEl[1]
				, startConfigEl[2]
				, startConfigEl[3]
				, multiplier);
		}

*/


		/* 2.1. find Human elbow angle, i.e.
			INPUT( <)UpperArm.LowerArm)	-> IK3 
				-> OUTPUT(J4) */
		for(size_t i = 0; i < nRows; i++) {
			Vector3d shoToElV = shoToElVs.col(i)
			, elToWrV = elToWrVs.col(i);
		
			double elTheta = ceres::acos(shoToElV.dot(elToWrV));
			bool retVal = true;	
			retVal =  retVal && (elTheta >= _jointMinAngles[3]);
			retVal =  retVal && (elTheta <= _jointMaxAngles[3]);
			if (retVal) 
				_JsWam(i,3) = M_PI - elTheta;
		}

	/*	IK CHAIN 2 IKfromWrToSh
			INPUT2( inPts))
			-> IK2 ->	OUTPUT2(J1,J2 & J3,J5) */
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> JsI
			=	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
				::Zero(1,7);
		for(size_t i = 0; i < nRows; i++) {
			double multiplier = double(i) * deltaSize; //0.01;

			JsI(0,0) = _JsWam(i,0);
			JsI(0,1) = _JsWam(i,1);
			JsI(0,2) = _JsWam(i,2);
			JsI(0,3) = _JsWam(i,3);

			//cf - Wr	with offset; i.e. Theta_3
			CostFunction *cfShToWr = IKfromWrToSh::
				CreateThirdOrderSevenParam (_fk,  JsI
				, inPtsRUAinBase.col(i), inPtsRELinBase.col(i)
				, inPtsRWinBase.col(i), multiplier
				, _jointMinAngles, _jointMaxAngles);
			problemLa.AddResidualBlock(cfShToWr
				, new ceres::CauchyLoss(0.5), startConfigUa
				, startConfigEl, startConfigLa);
		}		


		//set the joint limits into ceres for Theta_1, Theta_2
		ctr = 0;
		for(size_t joint = 0; joint < toDofUa; joint++) {
			problemLa.SetParameterLowerBound(startConfigUa, ctr
				, _jointMinAngles[joint]);
			problemLa.SetParameterUpperBound(startConfigUa, ctr
				, _jointMaxAngles[joint]);
			ctr++;
			for(size_t term = 0; term < polynomialOrder; term++) {
				problemLa.SetParameterLowerBound(
					startConfigUa, ctr, _jointMinAngles[joint]);
				problemLa.SetParameterUpperBound(
					startConfigUa, ctr, _jointMaxAngles[joint]);
				ctr++;
			}
		}	

		//set the joint limits into ceres for Theta_3 
		ctr = 0;
		problemLa.SetParameterLowerBound(startConfigEl, ctr
			, _jointMinAngles[2]);
		problemLa.SetParameterUpperBound(startConfigEl, ctr
			, _jointMaxAngles[2]);
		ctr++;
		for(size_t term = 0; term < polynomialOrder; term++) {
			problemLa.SetParameterLowerBound(
				startConfigEl, ctr, _jointMinAngles[2]);
			problemLa.SetParameterUpperBound(
				startConfigEl, ctr, _jointMaxAngles[2]);
		}

		//set the joint limits into ceres for Theta_5
		ctr = 0;
		problemLa.SetParameterLowerBound(startConfigLa, ctr
			, _jointMinAngles[4]);
		problemLa.SetParameterUpperBound(startConfigLa, ctr
			, _jointMaxAngles[4]);
		ctr++;
		for(size_t term = 0; term < polynomialOrder; term++) {
			problemLa.SetParameterLowerBound(
				startConfigLa, ctr, _jointMinAngles[4]);
			problemLa.SetParameterUpperBound(
				startConfigLa, ctr, _jointMaxAngles[4]);
			ctr++;
		}	
		ceres::Solve(options, &problemLa, &summaryLa);
		cout << summaryLa.FullReport() << endl;

		// Fit nth order polynomial to THETA5
		for(size_t i = 0; i < nRows; i++) {
			double multiplier = double(i) * deltaSize; 

			for(size_t j = 0; j < toDofUa; j++) {
				size_t jTimesFour = j * 4;
				_JsWam(i,j) = simpleThirdOrder(
					startConfigUa[jTimesFour]
					, startConfigUa[jTimesFour + 1]
					, startConfigUa[jTimesFour + 2]
					, startConfigUa[jTimesFour + 3]
					, multiplier);
			}	
			_JsWam(i,toDofUa) = simpleThirdOrder(
				startConfigEl[0]
				, startConfigEl[1]
				, startConfigEl[2]
				, startConfigEl[3]
				, multiplier);
			
			_JsWam(i,4) = simpleThirdOrder(
				startConfigLa[0]
				, startConfigLa[1]
				, startConfigLa[2]
				, startConfigLa[3]
				, multiplier);
		}


	}

/* 
IK CHAIN 1 
	INPUT1: (observedUaPt)		-> IK1 -> OUTPUT1: J1,J2,J3'	
*/

//template <typename T>
struct IKfromUAtoSh {
	IKfromUAtoSh( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedUaPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> origin
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles

		) :
		_fk(fk)
		, _observedUaPt(observedUaPt)
		, _origin(origin)
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		,_a(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _alpha(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _d(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1)) 
	{
		
		_a << T(0.0), T(0.0), T(0.045), T(-0.045), T(0.0)
				, T(0.0), T(0.0);
		_alpha << T(-M_PI/2.0), T(M_PI/2.0), T(-M_PI/2.0)
					 , T(M_PI/2.0), T(-M_PI/2.0), T(M_PI/2.0), T(0.0);
		_d << T(0.0), T(0.0), T(0.55), T(0.0), T(0.3)
			 , T(0.0), T(0.06);
	}

	template <typename U>
	bool operator()(const U* const candidateParamUa, U* residuals) const {
	int toDofUa = 2;
	U *thetas = new U[toDofUa];	
	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
		fkMatEl =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4), dhTransformEl(4, 4);

		//calculate theta1 & theta2
		for(size_t i = 0; i < toDofUa; i++) {
			size_t iTimesFour = i * 4;
	
			thetas[i] = candidateParamUa[iTimesFour] 
				+	U(_parametricTerm) * candidateParamUa[iTimesFour + 1] 
				+	U(_parametricTerm * _parametricTerm) 
					* candidateParamUa[iTimesFour + 2] 
				+	U(_parametricTerm * _parametricTerm * _parametricTerm) 
					* candidateParamUa[iTimesFour + 3];
	
			//if joints within limit, set to wam		
			bool retVal = true;
			retVal =  retVal && (thetas[i] >= _jointMinAngles[i]);
			retVal =  retVal && (thetas[i] <= _jointMaxAngles[i]);

			if (retVal) {
				//compute dh matrix (i-1)T(i)
				dhTransformEl(0,0) = U(ceres::cos(thetas[i]));
				dhTransformEl(0,1) = U(-ceres::cos(_alpha(i,0)) 
					* ceres::sin(thetas[i]));
				dhTransformEl(0,2) = U(ceres::sin(_alpha(i,0)) 
					* ceres::sin(thetas[i]));
				dhTransformEl(0,3) = U(_a(i,0) * ceres::cos(thetas[i]));

				dhTransformEl(1,0) = U(ceres::sin(thetas[i]));
				dhTransformEl(1,1) = U(ceres::cos(_alpha(i,0)) 
					* ceres:: cos(thetas[i]));
				dhTransformEl(1,2) = U(-ceres::cos(thetas[i]) 
					* ceres::sin(_alpha(i,0)));
				dhTransformEl(1,3) = U(_a(i,0) 
					* ceres::sin(thetas[i]));

				dhTransformEl(2,0) = U(0.0);
				dhTransformEl(2,1) = U(ceres::sin(_alpha(i,0)));
				dhTransformEl(2,2) = U(ceres::cos(_alpha(i,0)));
				dhTransformEl(2,3) = U(_d(i,0));

				dhTransformEl(3,0) = U(0.0);
				dhTransformEl(3,1) = U(0.0);
				dhTransformEl(3,2) = U(0.0);
				dhTransformEl(3,3) = U(1.0);
			}
			fkMatEl = fkMatEl * dhTransformEl;
		}

		//predicting each ELBOW-NO_OFFSET marker position wrt the base  
		//of the robot
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> wamEltoUaTrans =
			Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
				::Zero(4,1);
		wamEltoUaTrans << U(0.0), U(0.0), U(0.55), U(1.0);

		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
		 predElPt (3,1);
		predElPt(0,0) = U(fkMatEl(0,0)) * U(wamEltoUaTrans(0,0)) 
			+ U(fkMatEl(0,1)) * U(wamEltoUaTrans(1,0)) 
			+ U(fkMatEl(0,2)) * U(wamEltoUaTrans(2,0)) 
			+ U(fkMatEl(0,3)) * U(wamEltoUaTrans(3,0));//last element is 1 in hom coords
		predElPt(1,0) = U(fkMatEl(1,0)) * U(wamEltoUaTrans(0,0)) 
			+ U(fkMatEl(1,1)) * U(wamEltoUaTrans(1,0)) 
			+ U(fkMatEl(1,2)) * U(wamEltoUaTrans(2,0)) 
			+ U(fkMatEl(1,3)) * U(wamEltoUaTrans(3,0)); 
		predElPt(2,0) = U(fkMatEl(2,0)) * U(wamEltoUaTrans(0,0)) 
			+ U(fkMatEl(2,1)) * U(wamEltoUaTrans(1,0)) 
			+ U(fkMatEl(2,2)) * U(wamEltoUaTrans(2,0)) 
			+ U(fkMatEl(2,3)) * U(wamEltoUaTrans(3,0)); 

		residuals[0] = U(predElPt(0,0)) - U(_observedUaPt(0,0));
		residuals[1] = U(predElPt(1,0)) - U(_observedUaPt(1,0));
		residuals[2] = U(predElPt(2,0)) - U(_observedUaPt(2,0));

		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
		static ceres::CostFunction* CreateThirdOrderSevenParam(
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedElPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> origin
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		 ) {
 	return (new ceres::AutoDiffCostFunction<IKfromUAtoSh,3, 8>
		( new IKfromUAtoSh (fk, observedElPt, origin
				, parametricTerm, jointMinAngles, jointMaxAngles)));
	}

	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_observedUaPt, _origin;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	 	_a, _alpha, _d;
};
/* 
IK CHAIN 1 B
	INPUT1: (observedUaPt & observedElPt)		-> IK1 -> OUTPUT1: J1,J2,J3'	
*/

//template <typename T>
struct IKfromEltoSh {
	IKfromEltoSh( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
			observedUaPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
			observedElPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles

		) :
		_fk(fk)
		, _observedUaPt(observedUaPt)
		, _observedElPt(observedElPt)
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		,_a(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _alpha(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _d(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1)) 
	{
		
		_a << T(0.0), T(0.0), T(0.045), T(-0.045), T(0.0)
				, T(0.0), T(0.0);
		_alpha << T(-M_PI/2.0), T(M_PI/2.0), T(-M_PI/2.0)
					 , T(M_PI/2.0), T(-M_PI/2.0), T(M_PI/2.0), T(0.0);
		_d << T(0.0), T(0.0), T(0.55), T(0.0), T(0.3)
			 , T(0.0), T(0.06);
	}

	template <typename U>
	bool operator()(const U* const candidateParamUa
		, const U* const candidateParamEl, U* residuals) const {
	int toDofUa = 2;
	U *thetas = new U[3];	
	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
		fkMatEl =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4), dhTransformEl(4, 4);

		//calculate theta1 & theta2
		for(size_t i = 0; i < toDofUa; i++) {
			size_t iTimesFour = i * 4;
	
			thetas[i] = candidateParamUa[iTimesFour] 
				+	U(_parametricTerm) * candidateParamUa[iTimesFour + 1] 
				+	U(_parametricTerm * _parametricTerm) 
					* candidateParamUa[iTimesFour + 2] 
				+	U(_parametricTerm * _parametricTerm * _parametricTerm) 
					* candidateParamUa[iTimesFour + 3];
	
			//if joints within limit, set to wam		
			bool retVal = true;
			retVal =  retVal && (thetas[i] >= _jointMinAngles[i]);
			retVal =  retVal && (thetas[i] <= _jointMaxAngles[i]);

			if (retVal) {
				//compute dh matrix (i-1)T(i)
				dhTransformEl(0,0) = U(ceres::cos(thetas[i]));
				dhTransformEl(0,1) = U(-ceres::cos(_alpha(i,0)) 
					* ceres::sin(thetas[i]));
				dhTransformEl(0,2) = U(ceres::sin(_alpha(i,0)) 
					* ceres::sin(thetas[i]));
				dhTransformEl(0,3) = U(_a(i,0) * ceres::cos(thetas[i]));

				dhTransformEl(1,0) = U(ceres::sin(thetas[i]));
				dhTransformEl(1,1) = U(ceres::cos(_alpha(i,0)) 
					* ceres:: cos(thetas[i]));
				dhTransformEl(1,2) = U(-ceres::cos(thetas[i]) 
					* ceres::sin(_alpha(i,0)));
				dhTransformEl(1,3) = U(_a(i,0) 
					* ceres::sin(thetas[i]));

				dhTransformEl(2,0) = U(0.0);
				dhTransformEl(2,1) = U(ceres::sin(_alpha(i,0)));
				dhTransformEl(2,2) = U(ceres::cos(_alpha(i,0)));
				dhTransformEl(2,3) = U(_d(i,0));

				dhTransformEl(3,0) = U(0.0);
				dhTransformEl(3,1) = U(0.0);
				dhTransformEl(3,2) = U(0.0);
				dhTransformEl(3,3) = U(1.0);
			}
			fkMatEl = fkMatEl * dhTransformEl;
		}

		//predicting each ELBOW-NO_OFFSET marker position wrt the base  
		//of the robot
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> wamEltoUaTrans =
			Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
				::Zero(4,1);
		wamEltoUaTrans << U(0.0), U(0.0), U(0.55), U(1.0);

		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
		 predUaPt (3,1);
		predUaPt(0,0) = U(fkMatEl(0,0)) * U(wamEltoUaTrans(0,0)) 
			+ U(fkMatEl(0,1)) * U(wamEltoUaTrans(1,0)) 
			+ U(fkMatEl(0,2)) * U(wamEltoUaTrans(2,0)) 
			+ U(fkMatEl(0,3)) * U(wamEltoUaTrans(3,0));
		predUaPt(1,0) = U(fkMatEl(1,0)) * U(wamEltoUaTrans(0,0)) 
			+ U(fkMatEl(1,1)) * U(wamEltoUaTrans(1,0)) 
			+ U(fkMatEl(1,2)) * U(wamEltoUaTrans(2,0)) 
			+ U(fkMatEl(1,3)) * U(wamEltoUaTrans(3,0)); 
		predUaPt(2,0) = U(fkMatEl(2,0)) * U(wamEltoUaTrans(0,0)) 
			+ U(fkMatEl(2,1)) * U(wamEltoUaTrans(1,0)) 
			+ U(fkMatEl(2,2)) * U(wamEltoUaTrans(2,0)) 
			+ U(fkMatEl(2,3)) * U(wamEltoUaTrans(3,0)); 

		residuals[0] = U(predUaPt(0,0)) - U(_observedUaPt(0,0));
		residuals[1] = U(predUaPt(1,0)) - U(_observedUaPt(1,0));
		residuals[2] = U(predUaPt(2,0)) - U(_observedUaPt(2,0));

		//calculate theta3
			thetas[2] = candidateParamEl[0] 
				+	U(_parametricTerm) * candidateParamEl[1] 
				+	U(_parametricTerm * _parametricTerm) 
					* candidateParamEl[2] 
				+	U(_parametricTerm * _parametricTerm * _parametricTerm) 
					* candidateParamEl[3];
	
			//if joints within limit, set to wam		
			bool retVal = true;
			retVal =  retVal && (thetas[2] >= _jointMinAngles[2]);
			retVal =  retVal && (thetas[2] <= _jointMaxAngles[2]);

			if (retVal) {
				//compute dh matrix (i-1)T(i)
				dhTransformEl(0,0) = U(ceres::cos(thetas[2]));
				dhTransformEl(0,1) = U(-ceres::cos(_alpha(2,0)) 
					* ceres::sin(thetas[2]));
				dhTransformEl(0,2) = U(ceres::sin(_alpha(2,0)) 
					* ceres::sin(thetas[2]));
				dhTransformEl(0,3) = U(_a(2,0) * ceres::cos(thetas[2]));

				dhTransformEl(1,0) = U(ceres::sin(thetas[2]));
				dhTransformEl(1,1) = U(ceres::cos(_alpha(2,0)) 
					* ceres:: cos(thetas[2]));
				dhTransformEl(1,2) = U(-ceres::cos(thetas[2]) 
					* ceres::sin(_alpha(2,0)));
				dhTransformEl(1,3) = U(_a(2,0) 
					* ceres::sin(thetas[2]));

				dhTransformEl(2,0) = U(0.0);
				dhTransformEl(2,1) = U(ceres::sin(_alpha(2,0)));
				dhTransformEl(2,2) = U(ceres::cos(_alpha(2,0)));
				dhTransformEl(2,3) = U(_d(2,0));

				dhTransformEl(3,0) = U(0.0);
				dhTransformEl(3,1) = U(0.0);
				dhTransformEl(3,2) = U(0.0);
				dhTransformEl(3,3) = U(1.0);
			}
			fkMatEl = fkMatEl * dhTransformEl;

		//predicting each ELBOW WITH OFFSET marker position wrt the base  
		//of the robot
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> orgin =
			Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
				::Zero(4,1);
		orgin << U(0.0), U(0.0), U(0.0), U(1.0);

		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
		 predElPt (3,1);
		predElPt(0,0) = U(fkMatEl(0,0)) * U(orgin(0,0)) 
			+ U(fkMatEl(0,1)) * U(orgin(1,0)) 
			+ U(fkMatEl(0,2)) * U(orgin(2,0)) 
			+ U(fkMatEl(0,3)) * U(orgin(3,0));
		predElPt(1,0) = U(fkMatEl(1,0)) * U(orgin(0,0)) 
			+ U(fkMatEl(1,1)) * U(orgin(1,0)) 
			+ U(fkMatEl(1,2)) * U(orgin(2,0)) 
			+ U(fkMatEl(1,3)) * U(orgin(3,0)); 
		predElPt(2,0) = U(fkMatEl(2,0)) * U(orgin(0,0)) 
			+ U(fkMatEl(2,1)) * U(orgin(1,0)) 
			+ U(fkMatEl(2,2)) * U(orgin(2,0)) 
			+ U(fkMatEl(2,3)) * U(orgin(3,0)); 

		residuals[3] = U(predElPt(0,0)) - U(_observedElPt(0,0));
		residuals[4] = U(predElPt(1,0)) - U(_observedElPt(1,0));
		residuals[5] = U(predElPt(2,0)) - U(_observedElPt(2,0));


		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
		static ceres::CostFunction* CreateThirdOrderSevenParam(
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedUaPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedElPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		 ) {
 	return (new ceres::AutoDiffCostFunction<IKfromEltoSh,6,8,4>
		( new IKfromEltoSh (fk, observedUaPt, observedElPt
				, parametricTerm, jointMinAngles, jointMaxAngles)));
	}

	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_observedUaPt, _observedElPt;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	 	_a, _alpha, _d;
};



/*
IK CHAIN 1.2 IKfromUaToEl

	INPUT2: ((J1,J2,J4)	& observedElPtInBase)		
		-> IK2 -> OUTPUT2: J3 
*/
//template <typename T>
struct IKfromUaToEl {
	IKfromUaToEl( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
			 &JsWam
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
				observedElPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		) :
		_fk(fk)
		,_JsWam(JsWam)
		, _observedElPt(observedElPt)
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		,_a(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _alpha(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _d(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1)) 
	{
		_a << T(0.0), T(0.0), T(0.045), T(-0.045), T(0.0)
				, T(0.0), T(0.0);
		_alpha << T(-M_PI/2.0), T(M_PI/2.0), T(-M_PI/2.0)
					 , T(M_PI/2.0), T(-M_PI/2.0), T(M_PI/2.0), T(0.0);
		_d << T(0.0), T(0.0), T(0.55), T(0.0), T(0.3)
			 , T(0.0), T(0.06);
	}

/* IK CHAIN 2 
			INPUT1: (observedElPt)	-> IK2 
				-> OUTPUT1: J3		*/

	template <typename U>
	bool operator()(const U* const candidateParamJ3
		, U* residuals) const {
		U *thetas = new U[3];

		//calculate thetas 1 & 2 from UA (IK1)
		thetas[0] = U(_JsWam(0,0));
		thetas[1] = U(_JsWam(0,1));
		//theta 3 - to solve JUST BASED ON THE ELBOW OFFSET
		thetas[2] = candidateParamJ3[0] 
			+	U(_parametricTerm) * candidateParamJ3[1] 
			+	U(_parametricTerm * _parametricTerm) * candidateParamJ3[2] 
			+	U(_parametricTerm * _parametricTerm * _parametricTerm) 
				* candidateParamJ3[3];
		
		//compute dh matrix (i-1)T(i)
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> fkMatEl =
			Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
				::Identity(4,4);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			dhTransformEl(4, 4);
		for(size_t i = 0; i < 3; i++) {
			bool retVal = true;	
			retVal =  retVal && (thetas[i] >= _jointMinAngles[i]);
			retVal =  retVal && (thetas[i] <= _jointMaxAngles[i]);

			if (retVal) {
				dhTransformEl(0,0) = U(ceres::cos(thetas[i]));
				dhTransformEl(0,1) = U(-ceres::cos(_alpha(i,0)) 
					* ceres::sin(thetas[i]));
				dhTransformEl(0,2) = U(ceres::sin(_alpha(i,0)) 
					* ceres::sin(thetas[i]));
				dhTransformEl(0,3) = U(_a(i,0) * ceres::cos(thetas[i]));

				dhTransformEl(1,0) = U(ceres::sin(thetas[i]));
				dhTransformEl(1,1) = U(ceres::cos(_alpha(i,0)) 
					*ceres:: cos(thetas[i]));
				dhTransformEl(1,2) = U(-ceres::cos(thetas[i]) 
					* ceres::sin(_alpha(i,0)));
				dhTransformEl(1,3) = U(_a(i,0) * ceres::sin(thetas[i]));

				dhTransformEl(2,0) = U(0.0);
				dhTransformEl(2,1) = U(ceres::sin(_alpha(i,0)));
				dhTransformEl(2,2) = U(ceres::cos(_alpha(i,0)));
				dhTransformEl(2,3) = U(_d(i,0));

				dhTransformEl(3,0) = U(0.0);
				dhTransformEl(3,1) = U(0.0);
				dhTransformEl(3,2) = U(0.0);
				dhTransformEl(3,3) = U(1.0);

				fkMatEl = fkMatEl * dhTransformEl;
			}
		}

		//predict LA position of the marker on the wam - move in +Z4
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			orgin(4,1);
		orgin << U(0.0), U(0.0), U(0.0), U(1.0);

		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
		 predElPt (3,1);
		predElPt(0,0) = U(fkMatEl(0,0)) * U(orgin(0,0)) 
			+ U(fkMatEl(0,1)) * U(orgin(1,0)) 
			+ U(fkMatEl(0,2)) * U(orgin(2,0)) 
			+ U(fkMatEl(0,3)) * U(orgin(3,0));

		predElPt(1,0) = U(fkMatEl(1,0)) * U(orgin(0,0)) 
			+ U(fkMatEl(1,1)) * U(orgin(1,0)) 
			+ U(fkMatEl(1,2)) * U(orgin(2,0)) 
			+ U(fkMatEl(1,3)) * U(orgin(3,0)); 
		predElPt(2,0) = U(fkMatEl(2,0)) * U(orgin(0,0)) 
			+ U(fkMatEl(2,1)) * U(orgin(1,0)) 
			+ U(fkMatEl(2,2)) * U(orgin(2,0)) 
			+ U(fkMatEl(2,3)) * U(orgin(3,0)); 


		residuals[0] = U(predElPt(0,0)) - U(_observedElPt(0,0));
		residuals[1] = U(predElPt(1,0)) - U(_observedElPt(1,0));
		residuals[2] = U(predElPt(2,0)) - U(_observedElPt(2,0));
		
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
		static ceres::CostFunction* CreateThirdOrderSevenParam(
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &JsWam
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedElPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		 ) {
 	return (new ceres::AutoDiffCostFunction<IKfromUaToEl,3, 4>
		( new IKfromUaToEl (fk, JsWam, observedElPt, parametricTerm
			, jointMinAngles, jointMaxAngles)));
	}


	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &_JsWam;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	_observedElPt;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	 	_a, _alpha, _d;
};

/* IK CHAIN 2 
			INPUT1: (J1,J2,J3,J4)	-> IK2 
				-> OUTPUT1: J5		*/

//template <typename T>
struct IKfromElToLa {
	IKfromElToLa( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
			 &JsWam
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
				observedWrPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		) :
		_fk(fk)
		,_JsWam(JsWam)
		, _observedWrPt(observedWrPt)
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		,_a(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _alpha(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _d(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1)) 
	{
		_a << T(0.0), T(0.0), T(0.045), T(-0.045), T(0.0)
				, T(0.0), T(0.0);
		_alpha << T(-M_PI/2.0), T(M_PI/2.0), T(-M_PI/2.0)
					, T(M_PI/2.0), T(-M_PI/2.0), T(M_PI/2.0), T(0.0);
		_d << T(0.0), T(0.0), T(0.55), T(0.0), T(0.3)
			 , T(0.0), T(0.06);
	}


	template <typename U>
	bool operator()(const U* const candidateParamJ5 
		, U* residuals) const {
		U *thetas = new U[5];

		//calculate thetas 1 & 2 from UA (IK1)
		thetas[0] = U(_JsWam(0,0));
		thetas[1] = U(_JsWam(0,1));
		//theta 3 - from El
		thetas[2] = U(_JsWam(0,2));
		//calculate thetas 4 from human Elbow angle (IK2.1)
		thetas[3] = U(_JsWam(0,3));

		thetas[4] = candidateParamJ5[0] 
			+	U(_parametricTerm) * candidateParamJ5[1] 
			+	U(_parametricTerm * _parametricTerm) * candidateParamJ5[2] 
			+	U(_parametricTerm * _parametricTerm * _parametricTerm) 
				* candidateParamJ5[3];
	
		//compute dh matrix (i-1)T(i)
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> fkMatEl =
			Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
				::Identity(4,4);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			dhTransformEl(4, 4);
		for(size_t i = 0; i < 5; i++) {
			bool retVal = true;	
			retVal =  retVal && (thetas[i] >= _jointMinAngles[i]);
			retVal =  retVal && (thetas[i] <= _jointMaxAngles[i]);

			if (retVal) {
				dhTransformEl(0,0) = U(ceres::cos(thetas[i]));
				dhTransformEl(0,1) = U(-ceres::cos(_alpha(i,0)) 
					* ceres::sin(thetas[i]));
				dhTransformEl(0,2) = U(ceres::sin(_alpha(i,0)) 
					* ceres::sin(thetas[i]));
				dhTransformEl(0,3) = U(_a(i,0) * ceres::cos(thetas[i]));

				dhTransformEl(1,0) = U(ceres::sin(thetas[i]));
				dhTransformEl(1,1) = U(ceres::cos(_alpha(i,0)) 
					*ceres:: cos(thetas[i]));
				dhTransformEl(1,2) = U(-ceres::cos(thetas[i]) 
					* ceres::sin(_alpha(i,0)));
				dhTransformEl(1,3) = U(_a(i,0) * ceres::sin(thetas[i]));

				dhTransformEl(2,0) = U(0.0);
				dhTransformEl(2,1) = U(ceres::sin(_alpha(i,0)));
				dhTransformEl(2,2) = U(ceres::cos(_alpha(i,0)));
				dhTransformEl(2,3) = U(_d(i,0));

				dhTransformEl(3,0) = U(0.0);
				dhTransformEl(3,1) = U(0.0);
				dhTransformEl(3,2) = U(0.0);
				dhTransformEl(3,3) = U(1.0);

				fkMatEl = fkMatEl * dhTransformEl;
			}
		}

		//predict LA position of the marker on the wam - move in +Z4
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			j4toLaTransWam(4,1), orgin(4,1);
		j4toLaTransWam << U(0.0), U(0.0), U(0.3/2.0), U(1.0);
		orgin << U(0.0), U(0.0), U(0.0), U(1.0);

		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
		 predWrPt (3,1);
		predWrPt(0,0) = U(fkMatEl(0,0)) * U(orgin(0,0)) 
			+ U(fkMatEl(0,1)) * U(orgin(1,0)) 
			+ U(fkMatEl(0,2)) * U(orgin(2,0)) 
			+ U(fkMatEl(0,3)) * U(orgin(3,0));

		predWrPt(1,0) = U(fkMatEl(1,0)) * U(orgin(0,0)) 
			+ U(fkMatEl(1,1)) * U(orgin(1,0)) 
			+ U(fkMatEl(1,2)) * U(orgin(2,0)) 
			+ U(fkMatEl(1,3)) * U(orgin(3,0)); 
		predWrPt(2,0) = U(fkMatEl(2,0)) * U(orgin(0,0)) 
			+ U(fkMatEl(2,1)) * U(orgin(1,0)) 
			+ U(fkMatEl(2,2)) * U(orgin(2,0)) 
			+ U(fkMatEl(2,3)) * U(orgin(3,0)); 


		residuals[0] = U(predWrPt(0,0)) - U(_observedWrPt(0,0));
		residuals[1] = U(predWrPt(1,0)) - U(_observedWrPt(1,0));
		residuals[2] = U(predWrPt(2,0)) - U(_observedWrPt(2,0));
		
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
		static ceres::CostFunction* CreateThirdOrderSevenParam(
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &JsWam
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedWrPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		 ) {
 	return (new ceres::AutoDiffCostFunction<IKfromElToLa,3, 4>
		( new IKfromElToLa (fk, JsWam, observedWrPt, parametricTerm
			, jointMinAngles, jointMaxAngles)));
	}


	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &_JsWam;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	_observedWrPt;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	 	_a, _alpha, _d;
};




/*
IK CHAIN 2 IKfitFromElbToLa

	INPUT2: (OUTPUT1(J1,J2,J3')	& observedLaPtInBase & observed)		-> IK2 -> OUTPUT2: J3,J4,J5'
*/
//template <typename T>
struct IKfromUaToLa {
	IKfromUaToLa( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
			 &JsWam
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
				observedLaPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		) :
		_fk(fk)
		,_JsWam(JsWam)
		, _observedLaPt(observedLaPt)
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		,_a(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _alpha(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _d(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1)) 
	{
		_a << T(0.0), T(0.0), T(0.045), T(-0.045), T(0.0)
				, T(0.0), T(0.0);
		_alpha << T(-M_PI/2.0), T(M_PI/2.0), T(-M_PI/2.0)
					 , T(M_PI/2.0), T(-M_PI/2.0), T(M_PI/2.0), T(0.0);
		_d << T(0.0), T(0.0), T(0.55), T(0.0), T(0.3)
			 , T(0.0), T(0.06);
	}

/* IK CHAIN 2 
			INPUT1: (observedElPt)	-> IK2 
				-> OUTPUT1: J3		*/

	template <typename U>
	bool operator()(const U* const candidateParamJ3
		,const U* const candidateParamJ5 , U* residuals) const {
		U *thetas = new U[_fk.getNJoints()];

		//calculate thetas 1 & 2 from UA (IK1)
		thetas[0] = U(_JsWam(0,0));
		thetas[1] = U(_JsWam(0,1));
		//theta 3 - to solve
		thetas[2] = candidateParamJ3[0] 
			+	U(_parametricTerm) * candidateParamJ3[1] 
			+	U(_parametricTerm * _parametricTerm) * candidateParamJ3[2] 
			+	U(_parametricTerm * _parametricTerm * _parametricTerm) 
				* candidateParamJ3[3];
		//calculate thetas 4 from human Elbow angle (IK2.1)
		thetas[3] = U(_JsWam(0,3));

		thetas[4] = candidateParamJ5[0] 
			+	U(_parametricTerm) * candidateParamJ5[1] 
			+	U(_parametricTerm * _parametricTerm) * candidateParamJ5[2] 
			+	U(_parametricTerm * _parametricTerm * _parametricTerm) 
				* candidateParamJ5[3];
		thetas[5] = U(0.0);
		thetas[6] = U(0.0);
 
		//compute dh matrix (i-1)T(i)
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> fkMatEl =
			Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
				::Identity(4,4);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			dhTransformEl(4, 4);
		for(size_t i = 0; i < _fk.getNJoints(); i++) {
			bool retVal = true;	
			retVal =  retVal && (thetas[i] >= _jointMinAngles[i]);
			retVal =  retVal && (thetas[i] <= _jointMaxAngles[i]);

			if (retVal) {
				dhTransformEl(0,0) = U(ceres::cos(thetas[i]));
				dhTransformEl(0,1) = U(-ceres::cos(_alpha(i,0)) 
					* ceres::sin(thetas[i]));
				dhTransformEl(0,2) = U(ceres::sin(_alpha(i,0)) 
					* ceres::sin(thetas[i]));
				dhTransformEl(0,3) = U(_a(i,0) * ceres::cos(thetas[i]));

				dhTransformEl(1,0) = U(ceres::sin(thetas[i]));
				dhTransformEl(1,1) = U(ceres::cos(_alpha(i,0)) 
					*ceres:: cos(thetas[i]));
				dhTransformEl(1,2) = U(-ceres::cos(thetas[i]) 
					* ceres::sin(_alpha(i,0)));
				dhTransformEl(1,3) = U(_a(i,0) * ceres::sin(thetas[i]));

				dhTransformEl(2,0) = U(0.0);
				dhTransformEl(2,1) = U(ceres::sin(_alpha(i,0)));
				dhTransformEl(2,2) = U(ceres::cos(_alpha(i,0)));
				dhTransformEl(2,3) = U(_d(i,0));

				dhTransformEl(3,0) = U(0.0);
				dhTransformEl(3,1) = U(0.0);
				dhTransformEl(3,2) = U(0.0);
				dhTransformEl(3,3) = U(1.0);
			}
				fkMatEl = fkMatEl * dhTransformEl;
		}

		//predict LA position of the marker on the wam - move in +Z4
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			j4toLaTransWam(4,1), orgin(4,1);
		j4toLaTransWam << U(0.0), U(0.0), U(0.3/2.0), U(1.0);
		orgin << U(0.0), U(0.0), U(0.0), U(1.0);

		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
		 predLaPt (3,1);
		predLaPt(0,0) = U(fkMatEl(0,0)) * U(orgin(0,0)) 
			+ U(fkMatEl(0,1)) * U(orgin(1,0)) 
			+ U(fkMatEl(0,2)) * U(orgin(2,0)) 
			+ U(fkMatEl(0,3)) * U(orgin(3,0));

		predLaPt(1,0) = U(fkMatEl(1,0)) * U(orgin(0,0)) 
			+ U(fkMatEl(1,1)) * U(orgin(1,0)) 
			+ U(fkMatEl(1,2)) * U(orgin(2,0)) 
			+ U(fkMatEl(1,3)) * U(orgin(3,0)); 
		predLaPt(2,0) = U(fkMatEl(2,0)) * U(orgin(0,0)) 
			+ U(fkMatEl(2,1)) * U(orgin(1,0)) 
			+ U(fkMatEl(2,2)) * U(orgin(2,0)) 
			+ U(fkMatEl(2,3)) * U(orgin(3,0)); 


		residuals[0] = U(predLaPt(0,0)) - U(_observedLaPt(0,0));
		residuals[1] = U(predLaPt(1,0)) - U(_observedLaPt(1,0));
		residuals[2] = U(predLaPt(2,0)) - U(_observedLaPt(2,0));
		
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
		static ceres::CostFunction* CreateThirdOrderSevenParam(
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &JsWam
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedLaPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		 ) {
 	return (new ceres::AutoDiffCostFunction<IKfromUaToLa,3, 4,4>
		( new IKfromUaToLa (fk, JsWam, observedLaPt, parametricTerm
			, jointMinAngles, jointMaxAngles)));
	}


	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &_JsWam;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	_observedLaPt;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	 	_a, _alpha, _d;
};

/* 
IK CHAIN 1 B
	INPUT1: (observedUaPt & observedElPt)		-> IK1 -> OUTPUT1: J1,J2,J3'	
*/

//template <typename T>
struct IKfromWrToSh {
	IKfromWrToSh( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &JsWam
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
			observedUaPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
			observedElPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
			observedWrPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles

		) :
		_fk(fk)
		, _JsWam(JsWam)
		, _observedUaPt(observedUaPt)
		, _observedElPt(observedElPt)
		, _observedWrPt(observedWrPt)
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		,_a(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _alpha(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _d(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1)) 
	{
		
		_a << T(0.0), T(0.0), T(0.045), T(-0.045), T(0.0)
				, T(0.0), T(0.0);
		_alpha << T(-M_PI/2.0), T(M_PI/2.0), T(-M_PI/2.0)
					 , T(M_PI/2.0), T(-M_PI/2.0), T(M_PI/2.0), T(0.0);
		_d << T(0.0), T(0.0), T(0.55), T(0.0), T(0.3)
			 , T(0.0), T(0.06);
	}

	template <typename U>
	bool operator()(const U* const candidateParamUa
		, const U* const candidateParamEl
		, const U* const candidateParamWr
		, U* residuals) const {

	U *thetas = new U[_fk.getNJoints()];	
	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
		fkMat =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4), dhTransform(4, 4);

		//calculate theta1 & theta2
		for(size_t i = 0; i < 2; i++) {
			size_t iTimesFour = i * 4;
	
			thetas[i] = candidateParamUa[iTimesFour] 
				+	U(_parametricTerm) * candidateParamUa[iTimesFour + 1] 
				+	U(_parametricTerm * _parametricTerm) 
					* candidateParamUa[iTimesFour + 2] 
				+	U(_parametricTerm * _parametricTerm * _parametricTerm) 
					* candidateParamUa[iTimesFour + 3];
	
			//if joints within limit, set to wam		
			bool retVal = true;
			retVal =  retVal && (thetas[i] >= _jointMinAngles[i]);
			retVal =  retVal && (thetas[i] <= _jointMaxAngles[i]);

			if (retVal) {
				//compute dh matrix (i-1)T(i)
				dhTransform(0,0) = U(ceres::cos(thetas[i]));
				dhTransform(0,1) = U(-ceres::cos(_alpha(i,0)) 
					* ceres::sin(thetas[i]));
				dhTransform(0,2) = U(ceres::sin(_alpha(i,0)) 
					* ceres::sin(thetas[i]));
				dhTransform(0,3) = U(_a(i,0) * ceres::cos(thetas[i]));

				dhTransform(1,0) = U(ceres::sin(thetas[i]));
				dhTransform(1,1) = U(ceres::cos(_alpha(i,0)) 
					* ceres:: cos(thetas[i]));
				dhTransform(1,2) = U(-ceres::cos(thetas[i]) 
					* ceres::sin(_alpha(i,0)));
				dhTransform(1,3) = U(_a(i,0) 
					* ceres::sin(thetas[i]));

				dhTransform(2,0) = U(0.0);
				dhTransform(2,1) = U(ceres::sin(_alpha(i,0)));
				dhTransform(2,2) = U(ceres::cos(_alpha(i,0)));
				dhTransform(2,3) = U(_d(i,0));

				dhTransform(3,0) = U(0.0);
				dhTransform(3,1) = U(0.0);
				dhTransform(3,2) = U(0.0);
				dhTransform(3,3) = U(1.0);
			}
			fkMat = fkMat * dhTransform;
		}

		//predicting each ELBOW-NO_OFFSET marker position wrt the base  
		//of the robot
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> wamEltoUaTrans =
			Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
				::Zero(4,1);
		wamEltoUaTrans << U(0.0), U(0.0), U(0.55), U(1.0);

		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
		 predUaPt (3,1);
		predUaPt(0,0) = U(fkMat.row(0).dot(wamEltoUaTrans.col(0)));
		predUaPt(1,0) = U(fkMat.row(1).dot(wamEltoUaTrans.col(0)));
		predUaPt(2,0) = U(fkMat.row(2).dot(wamEltoUaTrans.col(0)));


		//calculate theta3
			thetas[2] = candidateParamEl[0] 
				+	U(_parametricTerm) * candidateParamEl[1] 
				+	U(_parametricTerm * _parametricTerm) 
					* candidateParamEl[2] 
				+	U(_parametricTerm * _parametricTerm * _parametricTerm) 
					* candidateParamEl[3];
	
			//if joints within limit, set to wam		
			bool retVal = true;
			retVal =  retVal && (thetas[2] >= _jointMinAngles[2]);
			retVal =  retVal && (thetas[2] <= _jointMaxAngles[2]);

			if (retVal) {
				//compute dh matrix (i-1)T(i)
				dhTransform(0,0) = U(ceres::cos(thetas[2]));
				dhTransform(0,1) = U(-ceres::cos(_alpha(2,0)) 
					* ceres::sin(thetas[2]));
				dhTransform(0,2) = U(ceres::sin(_alpha(2,0)) 
					* ceres::sin(thetas[2]));
				dhTransform(0,3) = U(_a(2,0) * ceres::cos(thetas[2]));

				dhTransform(1,0) = U(ceres::sin(thetas[2]));
				dhTransform(1,1) = U(ceres::cos(_alpha(2,0)) 
					* ceres:: cos(thetas[2]));
				dhTransform(1,2) = U(-ceres::cos(thetas[2]) 
					* ceres::sin(_alpha(2,0)));
				dhTransform(1,3) = U(_a(2,0) 
					* ceres::sin(thetas[2]));

				dhTransform(2,0) = U(0.0);
				dhTransform(2,1) = U(ceres::sin(_alpha(2,0)));
				dhTransform(2,2) = U(ceres::cos(_alpha(2,0)));
				dhTransform(2,3) = U(_d(2,0));

				dhTransform(3,0) = U(0.0);
				dhTransform(3,1) = U(0.0);
				dhTransform(3,2) = U(0.0);
				dhTransform(3,3) = U(1.0);
			}
			fkMat = fkMat * dhTransform;

		//predicting each ELBOW WITH OFFSET marker position wrt the base  
		//of the robot
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
		 orgin(4,1);
		orgin << U(0.0), U(0.0), U(0.0), U(1.0);

		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
		 predElPt (3,1);
		predElPt(0,0) = U(fkMat.row(0).dot(orgin.col(0)));
		predElPt(1,0) = U(fkMat.row(1).dot(orgin.col(0)));
		predElPt(2,0) = U(fkMat.row(2).dot(orgin.col(0)));



		/* Theta4 - find Human elbow angle, i.e.
			INPUT( <)UpperArm.LowerArm)	-> IK3 
				-> OUTPUT(J4) */	
		thetas[3] = U(_JsWam(0,3));

		/* calculate theta5 */
			thetas[4] = candidateParamWr[0] 
				+	U(_parametricTerm) * candidateParamWr[1] 
				+	U(_parametricTerm * _parametricTerm) 
					* candidateParamWr[2] 
				+	U(_parametricTerm * _parametricTerm * _parametricTerm) 
					* candidateParamWr[3];
			//Just to solve for Theta 5 with non-zero value
			thetas[5] = U(0.0);
			thetas[6] = U(0.0);
	
			for (size_t j = 3; j < _fk.getNJoints(); j++) {
				//if joints within limit, set to wam		
				bool retVal = true;
				retVal =  retVal && (thetas[j] >= _jointMinAngles[j]);
				retVal =  retVal && (thetas[j] <= _jointMaxAngles[j]);

				if (retVal) {
					//compute dh matrix (i-1)T(i)
					dhTransform(0,0) = U(ceres::cos(thetas[j]));
					dhTransform(0,1) = U(-ceres::cos(_alpha(j,0)) 
						* ceres::sin(thetas[j]));
					dhTransform(0,2) = U(ceres::sin(_alpha(j,0)) 
						* ceres::sin(thetas[j]));
					dhTransform(0,3) = U(_a(j,0) * ceres::cos(thetas[j]));

					dhTransform(1,0) = U(ceres::sin(thetas[j]));
					dhTransform(1,1) = U(ceres::cos(_alpha(j,0)) 
						* ceres:: cos(thetas[j]));
					dhTransform(1,2) = U(-ceres::cos(thetas[j]) 
						* ceres::sin(_alpha(j,0)));
					dhTransform(1,3) = U(_a(j,0) 
						* ceres::sin(thetas[j]));

					dhTransform(2,0) = U(0.0);
					dhTransform(2,1) = U(ceres::sin(_alpha(j,0)));
					dhTransform(2,2) = U(ceres::cos(_alpha(j,0)));
					dhTransform(2,3) = U(_d(j,0));

					dhTransform(3,0) = U(0.0);
					dhTransform(3,1) = U(0.0);
					dhTransform(3,2) = U(0.0);
					dhTransform(3,3) = U(1.0);
				}
				fkMat = fkMat * dhTransform;
			}

		//predicting each WRIST marker position wrt the base  
		//of the robot
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
		 predWrPt (3,1);
		predWrPt(0,0) = U(fkMat.row(0).dot(orgin.col(0)));
		predWrPt(1,0) = U(fkMat.row(1).dot(orgin.col(0)));
		predWrPt(2,0) = U(fkMat.row(2).dot(orgin.col(0)));


		residuals[0] = U(100.0) * U(predUaPt(0,0)) 
			-  U(100.0) * U(_observedUaPt(0,0));
		residuals[1] =  U(100.0) * U(predUaPt(1,0)) 
			-  U(100.0) * U(_observedUaPt(1,0));
		residuals[2] =  U(100.0) * U(predUaPt(2,0)) 
			-  U(100.0) * U(_observedUaPt(2,0));

		residuals[3] = U(100.0) * U(predElPt(0,0)) 
			- U(100.0) * U(_observedElPt(0,0));
		residuals[4] = U(100.0) * U(predElPt(1,0)) 
			-  U(100.0) * U(_observedElPt(1,0));
		residuals[5] = U(100.0) * U(predElPt(2,0)) 
			- U(100.0) * U(_observedElPt(2,0));

		residuals[6] = U(100.0) * U(predWrPt(0,0)) 
			- U(100.0) * U(_observedWrPt(0,0));
		residuals[7] = U(100.0) * U(predWrPt(1,0)) 
			- U(100.0) * U(_observedWrPt(1,0));
		residuals[8] = U(100.0) * U(predWrPt(2,0)) 
			- U(100.0) * U(_observedWrPt(2,0));


		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
		static ceres::CostFunction* CreateThirdOrderSevenParam(
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &JsWam
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedUaPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedElPt
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> observedWrPt
		, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles
		 ) {
 	return (new ceres::AutoDiffCostFunction<IKfromWrToSh,9,8,4,4>
		( new IKfromWrToSh (fk, JsWam, observedUaPt, observedElPt
			, observedWrPt, parametricTerm, jointMinAngles, jointMaxAngles)));
	}

	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &_JsWam;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_observedUaPt, _observedElPt, _observedWrPt;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	 	_a, _alpha, _d;
};

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_inPtsAlongRowsMB;
	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		_inPtsRUAinBase, _inPtsRELinBase, _inPtsRLAinBase
		, _inPtsRWinBase, _inPtsRTHinBase, _origin;
	T _parametricTerm;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &_JsWam;
};

#endif
