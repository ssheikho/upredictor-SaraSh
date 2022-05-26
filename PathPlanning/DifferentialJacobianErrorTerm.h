#ifndef DIFFERENTIAL_JACOBIAN_ERROR_TERM_H
#define DIFFERENTIAL_JACOBIAN_ERROR_TERM_H


#include "SpatialJacobian.h"
#include "DifferentialIKTypes.h"
#include "Eigen/Core"
#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include <Eigen/SVD>
#include <Eigen/Geometry>
#include <cmath>
#include "ceres/autodiff_cost_function.h"

using namespace Eigen;

namespace ceres {

struct DifferentialVelocityConstraint	{

	DifferentialVelocityConstraint (
		ForwardKin<double> &fk
		// observed scaled Vectors (Cartesian Positions)
		// on the wam at t_i
		, const double parametricTerm
		, std::vector<double> jointMinAngles
		, std::vector<double> jointMaxAngles
		//wrist center owr = o7 − d7.R.[0 0 1]^T			
		, const Eigen::Matrix<double,3,1> shToUaVi
		, const Eigen::Matrix<double,3,1> elOffsetVi 
		, const Eigen::Matrix<double,3,1> wrOffsetVi
		, const Eigen::Matrix<double,3,1> laToWrVi
		, const Eigen::Matrix<double,3,1> shoToWrVi
		//inverse of covariance matrix (x y z)
		, const Eigen::Matrix<double, 3, 3>& iMatShToUaWam
		, const Eigen::Matrix<double, 3, 3>& iMatElOffsetWam
		, const Eigen::Matrix<double, 3, 3>& iMatWrOffsetWam
		, const Eigen::Matrix<double, 3, 3>& iMatLaToWrWam
		, const Eigen::Matrix<double, 3, 3>& iMatShToWrWam
		) : _fk(fk)
		,_parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		,_a(_fk.getAvect(7))
		, _alpha(_fk.getAlphaVect(7))
		, _d(_fk.getDvect(7)) 
		,	_shToUaVi(shToUaVi)
		,	_elOffsetVi(elOffsetVi)
		,	_wrOffsetVi(wrOffsetVi)
		,	_laToWrVi(laToWrVi)
		,	_shoToWrVi(shoToWrVi)
		, _iMatShToUaWam(iMatShToUaWam)
		, _iMatElOffsetWam(iMatElOffsetWam)
		, _iMatWrOffsetWam(iMatWrOffsetWam)
		, _iMatLaToWrWam(iMatLaToWrWam)
		, _iMatShToWrWam(iMatShToWrWam)
	{}


		//oUA = o3 - a3.R3.[1 0 0]^T
	template <typename U>
	bool operator()( const U* const qBase_
		, const U* const candidateParamJs
		, U* residuals) const {
		
		//first 4 joints determine position constraints
		int toDof = 5;   
		
		// (0). model thetas as a cont. function of time 
		U *thetas_i = new U[toDof];	
		bool retVal = true;
		for(size_t i = 0; i < toDof; i++) {
			size_t iTimesFour = i * 4;
	
			thetas_i[i] = candidateParamJs[iTimesFour] 
				+	U(_parametricTerm) * candidateParamJs
					[iTimesFour + 1] +	U(_parametricTerm 
				* _parametricTerm) * candidateParamJs
					[iTimesFour + 2] +	U(_parametricTerm 
				* _parametricTerm * _parametricTerm) 
				* candidateParamJs[iTimesFour + 3];

			retVal =  retVal && (
				thetas_i[i] >= _jointMinAngles[i]);
			retVal =  retVal && (
				thetas_i[i] <= _jointMaxAngles[i]);
			
		}
		
		
		// (1). Form the homogeneous transformation matrices
		// 			A_i by substituting the parameters for each 
		//			link into DH-matrix of eqn. (3.10).

		//initialize the matrices to zero;
		//avoid calling external eigen::identity function
		Eigen::Matrix<U, 4, 4> tUa, tEl, tLa, tWr;
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {
				tUa(i,j) = U(0.0);	tEl(i,j) = U(0.0);
				tWr(i,j) = U(0.0); tLa(i,j) = U(0.0);
			}
			tUa(i,i) = U(1.0); tEl(i,i) = U(1.0);
				tWr(i,i) = U(1.0);	tLa(i,i) = U(1.0);
		}		

		Eigen::Matrix<U, 4, 4> A_j;
		for(int i = 0; i < toDof; i++) {
			//1st row 
			A_j(0,0) = 
				U(ceres::cos(thetas_i[i]));
			A_j(0,1) = U(-ceres::cos(_alpha(i,0)) 
				* ceres::sin(thetas_i[i]));
			A_j(0,2) = U(ceres::sin(_alpha(i,0)) 
				* ceres::sin(thetas_i[i]));
			A_j(0,3) = 
				U(_a(i,0) * ceres::cos(thetas_i[i]));
			//2nd row 
			A_j(1,0) = 
				U(ceres::sin(thetas_i[i]));
			A_j(1,1) = U(ceres::cos(_alpha(i,0)) 
				* ceres:: cos(thetas_i[i]));
			A_j(1,2) = U(-ceres::cos(thetas_i[i]) 
				* ceres::sin(_alpha(i,0)));
			A_j(1,3) = 
				U(_a(i,0) * ceres::sin(thetas_i[i]));
			//3rd row 
			A_j(2,0) = U(0.0);
			A_j(2,1) = U(ceres::sin(_alpha(i,0)));
			A_j(2,2) = U(ceres::cos(_alpha(i,0)));
			A_j(2,3) = U(_d(i,0));
			//4th row 
			A_j(3,0) = U(0.0);
			A_j(3,1) = U(0.0);
			A_j(3,2) = U(0.0);
			A_j(3,3) = U(1.0);

			
			// (2). Form T^n_0 = A 1 · · · A n . 
			//			This then gives the position (and
			//			orientation for wrist later) of the 
			//			tool frame expressed in base coordinates.
			//			

			//	"A1 * A2 * vUA": 
			//			-> Elbow center, pUa(i)				
			if (i < 2)
				tUa = tUa * A_j;
			//	"A1 * A2 * A3 * Origin": -> pEl(i)
			//	elbowOffset(i) = pUa(i)	- pEl(i)
			if (i < 3)
				tEl = tEl * A_j;
			if (i < 4)
				tLa = tLa * A_j;
			//	"A1 * A2 * ...* A5 * Origin": -> Wrist center
			//	pWr (center) = o5 = o6 OR (ALTERNATIVE)
			//							 = o4 + d5.R4.[0 0 1]^T			
			if (i < 5)
				tWr = tWr * A_j;
		}



    return true;
	}

	static ceres::CostFunction* Create(
		ForwardKin<double> &fk
		, const double parametricTerm
		, std::vector<double> jointMinAngles
		, std::vector<double> jointMaxAngles
		, const Eigen::Matrix<double,3,1> shToUaVi
		, const Eigen::Matrix<double,3,1> elOffsetVi 
		, const Eigen::Matrix<double,3,1> wrOffsetVi
		, const Eigen::Matrix<double,3,1> laToWrVi
		, const Eigen::Matrix<double,3,1> shoToWrVi
		//inverse of covariance matrix (x y z)
		, const Eigen::Matrix<double, 3, 3>& iMatShToUaWam
		, const Eigen::Matrix<double, 3, 3>& iMatElOffsetWam
		, const Eigen::Matrix<double, 3, 3>& iMatWrOffsetWam
		, const Eigen::Matrix<double, 3, 3>& iMatLaToWrWam
		, const Eigen::Matrix<double, 3, 3>& iMatShToWrWam
		)	
	{
		return (new ceres::AutoDiffCostFunction
			<DifferentialVelocityConstraint, 3 ,4, 12> (
				new DifferentialVelocityConstraint (
				fk , parametricTerm, jointMinAngles
				, jointMaxAngles, shToUaVi, elOffsetVi
				, wrOffsetVi,	laToWrVi, shoToWrVi
				,	iMatShToUaWam, iMatElOffsetWam
				, iMatWrOffsetWam, iMatLaToWrWam
				, iMatShToWrWam)));
	}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	ForwardKin<double> &_fk;
	const double _parametricTerm;
	std::vector<double> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<double, Eigen::Dynamic, 1> 
	 	_a, _alpha, _d;
	//measured/scaled Cartesian vectores on the wam at t_i
	const Eigen::Matrix<double, 3, 1> _shToUaVi, _elOffsetVi,
		_wrOffsetVi,	_laToWrVi, _shoToWrVi;
	//measurements covariance (x y z)
	const Eigen::Matrix<double, 3, 3> _iMatShToUaWam,
		_iMatElOffsetWam, _iMatWrOffsetWam,
	  _iMatLaToWrWam, _iMatShToWrWam;
};




class DifferentialJacobianErrorTerm {

public:
	DifferentialJacobianErrorTerm (
		ForwardKin<double> &fk
		//Observed Cartesian Positions t_(i)
		, Eigen::Matrix<double, Eigen::Dynamic, 1> inPi
		//measured linear velocity
		, Eigen::Matrix<double, Eigen::Dynamic, 1> 
				inLinearVi
		//only for the end-effector
		, Eigen::Matrix<double, 3,	1> inAngularVi
		//measurements covariance information (x y z φ θ ψ)
		, const Eigen::Matrix<double, 6, 6>& information
		) : _fk(fk)
		,_a(_fk.getAvect(7))
		, _alpha(_fk.getAlphaVect(7))
		, _d(_fk.getDvect(7)) 
		,	_pShoWam(inPi.block(Markers::RSh,0,3,1))
		,	_pUaWam(inPi.block(Markers::RUA,0,3,1))
		, _pElWam(inPi.block(Markers::REl,0,3,1))
		, _pLaWam(inPi.block(Markers::RLA,0,3,1))
		, _pWrWam(inPi.block(Markers::RWr,0,3,1))
		, _pThWam(inPi.block(Markers::RTh,0,3,1))
		, _pPiWam(inPi.block(Markers::RPi,0,3,1))
		,	_vShoWam(inLinearVi.block(Markers::RSh,0,3,1))
		,	_vUaWam(inLinearVi.block(Markers::RUA,0,3,1))
		, _vElWam(inLinearVi.block(Markers::REl,0,3,1))
		, _vLaWam(inLinearVi.block(Markers::RLA,0,3,1))
		, _vWrWam(inLinearVi.block(Markers::RWr,0,3,1))
		, _vThWam(inLinearVi.block(Markers::RTh,0,3,1))
		, _vPiWam(inLinearVi.block(Markers::RPi,0,3,1))
		//measured angular velocity of EE 
		, _omegaEE(inAngularVi)	
		//measurements covariance information (x y z φ θ ψ)
		, _information(information) {}

	template <typename U>
	bool operator()(const U* const EulAnglesBase_i
		, const U* const EulAnglesBase_iMinusOne
		, const U* const thetas_i
		, const U* const thetas_iMinusOne
		, U* residuals_ptr) const {
		int toDof = 7;
/*
		//(1).	compute Wrist dh matrix (B)T(n)	
		Eigen::Matrix	<U, 4, 4> fkMatToO7 =
			Eigen::Matrix<U, 4, 4>::Identity(4,4);
		DHTransformTo	(toDof, thetas_i, fkMatToO7);

		//(2).	compute Basic Jacobian w.r.t Base J(thetas)
		Eigen::Matrix<U, 6, Eigen::Dynamic> basicJ 
			=	Eigen::Matrix<U, 6, Eigen::Dynamic>
				::Zero(6, toDof);		
		basicJacobianMatTo	(toDof, thetas_i,	basicJ);
		
		//(3). compute Rot matrix Base w.r.t Camera 
		//		 from Euler Angles
		Eigen::Matrix<U, 3, 3>	RmatBase 
			=	Eigen::Matrix<U, 3, 3>::Identity(3,3);
		RmatFromEulAngles	(EulAnglesBase_i,	RmatBase);	

		//(4). Spatial Jacobian
		Eigen::Matrix<U, 6, Eigen::Dynamic>	spatialJMat 
			=	Eigen::Matrix<U, 6, Eigen::Dynamic>
				::Zero( 6,13);
		spatialJacobianMatTo (fkMatToO7, RmatBase, basicJ
			,  spatialJMat.data());		

		//(5). Generalized velocity vector u(t_i)
		Eigen::Matrix<U, 13, 1> uVelvec = 
			Eigen::Matrix<U, 13, 1>::Zero(13,1);

		//	(5.1).	MEASURED Linear Velocity Base
		// w.r.t Camera Frame
		uVelvec.template block<3, 1>(0, 0) = 
			_vShoWam.template cast<U>() ;
		//	(5.2).	Angular Velocity Base w.r.t Base-frame
		Eigen::Matrix<U, 3, 1> omegaBase = 
			Eigen::Matrix<U, 3, 1>::Zero(3,1);
		angularVel (EulAnglesBase_i, EulAnglesBase_iMinusOne
			, &omegaBase);
		//W.r.t base frame: B_ω = baseRmat^T * ω
		uVelvec.template block<3, 1>(3, 0) = 
			RmatBase.transpose() * omegaBase;
		//	(5.3). thetaDot Vector
		for (int i = 0; i < toDof; i++)	
			uVelvec(6+i, 0) = thetas_i[i] - thetas_iMinusOne[i];

		//(6). measured end-effector velocity
		Eigen::Matrix <U, 6, 1> measuredEeVel = 
			Eigen::Matrix <U, 6, 1>::Zero(6,1);
		measuredEeVel.template block<3, 1>(0, 0) = 
			_vWrWam.template cast<U>();
		measuredEeVel.template block<3, 1>(3, 0) = 
			_omegaEE.template cast<U>();

		//estimated wrist postion
		Eigen::Matrix	<U, 4, 4> fkMatToO5 =
			Eigen::Matrix<U, 4, 4>::Identity(4,4);
		DHTransformTo	(5, thetas_i, fkMatToO5);
		Eigen::Matrix	<U, 3, 1> pWrEstimatedInBase = 
			fkMatToO5.template block<3, 1>(0, 3);
		Eigen::Matrix	<U, 3, 1> pWrEstimated =
			RmatBase.transpose() * pWrEstimatedInBase
			+ _pShoWam.template cast<U>();
		//estimated wrist velocity
		Eigen::Matrix	<U, 3, 1> vWrEstimated =
			spatialJMat * uVelvec;

		//(7).	residuals
    Eigen::Map<Eigen::Matrix<U, 9, 1> >
			residuals(residuals_ptr);
		// position error
		residuals.template block<3, 1>(6, 0) =
			pWrEstimated - _pWrWam.template cast<U>();

		// velocity error
		residuals.template block<6, 1>(0, 0) = 
			measuredEeVel - (spatialJMat * uVelvec);
*/
    return true;
	} 


	static ceres::CostFunction* Create(ForwardKin<double> &fk
		,	const Eigen::Matrix<double , Eigen::Dynamic, 1> inPtsI
		,	const Eigen::Matrix<double , Eigen::Dynamic, 1> 
				inLinearVi
		,	const Eigen::Matrix<double , 3, 1> inAngularVWri
		//measurements covariance information (x y z φ θ ψ)
		, const Eigen::Matrix<double, 6, 6>& information)	{

		return (new ceres::AutoDiffCostFunction
			<DifferentialJacobianErrorTerm, 9, 3, 3, 7, 7> (
				new DifferentialJacobianErrorTerm (fk, inPtsI
					,inLinearVi ,	inAngularVWri, information)));
	}


//A.  solve at Velocity and delta=0
protected:
	ForwardKin<double> &_fk;
//	Eigen::Matrix<T, Eigen::Dynamic, 7> &_JsWam;
	Eigen::Matrix<double, Eigen::Dynamic, 1> 
	 	_a, _alpha, _d;
	double _nInPts;

	Eigen::Matrix<double, 3, 1> 
		_pShoWam, _pUaWam, _pElWam, _pLaWam
		, _pWrWam, _pThWam, _pPiWam;
	Eigen::Matrix<double, 3, 1> 
		_vShoWam, _vUaWam, _vElWam, _vLaWam
		, _vWrWam, _vThWam, _vPiWam;

	Eigen::Matrix<double, 3, 1> _omegaEE;

  // The square root of the measurement information matrix.
  const Eigen::Matrix<double, 6, 6> _information;

};
}
#endif
