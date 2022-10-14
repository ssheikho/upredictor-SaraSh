#ifndef DIFFERENTIAL_IK_ERROR_TERM_H
#define DIFFERENTIAL_IK_ERROR_TERM_H


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


/*	(2).	compute Basic Jacobian wrt BASE J(thetas)*/
template <typename U>
inline void basicJacobianMatTo	(int toDof
	, const U* const thetas_i
	,	Eigen::Matrix <U, 6, Eigen::Dynamic> &basicJ)	{

	Eigen::Matrix<U, 4, 4> fkMatToDof =
		Eigen::Matrix<U, 4, 4>::Identity(4,4);
	fkMatToDof = DHTransformTo	(toDof, thetas_i);
	Eigen::Matrix<U, 3, 1> On = fkMatToDof.block(0,3,3,1);

	Eigen::Matrix<U, 4, 4> fkMat 
			=	Eigen::Matrix<U, 4, 4>::Identity(4,4);
	for(int i = 0; i < toDof; i++) { 
		fkMat = DHTransformTo	(i, thetas_i);
		//	Jvi = Zi-1 x (On-Oi-1)
		Eigen::Matrix<U, 3, 1> ZiMinus1 =  fkMat.block(0,2,3,1)
			, OiMinus1 =  fkMat.block(0,3,3,1);
		basicJ.block(0,i,3,1) = ZiMinus1.cross(On- OiMinus1);
		//Jωi = Zi-1 
		basicJ.block(3,i,3,1) = ZiMinus1;
	}

}

//(3).	get RotMat from Euler Parameters 
template <typename U>
inline void RmatFromEulAngles	(const U* const EulAngles 
	,	Eigen::Matrix<U, 3, 3>	&Rmat)	{
		//1st col
		Rmat(0,0) = U(ceres::cos(EulAngles[1]) 
			* ceres::cos(EulAngles[0])
			* ceres::cos(EulAngles[2])) 
			- U(ceres::sin(EulAngles[0]) 
			* ceres::sin(EulAngles[2]));
		Rmat(1,0) = U(ceres::cos(EulAngles[0])
			*	ceres::sin(EulAngles[2])) 
			+ U(ceres::cos(EulAngles[1]) 
			* ceres::cos(EulAngles[2])
			* ceres::sin(EulAngles[0]));
		Rmat(2,0) = -U(ceres::cos(EulAngles[2])
			* ceres::sin(EulAngles[1]));
		//2nd col
		Rmat(0,1) = -U(ceres::cos(EulAngles[2]) 
			* ceres::sin(EulAngles[0])) 
			- U(ceres::cos(EulAngles[1]) 
			* ceres::cos(EulAngles[0]) 
			* ceres::sin(EulAngles[2]));
		Rmat(1,1) = U(ceres::cos(EulAngles[0]) 
			* ceres::cos(EulAngles[2])) 
			- U(ceres::cos(EulAngles[1])
			* ceres::sin(EulAngles[0])
			* ceres::sin(EulAngles[2]));
		Rmat(2,1) = U(ceres::sin(EulAngles[1])
			* ceres::sin(EulAngles[2]));
		//3rd col
		Rmat(0,2) = U(ceres::cos(EulAngles[0])
			* ceres::sin(EulAngles[1]));
		Rmat(1,2) = U(ceres::sin(EulAngles[0])
			* ceres::sin(EulAngles[1]));
		Rmat(2,2) = U(ceres::cos(EulAngles[1]));

}

//(3).	get Euler Parameters from RotMat
template <typename U>
inline void EulAnglesFromRmat	(U* EulAngles 
	,	Eigen::Matrix<U, 3, 3>	&Rmat)	{
	
	EulAngles[0] = ceres::atan2(Rmat(1,2), Rmat(0,2));

	EulAngles[1] = ceres::atan2(ceres::sqrt(Rmat(0,2) *
		Rmat(0,2) + Rmat(1,2) * Rmat(1,2)), Rmat(2,2));

	EulAngles[2] = ceres::atan2(Rmat(2,1), -Rmat(2,0));
}

/*(4). Spatial Jacobian*/
template <typename U>
inline void spatialJacobianMatTo (
	const Eigen::Matrix <U, 4, 4> &fkMat
	, const Eigen::Matrix<U, 3, 3>	&RmatBase
	,	const Eigen::Matrix <U, 6, Eigen::Dynamic> &basicJMat
	, U* spatialJptr )	{
  Eigen::Map<Eigen::Matrix<U, 6, 13> >
		spatialJmat(spatialJptr);

		//column 1
		spatialJmat.block(0,0,3,3) 
			= Eigen::Matrix<U, 3, 3>::Identity(3,3);
		spatialJmat.block(3,0,3,3) 
			= Eigen::Matrix<U, 3, 3>::Zero(3,3);
		//column 2
		//	skew-symmetric [(b)O(n-b)]x
		Eigen::Matrix<U, 3, 3> skewSymMatOrgN 
				=	Eigen::Matrix<U, 3, 3>::Zero(3,3);
		skewSymMatOrgN(1,0) = fkMat(2,3);
		skewSymMatOrgN(2,0) = -fkMat(1,3);
		skewSymMatOrgN(0,1) = -fkMat(2,3);
		skewSymMatOrgN(2,1) = fkMat(0,3);
		skewSymMatOrgN(0,2) = fkMat(1,3);
		skewSymMatOrgN(0,2) = -fkMat(0,3);

		spatialJmat.block(0,3,3,3) =	
			-RmatBase * skewSymMatOrgN;
		spatialJmat.block(3,3,3,3) =	RmatBase;
		//column 3
		spatialJmat.block(0,6,3,7) = 
			RmatBase * basicJMat.topRows(3);
		spatialJmat.block(3,6,3,7) = 
			RmatBase * basicJMat.bottomRows(3);

}


//(5). B(α) - map EulAngles to anglular velocity;	
//		 if α = [φ, θ, ψ]^T => ω = B(α)α̇
template <typename U>
inline void mapFncEulAnglesToAngularVel (
	const U* const EulAngles, U* mappingMat_Ptr)	{

  Eigen::Map<Eigen::Matrix<U, 3, 3> >
	 	mappingMat(mappingMat_Ptr);

	//	1st col
	mappingMat(0, 0) =
 		U(ceres::cos(EulAngles[2]) * ceres::sin(EulAngles[1]));
	mappingMat(1,0) = U(ceres::sin(EulAngles[2])
		*ceres::sin(EulAngles[1]));
	mappingMat(2,0) = 
		U(ceres::cos(EulAngles[1]));
	//			2nd col
	mappingMat(0,1) = 
		-U(ceres::sin(EulAngles[2]));
	mappingMat(1,1) = 
		U(ceres::cos(EulAngles[2]));
	//			3rd col
	mappingMat(2,2) = U(1.0);

}
// ω = B(α)α̇
template <typename U>
inline void angularVel (const U* const EulAnglesA
	, const U* const EulAnglesB
	, Eigen::Matrix<U, 3, 1> *omegaVel_ptr)	{

	//B(α)
	U* mappingMat_ptr = new U[9];
	mapFncEulAnglesToAngularVel (EulAnglesA, mappingMat_ptr);
  Eigen::Map< Eigen::Matrix<U, 3, 3> >
	 	mappingMat(mappingMat_ptr);

	//α̇
	Eigen::Matrix<U,3,1> EulAnglesDot = 
		Eigen::Matrix<U,3,1>::Identity(3,1);
	EulAnglesDot(0,0) = EulAnglesA[0] - EulAnglesB[0];
	EulAnglesDot(1,0) = EulAnglesA[1] - EulAnglesB[1];
	EulAnglesDot(2,0) = EulAnglesA[2] - EulAnglesB[2];
	//ω = B(α)α̇
//  Eigen::Map< Eigen::Matrix<U, 3, 1> >
//	 	omegaVel(omegaVel_ptr);
	omegaVel_ptr->template block<3, 1>(0, 0) =
		mappingMat * EulAnglesDot;

}
//(6). Linear Velocity
template <typename U>
inline void diffLinearVel (const Eigen::Matrix <U, 3, 1> &pA
	, const Eigen::Matrix <U, 3, 1> &pB, U timeFrame
	, U* linearVel_ptr)	{
	  Eigen::Map< Eigen::Matrix<U, 3, 1> >
			linearVel(linearVel_ptr);
		linearVel.template block<3, 1>(0, 0) =
			(pA - pB) / timeFrame;
}

//compute the (4*4) DH Transformation matrix
//	A. Eigen Version
template <typename U>
inline Eigen::Matrix <U, 4, 4> DHTransformTo	(
	int toDof, const U* const thetas_i)	{

	Eigen::Matrix<U, 7, 1> 
		_a(7,1), _alpha(7,1), _d(7,1);
	_a << U(0.0), U(0.0), U(0.045), U(-0.045)
		, U(0.0), U(0.0), U(0.0);
	_alpha << U(-M_PI/2.0), U(M_PI/2.0), U(-M_PI/2.0)
		, U(M_PI/2.0), U(-M_PI/2.0), U(M_PI/2.0), U(0.0);
	_d << U(0.0), U(0.0), U(0.55), U(0.0), U(0.3)
		, U(0.0), U(0.06);
	
	Eigen::Matrix<U, 4, 4> DHmat =
		Eigen::Matrix<U, 4, 4>::Identity(4,4);
	Eigen::Matrix	<U, 4, 4> fkMat =
		Eigen::Matrix<U, 4, 4>::Identity(4,4);
	for(int i = 0; i < toDof; i++) {
		//1st row 
		DHmat(0,0) = 
			U(ceres::cos(thetas_i[i]));
		DHmat(0,1) = U(-ceres::cos(_alpha(i,0)) 
			* ceres::sin(thetas_i[i]));
		DHmat(0,2) = U(ceres::sin(_alpha(i,0)) 
			* ceres::sin(thetas_i[i]));
		DHmat(0,3) = 
			U(_a(i,0) * ceres::cos(thetas_i[i]));
		//2nd row 
		DHmat(1,0) = 
			U(ceres::sin(thetas_i[i]));
		DHmat(1,1) = U(ceres::cos(_alpha(i,0)) 
			* ceres:: cos(thetas_i[i]));
		DHmat(1,2) = U(-ceres::cos(thetas_i[i]) 
			* ceres::sin(_alpha(i,0)));
		DHmat(1,3) = 
			U(_a(i,0) * ceres::sin(thetas_i[i]));
		//3rd row 
		DHmat(2,0) = U(0.0);
		DHmat(2,1) = U(ceres::sin(_alpha(i,0)));
		DHmat(2,2) = U(ceres::cos(_alpha(i,0)));
		DHmat(2,3) = U(_d(i,0));
		//4th row 
		DHmat(3,0) = U(0.0);
		DHmat(3,1) = U(0.0);
		DHmat(3,2) = U(0.0);
		DHmat(3,3) = U(1.0);

		fkMat = fkMat * DHmat;
	}

}

//	B. Ceres version
template <typename U>
inline void dhTransform	(const U a, const U alpha,
	const U d, const U theta, U result[16])	{
	//1st Row 
	result[0] = U(cos(theta));
	result[1] = U(-cos(alpha) * sin(theta));
	result[2] = U(sin(alpha) * sin(theta));
	result[3] = a* U(cos(theta));
	//2nd Row 
	result[4] = U(sin(theta));
	result[5] = U(cos(alpha) * cos(theta));
	result[6] = U(-cos(theta) * sin(alpha));
	result[7] = a * U(sin(theta));
	//3rd Row 
	result[8] = U(0.0);
	result[9] = U(sin(alpha));
	result[10] = U(cos(alpha));
	result[11] = d;
	//4th Row, i.e. homogonous row
	result[12] = U(0.0);
	result[13] = U(0.0);
	result[14] = U(0.0);
	result[15] = U(1.0);
}
/*
1.position constraint
2. wrist orientation constraint
3. linear velocity constraint 
4. angular velocity constraint
*/
/*
Step 1: Find θ1, θ2, θ3, θ4 such that the 
wrist center oc has coordinates given by:
oc = o − d7.R.[0 0 1]^T					(3.89)
*/

struct DifferentialPositionConstraintFixedBase {
	DifferentialPositionConstraintFixedBase (
		ForwardKin<double> &fk
		// observed scaled Vectors (Cartesian Positions)
		// on the wam at t_i
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
	bool operator()(const U* const thetas_i
		, U* residuals) const {

		// (1). Form the homogeneous transformation matrices
		// 			A_i by substituting the parameters for each 
		//			link into DH-matrix of eqn. (3.10).
		//first 5 joints determine position constraints
		int toDof = 5;   
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

			
			// (2). Form T^n_0 = A 1 · · · A n . This then gives 
			//			the position (and orientation for wrist later)
			//			of the tool frame expressed in base .
			//			coordinates

			//	"A1 * A2 * vUA": -> Elbow center, pUa(i)				
			if (i < 2)
				tUa = tUa * A_j;
			//	"A1 * A2 * A3 * Origin": -> pEl(i)
			//	elbowOffset(i) = pUa(i)	- pEl(i)
			if (i < 3)
				tEl = tEl * A_j;
			if (i < 4)
				tLa = tLa * A_j;
				//"A1 * A2 * ...* A5 * Origin": -> Wrist center
				//pWr (center) = o5 = o6 OR (ALTERNATIVE)
				// = o4 + d5.R4.[0 0 1]^T			
			if (i < 5)
				tWr = tWr * A_j;
		}


		//(3).	estimate marker poitions w.r.t. FIXED base
		//	"θ1, θ2": -> shoToUa(i)				
		Eigen::Matrix<U, 4, 1> shToUaTrans;
		shToUaTrans << U(0.0), U(0.0), U(0.55), U(1.0);
		Eigen::Matrix	<U, 3, 1> pUaEstimatedInBase = 
			(tUa * shToUaTrans
				).template block<3, 1>(0, 0);//.hnormalized();

		//	"θ1, θ2, θ3": -> elbowOffset(i)				
		Eigen::Matrix	<U, 3, 1> pElEstimatedInBase =
			tEl.template block<3, 1>(0, 3);
		//	"θ1, θ2, θ3, θ4": -> WrOffset(i)				
		Eigen::Matrix	<U, 3, 1> pLaEstimatedInBase =
			tLa.template block<3, 1>(0, 3);
		//	"θ1, θ2, θ3 θ4": -> wrist CENTER pWr(i)
		Eigen::Matrix	<U, 3, 1> pWrEstimatedInBase =
			tWr.template block<3, 1>(0, 3);

		
	    //(5).	Compute the residuals.
	    // [ position         ]   [ delta_p    ]
	    Eigen::Map<Eigen::Matrix<U, 12, 1> >
			residualsMat(residuals);

		//	"shToUa(i)"
		residualsMat.template block<3, 1>(0, 0) =
			pUaEstimatedInBase - _shToUaVi.template cast<U>();
		//elbowOffset(i)		
		residualsMat.template block<3, 1>(3, 0) =
			pElEstimatedInBase - pUaEstimatedInBase -
				_elOffsetVi.template cast<U>();

		//WristOffset(i)		
		residualsMat.template block<3, 1>(6, 0) =
			pLaEstimatedInBase - pElEstimatedInBase -
				_wrOffsetVi.template cast<U>();

		//	"LaToWr(i)"
		residualsMat.template block<3, 1>(9, 0) = 
			pWrEstimatedInBase - pLaEstimatedInBase - 
				_laToWrVi.template cast<U>();

		//(5). 	Scale the residuals by the
		//	measurement uncertainty.

		//	"iMatShToUa(i)"
    		residualsMat.template block<3, 1>(0, 0).applyOnTheLeft(
			_iMatShToUaWam.template cast<U>());

		//	"iMatElOffset(i)"
    		residualsMat.template block<3, 1>(3, 0).applyOnTheLeft(
			_iMatElOffsetWam.template cast<U>());

		//	"iMatWrOffset(i)"
    		residualsMat.template block<3, 1>(6, 0).applyOnTheLeft(
			_iMatWrOffsetWam.template cast<U>());

		//	"iMatLaToWr(i)"
		residualsMat.template block<3, 1>(9, 0).applyOnTheLeft(
			_iMatLaToWrWam.template cast<U>()); 

    		return true;
	}


	static ceres::CostFunction* Create(ForwardKin<double> &fk
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
		, const Eigen::Matrix<double, 3, 3>& iMatShToWrWam)	{

		return (new ceres::AutoDiffCostFunction
			<DifferentialPositionConstraintFixedBase, 12, 7> (
				new DifferentialPositionConstraintFixedBase (
					fk, shToUaVi, elOffsetVi, wrOffsetVi
					, laToWrVi, shoToWrVi, iMatShToUaWam
					, iMatElOffsetWam, iMatWrOffsetWam
					, iMatLaToWrWam, iMatShToWrWam)));
	}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	ForwardKin<double> &_fk;
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






struct DifferentialPositionConstraint	{
	DifferentialPositionConstraint (
		ForwardKin<double> &fk
		// observed scaled Vectors (Cartesian Positions)
		// on the wam at t_i
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
	bool operator()(const U* const qBase_i
		, const U* const thetas_i
		, U* residuals) const {

		// (1). Form the homogeneous transformation matrices
		// 			A_i by substituting the parameters for each 
		//			link into DH-matrix of eqn. (3.10).
		//first 5 joints determine position constraints
		int toDof = 5;   
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

			
			// (2). Form T^n_0 = A 1 · · · A n . This then gives 
			//			the position (and orientation for wrist later)
			//			of the tool frame expressed in base .
			//			coordinates

			//	"A1 * A2 * vUA": -> Elbow center, pUa(i)				
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


		//(3).	estimate marker poitions w.r.t. moving base
		//	"θ1, θ2": -> shoToUa(i)				
		Eigen::Matrix<U, 4, 1> shToUaTrans;
		shToUaTrans << U(0.0), U(0.0), U(0.55), U(1.0);
		Eigen::Matrix	<U, 3, 1> pUaEstimatedInBase = 
			(tUa * shToUaTrans
				).template block<3, 1>(0, 0);//.hnormalized();

		//	"θ1, θ2, θ3": -> elbowOffset(i)				
		Eigen::Matrix	<U, 3, 1> pElEstimatedInBase =
			tEl.template block<3, 1>(0, 3);
		//	"θ1, θ2, θ3, θ4": -> WrOffset(i)				
		Eigen::Matrix	<U, 3, 1> pLaEstimatedInBase =
			tLa.template block<3, 1>(0, 3);
		//	"θ1, θ2, θ3 θ4": -> wrist CENTER pWr(i)
		Eigen::Matrix	<U, 3, 1> pWrEstimatedInBase =
			tWr.template block<3, 1>(0, 3);

		//(4).	represent estimated marker poitions in
		//			the camera frame
		//	void QuaternionRotatePoint(
		//		const T q[4], const T pt[3], T result[3])
    // We use QuaternionRotatePoint as it does not assume 	
		// that the quaternion is normalized, since one of the 
    // ways to run the bundle adjuster is to let Ceres 
    // optimize all 4 quaternion parameters without a
    //  local parameterization.

		//	"shToUa(i)"
		Eigen::Matrix	<U, 3, 1> shToUaViEstimated;
    QuaternionRotatePoint(qBase_i,
			pUaEstimatedInBase.data(), shToUaViEstimated.data());
		//	"shToEl(i)"
		Eigen::Matrix	<U, 3, 1> shToElViEstimated;
    QuaternionRotatePoint(qBase_i,
			pElEstimatedInBase.data(), shToElViEstimated.data());
		//	"shToLa(i)"
		Eigen::Matrix	<U, 3, 1> shToLaViEstimated;
    QuaternionRotatePoint(qBase_i,
			pLaEstimatedInBase.data(), shToLaViEstimated.data());
		//	"shToWr(i)"
		Eigen::Matrix	<U, 3, 1> shToWrViEstimated;
    QuaternionRotatePoint(qBase_i,
			pWrEstimatedInBase.data(), shToWrViEstimated.data());

    //(5).	Compute the residuals.
    // [ position         ]   [ delta_p        ]
    Eigen::Map<Eigen::Matrix<U, 12, 1> >
			residualsMat(residuals);

		//	"shToUa(i)"
		residualsMat.template block<3, 1>(0, 0) =
			shToUaViEstimated - _shToUaVi.template cast<U>();
		//elbowOffset(i)		
		residualsMat.template block<3, 1>(3, 0) =
			shToElViEstimated - shToUaViEstimated -
					_elOffsetVi.template cast<U>();

		//WristOffset(i)		
		residualsMat.template block<3, 1>(6, 0) =
			shToLaViEstimated - shToElViEstimated -
					_wrOffsetVi.template cast<U>();

		//	"LaToWr(i)"
		residualsMat.template block<3, 1>(9, 0) = 
			shToWrViEstimated - shToLaViEstimated - 
				_laToWrVi.template cast<U>();

		//(5). 	Scale the residuals by the
		//			measurement uncertainty.

		//	"iMatShToUa(i)"
    residualsMat.template block<3, 1>(0, 0).applyOnTheLeft(
			_iMatShToUaWam.template cast<U>());

		//	"iMatElOffset(i)"
    residualsMat.template block<3, 1>(3, 0).applyOnTheLeft(
			_iMatElOffsetWam.template cast<U>());

		//	"iMatWrOffset(i)"
    residualsMat.template block<3, 1>(6, 0).applyOnTheLeft(
			_iMatWrOffsetWam.template cast<U>());

		//	"iMatLaToWr(i)"
		residualsMat.template block<3, 1>(9, 0).applyOnTheLeft(
			_iMatLaToWrWam.template cast<U>()); 

    return true;
	}

	static ceres::CostFunction* Create(ForwardKin<double> &fk
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
		, const Eigen::Matrix<double, 3, 3>& iMatShToWrWam)	{

		return (new ceres::AutoDiffCostFunction
			<DifferentialPositionConstraint, 12, 4, 7> (
				new DifferentialPositionConstraint (fk, shToUaVi,
					elOffsetVi, wrOffsetVi,	laToWrVi, shoToWrVi,
					iMatShToUaWam, iMatElOffsetWam, iMatWrOffsetWam, 						iMatLaToWrWam, iMatShToWrWam)));
	}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	ForwardKin<double> &_fk;
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



struct DifferentialOrientationConstraint	{
	DifferentialOrientationConstraint (
		ForwardKin<double> &fk
		, Eigen::Quaternion<double> qbase_i
		// The third column of R expresses the direction 
		//	of z7 with respect to the base frame.
		// Oee = Owr + d7.Ree.[0 0 1]^T			
		, const Eigen::Matrix<double,3,1> shoToEeVi
		, const Eigen::Matrix<double, 3, 3>& iMatShToEeWam
		) : 
		_fk(fk)
		, _qbase_i(qbase_i)
		, _AmatJ1toJ4(fk.getMat(4))
		,_a(_fk.getAvect(7))
		, _alpha(_fk.getAlphaVect(7))
		, _d(_fk.getDvect(7)) 
		,	_shoToEeVi(shoToEeVi)
		, _iMatShToEeWam(iMatShToEeWam)
	{}


		//oUA = o3 - a3.R3.[1 0 0]^T
	template <typename U>
	bool operator()( const U* const thetas_i
		, U* residuals) const {

		// (1). Form the homogeneous transformation matrices
		// 			A_i by substituting the parameters for each 
		//			link into DH-matrix of eqn. (3.10).
		//first 4 joints determine position constraints
		int toDof = 4,
			nJoints = 7;   
		// last 3 joints determine endeffector orientation.
		//initialize the matrices to zero;
		//avoid calling external eigen::identity function
		Eigen::Matrix<U, 4, 4> tEe;
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) 
				tEe(i,j) = U(0.0);	
			tEe(i,i) = U(1.0);
		}		

		Eigen::Matrix<U, 4, 4> A_j;
		for(int i = toDof; i < nJoints; i++) {
			//1st row 
			A_j(0,0) = 
				U(ceres::cos(thetas_i[i-toDof]));
			A_j(0,1) = U(-ceres::cos(_alpha(i,0)) 
				* ceres::sin(thetas_i[i-toDof]));
			A_j(0,2) = U(ceres::sin(_alpha(i,0)) 
				* ceres::sin(thetas_i[i-toDof]));
			A_j(0,3) = 
				U(_a(i,0) * ceres::cos(thetas_i[i-toDof]));
			//2nd row 
			A_j(1,0) = 
				U(ceres::sin(thetas_i[i-toDof]));
			A_j(1,1) = U(ceres::cos(_alpha(i,0)) 
				* ceres:: cos(thetas_i[i-toDof]));
			A_j(1,2) = U(-ceres::cos(thetas_i[i-toDof]) 
				* ceres::sin(_alpha(i,0)));
			A_j(1,3) = 
				U(_a(i,0) * ceres::sin(thetas_i[i-toDof]));
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

			
			/** (2). Form T^n_0 = A 1 · · · A n . 
								This then gives the position (and 
								orientation for wrist later) of the 
								tool frame expressed in base frame:
			**/
			tEe = tEe * A_j;
		}
		tEe = _AmatJ1toJ4.template cast<U>()
			 * tEe;

		//(3).	estimate Ee poition w.r.t. moving base
		//	"θ1, θ2, ...., θ7": -> endEffector pEe(i)		
		Eigen::Matrix	<U, 3, 1> pEeEstimatedInBase =
			tEe.template block<3, 1>(0, 3);


		//(4).	represent estimated endeffector poitions in
		//			the camera frame
		//	void QuaternionRotatePoint(
		//		const T q[4], const T pt[3], T result[3])
		Eigen::Matrix	<U, 3, 1> shToEeViEstimated;

		U q_b_quat[4];
		q_b_quat[0] = U(_qbase_i.x()); 
		q_b_quat[1] = U(_qbase_i.y()); 
		q_b_quat[2] = U(_qbase_i.z()); 
		q_b_quat[3] = U(_qbase_i.w()); 
    QuaternionRotatePoint(q_b_quat,
			pEeEstimatedInBase.data(), shToEeViEstimated.data());

    //(5).	Compute the residuals.
    // [ position         ]   [ delta_p        ]
    Eigen::Map<Eigen::Matrix<U, 3, 1> >
			residualsMat(residuals);
		residualsMat.template block<3, 1>(0, 0) =
			shToEeViEstimated - _shoToEeVi.template cast<U>();

		//(5). 	Scale the residuals by the
		//			measurement uncertainty.
    residualsMat.template block<3, 1>
    	(0, 0).applyOnTheLeft(
				_iMatShToEeWam.template cast<U>());

    return true;
	}

	static ceres::CostFunction* Create(ForwardKin<double> &fk
		, Eigen::Quaternion<double> qbase_i
		, const Eigen::Matrix<double,3,1> shoToEeVi
		, const Eigen::Matrix<double, 3, 3>& iMatShToEeWam
		)	{
		return (new ceres::AutoDiffCostFunction
			<DifferentialOrientationConstraint, 3, 3> (
				new DifferentialOrientationConstraint (fk, qbase_i,
					shoToEeVi, iMatShToEeWam)));
	}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	ForwardKin<double> &_fk;
//	const double* _qbase_i;
	Eigen::Quaternion<double> _qbase_i;
	Eigen::Matrix<double, 4,4> _AmatJ1toJ4;
	Eigen::Matrix<double, Eigen::Dynamic, 1> 
	 	_a, _alpha, _d;
	//measured/scaled Cartesian vectores on the wam at t_i
	const Eigen::Matrix<double, 3, 1> _shoToEeVi;
	//measurements covariance (x y z)
	const Eigen::Matrix<double, 3, 3> _iMatShToEeWam;
};





class DifferentialIKErrorTerm {

public:
	DifferentialIKErrorTerm (
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
			<DifferentialIKErrorTerm, 9, 3, 3, 7, 7> (
				new DifferentialIKErrorTerm (fk, inPtsI
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


