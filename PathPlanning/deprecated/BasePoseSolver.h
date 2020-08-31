/* 
IK CHAIN + CAMERA T
	Given inPts tracked camera markers, 
	the goal is to find cameraT (3 orient. 3 trans)
	that	minimize the re-Transformation error. 
*/
#ifndef BASE_POSE_SOLVER_H
#define BASE_POSE_SOLVER_H

#include "BasicFormulas.h"
#include "ForwardKin.h"
#include "WamIkProblem.h"
#include "WamIkPoseSolver.h"
#include <UBCUtil.h>

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "gflags/gflags.h"
#include "glog/logging.h"

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

namespace ceres {


template <typename T>
struct BasePoseProblem {
	
	BasePoseProblem( 
		ForwardKin<T> &fk
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &JsWam
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> shoToUaVsWam
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> shoToElVsWam
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> shoToLaVsWam
		, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> shoToWrVsWam
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles

		) :
		_fk(fk)
		, _JsWam(JsWam)
		, _shoToUaVsWam(shoToUaVsWam)
		, _shoToElVsWam(shoToElVsWam)
		, _shoToLaVsWam(shoToLaVsWam)
		, _shoToWrVsWam(shoToWrVsWam)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		,_a(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _alpha(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1))
		, _d(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(7,1)) 
		, _startConfigJoints(new T[getPtsBlockSize()])
		, _startConfigCamera(new T[getCameraBlockSize()])
	{
		
		_a << T(0.0), T(0.0), T(0.045), T(-0.045), T(0.0)
				, T(0.0), T(0.0);
		_alpha << T(-M_PI/2.0), T(M_PI/2.0), T(-M_PI/2.0)
					 , T(M_PI/2.0), T(-M_PI/2.0), T(M_PI/2.0), T(0.0);
		_d << T(0.0), T(0.0), T(0.55), T(0.0), T(0.3)
			 , T(0.0), T(0.06);
	}


const T* getPtI (
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPts
	, int index)	{
	T* pt = new T[3];
	pt[0] = inPts(0,index);
	pt[1] = inPts(1,index);
	pt[2] = inPts(2,index);
	return pt;
}
const T getNpts()	{
	return _shoToUaVsWam.cols();
}

const T getPtsBlockSize()	{
	const int  nJoints = _fk.getNJoints();
	return  nJoints * (_polynomialOrder + 1);
}

T* getStartConfigJoints()	{
	return _startConfigJoints;
}

const T getCameraBlockSize()	{
	//3 for rotation, 3 for trans
	const int cameraBlockSize = 6; 
	return cameraBlockSize;
}

T* getStartConfigCamera()	{
	return _startConfigCamera;
}

struct BasePoseReprojectionError {
	BasePoseReprojectionError (T* shoToElVi
		, T* elOffsetVi , T* wrOffsetVi
		, T* elToWrVi, T parametricTerm
		, std::vector<T> jointMinAngles
		, std::vector<T> jointMaxAngles)
		: _shoToElVi(shoToElVi)
		, _elOffsetVi(elOffsetVi)
		, _wrOffsetVi(wrOffsetVi)
		, _elToWrVi(elToWrVi)
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

	// Each Residual block takes 4 marker points and a camera as input
	// and outputs a ???? dimensional residual.
	template <typename U>
	bool operator()(const U* const camera
		, const U* const candidateParamJs
		, U* residuals) const {


		U *thetas = new U[7];	
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
		fkMat =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4), dhTransform(4, 4);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
		fkMatEl =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4);

		//calculate thetas
		//if joints within limit, set to wam		
		bool retVal = true;
		int toDof = 3; 
		for(size_t i = 0; i < toDof; i++) {
			size_t iTimesFour = i * 4;
	
			thetas[i] = candidateParamJs[iTimesFour] 
				+	U(_parametricTerm) * candidateParamJs[iTimesFour + 1] 
				+	U(_parametricTerm * _parametricTerm) 
					* candidateParamJs[iTimesFour + 2] 
				+	U(_parametricTerm * _parametricTerm * _parametricTerm) 
					* candidateParamJs[iTimesFour + 3];
	
				retVal =  retVal && (thetas[i] >= _jointMinAngles[i]);
				retVal =  retVal && (thetas[i] <= _jointMaxAngles[i]);
			}

			if (retVal) {
			
				for(size_t i = 0; i < toDof; i++) {
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
				
					if (i < 2)
						fkMat = fkMat * dhTransform;
				
					if (i < 3)
						fkMatEl = fkMatEl * dhTransform;
				}
			}

		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			origin(4,1);
		origin << U(0.0), U(0.0), U(0.0), U(1.0);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			wamEltoUaTrans(4,1);
		wamEltoUaTrans << U(0.0), U(0.0), U(0.55), U(1.0);
		
		//predicting each ELBOW-NO_OFFSET marker position wrt the base  
		//of the robot
		U *predShoToElViInBase = new U[3];	
		predShoToElViInBase[0] = 
			U(fkMat(0,0)) * U(wamEltoUaTrans(0,0)) 
			+ U(fkMat(0,1)) * U(wamEltoUaTrans(1,0)) 
			+ U(fkMat(0,2)) * U(wamEltoUaTrans(2,0)) 
			+ U(fkMat(0,3)) * U(wamEltoUaTrans(3,0));

		predShoToElViInBase[1] = 
			U(fkMat(1,0)) * U(wamEltoUaTrans(0,0)) 
			+ U(fkMat(1,1)) * U(wamEltoUaTrans(1,0)) 
			+ U(fkMat(1,2)) * U(wamEltoUaTrans(2,0)) 
			+ U(fkMat(1,3)) * U(wamEltoUaTrans(3,0)); 
		predShoToElViInBase[2] =
			U(fkMat(2,0)) * U(wamEltoUaTrans(0,0)) 
			+ U(fkMat(2,1)) * U(wamEltoUaTrans(1,0)) 
			+ U(fkMat(2,2)) * U(wamEltoUaTrans(2,0)) 
			+ U(fkMat(2,3)) * U(wamEltoUaTrans(3,0)); 

	//ROTATE (Vi IN BASE) AROUND CAMERA AXIS TO GET TO
	// (Vi) IN CAMERA FRAME
	// camera[0,1,2] are the angle-axis rotation.

		U *observedShoToElVi = new U[3];
		observedShoToElVi[0] = U(_shoToElVi[0]);
		observedShoToElVi[1] = U(_shoToElVi[1]);
		observedShoToElVi[2] = U(_shoToElVi[2]);

    U predShoToElVi[3];
    ceres::AngleAxisRotatePoint(camera, predShoToElViInBase
			, predShoToElVi);   

		residuals[0] = 
			U(_shoToElVi[0]) - U(predShoToElVi[0]);
		residuals[1] = 
			U(_shoToElVi[1]) - U(predShoToElVi[1]);
		residuals[2] = 
			U(_shoToElVi[2]) - U(predShoToElVi[2]);
/*
		residuals[3] = U(0.55 * 0.55) 
			- U(predShoToElVi[0] * predShoToElVi[0]
			+ predShoToElVi[1] * predShoToElVi[1]
			+ predShoToElVi[2] * predShoToElVi[2]);

	//predicting each ELBOW marker position wrt the base  
		//of the robot

		U *predElPt = new U[3];	
		predElPt[0] = U(fkMatEl(0,0)) * U(origin(0,0)) 
			+ U(fkMatEl(0,1)) * U(origin(1,0)) 
			+ U(fkMatEl(0,2)) * U(origin(2,0)) 
			+ U(fkMatEl(0,3)) * U(origin(3,0));

		predElPt[1] = U(fkMatEl(1,0)) * U(origin(0,0)) 
			+ U(fkMatEl(1,1)) * U(origin(1,0)) 
			+ U(fkMatEl(1,2)) * U(origin(2,0)) 
			+ U(fkMatEl(1,3)) * U(origin(3,0)); 
		predElPt[2] = U(fkMatEl(2,0)) * U(origin(0,0)) 
			+ U(fkMatEl(2,1)) * U(origin(1,0)) 
			+ U(fkMatEl(2,2)) * U(origin(2,0)) 
			+ U(fkMatEl(2,3)) * U(origin(3,0)); 

    U pEl[3];
		U *observedElPt = new U[3];
		observedElPt[0] = U(_observedEl[0]);
		observedElPt[1] = U(_observedEl[1]);
		observedElPt[2] = U(_observedEl[2]);
    ceres::AngleAxisRotatePoint(camera, observedElPt, pEl);
    // camera[3,4,5] are the translation.
//    pEl[0] += camera[3]; pEl[1] += camera[4]; pEl[2] += camera[5];


    U pWr[3];
		U *observedWrPt = new U[3];
		observedWrPt[0] = U(_observedWr[0]);
		observedWrPt[1] = U(_observedWr[1]);
		observedWrPt[2] = U(_observedWr[2]);
    ceres::AngleAxisRotatePoint(camera, observedWrPt, pWr);

		residuals[3] = U(predElPt[0]) - U(pEl[0]);
		residuals[4] = U(predElPt[1]) - U(pEl[1]);
		residuals[5] = U(predElPt[2]) - U(pEl[2]);





*/	


	return true;
	}


	// Factory to hide the construction of the CostFunction object from
	// the client code.
		static ceres::CostFunction* Create(const double* observedUa
			, const double* observedEl, const double* observedLa
			, const double* observedWr, const double parametricTerm
			, std::vector<double> jointMinAngles
			, std::vector<double> jointMaxAngles) {			
			return (new ceres::AutoDiffCostFunction
				<BasePoseReprojectionError,3,3,28> (
					new BasePoseReprojectionError ( observedUa, observedEl
						, observedLa, observedWr, parametricTerm
						, jointMinAngles, jointMaxAngles)));
		}

	T *_shoToElVi, *_elOffsetVi, *_wrOffsetVi, *_elToWrVi;
	T _parametricTerm;
	std::vector<double> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	 	_a, _alpha, _d;
};

/*
struct BasePoseReprojectionError {
	BasePoseReprojectionError (double* observedUa
		, double* observedEl , double* observedLa
		, double* observedWr, double parametricTerm
		, std::vector<double> jointMinAngles
		, std::vector<double> jointMaxAngles)
		: _observedUa(observedUa)
		, _observedEl(observedEl)
		, _observedLa(observedLa)
		, _observedWr(observedWr)
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

	// Each Residual block takes 4 marker points and a camera as input
	// and outputs a ???? dimensional residual.
	template <typename U>
	bool operator()(const U* const camera
		, const U* const candidateParamJs
		, U* residuals) const {


		U *thetas = new U[7];	
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
		fkMat =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4), dhTransform(4, 4);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
		fkMatEl =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4);

		//calculate thetas
		//if joints within limit, set to wam		
		bool retVal = true;
		for(size_t i = 0; i < 7; i++) {
			size_t iTimesFour = i * 4;
	
			thetas[i] = candidateParamJs[iTimesFour] 
				+	U(_parametricTerm) * candidateParamJs[iTimesFour + 1] 
				+	U(_parametricTerm * _parametricTerm) 
					* candidateParamJs[iTimesFour + 2] 
				+	U(_parametricTerm * _parametricTerm * _parametricTerm) 
					* candidateParamJs[iTimesFour + 3];
	
				retVal =  retVal && (thetas[i] >= _jointMinAngles[i]);
				retVal =  retVal && (thetas[i] <= _jointMaxAngles[i]);
			}

			if (retVal) {
			
				for(size_t i = 0; i < 7; i++) {
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
				
					if (i < 2)
						fkMat = fkMat * dhTransform;
				
					if (i < 3)
						fkMatEl = fkMatEl * dhTransform;
				}
			}

		//predicting each ELBOW-NO_OFFSET marker position wrt the base  
		//of the robot
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			wamEltoUaTrans(4,1);
		wamEltoUaTrans << U(0.0), U(0.0), U(0.55), U(1.0);

		U *predUaPt = new U[3];	
		predUaPt[0] = U(fkMat(0,0)) * U(wamEltoUaTrans(0,0)) 
			+ U(fkMat(0,1)) * U(wamEltoUaTrans(1,0)) 
			+ U(fkMat(0,2)) * U(wamEltoUaTrans(2,0)) 
			+ U(fkMat(0,3)) * U(wamEltoUaTrans(3,0));

		predUaPt[1] = U(fkMat(1,0)) * U(wamEltoUaTrans(0,0)) 
			+ U(fkMat(1,1)) * U(wamEltoUaTrans(1,0)) 
			+ U(fkMat(1,2)) * U(wamEltoUaTrans(2,0)) 
			+ U(fkMat(1,3)) * U(wamEltoUaTrans(3,0)); 
		predUaPt[2] = U(fkMat(2,0)) * U(wamEltoUaTrans(0,0)) 
			+ U(fkMat(2,1)) * U(wamEltoUaTrans(1,0)) 
			+ U(fkMat(2,2)) * U(wamEltoUaTrans(2,0)) 
			+ U(fkMat(2,3)) * U(wamEltoUaTrans(3,0)); 

	// camera[0,1,2] are the angle-axis rotation.
    U pUa[3];
		U *observedUaPt = new U[3];
		observedUaPt[0] = U(_observedUa[0]);
		observedUaPt[1] = U(_observedUa[1]);
		observedUaPt[2] = U(_observedUa[2]);
    ceres::AngleAxisRotatePoint(camera, observedUaPt, pUa);   
		// camera[3,4,5] are the translation.
//    pUa[0] += camera[3]; pUa[1] += camera[4]; pUa[2] += camera[5];

	//predicting each ELBOW marker position wrt the base  
		//of the robot
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			origin(4,1);
		origin << U(0.0), U(0.0), U(0.0), U(1.0);

		U *predElPt = new U[3];	
		predElPt[0] = U(fkMatEl(0,0)) * U(origin(0,0)) 
			+ U(fkMatEl(0,1)) * U(origin(1,0)) 
			+ U(fkMatEl(0,2)) * U(origin(2,0)) 
			+ U(fkMatEl(0,3)) * U(origin(3,0));

		predElPt[1] = U(fkMatEl(1,0)) * U(origin(0,0)) 
			+ U(fkMatEl(1,1)) * U(origin(1,0)) 
			+ U(fkMatEl(1,2)) * U(origin(2,0)) 
			+ U(fkMatEl(1,3)) * U(origin(3,0)); 
		predElPt[2] = U(fkMatEl(2,0)) * U(origin(0,0)) 
			+ U(fkMatEl(2,1)) * U(origin(1,0)) 
			+ U(fkMatEl(2,2)) * U(origin(2,0)) 
			+ U(fkMatEl(2,3)) * U(origin(3,0)); 

    U pEl[3];
		U *observedElPt = new U[3];
		observedElPt[0] = U(_observedEl[0]);
		observedElPt[1] = U(_observedEl[1]);
		observedElPt[2] = U(_observedEl[2]);
    ceres::AngleAxisRotatePoint(camera, observedElPt, pEl);
    // camera[3,4,5] are the translation.
//    pEl[0] += camera[3]; pEl[1] += camera[4]; pEl[2] += camera[5];

		residuals[0] = U(predUaPt[0]) - U(pUa[0]);
		residuals[1] = U(predUaPt[1]) - U(pUa[1]);
		residuals[2] = U(predUaPt[2]) - U(pUa[2]);

		residuals[3] = U(predElPt[0]) - U(pEl[0]);
		residuals[4] = U(predElPt[1]) - U(pEl[1]);
		residuals[5] = U(predElPt[2]) - U(pEl[2]);




	return true;
	}


	// Factory to hide the construction of the CostFunction object from
	// the client code.
		static ceres::CostFunction* Create(const double* observedUa
			, const double* observedEl, const double* observedLa
			, const double* observedWr, const double parametricTerm
			, std::vector<double> jointMinAngles
			, std::vector<double> jointMaxAngles) {			
			return (new ceres::AutoDiffCostFunction
				<BasePoseReprojectionError,6,3,28> (
					new BasePoseReprojectionError ( observedUa, observedEl
						, observedLa, observedWr, parametricTerm
						, jointMinAngles, jointMaxAngles)));
		}

	double *_observedUa, *_observedEl, *_observedLa, *_observedWr;
	double _parametricTerm;
	std::vector<double> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	 	_a, _alpha, _d;
};
*/
	void SetLinearSolver(Solver::Options* options) {
		CHECK(StringToLinearSolverType(FLAGS_linear_solver
			, &options->linear_solver_type));
		CHECK(StringToPreconditionerType(FLAGS_preconditioner
			, &options->preconditioner_type));
		CHECK(StringToVisibilityClusteringType(FLAGS_visibility_clustering
			, &options->visibility_clustering_type));
		CHECK(StringToSparseLinearAlgebraLibraryType(
			FLAGS_sparse_linear_algebra_library
			, &options->sparse_linear_algebra_library_type));
		CHECK(StringToDenseLinearAlgebraLibraryType(
		  FLAGS_dense_linear_algebra_library
		  , &options->dense_linear_algebra_library_type));
		options->use_explicit_schur_complement =
			FLAGS_explicit_schur_complement;
///		options->use_mixed_precision_solves = FLAGS_mixed_precision_solves;
//		options->max_num_refinement_iterations = FLAGS_max_num_refinement_iterations;
	}


	void SetMinimizerOptions(Solver::Options* options) {
		options->max_num_iterations = FLAGS_num_iterations;
		options->minimizer_progress_to_stdout = true;
		options->num_threads = FLAGS_num_threads;
		options->eta = FLAGS_eta;
		options->max_solver_time_in_seconds = FLAGS_max_solver_time;
		options->use_nonmonotonic_steps = FLAGS_nonmonotonic_steps;
		if (FLAGS_line_search) {
		  options->minimizer_type = ceres::LINE_SEARCH;
		}
		CHECK(StringToTrustRegionStrategyType	(
			FLAGS_trust_region_strategy
			, &options->trust_region_strategy_type));
		CHECK(StringToDoglegType(FLAGS_dogleg, &options->dogleg_type));
		options->use_inner_iterations = FLAGS_inner_iterations;
	}

	void SetSolverOptionsFromFlags(
		BasePoseProblem<double>* basePoseProblem
		, Solver::Options* options) {

		SetMinimizerOptions(options);
		SetLinearSolver(options);
		SetOrdering(basePoseProblem,options);
	}

	void SetOrdering(BasePoseProblem<double>* basePoseProblem
		 , Solver::Options* options) {
		const int nPts = basePoseProblem->getNpts();
		const int cameraBlockSize = 
			basePoseProblem->getCameraBlockSize();
		const int ptsBlockSize = 
			basePoseProblem->getPtsBlockSize();
		//double* startConfigJoints = new double[ptsBlockSize];
		//double* startConfigCamera = new double[cameraBlockSize];
		double* startConfigJoints = 
			basePoseProblem->getStartConfigJoints();
		double* startConfigCamera =
			basePoseProblem->getStartConfigCamera();

		if (options->use_inner_iterations) {
		  if (FLAGS_blocks_for_inner_iterations == "cameras") {
		    LOG(INFO) << "Camera blocks for inner iterations";
		    options->inner_iteration_ordering.reset(
					new ParameterBlockOrdering);
		    for (int i = 0; i < nPts; ++i) {
		      options->inner_iteration_ordering->AddElementToGroup
						(startConfigCamera + cameraBlockSize * i, 0);
		    }
		  }
			else if (FLAGS_blocks_for_inner_iterations == "points") {
				LOG(INFO) << "Point blocks for inner iterations";
				options->inner_iteration_ordering.reset(
					new ParameterBlockOrdering);
				for (int i = 0; i < nPts; ++i) {
					options->inner_iteration_ordering->AddElementToGroup
						(startConfigJoints + ptsBlockSize * i, 0);
				}
			}
			else if (FLAGS_blocks_for_inner_iterations =="cameras,points") 
			{
		  	LOG(INFO) << "Camera followed by point blocks for inner iterations";
		    options->inner_iteration_ordering.reset(
					new ParameterBlockOrdering);
		    for (int i = 0; i < nPts; ++i) {
		      options->inner_iteration_ordering->AddElementToGroup
						(startConfigCamera + cameraBlockSize * i, 0);
		    }
		    for (int i = 0; i < nPts; ++i) {
		      options->inner_iteration_ordering->AddElementToGroup
						(startConfigJoints + ptsBlockSize * i, 1);
		    }
		  } 
			else if (FLAGS_blocks_for_inner_iterations =="points,cameras") 
			{
		  	LOG(INFO) << "Point followed by camera blocks for inner iterations";
		    options->inner_iteration_ordering.reset(
					new ParameterBlockOrdering);
		    for (int i = 0; i < nPts; ++i) {
		      options->inner_iteration_ordering->AddElementToGroup
						(startConfigCamera + cameraBlockSize * i, 1);
		    }
		    for (int i = 0; i < nPts; ++i) {
		      options->inner_iteration_ordering->AddElementToGroup
						(startConfigJoints + ptsBlockSize * i, 0);
		    }
		  } 
			else if (FLAGS_blocks_for_inner_iterations == "automatic") {
		  	LOG(INFO) << "Choosing automatic blocks for inner iterations";
		  } 
			else {
		    LOG(FATAL) << "Unknown block type for inner iterations: "
		               << FLAGS_blocks_for_inner_iterations;
		  }
  
		}	
}
	void BuildProblem(BasePoseProblem<double>* basePoseProblem
		, Problem* problem) {

		// Observations for each marker is 3*num_observations long array 
		// observations = [u_1, u_2, ... , u_n], where each u_i is 3 
		//dimensional, the x, y and z positions of the mar
		const int nPts = basePoseProblem->getNpts();
		const int ptsBlockSize = 
			basePoseProblem->getPtsBlockSize();
		//double* startConfigJoints = new double[ptsBlockSize];
		//3 for rotation, 3 for trans
		const int cameraBlockSize = 
			basePoseProblem->getCameraBlockSize();
		//double* startConfigCamera = new double[cameraBlockSize];

		double* startConfigJoints = 
			basePoseProblem->getStartConfigJoints();
		double* startConfigCamera =
			basePoseProblem->getStartConfigCamera();

		const double deltaSize = 1.0/double(nPts);
		const int  nJoints = _fk.getNJoints();

		for(size_t i = 0; i < nPts; i++)	{
			const double parametricTerm = double(i) * deltaSize; 
			CostFunction* cost_function;
		  // Each Residual block takes 7 joint angles, which are 
			// identified by cubic spline coefficients joint angles,
			// and a cameraT as input  and outputs a ???? dimensional 
			// residual.
			const double* observedUaPt = 
				basePoseProblem->getPtI(_shoToUaVsWam, i);
			const double* observedElPt = 
				basePoseProblem->getPtI(_shoToElVsWam, i);
			const double* observedLaPt = 
				basePoseProblem->getPtI(_shoToLaVsWam, i);
			const double* observedWrPt = 
				basePoseProblem->getPtI(_shoToWrVsWam, i);

		  cost_function = BasePoseReprojectionError::Create (
				observedUaPt, observedElPt, observedLaPt
				, observedWrPt, parametricTerm, _jointMinAngles
				, _jointMaxAngles);


		  // If enabled use Huber's loss function.
		  LossFunction* loss_function = FLAGS_robustify ?
				new HuberLoss(1.0) : NULL;
		  // Each observation correponds to a pair of cameraT and 3/4
			//observed markers points which are identified by camera_index
			// and point_index()[i] respectively.
			problem->AddResidualBlock(cost_function
				, loss_function, startConfigCamera
				, startConfigJoints);

		}

		//set the joint limits into ceres 
		for (int joint = 0; joint < nJoints; joint++) 
 			setJointLimit(startConfigJoints, joint, problem);	

}

void vectorMapFromHtoR(
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsAlongRows)
 	{

}
	void SolveProblem(
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inPtsAlongRows)
 	{

		const int nRows = inPtsAlongRows.rows();
		//10 joints, 30 points
		_inPtsRW = inPtsAlongRows.leftCols(3).transpose(),
			_inPtsRTh = inPtsAlongRows.block(0,3,nRows,3).transpose(),
			_inPtsRPk = inPtsAlongRows.rightCols(3).transpose(),
			_inPtsRLA = inPtsAlongRows.block(0,6,nRows,3).transpose(),
			_inPtsREl = inPtsAlongRows.block(0,9,nRows,3).transpose(),
			_inPtsRUA = inPtsAlongRows.block(0,12,nRows,3).transpose(),
			_inPtsRSh = inPtsAlongRows.block(0,15,nRows,3).transpose(),
			_inPtsRCh = inPtsAlongRows.block(0,18,nRows,3).transpose(),
			_inPtsMCh = inPtsAlongRows.block(0,21,nRows,3).transpose(),
			_inPtsLCh = inPtsAlongRows.block(0,24,nRows,3).transpose(),
			_inPtsLSh = inPtsAlongRows.block(0,27,nRows,3).transpose();


		/**(0) Vector mapping from human motion data to to describe 
				the relative relationship between each joint
				instead of the location information.**/
	
		// (a) for human - sholder-to-elbow & -elbow-to-wrist 
		// 			vector in hom coordinates
		MatrixXd shoToUaVsH = MatrixXd::Zero(3,nRows)
			, shoToElVsH = MatrixXd::Zero(3,nRows)
			, elToLaVsH = MatrixXd::Zero(3,nRows)
			, elToWrVsH = MatrixXd::Zero(3,nRows)
			, WrToThVsH = MatrixXd::Zero(3,nRows);
		shoToUaVsH = (_inPtsRUA - _inPtsRSh).colwise().normalized();
		shoToElVsH = (_inPtsREl - _inPtsRSh).colwise().normalized();
		elToLaVsH = (_inPtsRLA - _inPtsREl).colwise().normalized();
		elToWrVsH = (_inPtsRW - _inPtsREl).colwise().normalized();
		WrToThVsH = (_inPtsRTh - _inPtsRW).colwise().normalized();

		double shoToElMeanLH = (_inPtsREl - _inPtsRSh)
			.colwise().norm().mean();
		double elToLaMeanLH = (_inPtsRLA - _inPtsREl)
			.colwise().norm().mean();
		double elToWrMeanLH = (_inPtsRW - _inPtsREl)
			.colwise().norm().mean();
		double WrToThMeanLH = (_inPtsRW - _inPtsRTh)
			.colwise().norm().mean();
		double WrToPkMeanLH = (_inPtsRW - _inPtsRPk)
			.colwise().norm().mean();

		// (b) for the WAM
		MatrixXd shoToUaVsWam = MatrixXd::Zero(3,nRows)
			,	shoToElVsWam = MatrixXd::Zero(3,nRows)
			, elToWrVsWam = MatrixXd::Zero(3,nRows)
			,	shoToLaVsWam = MatrixXd::Zero(3,nRows)
			,	shoToThVsWam = MatrixXd::Zero(3,nRows)
			,	shoToWrVsWam = MatrixXd::Zero(3,nRows);

		// WAM elb and Wr offset Vectors
		MatrixXd elOffsetVsWam = MatrixXd::Zero(3,nRows)
			, wOffsetVsWam = MatrixXd::Zero(3,nRows);
		double shoToElLwam = 0.55;
		double elToWrLwam = 0.3; 
		MatrixXd shoToElCrossEltoWrVsH = MatrixXd::Zero(nRows,3)
			, j4Wam = MatrixXd::Zero(nRows,1);

		//input vectors to CERES solver
		for (size_t i = 0; i < nRows; i++) {
			// (b.1.1) take sholder-to-elbow CROSS elbow-to-wrist
			Vector3d shoToElVh = shoToElVsH.col(i),
				elToWrVh = elToWrVsH.col(i),
				shoToElCrossEltoWrVsH = shoToElVh.cross(
					elToWrVh).normalized();
			// (b.1.2) find Human elbow angle, i.e. <)l0.l1
			double elThetaH = ceres::acos(shoToElVh.dot(elToWrVh));
			j4Wam(i,0) =  M_PI - elThetaH;
		
			// (b.2) Elbow offset Vector
			// x3 = -Z2.cross(Z3)
			elOffsetVsWam.col(i) = -shoToElVh.cross(
			shoToElCrossEltoWrVsH).normalized();

			shoToUaVsWam.col(i) = shoToElVsH.col(i) * (shoToElLwam);
			shoToElVsWam.col(i) = shoToElVsH.col(i) * shoToElLwam
				+ elOffsetVsWam.col(i) * 0.045;

			// (b.3) WRIST offset Vector
			// -x4 = -Y4{Z3}.cross(Z4)
			wOffsetVsWam.col(i) = 	
				-shoToElCrossEltoWrVsH.cross(elToWrVh).normalized();
			elToWrVsWam.col(i) = wOffsetVsWam.col(i) * 0.045
				+ elToWrVsH.col(i) * elToWrLwam;
			shoToLaVsWam.col(i) = shoToElVsWam.col(i) 
				+ wOffsetVsWam.col(i) * 0.045; 
		//		+ elToLaVsH.col(i) * (elToWrLwam/2.0);
			//i.e. to j7 approx MAxLLaToW =(0.3+0.06) when J6=0
			// and MinLLaToW =(0.3) when J6=0
			shoToWrVsWam.col(i) = shoToElVsWam.col(i) 
				+ elToWrVsWam.col(i); 
			shoToThVsWam.col(i) = shoToWrVsWam.col(i)
				+ WrToThVsH.col(i) * 0.06; 
		}
		double shoToElMeanLWam = shoToElVsWam.colwise().norm().mean();
		double elToWrMeanLWam = elToWrVsWam.colwise().norm().mean();
		double wrToThMeanLWam = 0.06;

		// human vs robot arm links 
		cout << "shoToElMeanLH: " << shoToElMeanLH << endl;
		cout << "shoToElMeanLWam: " << shoToElMeanLWam << endl;

		cout << "elToWrMeanLH: " << elToWrMeanLH << endl;
		cout << "elToWrMeanLWam: " << elToWrMeanLWam << endl;

		cout << "WrToThMeanLH: " << WrToThMeanLH << endl;
		cout << "wrToThMeanLWam: " << wrToThMeanLWam << endl;

		cout << "shoToElMeanLH / elToWrMeanLH: " 
				 << shoToElMeanLH / elToWrMeanLH << endl;
		cout << "shoToElMeanLWam / elToWrMeanLWam: " 
				 << shoToElMeanLWam / elToWrMeanLWam << endl;

		cout << "shoToElMeanLH / elToThMeanLH: " 
				 << shoToElMeanLH / (elToWrMeanLH + WrToThMeanLH)
				 << endl;
		cout << "shoToElMeanLWam / elToThMeanLWam: " 
				 << shoToElMeanLWam / (elToWrMeanLWam + wrToThMeanLWam)
				 << endl;


		/* (B) inpts Scaled to WAM */
		MatrixXd inPtsRShWam = MatrixXd::Zero(3,nRows)
			, inPtsRUaWam = MatrixXd::Zero(3,nRows)
			, inPtsRElWam = MatrixXd::Zero(3,nRows)
			, inPtsRLaWam = MatrixXd::Zero(3,nRows)
			,	inPtsRWrWam = MatrixXd::Zero(3,nRows);

		inPtsRShWam = _inPtsRSh;
		inPtsRUaWam = shoToUaVsWam + inPtsRShWam;
		inPtsRElWam = shoToElVsWam + inPtsRShWam;
		inPtsRLaWam = shoToLaVsWam + inPtsRShWam;
		inPtsRWrWam = shoToWrVsWam + inPtsRShWam;

	/** (2.)  given marker inPts, solve for IK of WAM  **/
	//generates the DH T matrix for all the joints	     
	ForwardKin<double> wam =
		ForwardKin<double>::generateBarrettWAM();
	const int  nJoints = wam.getNJoints();
  Problem problem;


/*	
	BasePoseProblem<double> basePoseProblem(wam, _JsWam
		, inPtsRUAinChestWam, inPtsRELinChestWam, inPtsRLAinChestWam
		, inPtsRWinChestWam, _jointMinAngles, _jointMaxAngles);
	BasePoseProblem<double> basePoseProblem(wam, _JsWam
		, shoToUaVsWam, shoToElVsWam, shoToLaVsWam
		, shoToWrVsWam, _jointMinAngles, _jointMaxAngles);
*/	

	BasePoseProblem<double> basePoseProblem(wam, _JsWam
		, shoToUaVsWam, elOffsetVsWam, elToWrVsH
		, shoToWrVsWam, _jointMinAngles, _jointMaxAngles);

 	BuildProblem(&basePoseProblem, &problem);
 	Solver::Options options;
	SetSolverOptionsFromFlags(&basePoseProblem, &options);
  options.gradient_tolerance = 1e-100;
  options.function_tolerance = 1e-100;
	options.parameter_tolerance = 1e-100;
	options.max_num_consecutive_invalid_steps = 100;
		
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  std::cout << summary.FullReport() << "\n";

	double* startConfigJoints = 
		basePoseProblem.getStartConfigJoints();
	double* startConfigCamera =
		basePoseProblem.getStartConfigCamera();

		const double deltaSize = 1.0/double(nRows);
	// Fit nth order polynomial to THETA5
		for(size_t i = 0; i < nRows; i++) {
			double multiplier = double(i) * deltaSize; 
			for(size_t j = 0; j < nJoints; j++) {
				size_t jTimesFour = j * 4;
				_JsWam(i,j) = simpleThirdOrder(
					startConfigJoints[jTimesFour]
					, startConfigJoints[jTimesFour + 1]
					, startConfigJoints[jTimesFour + 2]
					, startConfigJoints[jTimesFour + 3]
					, multiplier);
			}	
		}
/*
	// transformation from camera to base   
	MatrixXd transformationToBaseMat = MatrixXd::Identity(4,4);
	double rotToBaseMat[9];
		ceres::AngleAxisToRotationMatrix
			(startConfigCamera,rotToBaseMat);
		for(int i = 0; i < 3; i++)	{
			size_t iTimesThree = i * 3;
			for(int j = 0; j < 3; j++)
				transformationToBaseMat(i,j) = 
					rotToBaseMat[iTimesThree + j];
		}
	// Translation from base to camera 
//		transformationToBaseMat(0,3) = startConfigCamera[3];
//		transformationToBaseMat(1,3) = startConfigCamera[4];
//		transformationToBaseMat(2,3) = startConfigCamera[5];
		cout << "transformationToBaseMat: \n" 
				 << transformationToBaseMat << endl;
*/
	//OutPTS of WAM 
		MatrixXd origin = MatrixXd::Zero(4,1)
			, wamEltoUaTrans = MatrixXd ::Zero(4,1)
			, wamElToLaTrans = MatrixXd::Zero(4,1)
			, wamElToWrTrans =  MatrixXd::Zero(4,1)
			, zAxis = MatrixXd::Zero(4,1);
		origin << 0.0, 0.0, 0.0, 1.0;
		wamEltoUaTrans << (0.0), (0.0)
			, 0.55, (1.0);
		wamElToLaTrans << 0.0, 0.0, (0.3/2.0), 1.0;
		wamElToWrTrans << 0.0, 0.0, 0.3, 1.0;
		zAxis << 0.0, 0.0, 1.0, 1.0;

		MatrixXd outShoToUaVsWamInBase = MatrixXd::Zero(3,nRows);
		MatrixXd outShoToUaVsWam = MatrixXd::Zero(3,nRows);


		for(size_t i = 0; i < nRows; i++)	{
			wam.setThetaVect(_JsWam.row(i).transpose());

			outShoToUaVsWamInBase.col(i) = homToCart(
				(wam.getMat(2) * wamEltoUaTrans));

		//ROTATION from base to camera frame	
		  double outShoToUaViWam[3];
			double *outShoToUaViWamInBase = new double[3];
			outShoToUaViWamInBase[0] = outShoToUaVsWamInBase(0,i);
			outShoToUaViWamInBase[1] = outShoToUaVsWamInBase(1,i);
			outShoToUaViWamInBase[2] = outShoToUaVsWamInBase(2,i);

		  ceres::AngleAxisRotatePoint(startConfigCamera
				, outShoToUaViWamInBase, outShoToUaViWam);   

			outShoToUaVsWam(0,i) = outShoToUaViWam[0];
			outShoToUaVsWam(1,i) = outShoToUaViWam[1];
			outShoToUaVsWam(2,i) = outShoToUaViWam[2];

		}

	//Compute fit error 
		MatrixXd fitErrUAfJ2 = (shoToUaVsWam 
			- outShoToUaVsWam).colwise().norm();

		cout << "fitErrUAfJ2 Mean = " 
			 << fitErrUAfJ2.mean() << endl;

	/** outPts **/
	/* in Base */
		// (1) output Pts in BASE (i.e. SAME AS VECTORS)
	printEigenMathematica( shoToUaVsWam.transpose()
		, cout, "shoToUaVsWam");	
	printEigenMathematica( outShoToUaVsWam.transpose()
		, cout, "outShoToUaVsWam");	
	printEigenMathematica( outShoToUaVsWamInBase.transpose()
		, cout, "outShoToUaVsWamInBase");	
}

void setJointLimit(double* startConfigJ, int joint
	, Problem* problem)	{

		//set the joint limits into ceres for Theta_1, Theta_2
		int ctr = 0;
			problem->SetParameterLowerBound(startConfigJ, ctr
				, _jointMinAngles[joint]);
			problem->SetParameterUpperBound(startConfigJ, ctr
				, _jointMaxAngles[joint]);
			ctr++;
			for(size_t term = 0; term < _polynomialOrder; term++) {
				problem->SetParameterLowerBound(
					startConfigJ, ctr, _jointMinAngles[joint]);
				problem->SetParameterUpperBound(
					startConfigJ, ctr, _jointMaxAngles[joint]);
				ctr++;
			}

}
/*

	*/

	const int  _polynomialOrder = 3;
	ForwardKin<T> &_fk;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &_JsWam;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _inPtsRW,
			_inPtsRTh, _inPtsRPk, _inPtsRLA, _inPtsREl, _inPtsRUA
			, _inPtsRSh, _inPtsRCh, _inPtsMCh, _inPtsLCh, _inPtsLSh;

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _shoToUaVsWam
		, _shoToElVsWam, _shoToLaVsWam, _shoToWrVsWam;
	std::vector<T> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
	 	_a, _alpha, _d;
	T* _startConfigJoints;
	T* _startConfigCamera;
};

}

#endif
