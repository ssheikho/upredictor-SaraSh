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

DEFINE_string(trust_region_strategy, "levenberg_marquardt",
              "Options are: levenberg_marquardt, dogleg.");
DEFINE_string(dogleg, "traditional_dogleg", "Options are: traditional_dogleg,"
              "subspace_dogleg.");
DEFINE_bool(inner_iterations, false, "Use inner iterations to non-linearly "
            "refine each successful trust region step.");
DEFINE_string(blocks_for_inner_iterations, "automatic", "Options are: "
            "automatic, cameras, points, cameras,points, points,cameras");
DEFINE_string(linear_solver, "sparse_schur", "Options are: "
              "sparse_schur, dense_schur, iterative_schur, sparse_normal_cholesky, "
              "dense_qr, dense_normal_cholesky and cgnr.");
DEFINE_bool(explicit_schur_complement, false, "If using ITERATIVE_SCHUR "
            "then explicitly compute the Schur complement.");
DEFINE_string(preconditioner, "jacobi", "Options are: "
              "identity, jacobi, schur_jacobi, cluster_jacobi, "
              "cluster_tridiagonal.");
DEFINE_string(visibility_clustering, "canonical_views",
              "single_linkage, canonical_views");
DEFINE_string(sparse_linear_algebra_library, "suite_sparse",
              "Options are: suite_sparse and cx_sparse.");
DEFINE_string(dense_linear_algebra_library, "eigen",
              "Options are: eigen and lapack.");
DEFINE_string(ordering, "automatic", "Options are: automatic, user.");
DEFINE_bool(use_quaternions, false, "If true, uses quaternions to represent "
            "rotations. If false, angle axis is used.");
DEFINE_bool(use_local_parameterization, false, "For quaternions, use a local "
            "parameterization.");
DEFINE_bool(robustify, true, "Use a robust loss function.");
DEFINE_double(eta, 1e-6, "Default value for eta. Eta determines the "
             "accuracy of each linear solve of the truncated newton step. "
             "Changing this parameter can affect solve performance.");
DEFINE_int32(num_threads, 1, "Number of threads.");
DEFINE_int32(num_iterations, 50, "Number of iterations.");
DEFINE_double(max_solver_time, 1e32, "Maximum solve time in seconds.");
DEFINE_bool(nonmonotonic_steps, true, "Trust region algorithm can use"
            " nonmonotic steps.");
DEFINE_double(rotation_sigma, 0.0, "Standard deviation of camera rotation "
              "perturbation.");
DEFINE_double(translation_sigma, 0.0, "Standard deviation of the camera "
              "translation perturbation.");
DEFINE_double(point_sigma, 0.0, "Standard deviation of the point "
              "perturbation.");
DEFINE_int32(random_seed, 38401, "Random seed used to set the state "
             "of the pseudo random number generator used to generate "
             "the pertubations.");
DEFINE_bool(line_search, false, "Use a line search instead of trust region "
            "algorithm.");
DEFINE_bool(mixed_precision_solves, true, "Use mixed precision solves.");
DEFINE_int32(max_num_refinement_iterations, 10, "Iterative refinement iterations");

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

#endif
