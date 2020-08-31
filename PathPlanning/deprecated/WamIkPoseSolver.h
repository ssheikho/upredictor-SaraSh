#ifndef WAM_IK_POSE_SOLVER_H
#define WAM_IK_POSE_SOLVER_H

#include "BasicFormulas.h"
#include "ForwardKin.h"
#include "ForwardKinHand.h"

#include "WamIkProblem.h"
#include "DifferentialInverseKin.h"
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
DEFINE_double(eta, 1e-2, "Default value for eta. Eta determines the "
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

struct BasePoseReprojectionError {
	BasePoseReprojectionError (
		ForwardKin<double> &fk
		, const Eigen::Matrix<double,3,1> shoToUaVi
		, const Eigen::Matrix<double,3,1> elOffsetVi 
		, const Eigen::Matrix<double,3,1> wrOffsetVi
		, const Eigen::Matrix<double,3,1> laToWrVi
		, const Eigen::Matrix<double,3,1> shoToWrVi
		, const double parametricTerm
		, std::vector<double> jointMinAngles
		, std::vector<double> jointMaxAngles)
		: _fk(fk)
		, _shoToElVi(shoToUaVi)
		, _elOffsetVi(elOffsetVi)
		, _wrOffsetVi(wrOffsetVi)
		, _laToWrVi(laToWrVi)
		, _shoToWrVi(shoToWrVi)
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		,_a(_fk.getAvect(7))
		, _alpha(_fk.getAlphaVect(7))
		, _d(_fk.getDvect(7)) 
	{

 }

	// Each Residual block takes 4 marker points and a camera as input
	// and outputs a ???? dimensional residual.
	template <typename U>
	bool operator()(const U* const camera
		, const U* const candidateParamJs
		, U* residuals) const {
		//first 4 joints contribute to wrist position
		int toDof = 4;
		U *thetas = new U[toDof];	
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
		fkMat =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4), dhTransform(4, 4);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			fkMatEl =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			fkMatLa =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			fkMatWr =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			origin(4,1);
		origin << U(0.0), U(0.0), U(0.0), U(1.0);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			wamEltoUaTrans(4,1);
		wamEltoUaTrans << U(0.0), U(0.0), U(0.55), U(1.0);
	
		//calculate thetas
		bool retVal = true;
		for(size_t i = 0; i < toDof; i++) {
			size_t iTimesFour = i * 4;
	
			thetas[i] = candidateParamJs[iTimesFour] 
				+	U(_parametricTerm) * candidateParamJs[iTimesFour + 1] 
				+	U(_parametricTerm * _parametricTerm) 
					* candidateParamJs[iTimesFour + 2] 
				+	U(_parametricTerm * _parametricTerm 
					* _parametricTerm) * candidateParamJs[iTimesFour + 3];

			retVal =  retVal && (thetas[i] >= _jointMinAngles[i]);
			retVal =  retVal && (thetas[i] <= _jointMaxAngles[i]);
			}

		//if joints within limit, set to wam		
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
					if (i < 4)
						fkMatLa = fkMatLa * dhTransform;

					//fkMatWr = fkMatWr * dhTransform;
				
		}
		
			}
			//(1) "θ1 & θ2": predicting each ELBOW-NO_OFFSET marker
			// position wrt the base  of the robot
			U *predShoToUaViInBase = new U[3];	
			predShoToUaViInBase[0] = 
				U(fkMat(0,0)) * U(wamEltoUaTrans(0,0)) 
				+ U(fkMat(0,1)) * U(wamEltoUaTrans(1,0)) 
				+ U(fkMat(0,2)) * U(wamEltoUaTrans(2,0)) 
				+ U(fkMat(0,3)) * U(wamEltoUaTrans(3,0));
			predShoToUaViInBase[1] = 
				U(fkMat(1,0)) * U(wamEltoUaTrans(0,0)) 
				+ U(fkMat(1,1)) * U(wamEltoUaTrans(1,0)) 
				+ U(fkMat(1,2)) * U(wamEltoUaTrans(2,0)) 
				+ U(fkMat(1,3)) * U(wamEltoUaTrans(3,0)); 
			predShoToUaViInBase[2] =
				U(fkMat(2,0)) * U(wamEltoUaTrans(0,0)) 
				+ U(fkMat(2,1)) * U(wamEltoUaTrans(1,0)) 
				+ U(fkMat(2,2)) * U(wamEltoUaTrans(2,0)) 
				+ U(fkMat(2,3)) * U(wamEltoUaTrans(3,0)); 
			//ROTATE (Vi IN BASE) AROUND CAMERA AXIS TO GET TO
			// (Vi) IN CAMERA FRAME
		  U predShoToUaVi[3];
		  ceres::AngleAxisRotatePoint(camera, predShoToUaViInBase
				, predShoToUaVi); 
			//Compute Residuals 
			residuals[0] = U(_shoToElVi[0]) - U(predShoToUaVi[0]);
			residuals[1] = U(_shoToElVi[1]) - U(predShoToUaVi[1]);
			residuals[2] = U(_shoToElVi[2]) - U(predShoToUaVi[2]);
			residuals[3] = U(0.55 * 0.55) 
				- U(_shoToElVi[0] * U(predShoToUaVi[0])
					+ _shoToElVi[1] * U(predShoToUaVi[1])
					+ _shoToElVi[2] * U(predShoToUaVi[2]));

			//(2)."θ3": predicting each ELBOWoffset
			U *predShoToElViInBase = new U[3];	
			predShoToElViInBase[0] = 
				U(fkMatEl(0,0)) * U(origin(0,0)) 
				+ U(fkMatEl(0,1)) * U(origin(1,0)) 
				+ U(fkMatEl(0,2)) * U(origin(2,0)) 
				+ U(fkMatEl(0,3)) * U(origin(3,0));
			predShoToElViInBase[1] = 
				 U(fkMatEl(1,0)) * U(origin(0,0)) 
				+ U(fkMatEl(1,1)) * U(origin(1,0)) 
				+ U(fkMatEl(1,2)) * U(origin(2,0)) 
				+ U(fkMatEl(1,3)) * U(origin(3,0)); 
			predShoToElViInBase[2] =
			 U(fkMatEl(2,0)) * U(origin(0,0)) 
				+ U(fkMatEl(2,1)) * U(origin(1,0)) 
				+ U(fkMatEl(2,2)) * U(origin(2,0)) 
				+ U(fkMatEl(2,3)) * U(origin(3,0)); 
			//ROTATE (Vi IN BASE) AROUND CAMERA AXIS TO GET TO
			// (Vi) IN CAMERA FRAME
		  U predShoToElVi[3];
			ceres::AngleAxisRotatePoint(camera, predShoToElViInBase
				, predShoToElVi);    
			//Compute Residuals 
			residuals[4] = U(_elOffsetVi[0]) 
				- U(predShoToElVi[0] - predShoToUaVi[0]);
			residuals[5] = U(_elOffsetVi[1]) 
				- U(predShoToElVi[1] - predShoToUaVi[1]);
			residuals[6] = U(_elOffsetVi[2]) 
				- U(predShoToElVi[2] - predShoToUaVi[2]);
			//a.a = ||a||^2
			residuals[7] = U(0.045 * 0.045) - 
				U(_elOffsetVi[0]
				* U(predShoToElVi[0] - predShoToUaVi[0])
				+ _elOffsetVi[1]
				* U(predShoToElVi[1] - predShoToUaVi[1])
				+ _elOffsetVi[2]
				* U(predShoToElVi[2] - predShoToUaVi[2]));

			//(3)."θ4": predicting each WRISToffset
			U *predShoToLaViInBase = new U[3];	
			predShoToLaViInBase[0] = 
				U(fkMatLa(0,0)) * U(origin(0,0)) 
				+ U(fkMatLa(0,1)) * U(origin(1,0)) 
				+ U(fkMatLa(0,2)) * U(origin(2,0)) 
				+ U(fkMatLa(0,3)) * U(origin(3,0));
			predShoToLaViInBase[1] = 
				 U(fkMatLa(1,0)) * U(origin(0,0)) 
				+ U(fkMatLa(1,1)) * U(origin(1,0)) 
				+ U(fkMatLa(1,2)) * U(origin(2,0)) 
				+ U(fkMatLa(1,3)) * U(origin(3,0)); 
			predShoToLaViInBase[2] =
			 U(fkMatLa(2,0)) * U(origin(0,0)) 
				+ U(fkMatLa(2,1)) * U(origin(1,0)) 
				+ U(fkMatLa(2,2)) * U(origin(2,0)) 
				+ U(fkMatLa(2,3)) * U(origin(3,0)); 
			//ROTATE (Vi IN BASE) AROUND CAMERA AXIS TO GET TO
			// (Vi) IN CAMERA FRAME
		  U predShoToLaVi[3];
			ceres::AngleAxisRotatePoint(camera, predShoToLaViInBase
				, predShoToLaVi);    
			//Compute Residuals 
			residuals[8] = U(_wrOffsetVi[0]) 
				- U(predShoToLaVi[0] - predShoToElVi[0]);
			residuals[9] = U(_wrOffsetVi[1]) 
				- U(predShoToLaVi[1] - predShoToElVi[1]);
			residuals[10] = U(_wrOffsetVi[2]) 
				- U(predShoToLaVi[2] - predShoToElVi[2]);
			//a.a = ||a||^2
			residuals[11] = U(0.045 * 0.045) - 
				U(_wrOffsetVi[0]
				* U(predShoToLaVi[0] - predShoToElVi[0])
				+ _wrOffsetVi[1]
				* U(predShoToLaVi[1] - predShoToElVi[1])
				+ _wrOffsetVi[2]
				* U(predShoToLaVi[2] - predShoToElVi[2]));
/*
			//(4). predicting each La-to-Wrist
			U *predShoToWrViInBase = new U[3];	
			predShoToWrViInBase[0] = 
				U(fkMatWr(0,0)) * U(origin(0,0)) 
				+ U(fkMatWr(0,1)) * U(origin(1,0)) 
				+ U(fkMatWr(0,2)) * U(origin(2,0)) 
				+ U(fkMatWr(0,3)) * U(origin(3,0));
			predShoToWrViInBase[1] = 
				 U(fkMatWr(1,0)) * U(origin(0,0)) 
				+ U(fkMatWr(1,1)) * U(origin(1,0)) 
				+ U(fkMatWr(1,2)) * U(origin(2,0)) 
				+ U(fkMatWr(1,3)) * U(origin(3,0)); 
			predShoToWrViInBase[2] =
			 U(fkMatWr(2,0)) * U(origin(0,0)) 
				+ U(fkMatWr(2,1)) * U(origin(1,0)) 
				+ U(fkMatWr(2,2)) * U(origin(2,0)) 
				+ U(fkMatWr(2,3)) * U(origin(3,0)); 
			//ROTATE (Vi IN BASE) AROUND CAMERA AXIS TO GET TO
			// (Vi) IN CAMERA FRAME
		  U predShoToWrVi[3];
			ceres::AngleAxisRotatePoint(camera, predShoToWrViInBase
				, predShoToWrVi);    
			//Compute Residuals 
			residuals[12] = U(_laToWrVi[0]) 
				- U(predShoToWrVi[0] -predShoToLaVi[0]);
			residuals[13] =	U(_laToWrVi[1]) 
				- U(predShoToWrVi[1] -predShoToLaVi[1]);
			residuals[14] = U(_laToWrVi[2]) 
				- U(predShoToWrVi[2] -predShoToLaVi[2]);
			//a.a = ||a||^2
			residuals[15] = U(0.3 * 0.3) - 
				U(_laToWrVi[0]
				* U(predShoToWrVi[0] -predShoToLaVi[0])
				+ _laToWrVi[1]
				* U(predShoToWrVi[1] -predShoToLaVi[1])
				+ _laToWrVi[2]
				* U(predShoToWrVi[2] -predShoToLaVi[2]));
*/
	return true;
	}


		static ceres::CostFunction* Create(ForwardKin<double> &fk
		, const Eigen::Matrix<double,3,1> shoToUaVi
		, const Eigen::Matrix<double,3,1> elOffsetVi 
		, const Eigen::Matrix<double,3,1> wrOffsetVi
		, const Eigen::Matrix<double,3,1> laToWrVi
		, const Eigen::Matrix<double,3,1> shoToWrVi
		, const double parametricTerm
		, std::vector<double> jointMinAngles
		, std::vector<double> jointMaxAngles) {			
			return (new ceres::AutoDiffCostFunction
				<BasePoseReprojectionError,12,3,16> (
					new BasePoseReprojectionError (fk, shoToUaVi, elOffsetVi
						, wrOffsetVi, laToWrVi, shoToWrVi, parametricTerm
						, jointMinAngles, jointMaxAngles)));
		}

	ForwardKin<double> &_fk;
	const Eigen::Matrix<double,3,1> _shoToElVi, _elOffsetVi
		, _wrOffsetVi, _laToWrVi, _shoToWrVi;
	const double _parametricTerm;
	std::vector<double> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
	 	_a, _alpha, _d;
};

/*
	GOAL: to solve for theta_5, theta_6 and theta_7 for WAM
		based on input WRIST, THUMB, and PINKY markers
		Each marker inPt is referenced to J4 coordinate frame in 
		the center of the elbow of the wam.
*/

namespace ceres {

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


	void SetOrdering(WamIkProblem<double>* WamIkProblem
		 , Solver::Options* options) {
		const int nPts = WamIkProblem->getNpts();
		const int cameraBlockSize = 
			WamIkProblem->getCameraBlockSize();
		const int ptsBlockSize = 
			WamIkProblem->getPtsBlockSize();
		//double* startConfigJoints = new double[ptsBlockSize];
		//double* startConfigCamera = new double[cameraBlockSize];
		double* startConfigJoints = 
			WamIkProblem->getStartConfigJoints();
		double* startConfigCamera =
			WamIkProblem->getStartConfigCamera();

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
	void SetSolverOptionsFromFlags(
		WamIkProblem<double>* WamIkProblem
		, Solver::Options* options) {

		SetMinimizerOptions(options);
		SetLinearSolver(options);
		SetOrdering(WamIkProblem,options);
	}


	void BuildProblem(WamIkProblem<double>* WamIkProblem
		, Problem* problem) {

		// Observations = [u_1, u_2, ... , u_n], 3*num_observations
		// long array , where each u_i is 3 dimensional, 
		// the x, y and z positions of the mar
		const int nPts = WamIkProblem->getNpts();
		const int ptsBlockSize = 
			WamIkProblem->getPtsBlockSize();
		//double* startConfigJoints = new double[ptsBlockSize];
		//3 for rotation, 3 for trans
		const int cameraBlockSize = 
			WamIkProblem->getCameraBlockSize();
		//double* startConfigCamera = new double[cameraBlockSize];

		double* startConfigJoints = 
			WamIkProblem->getStartConfigJoints();
		double* startConfigCamera =
			WamIkProblem->getStartConfigCamera();

		const double deltaSize = 1.0/double(nPts);
		ForwardKin<double> &wam = WamIkProblem->getWam();
		size_t  nJoints = wam.getNJoints();

		for(size_t i = 0; i < nPts; i++)	{
		  double parametricTerm = double(i) * deltaSize; 
			CostFunction* cost_function;
		  // Each Residual block takes 7 joint angles, which are 
			// identified by cubic spline coefficients joint angles,
			// and a cameraT as input  and outputs a ???? dimensional 
			// residual.
			const Eigen::Matrix<double,3,1> shoToUaVi = 
				WamIkProblem->getShoToUaVsWamInChest().col(i)
				, elOffsetVi = WamIkProblem
					->getElOffsetVsWamInChest().col(i)
				, wrOffsetVi = WamIkProblem
					->getWrOffsetVsWamInChest().col(i)
				, laToWrVi = WamIkProblem->getLaToWrVsWamInChest().col(i)
				, shoToWrVi = WamIkProblem
					->getShoToWrVsWamInChest().col(i);

		  cost_function = BasePoseReprojectionError::Create (
				wam, shoToUaVi, elOffsetVi, wrOffsetVi
				, laToWrVi, shoToWrVi, parametricTerm, WamIkProblem
				->getJointMinAngles(), WamIkProblem->getJointMaxAngles());

		  // If enabled use Huber's loss function.
		  LossFunction* loss_function = FLAGS_robustify ?
				new HuberLoss(1.0) : NULL;
		  // Each observation correponds to a pair of cameraT and 3/4
			// observed markers points which are identified by 
			// camera_index and point_index()[i], respectively.
			problem->AddResidualBlock(cost_function, loss_function
				, WamIkProblem->getStartConfigCamera()
				, WamIkProblem->getStartConfigJoints());
		}

		//set the joint limits into ceres 
 		WamIkProblem->setJointLimits(problem
			, WamIkProblem->getStartConfigJoints(),4);	
}


//find J4 from Human elbow angle
void solveJ4FromHumanElbAngle(
	WamIkProblem<double>* WamIkProblem) {
	size_t nPoints = WamIkProblem->getNpts();

	for(size_t i = 0; i < nPoints; i++) {	
		//lengths of the upper and lower arm segments
		double uAL = WamIkProblem->getShoToElVsWam().col(i).norm();
		double lAL = WamIkProblem->getElToWrVsWam().col(i).norm();
		double j4_top = uAL*uAL + lAL*lAL 
			- pow(WamIkProblem->getShoToWrVsWam().col(i).norm(), 2.0);
		double j4_bottom = 2.0 * uAL * lAL;
		WamIkProblem->getJsWam()(i,3) = M_PI - (j4_top / j4_bottom );
	}
}


void computeFitErr(WamIkProblem<double>* WamIkProblem) {
	ForwardKin<double>&wam = WamIkProblem->getWam();
	size_t nPoints = WamIkProblem->getNpts();
	
	MatrixXd origin = MatrixXd::Zero(4,1)
		, wamShtoUaTrans = MatrixXd ::Zero(4,1)
		, wamElToLaTrans = MatrixXd::Zero(4,1)
		, wamElToWrTrans =  MatrixXd::Zero(4,1);
	origin << 0.0, 0.0, 0.0, 1.0;;
	wamShtoUaTrans << (0.0), (0.0)
		, (0.55), (1.0);
	wamElToLaTrans << 0.0, 0.0, (0.3/2.0), 1.0;
	wamElToWrTrans << 0.0, 0.0, 0.3, 1.0;


	//OutPTS of WAM 	
	//(a) outShoToUa
	MatrixXd outShoToUaVsWamInBase = MatrixXd::Zero(3,nPoints)
	, outShoToUaVsWamInChest = MatrixXd::Zero(3,nPoints)
	, outShoToUaVsWam = MatrixXd::Zero(3,nPoints);
	//(b) outShoToEl
	MatrixXd outShoToElVsWamInBase = MatrixXd::Zero(3,nPoints)
	, outShoToElVsWamInChest = MatrixXd::Zero(3,nPoints)
	, outShoToElVsWam = MatrixXd::Zero(3,nPoints);
	//(c) outElOffset
	MatrixXd outElOffsetVsWamInBase = MatrixXd::Zero(3,nPoints)
	, outElOffsetVsWamInChest = MatrixXd::Zero(3,nPoints)
	, outElOffsetVsWam = MatrixXd::Zero(3,nPoints);
	//(d) outShoToLa
	MatrixXd outShoToLaVsWamInBase = MatrixXd::Zero(3,nPoints)
	, outShoToLaVsWamInChest = MatrixXd::Zero(3,nPoints)
	, outShoToLaVsWam = MatrixXd::Zero(3,nPoints);
	//(e) outWrOffset
	MatrixXd outWrOffsetVsWamInBase = MatrixXd::Zero(3,nPoints)
	, outWrOffsetVsWamInChest = MatrixXd::Zero(3,nPoints)
	, outWrOffsetVsWam = MatrixXd::Zero(3,nPoints);
	//(f) outShoToWr
	MatrixXd outShoToWrVsWamInBase = MatrixXd::Zero(3,nPoints)
	, outShoToWrVsWamInChest = MatrixXd::Zero(3,nPoints)
	, outShoToWrVsWam = MatrixXd::Zero(3,nPoints);
	//(g) outShoToEe
	MatrixXd outShoToEeVsWamInBase = MatrixXd::Zero(3,nPoints)
	, outShoToEeVsWamInChest = MatrixXd::Zero(3,nPoints)
	, outShoToEeVsWam = MatrixXd::Zero(3,nPoints);


	for(size_t i = 0; i < nPoints; i++)	{
		wam.setThetaVect(WamIkProblem->getJsWam().row(i).transpose());

		//(a) ShoToUa
		outShoToUaVsWamInBase.col(i) = homToCart(
			(wam.getMat(2) * wamShtoUaTrans));
		//ROTATE (Vi IN BASE) AROUND CAMERA AXIS TO GET TO
		// (Vi) IN CHEST FRAME
	  double outShoToUaViWam[3];
		double *outShoToUaViWamInBase = new double[3];
		outShoToUaViWamInBase[0] = outShoToUaVsWamInBase(0,i);
		outShoToUaViWamInBase[1] = outShoToUaVsWamInBase(1,i);
		outShoToUaViWamInBase[2] = outShoToUaVsWamInBase(2,i);
	  ceres::AngleAxisRotatePoint(
			WamIkProblem->getStartConfigCamera()
			, outShoToUaViWamInBase, outShoToUaViWam);   
		outShoToUaVsWamInChest(0,i) = outShoToUaViWam[0];
		outShoToUaVsWamInChest(1,i) = outShoToUaViWam[1];
		outShoToUaVsWamInChest(2,i) = outShoToUaViWam[2];
		//ROTATION from CHEST to CAMERA frame
		Eigen::Vector3d RCh = WamIkProblem->getInPtsRCh().col(i)
			, LCh = WamIkProblem->getInPtsLCh().col(i)
			, MCh = WamIkProblem->getInPtsMCh().col(i);
		Eigen::Matrix<double, 3, 3> chestToCameraR =
		 buildRefFramefrom3Pts(RCh, MCh, LCh).transpose();
		outShoToUaVsWam.col(i) = chestToCameraR 
			* outShoToUaVsWamInChest.col(i); 

		//(b) ShoToEl
		outShoToElVsWamInBase.col(i) = homToCart(
			(wam.getMat(3) * origin));
		//ROTATION from base to CHEST frame	
	  double outShoToElViWam[3];
		double *outShoToElViWamInBase = new double[3];
		outShoToElViWamInBase[0] = outShoToElVsWamInBase(0,i);
		outShoToElViWamInBase[1] = outShoToElVsWamInBase(1,i);
		outShoToElViWamInBase[2] = outShoToElVsWamInBase(2,i);
	  ceres::AngleAxisRotatePoint(
			WamIkProblem->getStartConfigCamera()
			, outShoToElViWamInBase, outShoToElViWam);   
		outShoToElVsWamInChest(0,i) = outShoToElViWam[0];
		outShoToElVsWamInChest(1,i) = outShoToElViWam[1];
		outShoToElVsWamInChest(2,i) = outShoToElViWam[2];
		//ROTATION from CHEST to CAMERA frame
		outShoToElVsWam.col(i) = chestToCameraR 
			* outShoToElVsWamInChest.col(i); 
		
		//(c) elOffset		
		outElOffsetVsWamInBase.col(i) = outShoToElVsWamInBase.col(i)
			- outShoToUaVsWamInBase.col(i);
		//ROTATION from base to CHEST frame	
	  double outElOffsetViWam[3];
		double *outElOffsetViWamInBase = new double[3];
		outElOffsetViWamInBase[0] = outElOffsetVsWamInBase(0,i);
		outElOffsetViWamInBase[1] = outElOffsetVsWamInBase(1,i);
		outElOffsetViWamInBase[2] = outElOffsetVsWamInBase(2,i);
	  ceres::AngleAxisRotatePoint(
			WamIkProblem->getStartConfigCamera()
			, outElOffsetViWamInBase, outElOffsetViWam);   
		outElOffsetVsWamInChest(0,i) = outElOffsetViWam[0];
		outElOffsetVsWamInChest(1,i) = outElOffsetViWam[1];
		outElOffsetVsWamInChest(2,i) = outElOffsetViWam[2];
		//ROTATION from CHEST to CAMERA frame
		outElOffsetVsWam.col(i) = chestToCameraR 
			* outElOffsetVsWamInChest.col(i); 

		//(d) ShoToLa
		outShoToLaVsWamInBase.col(i) = homToCart(
			(wam.getMat(4) * origin));
		//ROTATION from base to CHEST frame	
	  double outShoToLaViWam[3];
		double *outShoToLaViWamInBase = new double[3];
		outShoToLaViWamInBase[0] = outShoToLaVsWamInBase(0,i);
		outShoToLaViWamInBase[1] = outShoToLaVsWamInBase(1,i);
		outShoToLaViWamInBase[2] = outShoToLaVsWamInBase(2,i);
	  ceres::AngleAxisRotatePoint(
			WamIkProblem->getStartConfigCamera()
			, outShoToLaViWamInBase, outShoToLaViWam);   
		outShoToLaVsWamInChest(0,i) = outShoToLaViWam[0];
		outShoToLaVsWamInChest(1,i) = outShoToLaViWam[1];
		outShoToLaVsWamInChest(2,i) = outShoToLaViWam[2];
		//ROTATION from CHEST to CAMERA frame
		outShoToLaVsWam.col(i) = chestToCameraR 
			* outShoToLaVsWamInChest.col(i); 

		//(e) Wrist offset
		outWrOffsetVsWamInBase.col(i) = 
			outShoToLaVsWamInBase.col(i)
			- outShoToElVsWamInBase.col(i);
		//ROTATION from base to CHEST frame	
	  double outWrOffsetViWam[3];
		double *outWrOffsetViWamInBase = new double[3];
		outWrOffsetViWamInBase[0] = outWrOffsetVsWamInBase(0,i);
		outWrOffsetViWamInBase[1] = outWrOffsetVsWamInBase(1,i);
		outWrOffsetViWamInBase[2] = outWrOffsetVsWamInBase(2,i);
	  ceres::AngleAxisRotatePoint(
			WamIkProblem->getStartConfigCamera()
			, outWrOffsetViWamInBase, outWrOffsetViWam);   
		outWrOffsetVsWamInChest(0,i) = outWrOffsetViWam[0];
		outWrOffsetVsWamInChest(1,i) = outWrOffsetViWam[1];
		outWrOffsetVsWamInChest(2,i) = outWrOffsetViWam[2];
		//ROTATION from CHEST to CAMERA frame
		outWrOffsetVsWam.col(i) = chestToCameraR 
			* outWrOffsetVsWamInChest.col(i); 

		//(4) Wr
		//Eigen::Matrix<double,3,1> elTo
		outShoToWrVsWamInBase.col(i) = homToCart(
			(wam.getMat(wam.getNJoints()-1) * origin));
		//ROTATION from base to CHEST frame	
	  double outShoToWrViWam[3];
		double *outShoToWrViWamInBase = new double[3];
		outShoToWrViWamInBase[0] = outShoToWrVsWamInBase(0,i);
		outShoToWrViWamInBase[1] = outShoToWrVsWamInBase(1,i);
		outShoToWrViWamInBase[2] = outShoToWrVsWamInBase(2,i);
	  ceres::AngleAxisRotatePoint(
			WamIkProblem->getStartConfigCamera()
			, outShoToWrViWamInBase, outShoToWrViWam);   
		outShoToWrVsWamInChest(0,i) = outShoToWrViWam[0];
		outShoToWrVsWamInChest(1,i) = outShoToWrViWam[1];
		outShoToWrVsWamInChest(2,i) = outShoToWrViWam[2];
		//ROTATION from CHEST to CAMERA frame
		outShoToWrVsWam.col(i) = chestToCameraR 
			* outShoToWrVsWamInChest.col(i); 

		//(4) EE
		//Eigen::Matrix<double,3,1> elTo
		outShoToEeVsWamInBase.col(i) = homToCart(
			(wam.getMat(wam.getNJoints()) * origin));
		//ROTATION from base to CHEST frame	
	  double outShoToEeViWam[3];
		double *outShoToEeViWamInBase = new double[3];
		outShoToEeViWamInBase[0] = outShoToEeVsWamInBase(0,i);
		outShoToEeViWamInBase[1] = outShoToEeVsWamInBase(1,i);
		outShoToEeViWamInBase[2] = outShoToEeVsWamInBase(2,i);
	  ceres::AngleAxisRotatePoint(
			WamIkProblem->getStartConfigCamera()
			, outShoToEeViWamInBase, outShoToEeViWam);   
		outShoToEeVsWamInChest(0,i) = outShoToEeViWam[0];
		outShoToEeVsWamInChest(1,i) = outShoToEeViWam[1];
		outShoToEeVsWamInChest(2,i) = outShoToEeViWam[2];
		//ROTATION from CHEST to CAMERA frame
		outShoToEeVsWam.col(i) = chestToCameraR 
			* outShoToEeVsWamInChest.col(i); 


	}

	//Compute fit error 
	MatrixXd fitErrUAfJ2InChest = (
		WamIkProblem->getShoToUaVsWamInChest() 
		- outShoToUaVsWamInChest).colwise().norm();
	MatrixXd fitErrUAfJ2 = (WamIkProblem->getShoToUaVsWam() 
		- outShoToUaVsWam).colwise().norm();
	MatrixXd fitErrElfJ3InChest = (
		WamIkProblem->getShoToElVsWamInChest()  
		- outShoToElVsWamInChest).colwise().norm();
	MatrixXd fitErrElfJ3 = (WamIkProblem->getShoToElVsWam()  
		- outShoToElVsWam).colwise().norm();
	MatrixXd fitErrLafJ4InChest = (
		WamIkProblem->getShoToLaVsWamInChest()  
		- outShoToLaVsWamInChest).colwise().norm();
	MatrixXd fitErrLafJ4 = (WamIkProblem->getShoToLaVsWam()  
		- outShoToLaVsWam).colwise().norm();
	MatrixXd fitErrWrfJ6InChest = (
		WamIkProblem->getShoToWrVsWamInChest()  
		- outShoToWrVsWamInChest).colwise().norm();
	MatrixXd fitErrWrfJ6 = (WamIkProblem->getShoToWrVsWam()  
		- outShoToWrVsWam).colwise().norm();
	MatrixXd fitErrEefJ7InChest = (
		WamIkProblem->getShoToEeVsWamInChest()  
		- outShoToEeVsWamInChest).colwise().norm();
	MatrixXd fitErrEefJ7 = (WamIkProblem->getShoToEeVsWam()  
		- outShoToEeVsWam).colwise().norm();

	cout << "fitErrUAfJ2InChest Mean = " 
		 << fitErrUAfJ2InChest.mean() << endl;
	cout << "fitErrUAfJ2 Mean = " 
		 << fitErrUAfJ2.mean() << endl;
	cout << "fitErrElfJ3InChest Mean = " 
		 << fitErrElfJ3InChest.mean() << endl;
	cout << "fitErrElfJ3 Mean = " 
		 << fitErrElfJ3.mean() << endl;
	cout << "fitErrLafJ4InChest Mean = " 
		 << fitErrLafJ4InChest.mean() << endl;
	cout << "fitErrLafJ4 Mean = " 
		 << fitErrLafJ4.mean() << endl;
	cout << "fitErrWrfJ6InChest Mean = " 
		 << fitErrWrfJ6InChest.mean() << endl;
	cout << "fitErrWrfJ6 Mean = " 
		 << fitErrWrfJ6.mean() << endl;
	cout << "fitErrEefJ7InChest Mean = " 
		 << fitErrEefJ7InChest.mean() << endl;
	cout << "fitErrEefJ7 Mean = " 
		 << fitErrEefJ7.mean() << endl;

	printEigenMathematica(fitErrUAfJ2, cout, "fitErrUAfJ2");
	printEigenMathematica(fitErrElfJ3, cout, "fitErrElfJ3");
	printEigenMathematica(fitErrLafJ4, cout, "fitErrLafJ4");
	printEigenMathematica(fitErrWrfJ6, cout, "fitErrWrfJ6");	
	printEigenMathematica(fitErrEefJ7, cout, "fitErrEefJ7");	
	/** outPts **/
	// (1) output Pts in BASE (i.e. SAME AS VECTORS)
	//(a) shoToUa
	printEigenMathematica(WamIkProblem->getShoToUaVsWam()
		.transpose(), cout, "shoToUaVsWam");	
	printEigenMathematica(WamIkProblem->getShoToUaVsWamInChest()
		.transpose(), cout, "shoToUaVsWamInChest");	
	printEigenMathematica(outShoToUaVsWamInChest.transpose()
		, cout, "outShoToUaVsWamInChest");	
	printEigenMathematica( outShoToUaVsWam.transpose()
		, cout, "outShoToUaVsWam");	
	printEigenMathematica( outShoToUaVsWamInBase.transpose()
		, cout, "outShoToUaVsWamInBase");	
	//(b) shoToEl
	printEigenMathematica(WamIkProblem->
		getShoToElVsWam().transpose(), cout, "shoToElVsWam");	
	printEigenMathematica(WamIkProblem->getShoToElVsWamInChest()
		.transpose(), cout, "shoToElVsWamInChest");	
	printEigenMathematica(outShoToElVsWamInChest.transpose()
		, cout, "outShoToElVsWamInChest");	
	printEigenMathematica( outShoToElVsWam.transpose()
		, cout, "outShoToElVsWam");	
	printEigenMathematica( outShoToElVsWamInBase.transpose()
		, cout, "outShoToElVsWamInBase");	
	//(c) elOffset
	printEigenMathematica(WamIkProblem->
		getElOffsetVsWam().transpose(), cout, "elOffsetVsWam");	
	printEigenMathematica(WamIkProblem->
		getElOffsetVsWamInChest().transpose(), cout
		, "elOffsetVsWamInChest");	
	printEigenMathematica(outElOffsetVsWamInChest.transpose()
		, cout, "outElOffsetVsWamInChest");	
	printEigenMathematica( outElOffsetVsWam.transpose()
		, cout, "outElOffsetVsWam");	
	printEigenMathematica( outElOffsetVsWamInBase.transpose()
		, cout, "outElOffsetVsWamInBase");	
	//(d) shoToLa
	printEigenMathematica( WamIkProblem->getShoToLaVsWam().transpose()
		, cout, "shoToLaVsWam");	
	printEigenMathematica(WamIkProblem->getShoToLaVsWamInChest()
		.transpose(), cout, "shoToLaVsWamInChest");	
	printEigenMathematica(outShoToLaVsWamInChest.transpose()
		, cout, "outShoToLaVsWamInChest");	
	printEigenMathematica( outShoToLaVsWam.transpose()
		, cout, "outShoToLaVsWam");	
	printEigenMathematica( outShoToLaVsWamInBase.transpose()
		, cout, "outShoToLaVsWamInBase");	
	//(e) WrOffset
	printEigenMathematica( WamIkProblem->
		getWrOffsetVsWam().transpose(), cout, "WrOffsetVsWam");	
	printEigenMathematica(WamIkProblem->
		getWrOffsetVsWamInChest().transpose()
		, cout, "wrOffsetVsWamInChest");	
	printEigenMathematica(outWrOffsetVsWamInChest.transpose()
		, cout, "outWrOffsetVsWamInChest");	
	printEigenMathematica( outWrOffsetVsWam.transpose()
		, cout, "outWrOffsetVsWam");	
	printEigenMathematica( outWrOffsetVsWamInBase.transpose()
		, cout, "outWrOffsetVsWamInBase");	
	//(f) Wr
	printEigenMathematica(WamIkProblem->
		getShoToWrVsWam().transpose(), cout, "shoToWrVsWam");	
	printEigenMathematica(WamIkProblem->
		getShoToWrVsWamInChest().transpose()
		, cout, "shoToWrVsWamInChest");	
	printEigenMathematica(outShoToWrVsWamInChest.transpose()
		, cout, "outShoToWrVsWamInChest");	
	printEigenMathematica( outShoToWrVsWam.transpose()
		, cout, "outShoToWrVsWam");	
	printEigenMathematica( outShoToWrVsWamInBase.transpose()
		, cout, "outShoToWrVsWamInBase");	
	//(g) Ee
	printEigenMathematica(WamIkProblem->
		getShoToEeVsWam().transpose(), cout, "shoToEeVsWam");	
	printEigenMathematica(WamIkProblem->
		getShoToEeVsWamInChest().transpose()
		, cout, "shoToEeVsWamInChest");	
	printEigenMathematica(outShoToEeVsWamInChest.transpose()
		, cout, "outShoToEeVsWamInChest");	
	printEigenMathematica( outShoToEeVsWam.transpose()
		, cout, "outShoToEeVsWam");	
	printEigenMathematica( outShoToEeVsWamInBase.transpose()
		, cout, "outShoToEeVsWamInBase");	



}

void SolveProblem(WamIkProblem<double>* WamIkProblem) {
	Problem problem;

 	BuildProblem(WamIkProblem, &problem);
 	Solver::Options options;
	SetSolverOptionsFromFlags(WamIkProblem, &options);
  options.gradient_tolerance = 1e-30;
  options.function_tolerance = 1e-30;
	options.parameter_tolerance = 1e-30;
	std::cout << "(* " << endl;
  Solver::Summary summary;
  Solve(options, &problem, &summary);
  std::cout << summary.FullReport() << "*) \n";

	double* startConfigJoints = 
		WamIkProblem->getStartConfigJoints();
	double* startConfigCamera =
		WamIkProblem->getStartConfigCamera();
	cout << "Angle/Axes Rot: " << endl;
	for(int i = 0; i < WamIkProblem->getCameraBlockSize(); i++)
		cout << startConfigCamera[i] << endl;

		size_t nPoints = WamIkProblem->getNpts();
		const double deltaSize = 1.0/double(nPoints);
		size_t  nJoints = WamIkProblem->getWam().getNJoints();

	// Fit nth order polynomial to THETA5
		for(size_t i = 0; i < nPoints; i++) {
			double multiplier = double(i) * deltaSize; 
			for(size_t j = 0; j < nJoints-3; j++) {
				size_t jTimesFour = j * 4;
				
				WamIkProblem->getJsWam()(i,j) = simpleThirdOrder(
					startConfigJoints[jTimesFour]
					, startConfigJoints[jTimesFour + 1]
					, startConfigJoints[jTimesFour + 2]
					, startConfigJoints[jTimesFour + 3]
					, multiplier);

			cout << "t = " << multiplier << endl;
			cout << startConfigJoints[jTimesFour]
					 <<	", " << startConfigJoints[jTimesFour + 1]
					 <<	", " << startConfigJoints[jTimesFour + 2]
					 <<	", " << startConfigJoints[jTimesFour + 3]
					 << endl;
			}		
		}
	computeFitErr(WamIkProblem);
	//solveJ4FromHumanElbAngle(WamIkProblem);
	//solveWristAngles(WamIkProblem);
	computeFitErr(WamIkProblem);

}

}


/** Solving for the Wrist Angles 
	T5.T6.T7 = (T1.T2.T3.T4)-1.Twrist 

The orientation of the EE frame with respect to the wrist frame. The FK (rotation matrix) is comprised purely of the  last three joint angles **/
/**
void solveWristAngles(WamIkProblem<double>* WamIkProblem)	{
**
	ForwardKin<double>&wam = WamIkProblem->getWam();
	size_t nPoints = WamIkProblem->getNpts();
	for(size_t i = 0; i < nPoints; i++) {
		wam.setThetaVect(
			WamIkProblem->getJsWam().row(i).transpose());
		//(A). Desired orientation of the EE 
		//depending on wether pinky marker is 
		//tracked or not
		Eigen::Matrix<double, 3, 3> EeOrientation
			= WamIkProblem->buildWrRotReferential(i);
		
		//(B). orientation of the Elbow(j4) solved 
		//     from (θ_1,θ_2,θ_3,θ_4)
		//(b.1). Rwr in Base
		Eigen::Matrix<double, 3, 3> fkRj4WamInBase
			= wam.getMat(4).block(0,0,3,3);
		//(b.2).Rwr in CHEST 	
		double rotFbaseTCh[9];
		ceres::AngleAxisToRotationMatrix (WamIkProblem
			->getStartConfigCamera(),rotFbaseTCh);
		Eigen::Matrix<double, 3, 3> rotMatFbaseTCh;
		for (size_t i = 0; i < 3; i++)	{
			for (size_t j = 0; j < 3; j++)	{
				size_t ithRow = i * 3;
				rotMatFbaseTCh(i,j) = rotFbaseTCh[ithRow + j];
			}
		}		
		Eigen::Matrix<double, 3, 3> fkRj4WamInChest
			= rotMatFbaseTCh * fkRj4WamInBase;
		//(b.2).Rwr in CAMERA 
		Eigen::Vector3d RCh = 
			WamIkProblem->getInPtsRCh().col(i)
			, LCh = WamIkProblem->getInPtsLCh().col(i)
			, MCh = WamIkProblem->getInPtsMCh().col(i);
		Eigen::Matrix<double, 3, 3> chestToCameraR =
		 buildRefFramefrom3Pts(RCh, MCh, LCh).transpose();
		Eigen::Matrix<double, 3, 3> fkRj4Wam =
			chestToCameraR * fkRj4WamInChest; 
		//(C). orientation of the EE frame with 
		//respect to the wrist frame.
		Eigen::Matrix<double, 3, 3> EeRotWrtWrFr 
			= fkRj4Wam.transpose() * EeOrientation;

		// CHECK - use negative sign ; i.e sin(θ_6) < 0
		//solve θ_6(t)
		double sqrtTerm = sqrt( 1.0 -
			pow(EeRotWrtWrFr(2,2),2.0));

		double testAj6 = 
			atan2(EeRotWrtWrFr(2,2),sqrtTerm);
		double testBj6 = 
			atan2(EeRotWrtWrFr(2,2),-sqrtTerm);
		cout << "testAj6: " << testAj6 << endl;
		cout << "testBj6: " << testBj6 << endl;

		WamIkProblem->getJsWam()(i,5) = 
			atan2(EeRotWrtWrFr(2,2),sqrtTerm);
		//solve θ_5(t)
		WamIkProblem->getJsWam()(i,4) = atan2(
			EeRotWrtWrFr(0,2),EeRotWrtWrFr(1,2)); 
		//solve θ_7(t)
		WamIkProblem->getJsWam()(i,6) = atan2(
			-EeRotWrtWrFr(2,0),EeRotWrtWrFr(2,1)); 

		}

**

	ForwardKin<double>&wam = WamIkProblem->getWam();
	//(FIX: MISSING GST
	Eigen::Matrix<double, 4, 4> gst 
		= Eigen::Matrix<double, 4, 4>::Identity(4,4);
	//gst_inv
	Eigen::Matrix<double, 4, 4> gst_inv 
		= InverseTmat(gst);

	size_t nPoints = WamIkProblem->getNpts();
	for(size_t i = 0; i < nPoints; i++) {
		wam.setThetaVect(
			WamIkProblem->getJsWam().row(i).transpose());

		//Eqn 13.A p0 == p7		
		//lengths of the upper and lower arm segments
		double uAL = WamIkProblem
			->getShoToElVsWam().col(i).norm();
		double lAL = WamIkProblem
			->getElToWrVsWam().col(i).norm();
		Eigen::Matrix<double, 3, 1> P7;
		P7 << 1.0, 0.0, (-uAL-lAL);
		//Eqn 13.B pd == (T1.T2.T3.T4)-1 .Twrist. gst_inv . P7
		//13.B.1 (T1.T2.T3.T4)-1
		Eigen::Matrix<double, 4, 4> fkMatToJ4 = wam.getMat(4);
		Eigen::Matrix<double, 4, 4> fkMatToJ4Inv 
			= InverseTmat(fkMatToJ4);
		//13.B.2 Twrist
		Eigen::Matrix<double, 4, 4> HwrInCameraF 
			= WamIkProblem->getHwrInCamera(i);
		//13.B.3 pd == (T1.T2.T3.T4)-1 .Twrist. gst_inv . P7
		Eigen::Matrix<double, 3, 1> Pd = homToCart(
			fkMatToJ4Inv * HwrInCameraF 
				* gst_inv * cartToHom(P7)) ;
		//13.B.4 Pr is the point where the WAM wrist
		// axes intersect; i.e inPtsRWrWam. 
		Eigen::Matrix<double, 3, 1> Pr = 
			WamIkProblem->getInPtsRSh().col(i)
			+	WamIkProblem->getShoToWrVsWam().col(i);
		//13.B.5. ω i and ω j: rotation axes of T i and T j
		//ω i i.e. Z5 which is || to inElToWrV
		Eigen::Matrix<double, 3, 1> w_i = 
			WamIkProblem->getLaToWrVsWam().col(i).normalized();
		//ω j i.e. Z6 == LAplane_normalVec X Z5
		const Eigen::Matrix<double, 3, 1> UaCrossLaVhN
			= WamIkProblem->getUaCrossLaVhN(i);
		Eigen::Matrix<double, 3, 1> w_j = 
			UaCrossLaVhN.cross(w_i).normalized();

		//Eqn 15 - α
		double alpha_top = w_i.transpose().dot(w_j) 	
			* w_j.transpose().dot(P7-Pr) 
			- w_i.transpose().dot(Pd-Pr);
		double alpha_bottom = w_i.transpose().dot(w_j) 
			* w_i.transpose().dot(w_j) - 1.0;
		double alpha = alpha_top / alpha_bottom;
		//Eqn 16 - β
		double beta_top = w_i.transpose().dot(w_j) 	
			* w_i.transpose().dot(Pd-Pr) 
			- w_j.transpose().dot(P7-Pr);
		double beta_bottom = alpha_bottom;
		double beta = beta_top / beta_bottom;
		//Eqn 17 - γ
		double gamma_top = (P7-Pr).norm() * (P7-Pr).norm()
			- alpha * alpha - beta * beta - 2.0 * alpha * beta 
			* w_i.transpose().dot(w_j);
		double gamma_bottom = (w_i.cross(w_j)).norm()
			* (w_i.cross(w_j)).norm();
		double gamma = gamma_top / gamma_bottom;
		//Eqn 14 - Pg
		Eigen::Matrix<double, 3, 1> Pg 
			= alpha * w_i + beta * w_j 
			- sqrt(gamma) * w_i.cross(w_j) + Pr;
		//if a real solution for Pg exists, then θ_i and θ_j 
		//can be found with subproblem one: 
		//(18). T_i(−θ_i).Pd = Pg && (19). T_j (θ_j)P0 = Pg .
		//Eqn 18 - T_5(−θ_5).Pd = Pg
		Eigen::Matrix<double, 3, 1> uVec = (Pd-Pr) 
			- w_i * w_i.transpose() * (Pd-Pr); 
		Eigen::Matrix<double, 3, 1> vVec = (Pg-Pr) 
			- w_i * w_i.transpose() * (Pg-Pr); 
		//ensure θ_5(t) is a monotonic function 
		if (i > 1) {
			double currJ5A = -atan2( w_i.transpose() *
				uVec.cross(vVec), uVec.transpose() * vVec);
			double currJ5B = atan2( w_i.transpose() *
				uVec.cross(vVec), uVec.transpose() * vVec);
			double distA = currJ5A 
				- WamIkProblem->getJsWam()(i-1,4);
			double distB = currJ5B 
				- WamIkProblem->getJsWam()(i-1,4);
		

			WamIkProblem->getJsWam()(i,4) = 
				abs(distA) > abs(distB) ? currJ5B : currJ5A; 
		}
		else {
			WamIkProblem->getJsWam()(i,4) = -atan2( 
				w_i.transpose() * uVec.cross(vVec)
				, uVec.transpose() * vVec);
		}

		//Eqn  (19). T_6 (θ_6)P7 = Pg
		Eigen::Matrix<double, 3, 1> uVec6 = (P7-Pr) 
			- w_j * w_j.transpose() * (P7-Pr); 
		Eigen::Matrix<double, 3, 1> vVec6 = (Pg-Pr) 
			- w_j * w_j.transpose() * (Pg-Pr); 
		WamIkProblem->getJsWam()(i,5) = atan2( 
			w_j.transpose() *	uVec6.cross(vVec6)
			, uVec6.transpose() * vVec6);

**
		//θ_5(t)
		double currJ5 = -atan2( w_i.transpose() *
			uVec.cross(vVec), uVec.transpose() * vVec);
		//ensure θ_5(t) is within joint limits
		while(currJ5 >= WamIkProblem->getJointMaxAngles())
			currJ5 -= M_PI;
		while(currJ5 <= WamIkProblem->getJointMinAngles())
			currJ5 += M_PI;

		if (i > 1) {

		double testA = currJ5 + T(2.0 * M_PI);
		double testB = currJ5 - T(2.0 * M_PI);
		
		double distNull = currJ5 - WamIkProblem->getJsWam()(i-1,4);
		T distA = testA - retVal(i - 1, 0);
		T distB = testB - retVal(i - 1, 0);

			Eigen::Matrix<double, Dynamic, 1> currJs5 
				= WamIkProblem->getJsWam().block(0,4,i+1,1);	
			printEigenMathematica(currJs5, cout, "currJs5");
			WamIkProblem->getJsWam().block(0,4,i+1,1)
				= fixThetas(currJs5);
		}


		cout << "alpha_top = " << alpha_top << ", "
			<< "alpha_bottom = " << alpha_bottom << ", "
			<< "alpha = " << alpha << endl;
		cout << "beta_top = " << beta_top << ", "
			<< "beta_bottom = " << beta_bottom << ", "
			<< "beta = " << beta << endl;
		cout << "gamma_top = " << gamma_top << ", "
			<< "gamma_bottom = " << gamma_bottom << ", "
			<< "gamma = " << gamma << endl;

		printEigenMathematica(P7, cout, "P7");
		printEigenMathematica(Pd, cout, "Pd");
		printEigenMathematica(Pr, cout, "Pr");
		printEigenMathematica(Pg, cout, "Pg");
**

	}


}
*/
#endif

