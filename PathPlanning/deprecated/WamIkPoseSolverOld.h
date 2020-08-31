#ifndef WAM_IK_POSE_SOLVER_H
#define WAM_IK_POSE_SOLVER_H

#include "BasicFormulas.h"
#include "ForwardKin.h"
#include "WamIkProblem.h"
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
DEFINE_double(eta, 1e-3, "Default value for eta. Eta determines the "
             "accuracy of each linear solve of the truncated newton step. "
             "Changing this parameter can affect solve performance.");
DEFINE_int32(num_threads, 1, "Number of threads.");
DEFINE_int32(num_iterations, 100, "Number of iterations.");
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
		ForwardKin<double> &fk, double* shoToElVi
		, double* elOffsetVi , double* wrOffsetVi
		, double* elToWrVi, double* shoToWrVi
		, double parametricTerm
		, std::vector<double> jointMinAngles
		, std::vector<double> jointMaxAngles)
		: _fk(fk)
		, _shoToElVi(shoToElVi)
		, _elOffsetVi(elOffsetVi)
		, _wrOffsetVi(wrOffsetVi)
		, _elToWrVi(elToWrVi)
		, _shoToWrVi(shoToWrVi)
		, _parametricTerm(parametricTerm)
		, _jointMinAngles(jointMinAngles)
		, _jointMaxAngles(jointMaxAngles)
		,_a(_fk.getAvect(7))
		, _alpha(_fk.getAlphaVect(7))
		, _d(_fk.getDvect(7)) 
	{
/*
	printEigenMathematica( _a, cout
		, "_a");	
	printEigenMathematica( _alpha, cout
		, "_alpha");	
	printEigenMathematica( _d, cout
		, "_d");	
*/	
 }

	// Each Residual block takes 4 marker points and a camera as input
	// and outputs a ???? dimensional residual.
	template <typename U>
	bool operator()(const U* const camera
		, const U* const candidateParamJs
		, U* residuals) const {

		int toDof = 3;
		U *thetas = new U[toDof];	
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
		fkMat =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4), dhTransform(4, 4);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			fkMatEl =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4);
		Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic> 
			fkMatWr =	Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>
			::Identity(4,4);

		//calculate thetas
		//if joints within limit, set to wam		

		bool retVal = true;
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

//					if (i < 5)
//					fkMatWr = fkMatWr * dhTransform;
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


		//predicting each ELBOW + offset
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

    U predShoToElVi[3];
		ceres::AngleAxisRotatePoint(camera, predShoToElViInBase
			, predShoToElVi);    
/*
		//predicting each Wrist
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

    U predShoToWrVi[3];
		ceres::AngleAxisRotatePoint(camera, predShoToWrViInBase
			, predShoToWrVi);    


*/
	
		residuals[0] = 
			U(10.0) * U(_shoToElVi[0]) 
			- U(10.0) * U(predShoToUaVi[0]);
		residuals[1] = 
			U(10.0) * U(_shoToElVi[1]) 
			- U(10.0) * U(predShoToUaVi[1]);
		residuals[2] = 
			U(10.0) * U(_shoToElVi[2]) 
			- U(10.0) * U(predShoToUaVi[2]);
/*
		residuals[3] = 
			U(100.0) *	U(_shoToElVi[0] + _elOffsetVi[0]) 
			- 	U(100.0) * U(predShoToElVi[0]);
		residuals[4] = 
				U(100.0) * U(_shoToElVi[1] + _elOffsetVi[1]) 
				-	U(100.0) *  U(predShoToElVi[1]);
		residuals[5] = 
				U(100.0) * U(_shoToElVi[2] + _elOffsetVi[2]) 
				-	U(100.0) *  U(predShoToElVi[2]);

		residuals[6] = 
			U(_shoToWrVi[0]) - U(predShoToWrVi[0]);
		residuals[7] = 
			U(_shoToWrVi[1]) - U(predShoToWrVi[1]);
		residuals[8] = 
			U(_shoToWrVi[2]) - U(predShoToWrVi[2]);
*/
	return true;
	}


	// Factory to hide the construction of the CostFunction object from
	// the client code.
		static ceres::CostFunction* Create(ForwardKin<double> &fk
		, const double* shoToUaVi, const double* elOffsetVi 
		, const double* wrOffsetVi, const double* elToWrVi
		, const double* shoToWrVi, const double parametricTerm
		, std::vector<double> jointMinAngles
		, std::vector<double> jointMaxAngles) {			
			return (new ceres::AutoDiffCostFunction
				<BasePoseReprojectionError,3,3,28> (
					new BasePoseReprojectionError (fk, shoToUaVi, elOffsetVi
						, wrOffsetVi, elToWrVi, shoToWrVi, parametricTerm
						, jointMinAngles, jointMaxAngles)));
		}

	ForwardKin<double> &_fk;
	double *_shoToElVi, *_elOffsetVi, *_wrOffsetVi
		, *_elToWrVi, *_shoToWrVi;
	double _parametricTerm;
	std::vector<double> _jointMinAngles, _jointMaxAngles;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
	 	_a, _alpha, _d;
};


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

		// Observations for each marker is 3*num_observations long array 
		// observations = [u_1, u_2, ... , u_n], where each u_i is 3 
		//dimensional, the x, y and z positions of the mar
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
			const double parametricTerm = double(i) * deltaSize; 
			CostFunction* cost_function;
		  // Each Residual block takes 7 joint angles, which are 
			// identified by cubic spline coefficients joint angles,
			// and a cameraT as input  and outputs a ???? dimensional 
			// residual.
			const double* shoToUaVi = 
				WamIkProblem->getShoToUaVi(i);
			const double* elOffsetVi = WamIkProblem->getElOffsetVi(i);
			const double* wrOffsetVi = WamIkProblem->getWrOffsetVi(i);
			const double* elToWrVi = WamIkProblem->getElToWrVi(i);
			const double* shoToWrVi = WamIkProblem->getShoToWrVi(i);

		  cost_function = BasePoseReprojectionError::Create (
				wam, shoToUaVi, elOffsetVi, wrOffsetVi
				, elToWrVi, shoToWrVi, parametricTerm, WamIkProblem
					->getJointMinAngles(), WamIkProblem->getJointMaxAngles());

		  // If enabled use Huber's loss function.
		  LossFunction* loss_function = FLAGS_robustify ?
				new HuberLoss(1.0) : NULL;
		  // Each observation correponds to a pair of cameraT and 3/4
			//observed markers points which are identified by camera_index
			// and point_index()[i] respectively.
			problem->AddResidualBlock(cost_function
				, loss_function, WamIkProblem->getStartConfigCamera()
				, WamIkProblem->getStartConfigJoints());

		}

		//set the joint limits into ceres 
 		WamIkProblem->setJointLimits(problem
			, WamIkProblem->getStartConfigJoints());	

}


void computeFitErr(WamIkProblem<double>* WamIkProblem) {
	ForwardKin<double>&wam = WamIkProblem->getWam();
	size_t nPoints = WamIkProblem->getNpts();
	
	MatrixXd origin = MatrixXd::Zero(4,1)
		, wamEltoUaTrans = MatrixXd ::Zero(4,1)
		, wamElToLaTrans = MatrixXd::Zero(4,1)
		, wamElToWrTrans =  MatrixXd::Zero(4,1);
	origin << 0.0, 0.0, 0.0, 1.0;;
	wamEltoUaTrans << (0.0), (0.0)
		, (0.55), (1.0);
	wamElToLaTrans << 0.0, 0.0, (0.3/2.0), 1.0;
	wamElToWrTrans << 0.0, 0.0, 0.3, 1.0;


	//OutPTS of WAM 
	MatrixXd outShoToUaVsWamInBase = MatrixXd::Zero(3,nPoints)
	, outShoToUaVsWam = MatrixXd::Zero(3,nPoints);

	MatrixXd outShoToElVsWamInBase = MatrixXd::Zero(3,nPoints)
	, outShoToElVsWam = MatrixXd::Zero(3,nPoints);

	MatrixXd outShoToWrVsWamInBase = MatrixXd::Zero(3,nPoints)
	, outShoToWrVsWam = MatrixXd::Zero(3,nPoints);


	printEigenMathematica( WamIkProblem->getJsWam()
		, cout, "JsWamInComputeFitErr");

	for(size_t i = 0; i < nPoints; i++)	{
		wam.setThetaVect(WamIkProblem->getJsWam().row(i).transpose());

		//(1) UA
		outShoToUaVsWamInBase.col(i) = homToCart(
			(wam.getMat(2) * wamEltoUaTrans));
		//ROTATION from base to camera frame	
	  double outShoToUaViWam[3];
		double *outShoToUaViWamInBase = new double[3];
		outShoToUaViWamInBase[0] = outShoToUaVsWamInBase(0,i);
		outShoToUaViWamInBase[1] = outShoToUaVsWamInBase(1,i);
		outShoToUaViWamInBase[2] = outShoToUaVsWamInBase(2,i);

	  ceres::AngleAxisRotatePoint(WamIkProblem->getStartConfigCamera()
			, outShoToUaViWamInBase, outShoToUaViWam);   

		outShoToUaVsWam(0,i) = outShoToUaViWam[0];
		outShoToUaVsWam(1,i) = outShoToUaViWam[1];
		outShoToUaVsWam(2,i) = outShoToUaViWam[2];

		//(2) El
		outShoToElVsWamInBase.col(i) = homToCart(
			(wam.getMat(3) * origin));

		//ROTATION from base to camera frame	
	  double outShoToElViWam[3];
		double *outShoToElViWamInBase = new double[3];
		outShoToElViWamInBase[0] = outShoToElVsWamInBase(0,i);
		outShoToElViWamInBase[1] = outShoToElVsWamInBase(1,i);
		outShoToElViWamInBase[2] = outShoToElVsWamInBase(2,i);

	  ceres::AngleAxisRotatePoint(WamIkProblem->getStartConfigCamera()
			, outShoToElViWamInBase, outShoToElViWam);   

		outShoToElVsWam(0,i) = outShoToElViWam[0];
		outShoToElVsWam(1,i) = outShoToElViWam[1];
		outShoToElVsWam(2,i) = outShoToElViWam[2];


		//(3) Wr
		outShoToWrVsWamInBase.col(i) = homToCart(
			(wam.getMat(7) * origin));

		//ROTATION from base to camera frame	
	  double outShoToWrViWam[3];
		double *outShoToWrViWamInBase = new double[3];
		outShoToWrViWamInBase[0] = outShoToWrVsWamInBase(0,i);
		outShoToWrViWamInBase[1] = outShoToWrVsWamInBase(1,i);
		outShoToWrViWamInBase[2] = outShoToWrVsWamInBase(2,i);

	  ceres::AngleAxisRotatePoint(WamIkProblem->getStartConfigCamera()
			, outShoToWrViWamInBase, outShoToWrViWam);   

		outShoToWrVsWam(0,i) = outShoToWrViWam[0];
		outShoToWrVsWam(1,i) = outShoToWrViWam[1];
		outShoToWrVsWam(2,i) = outShoToWrViWam[2];


	}

	//Compute fit error 
	MatrixXd fitErrUAfJ2 = (WamIkProblem->getShoToUaVsWam() 
		- outShoToUaVsWam).colwise().norm();
	MatrixXd fitErrElfJ3 = (WamIkProblem->getShoToElVsWam()  
		- outShoToElVsWam).colwise().norm();
	MatrixXd fitErrWrfJ7 = (WamIkProblem->getShoToWrVsWam()  
		- outShoToWrVsWam).colwise().norm();

	cout << "fitErrUAfJ2 Mean = " 
		 << fitErrUAfJ2.mean() << endl;
	cout << "fitErrElfJ3 Mean = " 
		 << fitErrElfJ3.mean() << endl;
	cout << "fitErrWrfJ7 Mean = " 
		 << fitErrWrfJ7.mean() << endl;

	/** outPts **/
		// (1) output Pts in BASE (i.e. SAME AS VECTORS)
	printEigenMathematica( WamIkProblem->getShoToUaVsWam().transpose()
		, cout, "shoToUaVsWam");	
	printEigenMathematica( outShoToUaVsWam.transpose()
		, cout, "outShoToUaVsWam");	
	printEigenMathematica( outShoToUaVsWamInBase.transpose()
		, cout, "outShoToUaVsWamInBase");	

	printEigenMathematica( WamIkProblem->getShoToElVsWam().transpose()
		, cout, "shoToElVsWam");	
	printEigenMathematica( outShoToElVsWam.transpose()
		, cout, "outShoToElVsWam");	
	printEigenMathematica( outShoToElVsWamInBase.transpose()
		, cout, "outShoToElVsWamInBase");	



}


void SolveProblem(WamIkProblem<double>* WamIkProblem) {
	Problem problem;

 	BuildProblem(WamIkProblem, &problem);
 	Solver::Options options;
	SetSolverOptionsFromFlags(WamIkProblem, &options);
  options.gradient_tolerance = 1e-10;
  options.function_tolerance = 1e-10;
	options.parameter_tolerance = 1e-10;
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
			for(size_t j = 0; j < nJoints; j++) {
				size_t jTimesFour = j * 4;
				WamIkProblem->getJsWam()(i,j) = simpleThirdOrder(
					startConfigJoints[jTimesFour]
					, startConfigJoints[jTimesFour + 1]
					, startConfigJoints[jTimesFour + 2]
					, startConfigJoints[jTimesFour + 3]
					, multiplier);
			}	
		}

	computeFitErr(WamIkProblem);

}

}

#endif

