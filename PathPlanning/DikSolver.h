#ifndef DIK_SOLVER_H
#define DIK_SOLVER_H


#include "BasicFormulas.h"
#include "SpatialJacobian.h"
#include "DikProblem.h"
#include "DifferentialIKErrorTerm.h"
#include "UBCUtil.h"

#include <Eigen/SVD>
#include <Eigen/Geometry> 
#include <cmath>
#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "gflags/gflags.h"
#include "glog/logging.h"

#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;
typedef Eigen::Matrix<double, 3, 3> Mat3;
typedef Eigen::Matrix<double, 7, 1> Vec7;
typedef Eigen::Matrix<double, 3, 1> Vec3;
typedef Eigen::Vector4d Vec4;
using std::vector;

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
DEFINE_bool(use_quaternions, true, "If true, uses quaternions to represent "
            "rotations. If false, angle axis is used.");
DEFINE_bool(use_local_parameterization, true, "For quaternions, use a local "
            "parameterization.");
DEFINE_bool(robustify, true, "Use a robust loss function.");
DEFINE_double(eta, 1e-10, "Default value for eta. Eta determines the "
             "accuracy of each linear solve of the truncated newton step. "
             "Changing this parameter can affect solve performance.");
DEFINE_int32(num_threads, 1, "Number of threads.");
DEFINE_int32(num_iterations, 500, "Number of iterations.");
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

//Base Orientation as:		
//(a) quaternion; Warning  -  Operations have
//	interpreting the quaternion as rotation undefined 
//	behavior if the quaternion is not normalized.
typedef std::vector<Eigen::Quaterniond , 		
	Eigen::aligned_allocator<Eigen::Quaterniond > >
		VectorOfQuaternions;
//(b) the roll, pitch, and yaw angles, described as
// a product of successive rotations about the 
// principal coordinate axes x0 , y0 , and z0
typedef std::vector<Mat3, 		
	Eigen::aligned_allocator<Eigen::Matrix3d > >
		VectorOfRmats;
//(c) Rotation Matrix
typedef std::vector<Vec3, 		
	Eigen::aligned_allocator<Eigen::Vector3d> >
		VectorOfXYZeulers;


namespace ceres {

	void SetLinearSolver(Solver::Options* options) {
		CHECK(StringToLinearSolverType(FLAGS_linear_solver
			, &options->linear_solver_type));
		CHECK(StringToPreconditionerType(FLAGS_preconditioner
			, &options->preconditioner_type));
		CHECK(StringToVisibilityClusteringType(
			FLAGS_visibility_clustering
			, &options->visibility_clustering_type));
		CHECK(StringToSparseLinearAlgebraLibraryType(
			FLAGS_sparse_linear_algebra_library
			, &options->sparse_linear_algebra_library_type));
		CHECK(StringToDenseLinearAlgebraLibraryType(
		  FLAGS_dense_linear_algebra_library
		  , &options->dense_linear_algebra_library_type));
		options->use_explicit_schur_complement =
			FLAGS_explicit_schur_complement;
	}

	void SetMinimizerOptions(Solver::Options* options) {
		options->max_num_iterations = FLAGS_num_iterations;
		options->minimizer_progress_to_stdout = true;
		options->num_threads = FLAGS_num_threads;
		options->eta = FLAGS_eta;
		options->max_solver_time_in_seconds = 
			FLAGS_max_solver_time;
		options->use_nonmonotonic_steps =
			FLAGS_nonmonotonic_steps;

//	  options->linear_solver_type = 	
//			ceres::SPARSE_NORMAL_CHOLESKY;

		options->function_tolerance = 1e-100;
		options->parameter_tolerance = 1e-100;
		if (FLAGS_line_search) 
		  options->minimizer_type = ceres::LINE_SEARCH;

		CHECK(StringToTrustRegionStrategyType	(
			FLAGS_trust_region_strategy
			, &options->trust_region_strategy_type));
		CHECK(StringToDoglegType(FLAGS_dogleg, 
			&options->dogleg_type));
		options->use_inner_iterations = FLAGS_inner_iterations;
		options->max_num_consecutive_invalid_steps = 20;
	}

	void SetSolverOptionsFromFlags(
		DikProblem* DIKProblem
		, Solver::Options* options) {

		SetMinimizerOptions(options);
		SetLinearSolver(options);
//		SetOrdering(DIKProblem,options);
	}


	void BuildProblem	(DikProblem* DIKProblem
		, Problem* problem) {		
		const int theta_block_size =
			DIKProblem->theta_block_size();
		const int qBase_block_size =
			DIKProblem->qBase_block_size();
		double* thetas = DIKProblem->mutable_thetas();
		double* qsBase = DIKProblem->mutable_qsBase();

		ForwardKin<double> &wam = DIKProblem->getWam();
		size_t  nJoints = wam.getNJoints();

		// Position measurements 
		
		// Observations = [u_1, u_2, ... , u_n],
		// where each u_i is nMarkers * 3 dimensional,
		// the x, y and z positions of each marker at t_i
		const int nObservations = DIKProblem->nObservations();
		Eigen::MatrixXd
			shToUaVs =	DIKProblem->shoToUaVsWam(),
			elOffsetVs = DIKProblem->elOffsetVsWam(),
			wrOffsetVs = DIKProblem->wrOffsetVsWam(),
			laToWrVs = DIKProblem->laToWrVsWam(),
			shoToWrVs = DIKProblem->shoToWrVsWam();
		//info mat -> inverse of covariance (x y z)T
		Eigen::Matrix<double, 3, 3> 
			iMatShToUa = varianceSolve(shToUaVs).inverse(),
			iMatElOffset = varianceSolve(elOffsetVs).inverse(),
			iMatWrOffset = varianceSolve(wrOffsetVs).inverse(),
			iMatLaToWr = varianceSolve(laToWrVs).inverse(),
			iMatShToWr = varianceSolve(shoToWrVs).inverse();

		for (int i = 0; i < nObservations; i++)	{
			// Each Residual block takes a theta and a qBase 
			// as input and outputs a 3*nMarkers dimensional 
    	// residual. 
   		CostFunction* cfPosition;
			cfPosition = 
		  	(FLAGS_use_quaternions)
    		? DifferentialPositionConstraint::Create (
						wam, 
						shToUaVs.col(i),
						elOffsetVs.col(i), 
						wrOffsetVs.col(i),	
						laToWrVs.col(i),
					  shoToWrVs.col(i), 
						iMatShToUa, 
						iMatElOffset, 
						iMatWrOffset, 
						iMatLaToWr, 
						iMatShToWr)
				//EDIT
				: DifferentialPositionConstraint::Create (
						wam, 
						shToUaVs.col(i),
						elOffsetVs.col(i), 
						wrOffsetVs.col(i),	
						laToWrVs.col(i),
					  shoToWrVs.col(i),
						iMatShToUa, 
						iMatElOffset, 
						iMatWrOffset, 
						iMatLaToWr, 
						iMatShToWr);

    	// If enabled use Huber's loss function.
    	LossFunction* loss_function = FLAGS_robustify 
				? new HuberLoss(0.5) : NULL;
    	// Each observation correponds to a pair of 
			// a qBase and a theta which are identified  
    	// by qBase_index()[i] and theta_index()[i]
    	// respectively.
   		double* qbase = qsBase + qBase_block_size * i;
   	 	double* theta = thetas + theta_block_size * i;

			problem->AddResidualBlock(cfPosition,
				loss_function, qbase, theta);
  	}
    // Set the limit for the theta parameter at position 
		// index in the thetas_i parameter block 

	std::vector<double> jointMinAngles = 
		DIKProblem->jointMinAngles();
	std::vector<double> jointMaxAngles = 
		DIKProblem->jointMaxAngles();

		for (int i = 0; i < nObservations; ++i) {
   	 	double* thetas_i = thetas + theta_block_size * i;
			for (int joint = 0; joint < nJoints &&
				joint < theta_block_size; joint++) {
				problem->SetParameterLowerBound(thetas_i, joint
					, jointMinAngles[joint]);
				problem->SetParameterUpperBound(thetas_i, joint
					, jointMaxAngles[joint]);
			}
		}

		//For quaternions, use a local parameterization
		if (FLAGS_use_quaternions && 		
			FLAGS_use_local_parameterization) {
		  
			LocalParameterization* qBase_parameterization =
		  	new QuaternionParameterization();
		  for (int i = 0; i < nObservations; ++i) {
		    problem->SetParameterization(qsBase +
					qBase_block_size * i,
		      	qBase_parameterization);
		  }
		}
	

	}


	void BuildProblemAt	(DikProblem* DIKProblem
		, Problem* problem, int index) {		
		const int theta_block_size =
			DIKProblem->theta_block_size();
		const int qBase_block_size =
			DIKProblem->qBase_block_size();
		double* thetas = DIKProblem->mutable_thetas();
		double* qsBase = DIKProblem->mutable_qsBase();

		ForwardKin<double> &wam = DIKProblem->getWam();
		size_t  nJoints = wam.getNJoints();

		// Position measurement at t(i)
		Eigen::MatrixXd
			shToUaVs =	DIKProblem->shoToUaVsWam(),
			elOffsetVs = DIKProblem->elOffsetVsWam(),
			wrOffsetVs = DIKProblem->wrOffsetVsWam(),
			laToWrVs = DIKProblem->laToWrVsWam(),
			shoToWrVs = DIKProblem->shoToWrVsWam();
		//info mat -> inverse of covariance (x y z)T
		Eigen::Matrix<double, 3, 3> 
			iMatShToUa = varianceSolve(shToUaVs).inverse(),
			iMatElOffset = varianceSolve(elOffsetVs).inverse(),
			iMatWrOffset = varianceSolve(wrOffsetVs).inverse(),
			iMatLaToWr = varianceSolve(laToWrVs).inverse(),
			iMatShToWr = varianceSolve(shoToWrVs).inverse();

   	CostFunction* cfPosition;
		cfPosition = 
	  	(FLAGS_use_quaternions)
  		? DifferentialPositionConstraint::Create (
					wam, 
					shToUaVs.col(index),
					elOffsetVs.col(index), 
					wrOffsetVs.col(index),	
					laToWrVs.col(index),
				  shoToWrVs.col(index), 
					iMatShToUa, 
					iMatElOffset, 
					iMatWrOffset, 
					iMatLaToWr, 
					iMatShToWr)
			//EDIT
			: DifferentialPositionConstraint::Create (
					wam, 
					shToUaVs.col(index),
					elOffsetVs.col(index), 
					wrOffsetVs.col(index),	
					laToWrVs.col(index),
				  shoToWrVs.col(index),
					iMatShToUa, 
					iMatElOffset, 
					iMatWrOffset, 
					iMatLaToWr, 
					iMatShToWr);

    	// If enabled use Huber's loss function.
    	LossFunction* loss_function = FLAGS_robustify 
				? new HuberLoss(0.5) : NULL;
    	// Each observation correponds to a pair of 
			// a qBase and a theta which are identified  
    	// by qBase_index()[i] and theta_index()[i]
    	// respectively.
   		double* qbase = qsBase + qBase_block_size * index;
   	 	double* theta = thetas + theta_block_size * index;

			problem->AddResidualBlock(cfPosition,
				loss_function, qbase, theta);
  	
   		// Set the limit for the theta parameter at position 
			// index in the thetas_i parameter block 
			std::vector<double> jointMinAngles = 
				DIKProblem->jointMinAngles();
			std::vector<double> jointMaxAngles = 
				DIKProblem->jointMaxAngles();
  		double* thetas_i = thetas + theta_block_size * index;
			for (int joint = 0; joint < nJoints &&
				joint < theta_block_size; joint++) {
				problem->SetParameterLowerBound(thetas_i, joint
					, jointMinAngles[joint]);
				problem->SetParameterUpperBound(thetas_i, joint
					, jointMaxAngles[joint]);
			}

		//For quaternions, use a local parameterization
		if (FLAGS_use_quaternions && 		
			FLAGS_use_local_parameterization) {
		  
			LocalParameterization* qBase_parameterization =
		  	new QuaternionParameterization();
		    problem->SetParameterization(qsBase +
					qBase_block_size * index,
		      	qBase_parameterization);
		}
	
		

	}


	void  BuildOrientationProblem	(
		DikProblem* DIKProblem, Problem* problem) {		
		const int theta_block_size =
			DIKProblem->theta_block_size();
		const int qBase_block_size =
			DIKProblem->qBase_block_size();
		double* thetas = DIKProblem->mutable_thetas();
		double* qsBase = DIKProblem->mutable_qsBase();
  
		ForwardKin<double> &wam = DIKProblem->getWam();
		size_t  nJoints = wam.getNJoints();

		// Position measurements 
		
		// Observations = [u_1, u_2, ... , u_n],
		// where each u_i is nMarkers * 3 dimensional,
		// the x, y and z positions of each marker at t_i
		const int nObservations = DIKProblem->nObservations();
		Eigen::MatrixXd
			shToEeVs =	DIKProblem->shoToEeVsWam();
		//info mat -> inverse of covariance (x y z)T
		Eigen::Matrix<double, 3, 3> 
			iMatShToEe = varianceSolve(shToEeVs).inverse();


		for (int i = 0; i < nObservations; i++)	{
			// Each Residual block takes wrist thetas as input
			// and outputs a 3 dimensional residual; 

  	 	double* qbase_i = qsBase + qBase_block_size * i;
			Eigen::Map< Eigen::Quaternion<double> > 
				q_b_quat(qbase_i);
   		CostFunction* cfOrientation;
			cfOrientation = DifferentialOrientationConstraint::
				Create (wam, q_b_quat,shToEeVs.col(i)
					, iMatShToEe);

    	// If enabled use Huber's loss function.
    	LossFunction* loss_function = FLAGS_robustify 
				? new HuberLoss(0.5) : NULL;
    	// Each observation correponds to 3 wrist
			// thetas which are identified by theta()[i]+ 4;
   	 	double* theta = thetas + theta_block_size * i + 4;
			problem->AddResidualBlock(cfOrientation,
				loss_function, theta);
  	}

		//	void Problem::SetParameterLowerBound(
		//		double *values, int index, double lower_bound)
		//    Set the lower bound for the parameter at position 
		//		index in the parameter block corresponding to values.	
		std::vector<double> jointMinAngles = 
			DIKProblem->jointMinAngles();
		std::vector<double> jointMaxAngles = 
			DIKProblem->jointMaxAngles();

		for (int i = 0; i < nObservations; ++i) {
			double* thetas_i = thetas + theta_block_size * i + 4;
			for (int joint = 4; joint < nJoints &&
				joint < theta_block_size; joint++) {
				problem->SetParameterLowerBound(thetas_i, joint-4
					, jointMinAngles[joint]);
				problem->SetParameterUpperBound(thetas_i, joint-4
					, jointMaxAngles[joint]);
			}
		}

	}

	void BuildOrientationProblemAt	(DikProblem* DIKProblem
		, Problem* problem, int index) {		
		const int theta_block_size =
			DIKProblem->theta_block_size();
		const int qBase_block_size =
			DIKProblem->qBase_block_size();
		double* thetas = DIKProblem->mutable_thetas();
		double* qsBase = DIKProblem->mutable_qsBase();
   	double* qbase_i = qsBase + qBase_block_size * index;
		Eigen::Map< Eigen::Quaternion<double> > 
			q_b_quat(qbase_i);

		ForwardKin<double> &wam = DIKProblem->getWam();
		size_t  nJoints = wam.getNJoints();

		// Position measurement at t(i)
		Eigen::MatrixXd
			shToEeVs =	DIKProblem->shoToEeVsWam();
		//info mat -> inverse of covariance (x y z)T
		Eigen::Matrix<double, 3, 3> 
			iMatShToEe = varianceSolve(shToEeVs).inverse();

   	CostFunction* cfOrientation;
		cfOrientation = DifferentialOrientationConstraint::
			Create (wam, q_b_quat,shToEeVs.col(index)
				, iMatShToEe);

    	// If enabled use Huber's loss function.
    	LossFunction* loss_function = FLAGS_robustify 
				? new HuberLoss(0.5) : NULL;
    	// Each observation correponds to a pair of 
			// a qBase and a theta which are identified  
    	// by qBase_index()[i] and theta_index()[i]
    	// respectively.
   	 	double* theta = thetas + theta_block_size * index + 4;

			problem->AddResidualBlock(cfOrientation,
				loss_function, theta);
  	
   		// Set the limit for the theta parameter at position 
			// index in the thetas_i parameter block 
			//	void Problem::SetParameterLowerBound(
			//		double *values, int index, double lower_bound)
			//    Set the lower bound for the parameter at position 
			//		index in the parameter block corresponding to values.
			std::vector<double> jointMinAngles = 
				DIKProblem->jointMinAngles();
			std::vector<double> jointMaxAngles = 
				DIKProblem->jointMaxAngles();
  		double* thetas_i = thetas + theta_block_size * index + 4;
			for (int joint = 4; joint < nJoints &&
				joint < theta_block_size; joint++) {
				problem->SetParameterLowerBound(thetas_i, joint-4
					, jointMinAngles[joint]);
				problem->SetParameterUpperBound(thetas_i, joint-4
					, jointMaxAngles[joint]);
			}
	
	}



	//Compute fit error at t(i)
	void computeFitErrPtsRAt(
		DikProblem* DIKProblem, int index
			, bool solveOrientation) {

		const int theta_block_size =
			DIKProblem->theta_block_size();
		const int qBase_block_size =
			DIKProblem->qBase_block_size();
		double* thetas = DIKProblem->mutable_thetas();
		double* qsBase = DIKProblem->mutable_qsBase();

		ForwardKin<double> &wam = DIKProblem->getWam();
		size_t  nJoints = wam.getNJoints();

		MatrixXd origin = MatrixXd::Zero(4,1)
			, wamShtoUaTrans = MatrixXd ::Zero(4,1)
			, wamElToLaTrans = MatrixXd::Zero(4,1)
			, wamElToWrTrans =  MatrixXd::Zero(4,1);
		origin << 0.0, 0.0, 0.0, 1.0;
		wamShtoUaTrans << (0.0), (0.0)
			, (0.55), (1.0);
		wamElToLaTrans << 0.0, 0.0, (0.3/2.0), 1.0;
		wamElToWrTrans << 0.0, 0.0, 0.3, 1.0;

 		double* qbase_i = qsBase + qBase_block_size * index;
 	 	double* theta_i = thetas + theta_block_size * index;
		Eigen::Map<const Eigen::Matrix<double, 7, 1> > 
			theta_Vi(theta_i);
		wam.setThetaVect(theta_Vi);
	//(1).	estimate marker poitions w.r.t. moving base
	//"A(θ1) * A(θ2) * vUA": -> pUa(i)				
	Vec3 outShoToUaViWamInBase = homToCart(
		wam.getMat(2) * wamShtoUaTrans);
	//"A(θ1) * A(θ2) * A(θ3) * Origin": -> pEl(i)		
	Vec3 outShoToElViWamInBase = homToCart(
		wam.getMat(3) * origin);
	//"A(θ1) * A(θ2) * A(θ3) * Origin": -> pLa(i)		
	Vec3 outShoToLaViWamInBase = homToCart(
		wam.getMat(4) * origin);
	//"A(θ1) * A(θ2) * ...* A(θ5) * Origin": -> pW(i)				
	Vec3 outShoToWrViWamInBase = homToCart(
		wam.getMat(5) * origin);
	//"A(θ1) * A(θ2) * ...* A(θ7) * Origin": -> pEe(i)				
	Vec3 outShoToEeViWamInBase = homToCart(
		wam.getMat(7) * origin);

	//(2). represent estimated marker poition in
	//			the camera frame
	//	"shToUa(i)"
	Vec3 outShoToUaViWam;
	QuaternionRotatePoint(qbase_i,
		outShoToUaViWamInBase.data(),
			outShoToUaViWam.data());
	//	"shToEl(i)"
	Vec3 outShoToElViWam;
	QuaternionRotatePoint(qbase_i,
		outShoToElViWamInBase.data(),
			outShoToElViWam.data());
	//	"shToLa(i)"
	Vec3 outShoToLaViWam;
	QuaternionRotatePoint(qbase_i,
		outShoToLaViWamInBase.data(),
			outShoToLaViWam.data());
	//	"shToWr(i)"
	Vec3 outShoToWrViWam;
	QuaternionRotatePoint(qbase_i,
		outShoToWrViWamInBase.data(),
			outShoToWrViWam.data());
	//	"shToEe(i)"
	Vec3 outShoToEeViWam;
	QuaternionRotatePoint(qbase_i,
		outShoToEeViWamInBase.data(),
			outShoToEeViWam.data());

	Eigen::Map<const Eigen::Quaternion<double> > 
		q_b_quat(qbase_i);

	Mat3 rMatsBase = 
		q_b_quat.normalized().toRotationMatrix();


	//Compute fit error 
	MatrixXd fitErrUA = (DIKProblem->
		shoToUaVsWam().col(index) - outShoToUaViWam);
	MatrixXd fitErrEl = (DIKProblem->
		shoToElVsWam().col(index) - outShoToElViWam);
	MatrixXd fitErrLa = (DIKProblem->
		shoToLaVsWam().col(index) - outShoToLaViWam);
	MatrixXd fitErrWr = (DIKProblem->
		shoToWrVsWam().col(index) - outShoToWrViWam);
	MatrixXd fitErrEe = (DIKProblem->
		shoToEeVsWam().col(index) - outShoToEeViWam);

	printEigenMathematica(fitErrUA.transpose()
		, cout, "fitErrUA"+to_string(index));	
	printEigenMathematica(fitErrEl.transpose()
		, cout, "fitErrEl"+to_string(index));	
	printEigenMathematica(fitErrLa.transpose()
		, cout, "fitErrLa"+to_string(index));	
	printEigenMathematica(fitErrWr.transpose()
		, cout, "fitErrWr"+to_string(index));	
	printEigenMathematica(fitErrEe.transpose()
		, cout, "fitErrEe"+to_string(index));	
	printEigenMathematica(rMatsBase
		, cout, "rMatsBase"+to_string(index));	
	printEigenMathematica(theta_Vi.transpose()
		, cout, "outThetasWam"+to_string(index));	
	if (!solveOrientation)	{
		//initialize next observation
		if (fitErrUA.norm()>0.01)	{
			cout << "fitErrUA.norm()" << fitErrUA.norm() <<endl;
					double* qbaseIpOne = qsBase + 
						qBase_block_size * 	(index-1);
	 		 		double* theta_ipone = thetas + 
						theta_block_size * (index-1);
				for (int j = 0; j < qBase_block_size; j++)	
					*qbase_i++ = *qbaseIpOne++;
				for (int j = 0; j < theta_block_size; j++)	
					*theta_i++ = *theta_ipone++;
			Problem problem;
			BuildProblemAt(DIKProblem, &problem, index);
	 		Solver::Options options;
			SetSolverOptionsFromFlags(DIKProblem, &options);
			std::cout << "(* " << std::endl;
	 	 	Solver::Summary summary;
			Solve(options, &problem, &summary);
			std::cout << summary.FullReport() << std::endl;
			std::cout << "*) " << std::endl;
				}
		}
/*
		if (solveOrientation)	{
			//initialize next observation
			if (fitErrEe.norm()>0.01)	{
				cout << "fitErrEe.norm()" << fitErrEe.norm() <<endl;
		  	// Each observation correponds to 3 wrist
				// thetas which are identified by theta()[i]+ 4;
	 		 		double* theta_ipone = thetas + 
						theta_block_size * (index-1) + 4;
				for (int j = 0; j < theta_block_size-4; j++)	
					*theta_i++ = *theta_ipone++;
			Problem problem;
			BuildOrientationProblemAt(DIKProblem, &problem, index);
	 		Solver::Options options;
			SetSolverOptionsFromFlags(DIKProblem, &options);
			std::cout << "(* " << std::endl;
	 	 	Solver::Summary summary;
			Solve(options, &problem, &summary);
			std::cout << summary.FullReport() << std::endl;
			std::cout << "*) " << std::endl;
				}
		}
*/
	}

	
	void computeFitErrPtsR (
		DikProblem* DIKProblem) {

		const int theta_block_size =
			DIKProblem->theta_block_size();
		const int qBase_block_size =
			DIKProblem->qBase_block_size();
		double* thetas = DIKProblem->mutable_thetas();
		double* qsBase = DIKProblem->mutable_qsBase();

		ForwardKin<double> &wam = DIKProblem->getWam();
		size_t  nJoints = wam.getNJoints();
		const int nObservations = DIKProblem->nObservations();

		MatrixXd origin = MatrixXd::Zero(4,1)
			, wamShtoUaTrans = MatrixXd ::Zero(4,1)
			, wamElToLaTrans = MatrixXd::Zero(4,1)
			, wamElToWrTrans =  MatrixXd::Zero(4,1);
		origin << 0.0, 0.0, 0.0, 1.0;
		wamShtoUaTrans << (0.0), (0.0)
			, (0.55), (1.0);
		wamElToLaTrans << 0.0, 0.0, (0.3/2.0), 1.0;
		wamElToWrTrans << 0.0, 0.0, 0.3, 1.0;

		//OutPTS of WAM 
		Eigen::MatrixXd 
			outPShWam = MatrixXd::Zero(3,nObservations), 
			outPUaWam = MatrixXd::Zero(3,nObservations),
			outPElWam = MatrixXd::Zero(3,nObservations),
			outPLaWam = MatrixXd::Zero(3,nObservations), 
			outPWrWam = MatrixXd::Zero(3,nObservations), 
			outPThWam = MatrixXd::Zero(3,nObservations),
		 	outPPiWam = MatrixXd::Zero(3,nObservations),
			outPEeWam = MatrixXd::Zero(3,nObservations);
		//OutVecs of WAM
		//(a) outShoToUa
		MatrixXd outShoToUaVsWamInBase = 
				MatrixXd::Zero(3,nObservations)
		, outShoToUaVsWamInChest = 
				MatrixXd::Zero(3,nObservations)
		, outShoToUaVsWam = MatrixXd::Zero(3,nObservations),
		//(b) outShoToEl
		outShoToElVsWamInBase = MatrixXd::Zero(3,nObservations)
		, outShoToElVsWamInChest =
				MatrixXd::Zero(3,nObservations)
		, outShoToElVsWam = MatrixXd::Zero(3,nObservations),
		//(c) outElOffset
		outElOffsetVsWamInBase = MatrixXd::Zero(3,nObservations)
		, outElOffsetVsWamInChest = 
				MatrixXd::Zero(3,nObservations)
		, outElOffsetVsWam = MatrixXd::Zero(3,nObservations),
		//(d) outShoToLa
		outShoToLaVsWamInBase = MatrixXd::Zero(3,nObservations)
		, outShoToLaVsWamInChest =
				MatrixXd::Zero(3,nObservations)
		, outShoToLaVsWam = MatrixXd::Zero(3,nObservations),
		//(e) outWrOffset
		outWrOffsetVsWamInBase = MatrixXd::Zero(3,nObservations)
		, outWrOffsetVsWamInChest = 
				MatrixXd::Zero(3,nObservations)
		, outWrOffsetVsWam = MatrixXd::Zero(3,nObservations),
		//(f) outShoToWr
		outShoToWrVsWamInBase = MatrixXd::Zero(3,nObservations)
		, outShoToWrVsWamInChest = 
				MatrixXd::Zero(3,nObservations)
		, outShoToWrVsWam = MatrixXd::Zero(3,nObservations),
		//(g) outShoToEe
		outShoToEeVsWamInBase = MatrixXd::Zero(3,nObservations)
		, outShoToEeVsWamInChest = 
				MatrixXd::Zero(3,nObservations)
		, outShoToEeVsWam = MatrixXd::Zero(3,nObservations);
		//outThetas
		MatrixXd outThetasWam =
			Eigen::MatrixXd::Zero(7,nObservations);

		//Base Orientation as:
		//mathematica matrices
		MatrixXd 
			quaternionsBase = MatrixXd::Zero(4,nObservations),
			rMatsBase = MatrixXd::Zero(3,3*nObservations),
			xyzEulersBase = MatrixXd::Zero(3,nObservations),
			rMatsBaseInMaya = 
				MatrixXd::Zero(3,3*nObservations),
			xyzEulersBaseInMaya = 
				MatrixXd::Zero(3,nObservations),		
			//at time t_i	
			rMatsBaseInChest = 
				MatrixXd::Zero(3,3*nObservations),
			xyzEulersBaseInChest = 
				MatrixXd::Zero(3,nObservations);
		//vector form
		VectorOfQuaternions VectorOfQuaternionsBase;
		VectorOfRmats VectorOfRmatsBase;
		VectorOfXYZeulers VectorOfXYZeulersBase;

		for (int i = 0; i < nObservations; i++)	{
   		double* qbase_i = qsBase + qBase_block_size * i;
   	 	double* theta_i = thetas + theta_block_size * i;
    	Eigen::Map<const Eigen::Matrix<double, 7, 1> > 
				theta_Vi(theta_i);
			outThetasWam.col(i) = theta_Vi;
			wam.setThetaVect(theta_Vi);
		//(1).	estimate marker poitions w.r.t. moving base
		//"A(θ1) * A(θ2) * vUA": -> pUa(i)				
		outShoToUaVsWamInBase.col(i) = homToCart(
			wam.getMat(2) * wamShtoUaTrans);
		//"A(θ1) * A(θ2) * A(θ3) * Origin": -> pEl(i)		
		outShoToElVsWamInBase.col(i) = homToCart(
			wam.getMat(3) * origin);
		//"A(θ1) * A(θ2) * A(θ3) * A(θ4) * Origin": 
		//																	-> pLa(i)		
		outShoToLaVsWamInBase.col(i) = homToCart(
			wam.getMat(4) * origin);
		//"A(θ1) * A(θ2) * ...* A(θ5) * Origin": -> pW(i)				
		outShoToWrVsWamInBase.col(i) = homToCart(
			wam.getMat(5) * origin);
		//"A(θ1) * A(θ2) * ...* A(θ7) * Origin": -> pEE(i)				
		outShoToEeVsWamInBase.col(i) = homToCart(
			wam.getMat(7) * origin);
		/*outShoToEeVsWamInBase.col(i) = outShoToWrVsWamInBase.col(i)
			+ 0.06 * wam.getMat(7).block(0,2,3,1);*/
/*		
		//(2). represent estimated marker poitions in
		//			the camera frame

		//Base Orientation as:
		//(a) quaternion- normalized
		Eigen::Map< Eigen::Quaternion<double> > 
			q_b_quat(qbase_i);
		VectorOfQuaternionsBase.push_back(
			q_b_quat.normalized());
		quaternionsBase(0,i) = q_b_quat.x(); 
		quaternionsBase(1,i) = q_b_quat.y(); 
		quaternionsBase(2,i) = q_b_quat.z(); 
		quaternionsBase(3,i) = q_b_quat.w(); 

		//(b) Rotation Matrix		
		Mat3 rmBase_i = 
			q_b_quat.normalized().toRotationMatrix();
		VectorOfRmatsBase.push_back(rmBase_i);
		rMatsBase.block(0,3*i,3,3) = rmBase_i;
		//(c) the roll, pitch, and yaw angles
		//		 x0 , y0 , and z0
		//in Vicon frame
		Vec3 eaBase_i = rmBase_i.eulerAngles(0, 1, 2); 
		VectorOfXYZeulersBase.push_back(eaBase_i);
		xyzEulersBase.col(i) = eaBase_i;	

		//in Maya frame
		//get vicon frame in Maya local base frame
		Mat3 rMatViconWrtBaseInMaya = 
			Eigen::MatrixXd::Zero(3,3);
		rMatViconWrtBaseInMaya(0,2) = 1.0;
		rMatViconWrtBaseInMaya(1,0) = 1.0;
		rMatViconWrtBaseInMaya(2,1) = 1.0;
		//Base Orientation w.r.t. Maya loca frame of the base:
		Mat3 rmBaseInMaya_i = 
			rMatViconWrtBaseInMaya * rmBase_i;
		rMatsBaseInMaya.block(0,3*i,3,3) =
			rmBaseInMaya_i;
		xyzEulersBaseInMaya.col(i) =
			rmBaseInMaya_i.eulerAngles(0, 1, 2); 


		//Base rotation in chest frame at time t_i
		Mat3 rmBaseInChest_i = DIKProblem->rMatsCh().block
			(0,3*i,3,3).transpose() * rmBase_i;
		rMatsBaseInChest.block(0,3*i,3,3) =
			rmBaseInChest_i;
		xyzEulersBaseInChest.col(i) =
			rmBaseInChest_i.eulerAngles(0, 1, 2); 
*/
		//	"shToUa(i)"
		QuaternionRotatePoint(qbase_i,
			outShoToUaVsWamInBase.col(i).data(),
				outShoToUaVsWam.col(i).data());
		//	"shToEl(i)"
		QuaternionRotatePoint(qbase_i,
			outShoToElVsWamInBase.col(i).data(),
				outShoToElVsWam.col(i).data());
		//	"shToLa(i)"
		QuaternionRotatePoint(qbase_i,
			outShoToLaVsWamInBase.col(i).data(),
				outShoToLaVsWam.col(i).data());
		//	"shToWr(i)"
		QuaternionRotatePoint(qbase_i,
			outShoToWrVsWamInBase.col(i).data(),
				outShoToWrVsWam.col(i).data());
		//	"shToEe(i)"
		QuaternionRotatePoint(qbase_i,
			outShoToEeVsWamInBase.col(i).data(),
				outShoToEeVsWam.col(i).data());

		//(3). add shoulder position
		outPUaWam.col(i) = DIKProblem->pShWam().col(i) +
			outShoToUaVsWam.col(i);
		outPElWam.col(i) = DIKProblem->pShWam().col(i) +
			outShoToElVsWam.col(i);
		outPLaWam.col(i) = DIKProblem->pShWam().col(i) +
			outShoToLaVsWam.col(i);
		outPWrWam.col(i) =  DIKProblem->pShWam().col(i) +
			outShoToWrVsWam.col(i);

		outPEeWam.col(i) =  DIKProblem->pShWam().col(i) +
			outShoToEeVsWam.col(i);

	}


	//Compute fit error 
	MatrixXd fitErrUA = (DIKProblem->
		shoToUaVsWam() - outShoToUaVsWam);
	MatrixXd fitErrEl = (DIKProblem->
		shoToElVsWam() - outShoToElVsWam);
	MatrixXd fitErrLa = (DIKProblem->
		shoToLaVsWam() - outShoToLaVsWam);
	MatrixXd fitErrWr = (DIKProblem->
		shoToWrVsWam() - outShoToWrVsWam);
	MatrixXd fitErrEe = (DIKProblem->
		shoToEeVsWam() - outShoToEeVsWam);

	cout << "fitErrUA Mean = " 
		 << fitErrUA.colwise().norm().mean() << endl;
	cout << "fitErrEl Mean = " 
		 << fitErrEl.colwise().norm().mean() << endl;
	cout << "fitErrLa Mean = " 
		 << fitErrLa.colwise().norm().mean() << endl;
	cout << "fitErrWr Mean = " 
		 << fitErrWr.colwise().norm().mean() << endl;
	cout << "fitErrEe Mean = " 
		 << fitErrEe.colwise().norm().mean() << endl;

	printEigenMathematica(fitErrUA.transpose()
		, cout, "fitErrUA");	
	printEigenMathematica(fitErrEl.transpose()
		, cout, "fitErrEl");	
	printEigenMathematica(fitErrLa.transpose()
		, cout, "fitErrLa");	
	printEigenMathematica(fitErrWr.transpose()
		, cout, "fitErrWr");	
	printEigenMathematica(fitErrEe.transpose()
		, cout, "fitErrEe");	

	printEigenMathematica(outShoToUaVsWam.transpose()
		, cout, "outShoToUaVsWam");	
	printEigenMathematica(outShoToElVsWam.transpose()
		, cout, "outShoToElVsWam");	
	printEigenMathematica(outShoToLaVsWam.transpose()
		, cout, "outShoToLaVsWam");	
	printEigenMathematica(outShoToWrVsWam.transpose()
		, cout, "outShoToWrVsWam");	
	printEigenMathematica(outShoToEeVsWam.transpose()
		, cout, "outShoToEeVsWam");	

	printEigenMathematica(outPUaWam.transpose()
		, cout, "outPUaWam");	
	printEigenMathematica(outPElWam.transpose()
		, cout, "outPElWam");	
	printEigenMathematica(outPLaWam.transpose()
		, cout, "outPLaWam");	
	printEigenMathematica(outPWrWam.transpose()
		, cout, "outPWrWam");	
	printEigenMathematica(outPEeWam.transpose()
		, cout, "outPEeWam");	

	//outJoint angles
	printEigenMathematica(outThetasWam.transpose()
		, cout, "outThetasWam");	
/*
  for (VectorOfQuaternions::const_iterator 
		QuaternionsBase_iter = VectorOfQuaternionsBase.begin();
    	QuaternionsBase_iter != 
				VectorOfQuaternionsBase.end(); 
					++QuaternionsBase_iter) {

	   	const Eigen::Quaterniond& q_b_quat = 
				*QuaternionsBase_iter;
			cout << q_b_quat.x() << " " 
				<< q_b_quat.y() << " " 
				<< q_b_quat.z() << " "   
				<< q_b_quat.w() << " " <<  endl;
	}
*/
/*
	//out Base orientation
	printEigenMathematica(rMatsBase
		, cout, "rMatsBase");	
	printEigenMathematica(xyzEulersBase.transpose()
		, cout, "xyzEulersBase");	
	printEigenMathematica(quaternionsBase.transpose()
		, cout, "quaternionsBase");	

	printEigenMathematica(rMatsBaseInMaya
		, cout, "rMatsBaseInMaya");	
	printEigenMathematica(xyzEulersBaseInMaya.transpose()
		, cout, "xyzEulersBaseInMaya");	

	printEigenMathematica(rMatsBaseInChest
		, cout, "rMatsBaseInChest");	
	printEigenMathematica(xyzEulersBaseInChest.transpose()
		, cout, "xyzEulersBaseInChest");				
*/
	}

	void SolveOrientationProblemAt(DikProblem* DIKProblem
		, int index) {
		Problem problem;
		BuildOrientationProblemAt(DIKProblem, &problem, index);
 		Solver::Options options;
		SetSolverOptionsFromFlags(DIKProblem, &options);
		std::cout << "(* " << std::endl;
		options.minimizer_progress_to_stdout = false;
 	 	Solver::Summary summary;
  	Solve(options, &problem, &summary);
  	std::cout << summary.FullReport() << std::endl;
		std::cout << "*) " << std::endl;
		
		computeFitErrPtsRAt(DIKProblem, index, true);

	}



	void SolveOrientationProblem(
		DikProblem* DIKProblem) {

		computeFitErrPtsR(DIKProblem);

		Problem problem;
		BuildOrientationProblem(DIKProblem, &problem);
 		Solver::Options options;
		SetSolverOptionsFromFlags(DIKProblem, &options);
		std::cout << "(* " << std::endl;
 	 	Solver::Summary summary;
  	Solve(options, &problem, &summary);
  	std::cout << summary.FullReport() << std::endl;
		std::cout << "*) " << std::endl;

		computeFitErrPtsR(DIKProblem);

		}
	


	void SolveProblemAt(DikProblem* DIKProblem, int index) {
		Problem problem;
		BuildProblemAt(DIKProblem, &problem, index);
 		Solver::Options options;
		SetSolverOptionsFromFlags(DIKProblem, &options);
		std::cout << "(* " << std::endl;
		options.minimizer_progress_to_stdout = false;
 	 	Solver::Summary summary;
  	Solve(options, &problem, &summary);
  	std::cout << summary.FullReport() << std::endl;
		std::cout << "*) " << std::endl;
		computeFitErrPtsRAt(DIKProblem, index, false);

	}


	void SolveProblem(DikProblem* DIKProblem) {
		Problem problem;
		BuildProblem(DIKProblem, &problem);
 		Solver::Options options;
		SetSolverOptionsFromFlags(DIKProblem, &options);
		std::cout << "(* " << std::endl;
 	 	Solver::Summary summary;
  	Solve(options, &problem, &summary);
  	std::cout << summary.FullReport() << std::endl;
		std::cout << "*) " << std::endl;
		computeFitErrPtsR(DIKProblem);

	}


Eigen::MatrixXd FitEllipseModel(
	Eigen::MatrixXd inPts) {

//	Ellipse3D e3dUA(cartToHom(inPts));

}

}
/*

	// if s θ > 0
	void rotMatToEulerPos(Mat3 inRotMat, double* theta_i)	{
		Eigen::Map<Eigen::Matrix<double, 7, 1> > 
		theta_Vi(theta_i);

		theta_Vi(5,0) =  atan2(inRotMat(2,2)
			, sqrt(1.0 -(inRotMat(2,2)*inRotMat(2,2))));
		theta_Vi(4,0) = atan2( inRotMat(0,2)
			, inRotMat(1,2));
		theta_Vi(6,0) = atan2(
			-inRotMat(2,0), inRotMat(2,1));
	}

	// if s θ > 0
	void rotMatToEulerNeg(Mat3 inRotMat, double* theta_i)	{
		Eigen::Map<Eigen::Matrix<double, 7, 1> > 
		theta_Vi(theta_i);

		theta_Vi(5,0) =  atan2(inRotMat(2,2), -sqrt(1.0 -
			(inRotMat(2,2)*inRotMat(2,2))));
		theta_Vi(4,0) = atan2(
			-inRotMat(0,2), -inRotMat(1,2));
		theta_Vi(6,0) = atan2(
			inRotMat(2,0), -inRotMat(2,1));
	}

	void solveInverseEeOrientationAt (
		DikProblem* DIKProblem, int index) {
		const int theta_block_size =
			DIKProblem->theta_block_size();
		const int qBase_block_size =
			DIKProblem->qBase_block_size();
		double* thetas = DIKProblem->mutable_thetas();
		double* qsBase = DIKProblem->mutable_qsBase();
 		double* qbase_i = qsBase + qBase_block_size * index;
		Eigen::Map< Eigen::Quaternion<double> > 
			q_b_quat(qbase_i);
 	 	double* theta_i = thetas + theta_block_size * index;
		Eigen::Map<const Eigen::Matrix<double, 7, 1> > 
			theta_Vi(theta_i);

		ForwardKin<double> &wam = DIKProblem->getWam();
		size_t  nJoints = wam.getNJoints();
		wam.setThetaVect(theta_Vi);

		// Step 1: Using the joint variables θ1, θ2, θ3 θ4
		// 	 determined from inverse position, 
		//		evaluate R_4^base.		
			Mat3 rotMatj4InBase = wam.getMat(4).block(0,0,3,3);
		//		evaluate R_4^Vicon.		
		Mat3 rotMatBaseInVicon = 
			q_b_quat.normalized().toRotationMatrix();
		Mat3 rotMatj4InVicon = rotMatBaseInVicon
			* rotMatj4InBase;
	
		// Step 2: Find R_end-effector^Vicon
		// determined from vicon markers of the hand
			Mat3 rotMatEeInVicon = DIKProblem
				->buildWrRotReferential(index);
			Mat3 rotMatEeInBase = rotMatBaseInVicon.transpose()
				* rotMatEeInVicon;
		// Step 3: Find R_7^4
		// 	  =(R_4^Vicon)T R_ee^Vicon
  	Mat3 rotMatj7InJ4 = rotMatj4InBase.transpose() 
			* rotMatEeInBase;

//	printEigenMathematica(rotMatj7InJ4
//		, cout, "rotMatj7InJ4");	

		// Step 3: Find a set of Euler angles corresponding 
		// 	to the rotation matrix R_7^4 , and then use the 
		//	mapping θ5, θ6, θ7 ==  φ, θ, ψ
		// if s θ > 0
		rotMatToEulerPos(rotMatj7InJ4, theta_i);
		wam.setThetaVect(theta_Vi);
		// Step 4: fit error Ee
		//"A(θ1) * A(θ2) * ...* A(θ7) * Origin": -> pEe(i)				


		//Pee = Pwr + d6*Ree.[0 0 1]T
		_wrToEeVsWam.col(i) = l5wam_wrToEe * 
			buildWrRotReferential(i).col(2);
		_shoToEeVsWam.col(i) = _shoToWrVsWam.col(i) 
			+ _wrToEeVsWam.col(i);

		MatrixXd origin = MatrixXd::Zero(4,1);
		origin << 0.0, 0.0, 0.0, 1.0;
		Vec3 outShoToEeViWamInBase = homToCart(
			wam.getMat(7) * origin);

		Vec3 outShoToEeViWam;
		QuaternionRotatePoint(qbase_i,
			outShoToEeViWamInBase.data(),
				outShoToEeViWam.data());


//		Vec3 fitErrEe = (DIKProblem->
//			shoToEeVsWam().col(index) - outShoToEeViWam);

//		printEigenMathematica(fitErrEe.transpose()
//			, cout, "fitErrEe");	

	}
*/

/*
	
			if (theta_Vi(5,0) > DIKProblem->jointMinAngles()[5]
			 &&
				theta_Vi(5,0) < DIKProblem->jointMaxAngles()[5]) {
		
				if (theta_Vi(4,0) > 
					DIKProblem->jointMaxAngles()[4])
					theta_Vi(4,0) -= 2.0 * M_PI;		
					theta_Vi(6,0) = atan2(
					-rotMatj7InJ4(2,0), rotMatj7InJ4(2,1));
			}

			else {
				theta_Vi(4,0) = atan2(-rotMatj7InJ4(0,2)
					, -rotMatj7InJ4(1,2));
				theta_Vi(5,0) = atan2(rotMatj7InJ4(2,2)
					, -sqrt(1.0 -
					(rotMatj7InJ4(2,2)*rotMatj7InJ4(2,2))));
				theta_Vi(6,0) = atan2(rotMatj7InJ4(2,0)
					, -rotMatj7InJ4(2,1));
			}
*/

#endif

