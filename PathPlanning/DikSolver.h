#ifndef DIK_SOLVER_H
#define DIK_SOLVER_H


#include "BasicFormulas.h"
#include "SpatialJacobian.h"
#include "DikProblem.h"
#include "DifferentialIKErrorTerm.h"
#include "qBaseErrorTerm.h"
#include "UBCUtil.h"
#include "MayaAnimation.h"

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

DEFINE_double(fitErrNtolerance,  0.00001 , "fit error toleramce from DikProb m.");
 
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
		CHECK(StringToPreconditionerType
			( FLAGS_preconditioner
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

		options->function_tolerance = .000000000001e-90;
		options->parameter_tolerance = .0000000000001e-90;
		if (FLAGS_line_search) 
		  options->minimizer_type = ceres::LINE_SEARCH;

		CHECK(StringToTrustRegionStrategyType	(
			FLAGS_trust_region_strategy
			, &options->trust_region_strategy_type));
		CHECK(StringToDoglegType(FLAGS_dogleg, 
			&options->dogleg_type));
		options->use_inner_iterations 
			= FLAGS_inner_iterations;
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
		const int nObservations = 
			DIKProblem->nObservations();
		Eigen::MatrixXd
			shToUaVs =	DIKProblem->shoToUaVsWam(),
			elOffsetVs = DIKProblem->elOffsetVsWam(),
			wrOffsetVs = DIKProblem->wrOffsetVsWam(),
			laToWrVs = DIKProblem->laToWrVsWam(),
			shoToWrVs = DIKProblem->shoToWrVsWam();
		//info mat -> inverse of covariance (x y z)T
		Eigen::Matrix<double, 3, 3> 
			iMatShToUa =
				varianceSolve(shToUaVs).inverse(),
			iMatElOffset = 
				varianceSolve(elOffsetVs).inverse(),
			iMatWrOffset = 
				varianceSolve(wrOffsetVs).inverse(),
			iMatLaToWr = varianceSolve(laToWrVs).inverse(),
			iMatShToWr = varianceSolve(shoToWrVs).inverse();

		for (int i = 0; i < nObservations; i++)	{
    	// Each observation correponds to a pair of 
			// a qBase and a theta which are identified  
    	// by qBase_index()[i] and theta_index()[i]
    	// respectively.
   		double* qbase = qsBase + qBase_block_size * i;
   	 	double* theta = thetas + theta_block_size * i;
  
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
			problem->AddResidualBlock(cfPosition,
				loss_function, qbase, theta);
		

  	}
    // Set the limit for the theta parameter at  
		// position index in the thetas_i parameter block 

		std::vector<double> jointMinAngles = 
			DIKProblem->jointMinAngles();
		std::vector<double> jointMaxAngles = 
			DIKProblem->jointMaxAngles();

		for (int i = 0; i < nObservations; ++i) {
   	 	double* theta = thetas + theta_block_size * i;
			for (int joint = 0; joint < nJoints &&
				joint < theta_block_size; joint++) {
				problem->SetParameterLowerBound(theta, joint
					, jointMinAngles[joint]);
				problem->SetParameterUpperBound(theta, joint
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
			iMatElOffset = 
				varianceSolve(elOffsetVs).inverse()
			, iMatWrOffset =
				varianceSolve(wrOffsetVs).inverse()
			, iMatLaToWr = varianceSolve(laToWrVs).inverse()
			, iMatShToWr = 
				varianceSolve(shoToWrVs).inverse();

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
   		double *qbase = qsBase + qBase_block_size * index
				, *theta = thetas + theta_block_size * index
				, *qbasePrev = 
					qsBase + qBase_block_size * (index-1)
				,	*thetaPrev = 
					thetas + theta_block_size * (index-1);

//				for (int j = 0; j < qBase_block_size; j++)	
//					*qbase++ = *qbasePrev++;
//				for (int j = 0; j < theta_block_size; j++)	
//					*theta++ = *thetaPrev++;


			problem->AddResidualBlock(cfPosition,
				loss_function, qbase, theta);
  	
   		// Set the limit for the theta parameter  
			//  at position index in the thetas_i 
			// parameter block 
			std::vector<double> jointMinAngles = 
				DIKProblem->jointMinAngles();
			std::vector<double> jointMaxAngles = 
				DIKProblem->jointMaxAngles();
  		double* thetas_i = 
  			thetas + theta_block_size * index;
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
  		double* thetas_i = thetas + 
				theta_block_size * index + 4;
			for (int joint = 4; joint < nJoints &&
				joint < theta_block_size; joint++) {
				problem->SetParameterLowerBound(
					thetas_i, joint-4, jointMinAngles[joint]);
				problem->SetParameterUpperBound(
					thetas_i, joint-4, jointMaxAngles[joint]);
			}
	}

	MatrixXd outPtsWamAt (
		DikProblem* DIKProblem, int index) {

		//(0). initialize Optimizer parameters;
		const int nMarkers = 
			DIKProblem->nMarkers();

		//(0.1).qs - BASE ORIENTATION
		const int qBase_block_size =
			DIKProblem->qBase_block_size();
		double* qsBase = DIKProblem->mutable_qsBase();
 		double* qbase_i = qsBase + 
			qBase_block_size * index;
		Eigen::Map<const Eigen::Quaternion<double> > 
			q_b_quat(qbase_i);
		Mat3 rMatsBase = 
			q_b_quat.normalized().toRotationMatrix();

		//(0.2).θs -  WAM JOINTS ORIENTATIONS
		//			A(θ1), A(θ2), ..., A(θ_Njoints) 
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

		const int theta_block_size =
			DIKProblem->theta_block_size();
		double* thetas = DIKProblem->mutable_thetas();
	 	double* theta_i = thetas 
	 		+ theta_block_size * index;
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
		//"A(θ1) * A(θ2) * ...* A(θ5) * Origin"
		//				: -> pW(i)				
		Vec3 outShoToWrViWamInBase = homToCart(
			wam.getMat(5) * origin);
		//"A(θ1)*A(θ2)* ...* A(θ7) * Origin"
		//				: -> pEe(i)				
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

		Vec3 outShWam = DIKProblem->pShWam().col(index);

		MatrixXd outPtsWam = MatrixXd::Zero(nMarkers*3,1);

		outPtsWam.block(Markers::RSh,0,3,1)
			= outShWam;
		outPtsWam.block(Markers::RUA,0,3,1)
			=	outShWam + outShoToUaViWam;
		outPtsWam.block(Markers::REl,0,3,1)
			= outShWam + outShoToElViWam;
		outPtsWam.block(Markers::RLA,0,3,1)
			=	outShWam + outShoToLaViWam;
		outPtsWam.block(Markers::RWr,0,3,1)
			=	outShWam + outShoToWrViWam;
		outPtsWam.block(Markers::RTh,0,3,1)
			=	outShWam + outShoToEeViWam;
		outPtsWam.block(Markers::RPi,0,3,1)
			=	outShWam + outShoToEeViWam;

/*
	//outJoint angle
	printEigenMathematica(theta_Vi, cout, "outThetaWam" 
	+ to_string(index)	+ DIKProblem->phase());	
*/
		return outPtsWam;
	}
	
	void printOutPtsWam (DikProblem* DIKProblem,	
		MatrixXd rMatsCh ) {
		const int nObservations = 
			DIKProblem->nObservations();
		const int nMarkers = 
			DIKProblem->nMarkers();
		string phase = DIKProblem->phase();
		MatrixXd outShWam = DIKProblem->pShWam()
			, outUaVsWam = MatrixXd::Zero(3,nObservations)
			, outElVsWam = MatrixXd::Zero(3,nObservations)
			, outLaVsWam = MatrixXd::Zero(3,nObservations)
			, outWrVsWam = MatrixXd::Zero(3,nObservations)
			, outEeVsWam = MatrixXd::Zero(3,nObservations);


		MatrixXd outPtsWam_i 
			= MatrixXd::Zero(nMarkers*3,1);
		for (int i = 0; i < nObservations; i++)	{
			Mat3 rMatsCh_i = rMatsCh.block(0,i*3,3,3);
	
			outPtsWam_i = outPtsWamAt (DIKProblem,i);		
			outUaVsWam.col(i) = rMatsCh_i *
				outPtsWam_i.block(Markers::RUA,0,3,1);
			outElVsWam.col(i) = rMatsCh_i *
				outPtsWam_i.block(Markers::REl,0,3,1);
			outLaVsWam.col(i) = rMatsCh_i *
				outPtsWam_i.block(Markers::RLA,0,3,1);
			outWrVsWam.col(i) = rMatsCh_i *
				outPtsWam_i.block(Markers::RWr,0,3,1);
			outEeVsWam.col(i) = rMatsCh_i *
				outPtsWam_i.block(Markers::RTh,0,3,1);
		}

		printEigenMathematica( outEeVsWam.transpose()
			, cout, "outEeVsWam" + phase);	
		printEigenMathematica( outWrVsWam.transpose()
			, cout, "outWrVsWam" + phase);	
		printEigenMathematica( outLaVsWam.transpose()
			, cout, "outLaVsWam" + phase);	
		printEigenMathematica( outElVsWam.transpose()
			, cout, "outElVsWam" + phase);	
		printEigenMathematica( outUaVsWam.transpose()
			, cout, "outUaVsWam" + phase);			
	}	 
		 
	//Compute fit error at t(i)
	MatrixXd computeFitErrPtsRAt(
		DikProblem* DIKProblem, int index
			, bool solveOrientation) {

		//  outVecsWam predictions
		const int nMarkers = DIKProblem->nMarkers();
		MatrixXd outMarkerPtsWam_i =
		 MatrixXd ::Zero(nMarkers*3,1);
		outMarkerPtsWam_i = 
			outPtsWamAt(DIKProblem, index);

		Vec3 
			outShWam = DIKProblem->pShWam().col(index),
			outShoToUaViWam = outMarkerPtsWam_i.block(
				Markers::RUA,0,3,1) - outShWam,
			outShoToElViWam = outMarkerPtsWam_i.block(
				Markers::REl,0,3,1) - outShWam,
			outShoToLaViWam = outMarkerPtsWam_i.block(
				Markers::RLA,0,3,1) - outShWam,
			outShoToWrViWam = outMarkerPtsWam_i.block(
				Markers::RWr,0,3,1)- outShWam,
			outShoToEeViWam = outMarkerPtsWam_i.block(
				Markers::RTh,0,3,1)- outShWam;

		//Compute fit Error
		MatrixXd fitErrMarkerPts = 
			MatrixXd::Zero(nMarkers*3,1);

		Vec3 fitErrUA = DIKProblem->shoToUaVsWam()
				.col(index) - outShoToUaViWam;
		fitErrMarkerPts.block(
			Markers::RUA,0,3,1) =	fitErrUA;

		Vec3 fitErrEl = (DIKProblem->shoToElVsWam()
			 .col(index) - outShoToElViWam);
		fitErrMarkerPts.block(
			Markers::REl,0,3,1) =	fitErrEl;

		Vec3 fitErrLA = (DIKProblem->shoToLaVsWam()
				.col(index) - outShoToLaViWam);
		fitErrMarkerPts.block(
			Markers::RLA,0,3,1) =	fitErrLA;

		Vec3 fitErrWr = (DIKProblem->shoToWrVsWam()
				.col(index) - outShoToWrViWam);
		fitErrMarkerPts.block(
			Markers::RWr,0,3,1) =	fitErrWr;

		Vec3 fitErrEe = (DIKProblem->shoToEeVsWam()
				.col(index) - outShoToEeViWam);
		fitErrMarkerPts.block(
			Markers::RTh,0,3,1) =	fitErrEe;

/*
		printEigenMathematica(fitErrUA.transpose()
			, cout, "fitErrUA"+to_string(index));	
		printEigenMathematica(fitErrEl.transpose()
			, cout, "fitErrEl"+to_string(index));	
		printEigenMathematica(fitErrLA.transpose()
			, cout, "fitErrLA"+to_string(index));	
		printEigenMathematica(fitErrWr.transpose()
			, cout, "fitErrWr"+to_string(index));	
		printEigenMathematica(fitErrEe.transpose()
			, cout, "fitErrEe"+to_string(index));	
*/
		return fitErrMarkerPts;
	}



	void computeFitErrPtsR (
		DikProblem* DIKProblem) {

		const int nObservations = 
			DIKProblem->nObservations();
		const int nMarkers = DIKProblem->nMarkers();

		//Fit Errors between scaled inVecs & OutVecs of WAM
		MatrixXd fitErrMarkerPts =
 			MatrixXd::Zero(3*nMarkers,nObservations);
		MatrixXd fitErrMarkerPts_i =
 			MatrixXd::Zero(3*nMarkers,1);

		for (int i = 0; i < nObservations; i++)	{
			fitErrMarkerPts_i = 
				computeFitErrPtsRAt(DIKProblem, i, false);

		fitErrMarkerPts.col(i) = fitErrMarkerPts_i;

		}
		//Compute fit error 
		MatrixXd 
			fitErrUA = fitErrMarkerPts.block(
				Markers::RUA,0,3,nObservations),
			fitErrEl = fitErrMarkerPts.block(
				Markers::REl,0,3,nObservations),
			fitErrLA = fitErrMarkerPts.block(
				Markers::RLA,0,3,nObservations),
			fitErrWr = fitErrMarkerPts.block(
				Markers::RWr,0,3,nObservations),
			fitErrEe = fitErrMarkerPts.block(
				Markers::RTh,0,3,nObservations);


		cout << "fitErrUAMean = " 
			 << fitErrUA.colwise().norm().mean() << endl;
		cout << "fitErrElMean = " 
			 << fitErrEl.colwise().norm().mean() << endl;
		cout << "fitErrLAMean = " 
			 << fitErrLA.colwise().norm().mean() << endl;
		cout << "fitErrWrMean = " 
			 << fitErrWr.colwise().norm().mean() << endl;
		cout << "fitErrEeMean = " 
			 << fitErrEe.colwise().norm().mean() << endl;

		printEigenMathematica(fitErrUA.transpose()
			, cout, "fitErrUA" + DIKProblem->phase());	
		printEigenMathematica(fitErrEl.transpose()
			, cout, "fitErrEl"+ DIKProblem->phase());	
		printEigenMathematica(fitErrLA.transpose()
			, cout, "fitErrLA"+ DIKProblem->phase());	
		printEigenMathematica(fitErrWr.transpose()
			, cout, "fitErrWr"+ DIKProblem->phase());	
		printEigenMathematica(fitErrEe.transpose()
			, cout, "fitErrEe"+ DIKProblem->phase());	

	}

	
	void SolveProblemFPrevPt(
		DikProblem* DIKProblem) {

		const int theta_block_size =
			DIKProblem->theta_block_size();
		double* thetas = DIKProblem->mutable_thetas();

		const int qBase_block_size =
			DIKProblem->qBase_block_size();
		double* qsBase = DIKProblem->mutable_qsBase();

		const int nObservations = 
			DIKProblem->nObservations();
		const int nMarkers = 
			DIKProblem->nMarkers();
		MatrixXd fitErrMarkerPts_i = 
			MatrixXd::Zero(3*nMarkers,1);
		double fitErrUANorm_i = 0.0;

		//Max optimizer step size 
		int startErrPt = 0;
		while (fitErrUANorm_i < 0.00001 
			&& startErrPt < nObservations-5)	{
			startErrPt++;

			fitErrMarkerPts_i = computeFitErrPtsRAt(
				DIKProblem, startErrPt, false);
			fitErrUANorm_i = fitErrMarkerPts_i.block(
				Markers::RUA,0,3,1).norm();
			}


		//Min optimizer step size 
		fitErrUANorm_i	=	0.0;
		int finalErrPt = nObservations-1;
		while (fitErrUANorm_i < 0.00001 
			&& finalErrPt > 0)	{
		
			fitErrMarkerPts_i = computeFitErrPtsRAt(
				DIKProblem, finalErrPt, false);
			fitErrUANorm_i = fitErrMarkerPts_i.block(
				Markers::RUA,0,3,1).norm();
			finalErrPt--;
		}


cout << "startErrPt: " << startErrPt << endl;
cout << "finalErrPt: " << finalErrPt << endl;

		for (int i = 2; i < finalErrPt; i++)	{	
	 		 	double* theta_i = 
	 		 		thetas + theta_block_size * i;
				double* qbase_i = 
					qsBase + qBase_block_size * i;	

				fitErrMarkerPts_i = 
					computeFitErrPtsRAt(DIKProblem, i,  false);
				fitErrUANorm_i = fitErrMarkerPts_i.block(
				Markers::RUA,0,3,1).norm();

				//A- SOLVE Pt_i FROM Pt_(i-1)	
				//Max optimizer step size 
				int stepSize = i;
				while (fitErrUANorm_i > 0.001 
					&& stepSize > i-5 && stepSize > -1)	{		
					stepSize--;

					double* qbaseIMOne = qsBase + 
						qBase_block_size * 	(stepSize);
				 	double* theta_iMone = thetas + 
						theta_block_size * (stepSize);
					for (int j = 0; j < qBase_block_size; j++)	
						*qbase_i++ = *qbaseIMOne++;
					for (int j = 0; j < theta_block_size; j++)	
						*theta_i++ = *theta_iMone++;
				
					Problem problem;
					BuildProblemAt(DIKProblem, &problem, i);
			 		Solver::Options options;
					SetSolverOptionsFromFlags(
						DIKProblem, &options);
					options.minimizer_progress_to_stdout = false;
			 	 	Solver::Summary summary;
					Solve(options, &problem, &summary);

					fitErrMarkerPts_i = 
						computeFitErrPtsRAt(DIKProblem, i, false);
					fitErrUANorm_i = fitErrMarkerPts_i.block(
					Markers::RUA,0,3,1).norm();
				}
		}

		for (int i = nObservations-1; 
			i > startErrPt; i--)	{	
			double* theta_i = 
			 	thetas + theta_block_size * i;
			double* qbase_i = 
				qsBase + qBase_block_size * i;	

			fitErrMarkerPts_i = 
				computeFitErrPtsRAt(DIKProblem, i,  false);
			fitErrUANorm_i = fitErrMarkerPts_i.block(
				Markers::RUA,0,3,1).norm();

		  int stepSize = i;
			while (fitErrUANorm_i > 0.01 
				&& stepSize < i+5
				&& stepSize < nObservations) {		
														
				stepSize++;
				double* qbaseIPone = qsBase + 
					qBase_block_size * 	(stepSize);
			 	double* theta_iPone = thetas + 
					theta_block_size * (stepSize);
				for (int j = 0; j < qBase_block_size; j++)	
					*qbase_i++ = *qbaseIPone++;
				for (int j = 0; j < theta_block_size; j++)	
					*theta_i++ = *theta_iPone++;

				Problem problem;
				BuildProblemAt(DIKProblem, &problem, i);
		 		Solver::Options options;
				SetSolverOptionsFromFlags(
					DIKProblem, &options);
				options.minimizer_progress_to_stdout = false;
		 	 	Solver::Summary summary;
				Solve(options, &problem, &summary);
//				computeFitErrPtsRAt(DIKProblem, i, true);


				fitErrMarkerPts_i = 
					computeFitErrPtsRAt(DIKProblem, i, false);
				fitErrUANorm_i = fitErrMarkerPts_i.block(
				Markers::RUA,0,3,1).norm();
			}
	}
}







	void BuildProbFPrevPtWqBase	(
		DikProblem* DIKProblem
		, Problem* problem) {		

		const int theta_block_size =
			DIKProblem->theta_block_size();
		double* thetas = DIKProblem->mutable_thetas();

		const int qBase_block_size =
			DIKProblem->qBase_block_size();
		double* qsBase = DIKProblem->mutable_qsBase();

		const int nObservations = 
			DIKProblem->nObservations();
		const int nMarkers = 
			DIKProblem->nMarkers();
		MatrixXd fitErrMarkerPts_i = 
			MatrixXd::Zero(3*nMarkers,1);
		double fitErrUANorm_i = 0.0;

		//Max optimizer step size 
		int firstPtWzeroErr = 0;
		fitErrMarkerPts_i = computeFitErrPtsRAt(
			DIKProblem, firstPtWzeroErr, false);
		fitErrUANorm_i = fitErrMarkerPts_i.block(
			Markers::RUA,0,3,1).norm();
		while (fitErrUANorm_i > FLAGS_fitErrNtolerance
			&& firstPtWzeroErr < nObservations)	{
			firstPtWzeroErr++;

			fitErrMarkerPts_i = computeFitErrPtsRAt(
				DIKProblem, firstPtWzeroErr, false);
			fitErrUANorm_i = fitErrMarkerPts_i.block(
				Markers::RUA,0,3,1).norm();
			}
	
		//Min optimizer step size 
		int LastPtWzeroErr = firstPtWzeroErr; 
		fitErrUANorm_i	=	0.0;
		fitErrMarkerPts_i = computeFitErrPtsRAt(
			DIKProblem, LastPtWzeroErr, false);
		fitErrUANorm_i = fitErrMarkerPts_i.block(
			Markers::RUA,0,3,1).norm();
	while (fitErrUANorm_i < FLAGS_fitErrNtolerance
			&& LastPtWzeroErr < nObservations)	{

			LastPtWzeroErr++;
			fitErrMarkerPts_i = computeFitErrPtsRAt(
				DIKProblem, LastPtWzeroErr, false);
			fitErrUANorm_i = fitErrMarkerPts_i.block(
				Markers::RUA,0,3,1).norm();
		}


cout << "firstPtWzeroErr: " << firstPtWzeroErr << endl;
cout << "LastPtWzeroErr: " << LastPtWzeroErr << endl;
/*
	for (int i = firstPtWzeroErr + 1;
		i > 0; 	i--)	{	
	 		 	double* theta_i = 
	 		 		thetas + theta_block_size * i;
				double* qbase_i = 
					qsBase + qBase_block_size * i;	

				fitErrMarkerPts_i = 
					computeFitErrPtsRAt(DIKProblem, i,  false);
				fitErrUANorm_i = fitErrMarkerPts_i.block(
				Markers::RUA,0,3,1).norm();

			//A- SOLVE Pt_i FROM Pt_(i-1)	
				//Max optimizer step size 
				int stepSize = i;
				while (fitErrUANorm_i > 0.001 
					&& stepSize > i-5 && stepSize > -1)	{		
					stepSize--;

					double* qbaseIMOne = qsBase + 
						qBase_block_size * 	(stepSize);
				 	double* theta_iMone = thetas + 
						theta_block_size * (stepSize);
					for (int j = 0; j < qBase_block_size; j++)	
						*qbase_i++ = *qbaseIMOne++;
					for (int j = 0; j < theta_block_size; j++)	
						*theta_i++ = *theta_iMone++;
				
					Problem problem;
					BuildProblemAt(DIKProblem, &problem, i);
			 		Solver::Options options;
					SetSolverOptionsFromFlags(
						DIKProblem, &options);
					options.minimizer_progress_to_stdout = false;
			 	 	Solver::Summary summary;
					Solve(options, &problem, &summary);

					fitErrMarkerPts_i = 
						computeFitErrPtsRAt(DIKProblem, i, false);
					fitErrUANorm_i = fitErrMarkerPts_i.block(
					Markers::RUA,0,3,1).norm();
				}
		}

		for (int i = nObservations-1; 
			i > firstPtWzeroErr; i--)	
		{	

			 	double* theta_i = 
			 		thetas + theta_block_size * i;
				double* qbase_i = qsBase 
					+ qBase_block_size * i;	

				fitErrMarkerPts_i = 
					computeFitErrPtsRAt(DIKProblem, i,  false);
				fitErrUANorm_i = fitErrMarkerPts_i.block(
				Markers::RUA,0,3,1).norm();

		  int stepSize = i;
			while (fitErrUANorm_i > 0.01 
				&& stepSize < i+5 
				&& stepSize < nObservations)	
				{		
														
				stepSize++;

				double* qbaseIPone = qsBase + 
					qBase_block_size * 	(stepSize);
			 	double* theta_iPone = thetas + 
					theta_block_size * (stepSize);
				for (int j = 0; j < qBase_block_size; j++)	
					*qbase_i++ = *qbaseIPone++;
				for (int j = 0; j < theta_block_size; j++)	
					*theta_i++ = *theta_iPone++;

				Problem problem;
				BuildProblemAt(DIKProblem, &problem, i);
		 		Solver::Options options;
				SetSolverOptionsFromFlags(
					DIKProblem, &options);
				options.minimizer_progress_to_stdout = false;
		 	 	Solver::Summary summary;
				Solve(options, &problem, &summary);
//				computeFitErrPtsRAt(DIKProblem, i, true);


				fitErrMarkerPts_i = 
					computeFitErrPtsRAt(DIKProblem, i, false);
				fitErrUANorm_i = fitErrMarkerPts_i.block(
				Markers::RUA,0,3,1).norm();
			}
	}
	*/
}





		 	 	
   	 	// ADD RESIDUAL for interpolation between
   	 	// series of BASE QUATERNION rotations.
/*
   	 	
			if (i < nObservations-2)	{
			
		 		double* qbaseNext = qsBase 
   	 			+ qBase_block_size * (i+1);
   	 				
*//*
   			CostFunction* cfQbase;
   		
   			double* qbasePrev = qsBase 
   	 			+ qBase_block_size * (i-1);
   

   	 			
				Eigen::Map< Eigen::Quaternion<double> > 
					q_b_quat(qbase);
   	 		
   			cfQbase = qBaseErrorTerm::Create (q_b_quat);
   			problem->AddResidualBlock(cfQbase, NULL
					, qbasePrev, qbase, qbaseNext);
*//*
				problem->AddResidualBlock( new 
					AutoDiffCostFunction<qBaseDeltaErrorTerm, 3
						, 4, 4>(new qBaseDeltaErrorTerm), nullptr
							,qbase, qbaseNext);
			}
*/






/*

	void SolveProblemFPrevPt(
		DikProblem* DIKProblem, int startErrPt) {

		const int nObservations = DIKProblem->nObservations();
		const int theta_block_size =
			DIKProblem->theta_block_size();
		double* thetas = DIKProblem->mutable_thetas();

		const int qBase_block_size =
			DIKProblem->qBase_block_size();
		double* qsBase = DIKProblem->mutable_qsBase();
		ForwardKin<double> &wam = DIKProblem->getWam();
		size_t  nJoints = wam.getNJoints();

		MatrixXd wamShtoUaTrans =  MatrixXd::Zero(4,1)
		, outShoToUaVWam = MatrixXd::Zero(3,1);
		wamShtoUaTrans << 
			(0.0), (0.0), (0.55), (1.0);


		for (int i = startErrPt; i > -1; i--)	{	

	 	 	double* theta_i = thetas + theta_block_size * i;
			double* qbase_i = qsBase + qBase_block_size * i;	

			MatrixXd outShoToUaVWamInBase =
				homToCart(wam.getMat(2) * wamShtoUaTrans);

			QuaternionRotatePoint(qbase_i,
				outShoToUaVWamInBase.data(),
					outShoToUaVWam.data());

			double fitErrUANorm_i = 
				(DIKProblem->shoToUaVsWam().col(i)
					- outShoToUaVWam).norm();

			//Max optimizer step size 
			int stepSize = i;
			while (fitErrUANorm_i > 0.01 && stepSize < 5)	{		
														
				double* qbaseIMOne = qsBase + 
					qBase_block_size * 	(stepSize+1);
			 	double* theta_iMone = thetas + 
					theta_block_size * (stepSize+1);
				for (int j = 0; j < qBase_block_size; j++)	
					*qbase_i++ = *qbaseIMOne++;
				for (int j = 0; j < theta_block_size; j++)	
					*theta_i++ = *theta_iMone++;

				Problem problem;
				BuildProblemAt(DIKProblem, &problem, i);
		 		Solver::Options options;
				SetSolverOptionsFromFlags(DIKProblem, &options);
//				std::cout << "(* " << std::endl;
				options.minimizer_progress_to_stdout = false;
		 	 	Solver::Summary summary;
				Solve(options, &problem, &summary);
//				std::cout << summary.FullReport() << std::endl;
//				std::cout << "*) " << std::endl;
//				computeFitErrPtsRAt(DIKProblem, i, true);
					fitErrUANorm_i = 
						(DIKProblem->shoToUaVsWam().col(i)
					- outShoToUaVWam).norm();
					stepSize++;
			}
		}

	}




	void SolveProblemFNextPt(
		DikProblem* DIKProblem, int finalErrPt) {

		const int nObservations = DIKProblem->nObservations();
		const int theta_block_size =
			DIKProblem->theta_block_size();
		double* thetas = DIKProblem->mutable_thetas();

		const int qBase_block_size =
			DIKProblem->qBase_block_size();
		double* qsBase = DIKProblem->mutable_qsBase();
		ForwardKin<double> &wam = DIKProblem->getWam();
		size_t  nJoints = wam.getNJoints();

		MatrixXd wamShtoUaTrans =  MatrixXd::Zero(4,1)
		, outShoToUaVWam = MatrixXd::Zero(3,1);
		wamShtoUaTrans << 
			(0.0), (0.0), (0.55), (1.0);


		for (int i = finalErrPt; i < nObservations; i++)	{	

	 	 	double* theta_i = thetas + theta_block_size * i;
			double* qbase_i = qsBase + qBase_block_size * i;	
			MatrixXd outShoToUaVWamInBase =
				homToCart(wam.getMat(2) * wamShtoUaTrans);
			QuaternionRotatePoint(qbase_i,
				outShoToUaVWamInBase.data(),
					outShoToUaVWam.data());

			double fitErrUANorm_i = 
				(DIKProblem->shoToUaVsWam().col(i)
					- outShoToUaVWam).norm();
FN
			if (fitErrUANorm_i > 0.01)	{		
														
				double* qbaseIpOne = qsBase + 
					qBase_block_size * 	(i-1);
			 	double* theta_ipone = thetas + 
					theta_block_size * (i-1);
				for (int j = 0; j < qBase_block_size; j++)	
					*qbase_i++ = *qbaseIpOne++;
				for (int j = 0; j < theta_block_size; j++)	
					*theta_i++ = *theta_ipone++;

				Problem problem;
				BuildProblemAt(DIKProblem, &problem, i);
		 		Solver::Options options;
				SetSolverOptionsFromFlags(DIKProblem, &options);
//				std::cout << "(* " << std::endl;
				options.minimizer_progress_to_stdout = false;		
		 	 	Solver::Summary summary;
				Solve(options, &problem, &summary);
//				std::cout << summary.FullReport() << std::endl;
//				std::cout << "*) " << std::endl;
//				computeFitErrPtsRAt(DIKProblem, i, true);
			}
		}

	}


*/

	// Output the BaseQuats with format: 
	// q_x q_y q_z q_w.
	Eigen::MatrixXd getOutQuatsBase(
		DikProblem* DIKProblem) {

		const int nObservations = 
			DIKProblem->nObservations();
		Eigen::MatrixXd outQuatsBase = 
			Eigen::MatrixXd::Zero(4,nObservations);

		for (int i = 0; i < nObservations; i++)	{
			Eigen::Quaternion<double> qBaseI 
				= DIKProblem->qbaseQuatAt(i);
				outQuatsBase(0,i) = qBaseI.x();
				outQuatsBase(1,i) = qBaseI.y();
				outQuatsBase(2,i) = qBaseI.z();
				outQuatsBase(3,i) = qBaseI.w();
		}
		return 	outQuatsBase;
	}

	Eigen::MatrixXd getOutThetasWam(
		DikProblem* DIKProblem) {

		const int nObservations = 
			DIKProblem->nObservations();
		Eigen::MatrixXd outThetasWam = 
			Eigen::MatrixXd::Zero(7,nObservations);

		ForwardKin<double> &wam = DIKProblem->getWam();

		const int theta_block_size =
			DIKProblem->theta_block_size();
		double* thetas = DIKProblem->mutable_thetas();
		size_t  nJoints = wam.getNJoints();

		for (int i = 0; i < nObservations; i++)	{
   	 	double* theta_i = thetas + theta_block_size * i;
    	Eigen::Map<const Eigen::Matrix<double, 7, 1> > 
				theta_Vi(theta_i);
			outThetasWam.col(i) = theta_Vi;
		}
		return 	outThetasWam;
	}

	void SolveOrientationProblemAt(DikProblem* DIKProblem
		, int index) {
		Problem problem;
		BuildOrientationProblemAt(
			DIKProblem, &problem, index);
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


	void SolveProblemAt(
		DikProblem* DIKProblem, int index) {
		
		Problem problem;
		BuildProblemAt(DIKProblem, &problem, index);
 		Solver::Options options;
		SetSolverOptionsFromFlags(DIKProblem, &options);
//		std::cout << "(* " << std::endl;
		options.minimizer_progress_to_stdout = false;
 	 	Solver::Summary summary;
  	Solve(options, &problem, &summary);
//  	std::cout << summary.FullReport() << std::endl;
//		std::cout << "*) " << std::endl;
//		computeFitErrPtsRAt(DIKProblem, index, false);

	}


	void SolveProblem(DikProblem* DIKProblem) {
		Problem problem;
		BuildProblem(DIKProblem, &problem);
 		Solver::Options options;
		SetSolverOptionsFromFlags(DIKProblem, &options);
//		std::cout << "(* " << std::endl;
		options.minimizer_progress_to_stdout = false;
 	 	Solver::Summary summary;
  	Solve(options, &problem, &summary);
//  	std::cout << summary.FullReport() << std::endl;
//		std::cout << "*) " << std::endl;
		computeFitErrPtsR(DIKProblem);

	}


Eigen::MatrixXd FitEllipseModel(
	Eigen::MatrixXd inPts) {

//	Ellipse3D e3dUA(cartToHom(inPts));

}

}

/*


	void SolveProblemFNextPt(
		DikProblem* DIKProblem) {
		const int theta_block_size =
			DIKProblem->theta_block_size();
		double* thetas = DIKProblem->mutable_thetas();

		const int qBase_block_size =
			DIKProblem->qBase_block_size();
		double* qsBase = DIKProblem->mutable_qsBase();

		const int nObservations = DIKProblem->nObservations();
		const int nMarkers = DIKProblem->nMarkers();
		MatrixXd fitErrMarkerPts_i = MatrixXd::Zero(3*nMarkers,1);
		double fitErrMarkersNorm_i = 0.0;


		//Max optimizer step size 
		int finalErrPt = 0;
		while (fitErrMarkersNorm_i < 0.00001 
			&& finalErrPt > 0)	{
		
			fitErrMarkerPts_i = 
				computeFitErrPtsRAt(DIKProblem, finalErrPt, false);
			double fitErrMarkersNorm_i = 
				fitErrMarkerPts_i.norm();
			finalErrPt--;

		}

	}
*/
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

