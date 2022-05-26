#ifndef Q_BASE_ERROR_TERM_H_
#define Q_BASE_ERROR_TERM_H_
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

// Computes the error term for two base quaternion 
// orientations that have a relative pose measurement
// between them. 
//  Let the hat variables be the measurement. 
//  We have two Qs x_a and x_b. 
//  Through DIKsolver we can measure the quat orientation
//  of frame B w.r.t frame A denoted as q_ab_hat. 
//  We can compute an erro metric between the current 
//  estimate of the qBase and the measurement.
//
// In this formulation, we have chosen to represent the 
// base orientation as a Hamiltonian quaternion, q, and 
// The quaternion ordering is [x, y, z, w].
// 
// The estimated measurement is:
//             [ q_ab ]    [ q_a^{-1] * q_b         ]
//
// where ^{-1} denotes the inverse and R(q) is the 
// rotation matrix for the quaternion.
//  Now we can compute an error metric between the 
//  estimated and measurement transformation. For the 
//  orientation error, we will use the standard  
//  multiplicative error resulting in:
//
//   error = [ 2.0 * Vec(q_ab * \hat{q}_ab^{-1}) ]
//
// where Vec(*) returns the vector (imaginary) part of 
// the quaternion. 

// Since the measurement has an uncertainty associated 
// with how accurate it is, we will weight the errors 
// by the square root of the measurement information
// matrix:
//
//   residuals = I^{1/2) * error
//
// where I is the information matrix which is the 
// inverse of the covariance matrix for the delta 
// orientation measurements.


struct qBaseErrorTerm {
 public:
  qBaseErrorTerm ( 
  	Eigen::Quaternion<double> deta_qBase_measured
//  	, const Eigen::Matrix<double, 7, 7>& 
//  			sqrt_information
  			) : 
  			
  	_deta_qBase_measured(deta_qBase_measured)
// 	, _sqrt_information(sqrt_information) 
			{}



  template <typename T> bool operator()(
  	const T* const q_Iminus1_ptr
  	, const T* const q_I_ptr
  	, const T* const q_Iplus1_ptr
    , T* residuals_ptr) const {
                  
    Eigen::Map<const Eigen::Quaternion<T>> 
    	q_Iminus1(q_Iminus1_ptr), q_I(q_I_ptr)
    	, q_Iplus1(q_Iplus1_ptr);

	 	// F = l(qi).||k(qi)||^2 + c[g(qi)]^2
		//			(A) .   (B)      +    (C)
   	// the Energy function that is to be minimized
	
    
    // (A). l_qi = f(q_Iminus1, q_i, q_Iplus1):
    // The "parameter width" of the interval we 
    // integrate over is termed l ( q i ). It can be
    // expressed as a centered average of the parameter
    // width in the intervals immediately before and 
    // after the i 'th quaternion.
   
    // returns the angle (in radian) 
		// 	between two rotations
		T angularDistanceFprev 
			= q_I.angularDistance(q_Iminus1);
		T angularDistanceFnext 
			= q_I.angularDistance(q_Iplus1);
		T l_qi = (angularDistanceFprev 
			+ angularDistanceFprev)/T(2.0);
			
   
   
   	// (B). k_DiscLocCurvature = f(q_i,q_i_dotDot):
   	// The discrete version of the local curvature:
   	
		// q_i_dotDot = f(q_Iminus1, q_i, q_Iplus1):
    // The second derivative of the discrete  
    // approximation of the interpolation curve
   	Eigen::Matrix<T, 4, 1> qIpp = (q_Iminus1.coeffs() 
   		- T(2.0)*q_I.coeffs() + q_Iplus1.coeffs())
   			/(l_qi*l_qi);
		//dot products
   	T qIpp_dot_qI 
   		= qIpp(0,0) * q_I_ptr[0] 
   		+ qIpp(1,0) * q_I_ptr[1] 
   		+ qIpp(2,0) * q_I_ptr[2] 
   		+ qIpp(3,0) * q_I_ptr[3];
 		T qI_dot_qI  
 			= q_I_ptr[0] * q_I_ptr[0] 
   		+ q_I_ptr[1] * q_I_ptr[1] 
   		+ q_I_ptr[2] * q_I_ptr[2] 
   		+ q_I_ptr[3] * q_I_ptr[3];
 		// OPTION B. A.A
		//	T qIDotqI = q_I.norm() * q_I.norm(); 

 		
   // k_DiscLocCurvature(4x1) = f(q_i,q_i_dotDot):
   Eigen::Matrix<T, 4, 1> k_qi = qIpp 
   	- (qIpp_dot_qI/qI_dot_qI) * q_I.coeffs();
   
   	// (C). c[g(qi)]^2 where
   	// g(q) = qq-1 (4x1) 
   	// The measure for determining if the  
   	// quaternions are unit quaternions
   	// Assuming c is suitably large, the energy 
   	// function F will have a minimum approximately
   	// where E has a minimum, and where all q i are
   	// approximately unit quaternions
   	
   	Eigen::Quaternion<T> q_I_inverse 
			= q_I.conjugate();
		Eigen::Quaternion<T> g_qi = q_I * q_I_inverse;

    // Compute the residuals.
    // [ orientation (3x1)] = [ 2 * delta_q(0:2) ]
    Eigen::Map<Eigen::Matrix<T, 3, 1>> 
    	residuals(residuals_ptr);
    	residuals.template block<3, 1>(0, 0) 
    		= T(2.0) * g_qi.vec() + l_qi  
    		* k_qi.template block<3, 1>(0, 0);
//    residuals.template block<4, 1>(0, 0) 
//   	= l_qi * k_qi;
//    residuals.template block<3, 1>(4, 0) 
//    	= T(2.0) * g_qi.vec();
    	

    // Scale residuals by the measurement uncertainty.
//   residuals.applyOnTheLeft(
//    	_sqrt_information.template cast<T>());
   
    return true;
  }
  
  
  static ceres::CostFunction* Create(
  	Eigen::Quaternion<double> deta_qBase_measured
//  	, const Eigen::Matrix<double,  7, 7>& 
//  			sqrt_information
  			) {
    
    return new ceres::AutoDiffCostFunction<
    	qBaseErrorTerm, 7, 4, 4, 4>(
        new qBaseErrorTerm(deta_qBase_measured
//        	, sqrt_information
        	));
  }
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
 private:
  // The measurement for the change in qBase at 
  // frame B relative to A expressed in the A frame.
  Eigen::Quaternion<double> _deta_qBase_measured;
//  const Eigen::Matrix<double,  7, 7> _sqrt_information;
  
 };
 
 
	struct qBaseDeltaErrorTerm {


  template <typename T> bool operator()(
  	const T* const q_Iminus1_ptr
  	, const T* const q_I_ptr
    , T* residuals_ptr) const {
                           
    Eigen::Map<const Eigen::Quaternion<T>> 
    	q_Iminus1(q_Iminus1_ptr), q_I(q_I_ptr);
/*
    	, q_Iplus1(q_Iplus1_ptr);
    // (A). returns the angle (in radian) 
		// 	between two rotations
		T angularDistanceFprev 
			= q_I.angularDistance(q_Iminus1);
		T angularDistanceFnext 
			= q_I.angularDistance(q_Iplus1);
*/			
   	// (B). Compute the relative transformation
   	// between the two frames.
   	Eigen::Quaternion<T> q_Iminus1_inverse 
   		= q_Iminus1.conjugate();
    Eigen::Quaternion<T> delta_q   
    	= q_Iminus1_inverse * q_I;
    // Compute the residuals.
    // [ orientation (3x1)] = [ 2 * delta_q(0:2) ]
    Eigen::Map<Eigen::Matrix<T, 3, 1>> 
    	residuals(residuals_ptr);
    residuals.template block<3, 1>(0, 0) =
    	T(2.0) * delta_q.vec();
	}
	
	};
 }
#endif  
