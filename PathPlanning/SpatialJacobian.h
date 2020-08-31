#ifndef SPATIAL_JACOBIAN_H
#define SPATIAL_JACOBIAN_H

#include "BasicJacobian.h"

#include <cmath>

#include <Eigen/SVD>

#include <cmath>

using namespace Eigen;
/* The mapping from generalized velocities u to the 	
	 	operational space twist [c_V_e c_ω_ce] of frame Q is
	  realized by the spatial Jacobian:
*/
template<typename T>
class SpatialJacobian {

public:
	SpatialJacobian(BasicJacobian<T> &basicJMat
		, Eigen::Matrix<T,3,1> inShoToWrInBaseF
		, Eigen::Matrix<T, 3, 3> rotFromBaseToCameraF) 
		: _basicJMat(basicJMat)
		, _inShoToWrInBaseF(inShoToWrInBaseF) 
		, _rotFromBaseToCameraF(rotFromBaseToCameraF){}

	~SpatialJacobian() {
	}

	Eigen::Matrix<T,  6, 13> getJacobianMat()	{
		return _spatialJMat;
	}

//the unit vectors of B expressed in camera frame 
protected:
	void computeJacobian()	{
		//column 1
		_spatialJMat.block(0,0,3,3) 
			= Eigen::Matrix<T, 3, 3>::Identity(3,3);
		_spatialJMat.block(3,0,3,3) 
			= Eigen::Matrix<T, 3, 3>::Zero(3,3);
		//column 2
		Eigen::Matrix<T, 3, 3> inShoToWrInBaseFssMat 
			= makeSkewSymmetricMatrix(_inShoToWrInBaseF);
		_spatialJMat.block(0,3,3,3) 
			= -_rotFromBaseToCameraF * inShoToWrInBaseFssMat;
		_spatialJMat.block(3,3,3,3) 
			= _rotFromBaseToCameraF;
		//column 3
		_spatialJMat.block(0,6,6,7) 
			= _rotFromBaseToCameraF * _basicJMat.getJacobianMat();
	}

	BasicJacobian<T> &_basicJMat;
	Eigen::Matrix<T, 3, 1> _inShoToWrInBaseF;
	// Rotation matrix that maps robots base frame 
	// '0' to camera frame 
	Eigen::Matrix<T, 3, 3> _rotFromBaseToCameraF;
	// 6 rows: 3 for J_p (position), and 3 for J_R (rotation)
	// 13 columns: 7(joints) + 3(c_V_b) + 3(b_ω_cb)
	Eigen::Matrix<T, 6, 13> _spatialJMat;


};
#endif
