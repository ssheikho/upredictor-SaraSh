#ifndef BASIC_JACOBIAN_H
#define BASIC_JACOBIAN_H

#include "DH.h"
#include "ForwardKin.h"

#include <cmath>

#include <Eigen/SVD>

#include <cmath>

using namespace Eigen;
/* The coordinates for z_i with respect to the base frame
		are given by the first three elements in the third column
		of ^0T_i .
	while o_i is given by the first three elements of the
		fourth column of ^0T_i . 
*/
template<typename T>
class BasicJacobian {

public:
	BasicJacobian(ForwardKin<T> &fk) : _fk(fk) {
		computeJacobian();
	}

	~BasicJacobian() {
	}

	Eigen::Matrix<T, 6, 7> getJacobianMat()	{
		return _jacobianMat;
	}

protected:
	void computeJacobian()	{
		Eigen::Matrix<T, 3, 1> o_n = 
			_fk.getMat(_fk.getNJoint()).block(0,3,3,1);
		for (size_t i = 0; i < _fk.getNJoints(); i++)	{
			// (.1) rotation axis z_i
			Eigen::Matrix<T, 3, 1> z_iMinus1 =
				_fk.getMat(i).block(0,2,3,1);
			// (.2) position vector
			Eigen::Matrix<T, 3, 1> o_iMinus1 
				= _fk.getMat(i).block(0,3,3,1);

			_jacobianMat.block(0, i, 3, 1) 
				= z_iMinus1.cross(o_n-o_iMinus1);
			_jacobianMat.block(3,i,3,1) = z_iMinus1;
		}
	}

	ForwardKin<T> &_fk;
	Eigen::Matrix<T, 6, 7> _jacobianMat;

};
#endif
