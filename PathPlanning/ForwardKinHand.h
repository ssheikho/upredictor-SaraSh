#ifndef FORWARD_KIN_HAND_H
#define FORWARD_KIN_HAND_H

#include "DhHand.h"

#include <cmath>

#include <Eigen/SVD>

#include <cmath>

using namespace Eigen;
/** Each finger is considered its own manipulator and is referenced to a wrist coordinate frame in the center of the palm. Use the forward kinematics calculated in this section to determine fingertip position and orientation with respect to the palm.
**/ 
template<typename T>
class ForwardKinHand {

public:
	ForwardKinHand(size_t nJoints, size_t nFingers) 
		: _nJoints(nJoints)
		, _nFingers(nFingers)
		, _joint(new DhHand<T>[_nFingers * _nJoints]) {
	}

	~ForwardKinHand() {
		delete[] _joint;
	}

	void setThetaVect(MatrixXd thetaVect) {
		for(size_t i = 0; i < _nJoints; i++)
			_joint[i].setTheta(thetaVect(i,0));
	}

 	//“r” is either [-1,1,0] for [F1,F2,F3] respectively.	
 	//“j” is either [1,1,-1] for [F1,F2,F3] respectively. 
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		getAvectFinger(size_t fingerIndex) {
		Eigen::Matrix<T, 4, 1> a;
		T r = T(0.0);
		switch (fingerIndex) {
			case 0:
				r = T(-1.0);
				break;
			case 1:
				r = T(1.0);
				break;
			case 2:
				r = T(0.0);
				break;
		}
		a << T(r * 0.025), T(0.05), T(0.07), T(0.05);
		return a;
	}

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		getAlphaVect() {
		Eigen::Matrix<T, 4, 1> alpha;
		alpha << T(0.0), T(M_PI/2.0), T(0.0), T(-M_PI/2.0);
		return alpha;
	}

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		getDvect() {
		Eigen::Matrix<T, 4, 1> d;
		d << T(0.084), T(0.0), T(0.0), T(0.0095);
		return d;
	}

	void setDhHandParameters(
			Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> a
			, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> alpha
			, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> d
			, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> theta
			) {
		for(size_t i = 0; i < _nJoints; i++)
			_joint[i].setValues(a(i,0), alpha(i,0)
				, d(i,0), theta(i,0));
	}

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		getMat(size_t toDOF) {

		Eigen::Matrix<T, 4,4> retVal =
			Eigen::Matrix<T, 4, 4>::Identity(4,4);

		for(size_t i = 0; (i < _nJoints) && (i < toDOF); i++)
			retVal = retVal * _joint[i].getMat();
	
		return retVal;
	}

	DhHand<T> &getDhHand(size_t index) {
		return _joint[index];
	}

	size_t getNJoints() {
		return _nJoints;
	}

	static ForwardKinHand generateBarrettHand() {
		ForwardKinHand retVal(12);
		Eigen::Matrix<T, 12, 1> a, alpha, d, theta;
		T r = T(0.0);
		T j = T(0.0);
		for (int f = 0; f < 3; f++) {
			switch (f) {
				case 0:
					r = T(-1.0);
					j = T(1.0);
					break;
				case 1:
					r = T(1.0);
					j = T(1.0);
					break;
				case 2:
					r = T(0.0);
					j = T(-1.0);
					break;
			}
			a << T(r * 0.025), T(0.05), T(0.07), T(0.05);
			alpha << T(0.0), T(M_PI/2.0), T(0.0), T(-M_PI/2.0);
			d << T(0.084), T(0.0), T(0.0), T(0.0095);
			theta << T(r * T(0.0) - T(M_PI/2.0) * j)
				, T(0.0), T(0.0) + T(42.0 * M_PI/180.0);
				                                                                                                                                            
			//compute the (4*4) DH Transformation matrix
			retVal.setDhHandParameters(a, alpha, d, theta);
		}
		return retVal;
	}


protected:
	size_t _nFingers, _nJoints;
	DhHand<T> *_joint;


};
#endif
