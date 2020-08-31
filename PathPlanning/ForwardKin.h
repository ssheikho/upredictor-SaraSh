#ifndef FORWARD_KIN_H
#define FORWARD_KIN_H

#include "DH.h"

#include <cmath>

#include <Eigen/SVD>

#include <cmath>

using namespace Eigen;

template<typename T>
class ForwardKin {

public:
	ForwardKin(size_t nJoints) : _nJoints(nJoints)
		, _joint(new DH<T>[_nJoints]) {
		//_joint[0]
		//setDHParameters();
	}

	~ForwardKin() {
		delete[] _joint;
	}

	void setThetaVect(MatrixXd thetaVect) {
		for(size_t i = 0; i < _nJoints; i++)
			_joint[i].setTheta(thetaVect(i,0));
	}

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		getAvect(size_t toDOF) {
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> a(7,1);
		a << T(0.0), T(0.0), T(0.045), T(-0.045)
			, T(0.0), T(0.0), T(0.0);
		return a.topRows(toDOF);
	}

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		getAlphaVect(size_t toDOF) {
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
			alpha(7,1);
		alpha << T(-M_PI/2.0), T(M_PI/2.0), T(-M_PI/2.0)
			, T(M_PI/2.0), T(-M_PI/2.0), T(M_PI/2.0), T(0.0);
		return alpha.topRows(toDOF);
	}

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
		getDvect(size_t toDOF) {
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
		 d(7,1);
		d << T(0.0), T(0.0), T(0.55), T(0.0)
			, T(0.3), T(0.0), T(0.06);
		return d.topRows(toDOF);
	}

	void setDHParameters(
			Matrix<T, Eigen::Dynamic, Eigen::Dynamic> a
			, Matrix<T, Eigen::Dynamic, Eigen::Dynamic> alpha
			, Matrix<T, Eigen::Dynamic, Eigen::Dynamic> d) {
		for(size_t i = 0; i < _nJoints; i++)
			_joint[i].setValues(a(i,0), alpha(i,0)
				, d(i,0), T(0.0));
	}

	Eigen::Matrix<T, 4, 4> getMat(size_t toDOF) {
		Eigen::Matrix<T, 4, 4> retVal =
			Eigen::Matrix<T, 4, 4>::Identity(4,4);

		for(size_t i = 0; (i < _nJoints) && (i < toDOF); i++)
			retVal = retVal * _joint[i].getMat();
		return retVal;
	}

	DH<T> &getDH(size_t index) {
		return _joint[index];
	}

	size_t getNJoints() {
		return _nJoints;
	}

	static ForwardKin generateBarrettWAM() {
		ForwardKin retVal(7);

		Eigen::Matrix<T, 7, 1> 
			a(7,1), alpha(7,1), d(7,1);
	
		a << T(0.0), T(0.0), T(0.045), T(-0.045)
			, T(0.0), T(0.0), T(0.0);
		alpha << T(-M_PI/2.0), T(M_PI/2.0), T(-M_PI/2.0)
			, T(M_PI/2.0), T(-M_PI/2.0), T(M_PI/2.0), T(0.0);
		d << T(0.0), T(0.0), T(0.55), T(0.0), T(0.3)
			, T(0.0), T(0.06);
		
		//compute the (4*4) DH Transformation matrix
		retVal.setDHParameters(a, alpha, d);

		return retVal;
	}


protected:
	size_t _nJoints;
	DH<T> *_joint;


};
#endif
