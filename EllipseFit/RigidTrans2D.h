#ifndef RIGID_TRANS_2D_H
#define RIGID_TRANS_2D_H

#include <Eigen/Dense>

#include <iostream>

using namespace std;

//template<typename T>
class RigidTrans2D {
public:
	RigidTrans2D(double theta, double cX, double cY) :
		_theta(theta), _cX(cX), _cY(cY), _rt(3,3) {
		setRT(_theta, _cX, _cY);
	}

	~RigidTrans2D() {}

	Eigen::Matrix<double, 3, 3> getRT() {
		return _rt;
	}

	void setRT(double theta, double cX, double cY) {
		_theta = theta;
		_cX = cX;
		_cY = cY;

		//cout << _rt << endl;

		_rt(0,0) = cos(_theta);
		_rt(0,1) = -sin(_theta);
		_rt(0,2) = _cX;

		_rt(1,0) = sin(_theta);
		_rt(1,1) = cos(_theta);
		_rt(1,2) = _cY;

		_rt(2,0) = 0.0;
		_rt(2,1) = 0.0;
		_rt(2,2) = 1.0;
	}

	double getTheta() {
		return _theta;
	}

	double getCX() {
		return _cX;
	}

	double getCY() {
		return _cY;
	}

protected:
	double _theta, _cX, _cY;
	Eigen::Matrix<double, 3, 3> _rt;
};

#endif
