#ifndef DH_HAND_H
#define DH_HAND_H

#include <Eigen/Core>

#include <iostream>

template<typename T>
class DhHand {
public:
	DhHand(T a = T(0.0), T alpha = T(0.0), T d = T(0.0), T theta = T(0.0)) :
		_a(a)
		, _alpha(alpha)
		, _d(d)
		, _theta(theta) {
		computeDH();
	}

	~DhHand() {}

	T getA() { return _a; }
	T getAlpha() { return _alpha; }
	T getD() { return _d; }
	T getTheta() { return _theta; }

	void setValues(T a, T alpha, T d, T theta) {
		_a = a;
		_alpha = alpha;
		_d = d;
		_theta = theta;
		computeDH();
	}

	void setTheta(T theta) {
		_theta = theta;
		computeDH();
	}

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> getMat() {
		return _mat;
	}

	void prettyPrint() {
		std::cout << "DH:	" << _a << "	" << _alpha << "	" << _d <<
			"	" <<
			_theta << std::endl;
		std::cout << _mat << std::endl;
	}

	static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
		dhHandTransform(T a, T alpha, T d, T theta) {
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> retVals(4, 4);

		retVals(0,0) = cos(theta);
		retVals(0,1) = -sin(theta);
		retVals(0,2) = T(0.0);
		retVals(0,3) = a;


		retVals(1,0) = sin(theta) * cos(alpha);
		retVals(1,1) = cos(theta) * cos(alpha);
		retVals(1,2) = -sin(alpha);
		retVals(1,3) = -sin(alpha) * d;

		retVals(2,0) = sin(theta) * sin(alpha);
		retVals(2,1) = cos(theta) * sin(alpha);
		retVals(2,2) = cos(alpha);
		retVals(2,3) = cos(alpha) * d;

		retVals(3,0) = T(0.0);
		retVals(3,1) = T(0.0);
		retVals(3,2) = T(0.0);
		retVals(3,3) = T(1.0);

		return retVals;
	}

protected:
	void computeDhHand() {
		_mat = dhTransform(_a, _alpha, _d, _theta);
	}

	T _a, _alpha, _d, _theta;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _mat;
};

#endif
