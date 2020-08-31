#ifndef FIND_MIN_H
#define FIND_MIN_H

#include "UBCUtil.h"

#include "ceres/ceres.h"

#include <Eigen/Core>
#include <string>

using namespace Eigen;
using namespace std;


//template <typename T>
struct FindMin {
	FindMin(MatrixXd &X, MatrixXd &y) : _X(X), _y(y) {}

	template <typename U>
	//bool operator()(const U* const candidatePt, U* residuals) const {
	bool operator()(const U* const candidateW, U* residuals) const {
		
		U Xw = _X * candidateW;
		U funObj = (1/2.0)*sum((Xw - _y)^2);

		residuals[0] = funObj;

		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	//static ceres::CostFunction* Create(
	//	const double cX, const double cY, const double alpha
	//	, const double a, const double b
	//	, const double x, const double y) {
	//		return (
	//			new ceres::AutoDiffCostFunction<EllipseThetaMin, 1, 1>(
		//		new EllipseThetaMin(cX, cY, alpha, a, b, x, y)));
	//}

	MatrixXd _X, _y;
};

#endif 
