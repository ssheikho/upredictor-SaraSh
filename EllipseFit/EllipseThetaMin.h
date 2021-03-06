#ifndef ELLIPSE_THETA_MIN_H
#define ELLIPSE_THETA_MIN_H

#include "Ellipse.h"
#include "UBCUtil.h"

#include "ceres/ceres.h"

//template <typename T>
struct EllipseThetaMin {
	EllipseThetaMin(
		double cX, double cY, double alpha, double a, double b
		, double x, double y) :
		_cX(cX), _cY(cY), _alpha(alpha), _a(a), _b(b), _x(x), _y(y) {}

	template <typename U>
	//bool operator()(const U* const candidatePt, U* residuals) const {
	bool operator()(const U* const candidateTheta, U* residuals) const {
		//_e.getA();//.getTheta();

		U x =  U(_cX) +  U(_a) * cos(U(_alpha)) * cos(candidateTheta[0]) -
			 U(_b) * sin(U(_alpha)) * sin(candidateTheta[0]);
		U y =  U(_cY) +  U(_a) * cos(candidateTheta[0]) * sin( U(_alpha)) +
			 U(_b) * cos( U(_alpha)) * sin(candidateTheta[0]);

	/* The error is the difference between 
		the predicted and observed position. */

    x = x - U(_x);
    y = y - U(_y);
   
		residuals[0] = sqrt(x * x + y * y);
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(
		const double cX, const double cY, const double alpha
		, const double a, const double b
		, const double x, const double y) {
			return (
				new ceres::AutoDiffCostFunction<EllipseThetaMin, 1, 1>(
				new EllipseThetaMin(cX, cY, alpha, a, b, x, y)));
	}

	double _cX, _cY, _alpha, _a, _b, _x, _y;
};

#endif 
