#ifndef ELLIPSE_CONIC_CONSTRAINTS_H
#define ELLIPSE_CONIC_CONSTRAINTS_H

#include "FitFunctions.h"
#include "Plane.h"
#include "Ellipse.h"
#include "UBCUtil.h"

#include <Eigen/Core>

#include "ceres/ceres.h"
#include "gflags/gflags.h"
#include "glog/logging.h"

using ceres::AutoDiffCostFunction;
using ceres::CauchyLoss;
using ceres::CostFunction;
using ceres::LossFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;

struct EllipseConicConstraints {
	EllipseConicConstraints(double x, double y) :
		_x(x), _y(y) {}

	template <typename U>
	bool operator()(const U* const conicConstraints, U* residuals) const {
		
		/* Each residual depends on the 6 conic constraints */
    residuals[0] = U (_x * _x) * conicConstraints[0] 
			+  U (_x * _y) * conicConstraints[1] 
			+  U (_y * _y) * conicConstraints[2] 
			+  U (_x) * conicConstraints[3] 
			+  U (_y) * conicConstraints[4] 
			+  conicConstraints[5];

		// 4ac - b^2 = 1.
    residuals[1] = U (4.0) * conicConstraints[0] 
			* conicConstraints[2]	- conicConstraints[1] 
			* conicConstraints[1] - U (1.0);
   
		return true;
	}
	/* creating a CostFunction for every observation. */
	static ceres::CostFunction* Create(double x, double y) {
			return (
				new ceres::AutoDiffCostFunction
					<EllipseConicConstraints, 2, 6> (
				new EllipseConicConstraints(x, y)) 
		);
	}

	double _x, _y;
};

#endif 
