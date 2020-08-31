#ifndef MLE_H
#define MLE_H

#include "BaysFitFunctions.h"
#include "UBCUtil.h"

#include <Eigen/Dense>

#include "ParseMathematica.h"
#include "ceres/ceres.h"
#include "ceres/dynamic_autodiff_cost_function.h"


using namespace Eigen;
using ceres::AutoDiffCostFunction;
using ceres::DynamicAutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;
using ceres::CauchyLoss;

//static const int _m = 1;

// The stride length of the dynamic_autodiff_cost_function evaluator.
static const int kStride = 10;

//template <typename T>
struct MLE {

	//typedef AutoDiffCostFunction<MLE, 1, static const int(_m> MLECostFunction;

	// The DynamicAutoDiffCostFunction is meant to be used in cases
	// where the number of parameter blocks or the sizes are not
	// known at compile time.
	typedef DynamicAutoDiffCostFunction<MLE, kStride> MLECostFunction;

	MLE(Eigen::MatrixXd x, double y, double stddev) :
		_x(x), _y(y), _stddev(stddev), _polyOrder(x.rows()-1) {}

	template <typename T>
	bool operator()(T const* const* A, T* residuals) const {
		
		T Ax_i(0);
		for (int i=0; i <= _polyOrder; i++) Ax_i += A[i][0]*T(_x(i,0));

		residuals[0] = (_y - Ax_i)/T(_stddev);
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static MLECostFunction* Create (Eigen::MatrixXd x
																	, double y
																	, double stddev
																	, vector<double>* Amat
																	, vector<double*>* parameterBlocks) {

		MLE* mle = new MLE(x, y, stddev);
		MLECostFunction* cf = new MLECostFunction(mle);
		
		// Add all the parameter blocks that affect this constraint.
    parameterBlocks->clear();
		int polyOrder = x.rows()-1;
    for (int i = 0; i <= polyOrder; ++i) {
    	parameterBlocks->push_back(&((*Amat)[i]));
    	cf->AddParameterBlock(1);
    }

    cf->SetNumResiduals(1);
		return (cf);
	}

	Eigen::MatrixXd _x;
	double _y, _stddev;
	int _polyOrder;
};

#endif 
