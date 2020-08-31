#ifndef MEE_H
#define MEE_H

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
static const int evStride = 20;

//template <typename T>
struct MEE {

	//typedef AutoDiffCostFunction<MEE, 1, static const int(_m> MEECostFunction;

	// The DynamicAutoDiffCostFunction is meant to be used in cases
	// where the number of parameter blocks or the sizes are not
	// known at compile time.
	typedef DynamicAutoDiffCostFunction<MEE, evStride> MEECostFunction;

	MEE(Eigen::MatrixXd Xmat, Eigen::MatrixXd y, double stddev, double alpha) :
		_x(Xmat), _y(y), _alpha(alpha), _stddev(stddev), _polyOrder(_x.rows()-1) {}

	template <typename T>
	bool operator()(T const* const* A, T* residuals) const {

		Eigen::MatrixXd Sxx = (_alpha+1.0) * (_x*_x.transpose());			
		Eigen::MatrixXd Syx = _y*_x.transpose();
		Eigen::MatrixXd RS = Syx * Sxx.inverse() * Syx.transpose();
		Eigen::MatrixXd Syy = _y*_y.transpose();		
		T Sy_x = T( Syy(0,0) - RS(0,0) );		


		T ASxx[_polyOrder+1];		
		for (int c=0; c <= _x.cols(); c++) {
			ASxx[c] = T(0);
			for (int i=0; i <= _polyOrder; i++) 
				ASxx[i] += A[i][0] * T(Sxx(i,c)); //(d.m) -> row vector
		}

		T ASxxAtrans(0);		
			for (int i=0; i <= _polyOrder; i++) 
				ASxxAtrans += ASxx[i] * A[i][0];

		T SyxAtrans(0);		
			for (int i=0; i <= _polyOrder; i++) 
				SyxAtrans += T(Syx(0,i)) * A[i][0];

		residuals[0] = T(1.0/_stddev) * 
			( ASxxAtrans - T(2.0)*SyxAtrans + T(Syy(0,0)) );

		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static MEECostFunction* Create (Eigen::MatrixXd Xmat
																	, Eigen::MatrixXd Ymat
																	, double stddev
																	, double alpha
																	, vector<double>* Amat
																	, vector<double*>* parameterBlocks) {

		MEE* mee = new MEE(Xmat, Ymat, stddev, alpha);
		MEECostFunction* cf = new MEECostFunction(mee);
		
		// Add all the parameter blocks that affect this constraint.
    parameterBlocks->clear();
		int polyOrder = Xmat.rows()-1;
    for (int i = 0; i <= polyOrder; ++i) {
    	parameterBlocks->push_back(&((*Amat)[i]));
    	cf->AddParameterBlock(1);
    }

    cf->SetNumResiduals(1);
		return (cf);
	}

	Eigen::MatrixXd _x, _y;
	double _alpha, _stddev;
	int _polyOrder;
};

#endif 
