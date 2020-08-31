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
static const int kStride = 50;

//template <typename T>
struct MLE {

	//typedef AutoDiffCostFunction<MLE, 1, static const int(_m> MLECostFunction;

	// The DynamicAutoDiffCostFunction is meant to be used in cases
	// where the number of parameter blocks or the sizes are not
	// known at compile time.
	typedef DynamicAutoDiffCostFunction<MLE, kStride> MLECostFunction;

	MLE(MatrixXd inXi, MatrixXd inYi, MatrixXd stddev) :
		_inXi(inXi)
		, _inYi(inYi)
		, _stddev(stddev)
		, _m(inXi.rows())
		, _d(inYi.rows()) {}

	template <typename T>
	bool operator()(T const* const* A, T* residuals) const {
	
		//d.m x m.1 = d.1
		T Ax_i[_d];
		for (int d = 0; d < _d; d++) {
			Ax_i[d] = T(0.0);
			for (int m = 0; m < _m; m++) 
				Ax_i[d] += A[d*_m+m][0] * T(_inXi(m,0));

			residuals[d] = (T(_inYi(d,0)) - Ax_i[d])/T(_stddev(d,d));
		
		}
/*
		T Alen(0);
		for (int i=0; i <= _polyOrder; i++) {
			//REGULIZER
			Alen += A[i][0];
		}
*/
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static MLECostFunction* Create (Eigen::MatrixXd inX_i
																	, Eigen::MatrixXd inY_i
																	, Eigen::MatrixXd stddev
																	, vector<double>* Amat
																	, vector<double*>* 
																		parameterBlocks) {

		MLE* mle = new MLE(inX_i, inY_i, stddev);
		MLECostFunction* cf = new MLECostFunction(mle);
		
		// Add all the parameter blocks that affect this constraint.
    parameterBlocks->clear();
		int inM = inX_i.rows();
		int inD = inY_i.rows();
		for (int d = 0; d < inD; d++) {
			for (int m = 0; m < inM; m++) {
				//d.m
    		parameterBlocks->push_back(&((*Amat)[d*inM + m]));
    		cf->AddParameterBlock(1);
			}
		}
    cf->SetNumResiduals(inD);
		return (cf);
}


	Eigen::MatrixXd _inXi, _inYi;
	Eigen::MatrixXd _stddev;
	int _m, _d;
};

#endif 
