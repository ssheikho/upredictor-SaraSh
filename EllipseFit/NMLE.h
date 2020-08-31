#ifndef NMLE_H
#define NMLE_H

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

//template <typename T>
struct NMLE {

/*
selecting number of parts instd of features
z_i = Wx_i -> W_(k.d) 
y_i = w^T z_i -> w_(k.1)
*/
static const int kStride = 10;
typedef DynamicAutoDiffCostFunction<NMLE, kStride> 	
	NMLECostFunction;

	NMLE(Eigen::MatrixXd x, double y, int k) :
		_x(x), _y(y), _k(k), _d(x.rows()) {}

/*
 	keep constraints as one long vector.
	W(k*d) followed by w(k) with a totak of k*(d+1) 
	parameters to solve for
*/
	template <typename T>

	// Use parameters[i] to access the i'th parameter block.
	bool operator()(T const* const* W
		,T* residuals) const {
/*
	bool operator()(const T* const W, const T* const w
		,T* residuals) const {
*/
		T Z_ic(0);
		T hZ_ic(0);
		T predicted_y(0);		
		for (int k=0; k < _k; k++) {
			for (int i=0; i < _d; i++) 
					Z_ic += W[i + k*_d][0]*T(_x(i,0)); //k.1
		hZ_ic = (T(1.0)/(T(1.0) + exp(-Z_ic)));
		predicted_y += hZ_ic * W[k + _k*_d][0];
		}
		residuals[0] = predicted_y - T(_y);

		return true;
}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static NMLECostFunction* Create (Eigen::MatrixXd x
																	, double y
																	, int k
																	, vector<double>* Wmat
																	, vector<double*>* 
																		parameterBlocks ) 
{

		NMLE* nmle = new NMLE(x, y, k);
		NMLECostFunction* cf = new NMLECostFunction(nmle);
		
		// Add all the parameter blocks that affect this constraint.
    parameterBlocks->clear();

		for (int j = 0; j < (k + k*x.rows()); j++) {
		  	parameterBlocks->push_back(&((*Wmat)[j]));
		  	cf->AddParameterBlock(1);
		  }

    cf->SetNumResiduals(1);
		return (cf);
	}

	Eigen::MatrixXd _x;
	double _y;
	int _k, _d;
};

#endif 
