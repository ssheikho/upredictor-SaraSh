#ifndef BAYS_FIT_FUNCTIONS_H
#define BAYS_FIT_FUNCTIONS_H

#include <Eigen/Dense>

#include "EllipseImplicitFit.h"
#include "MLE.h"
#include "MEE.h"
#include "Ellipse.h"
#include "Ellipse3D.h"

using namespace Eigen;


//m -> x-features (d)

// The scenario is that we are given a data set of exchangeable pairs 
// D = {(y 1 , x 1 ), ..., (y N , x N )}.
// Collect Y = [y 1 · · · y N ] and X = [x 1 · · · x N ].

//-----The evidence for linearity (p(Y|X, V)) 

//useful for selecting among different linear models, namely models with different inputs. 
//The different inputs might be different nonlinear transformations of the measurements. If we consider the different inputs as separate models with separate priors, then we compute (35) and (30) for each model and see which is largest. 
	//(35) α 
	//(30) p(Y|X, V) (The invariant prior reduces this to) -> p(Y|X, V, α)

// ----select polynomial order. 
//The data is synthetic with N = 50 and known variance V = 10. 
//For order k, the input vector is x = [1 x x 2 · · · x k ] T . 

Eigen::MatrixXd changeBasis(double inX, int k);

Eigen::MatrixXd getPtsXYAlongCols(MatrixXd inPtsAlongColsHom);

Eigen::MatrixXd buildXMatIn(Eigen::MatrixXd inXvec, int k);
double computeAlpha(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat);
Eigen::MatrixXd computeSy_x(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat);
double computePy_xvalpha(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat);
vector<double> findAEv(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat);
vector<double> findALi(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat);
double computePy_xav(Eigen::MatrixXd inXmat, Eigen::MatrixXd inYmat);


Eigen::MatrixXd predictY(Eigen::MatrixXd inXmat, Eigen::MatrixXd inA);


//Ellipse3D* getE3D(Eigen::MatrixXd inPtsAlongColsHom);

#endif
