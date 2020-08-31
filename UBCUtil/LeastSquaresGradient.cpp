#include "LeastSquaresGradient.h"

#include "FindMin.h"
void LeastSquaresGradient(MatrixXd &X, MatrixXd &y) {

	int n = X.rows();
	int d = X.cols();

	// Initial guess
	MatrixXd w = MatrixXd::Zero(d,1);

	// Function we're going to minimize (and that computes gradient)
	//funObj(w) = leastSquaresObj(w,X,y)

	// This is how you compute the function and gradient:
	//(f,g) = funObj(w)

	// Derivative check that the gradient code is correct:
	//g2 = numGrad(funObj,w)


	// Solve least squares problem
	//w = findMin(funObj,w)

	// Make linear prediction function
	//predict(Xhat) = Xhat*w
}
