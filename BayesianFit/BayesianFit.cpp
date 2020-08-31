
#include "Ellipse.h"
#include "FitFunctions.h"
#include "LinearAlgebra.h"
#include "EllipseImplicitFit.h"
#include "Plane.h"

#include "Ellipse3D.h"
#include "BaysFitFunctions.h"
#include "ParseMathematica.h"
#include "ParseCSV.h"
#include "UBCUtil.h"


#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include "RigidTrans2D.h"


//#include <cmath>
#include <iostream>

using namespace Eigen;
using namespace std;

int main(int argc, char **argv) {
	srand(time(NULL));
	//MatrixXd inPts2dAlongColsHom = parseMathematica(getWholeFile(argv[1]));
	
	int colNo = 0;
	int nRows = countRowsFColCSV(argv[1]);
	// 3xn
	MatrixXd inWPtsAlongCols(3,nRows);
	fillMatWholeCSV(colNo, argv[1], inWPtsAlongCols);

	//Ellipse3D::Ellipse3D e3d(cartToHom(inWPtsAlongCols));
	MatrixXd inPts2dAlongColsHom = 
		getPtsXYAlongCols(cartToHom(inWPtsAlongCols));


	int n = inPts2dAlongColsHom.cols();
	MatrixXd xVec = (inPts2dAlongColsHom.block(0,0,1,n)).transpose();
	MatrixXd yVec = (inPts2dAlongColsHom.block(1,0,1,n)).transpose();	
	//cout << "xVec:" << xVec << endl;
	//cout << "---------end---------" <<endl;

	int maxPolyOrder = 5;
	MatrixXd evidenceMat(maxPolyOrder+1,1);		
	MatrixXd AMetaLi = MatrixXd::Zero(maxPolyOrder+1,maxPolyOrder+1);	
	MatrixXd AMetaEv = MatrixXd::Zero(maxPolyOrder+1,maxPolyOrder+1);
	MatrixXd likelihoodMat(maxPolyOrder+1,1);
	Eigen::MatrixXd yVecHatMetaLi	= 
		MatrixXd::Zero(maxPolyOrder+1,n);
	Eigen::MatrixXd yVecHatMetaEv	= 
		MatrixXd::Zero(maxPolyOrder+1,n);

	for(size_t k = 0; k < maxPolyOrder; k++) {
		Eigen::MatrixXd Xmat = MatrixXd::Zero(k+1,n);
		Xmat = buildXMatIn(xVec, k);
	  //cout << "Xmat:\n" << Xmat << endl;

		//A -> fit row per model order m.1 
		//where "m = maxPolyOrder+1"
		vector<double> AvecLi = findALi(Xmat, yVec.transpose());	

//		vector<double> AvecEv = findAEv(Xmat, yVec.transpose());
//		cout << "AvecEv:" << AvecEv.size() << endl;	


//		AMetaLi.block(k,0,1,AvecLi.size()) = AvecLi;
//		AMetaEv.block(k,0,1,AvecEv.size()) = AvecEv;
		for(size_t j = 0; j < AvecLi.size(); j++) {
			AMetaLi(k,j) = AvecLi[j];
//			AMetaEv(k,j) = AvecEv[j];
		}

		//yVecHatLi 1.n
		Eigen::MatrixXd yVecHatLi = 
			predictY(Xmat, AMetaLi.block(k,0,1,k+1));
		yVecHatMetaLi.block(k,0,1,n) = yVecHatLi;
//		printEigenMathematica(yVecHatLi, cout, "yVecHatLi");	

		//Eigen::MatrixXd yVecHatEv = 
		//	predictY(Xmat, AMetaEv.block(k,0,1,k+1));
		//yVecHatMetaEv.block(k,0,1,n) = yVecHatEv;
		//printEigenMathematica(yVecHatEv, cout, "yVecHatEv");	

		evidenceMat(k,0) = 
			computePy_xvalpha(Xmat, yVec.transpose());
		likelihoodMat(k,0) = 
			computePy_xav(Xmat, yVec.transpose());
	}


	cout << /*"evidence for " << argv[1] <<":\n" <<*/ evidenceMat(0,0) /*<< "\n"*/ << endl;	
	//cout << /*"likelihood for " << argv[1] <<":\n" <<*/ likelihoodMat << "\n" << endl;
	//cout << "AMetaLi: \n" << AMetaLi << endl;
	
	//printEigenMathematica(AMetaLi, cout, "AMetaLi");
	//printEigenMathematica(AMetaEv, cout, "AMetaEv");
	//printEigenMathematica(evidenceMat, cout, "evidenceMat");
	//printEigenMathematica(likelihoodMat, cout, "likelihoodMat");

	
	return 0;
}









