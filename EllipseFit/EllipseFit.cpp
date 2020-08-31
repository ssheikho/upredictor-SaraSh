#include "ParseMathematica.h"
#include "ParseCSV.h"
#include "UBCUtil.h"
#include "FitFunctions.h"
#include "Ellipse3D.h"
#include "Plane.h"
#include "LinearAlgebra.cpp"
#include "RigidTrans2D.h"

#include "BaysFitFunctions.h"
#include "EllipseConicConstraints.h"

#include <math.h> 
#include <Eigen/Core>
#include <Eigen/Dense>
//#include <Eigen/SVD>

//#include <cmath>
#include <iostream>

using namespace Eigen;
using namespace std;
//ensure theta is within [-pi,pi], 
//	&& theta(t) is a monotonic function 
MatrixXd fixThetas(MatrixXd a) {
	MatrixXd retVal(a.rows(), a.cols());

	double curTheta = a(0,0);
	while(curTheta >= -M_PI) curTheta -= M_PI;
	while(curTheta <= M_PI) curTheta += M_PI;
	retVal(0,0) = curTheta;

	for(size_t i = 1; i < a.cols(); i++) {
		double curTheta = a(0,i);
		while(curTheta >= -M_PI) curTheta -= 2.0 * M_PI;
		while(curTheta <= M_PI) curTheta += 2.0 * M_PI;

		double testA = curTheta + 2.0 * M_PI;
		double testB = curTheta - 2.0 * M_PI;
		
		double distNull = curTheta - retVal(0, i - 1);
		double distA = testA - retVal(0, i - 1);
		double distB = testB - retVal(0, i - 1);

		retVal(0, i) = distNull < distA ? curTheta : testA;
		retVal(0, i) = distNull < distB ? curTheta : testB;
	}

	return retVal;
}

int main(int argc, char **argv) {

	string fileName = argv[1];
	cout << "(* " << fileName << "*)" << endl;
/*************************************/
/* H1.A. Ellipse Fitting Algorith */
/**********************************************/

	/*** Reading input CSV Motion Files ***/
	std::vector<std::string> wristMarkers =
 		{"x", "y", "z"};
	// Copy input marker data into Eigen Matrix 
	std::vector<std::string> wSMarkers =	
	{	
   "RWristX",	"RWristY",	"RWristZ" 
	,	"RThumbX",	"RThumbY",	"RThumbZ" 
	, "RLArmX",	"RLArmY",	"RLArmZ"		
	, "RElbX",	"RElbY",	"RElbZ"				
	, "RUArmX",	"RUArmY",	"RUArmZ"		
	, "RShoX",	"RShoY",	"RShoZ"
	, "RChestX",	"RChestY",	"RChestZ"	
	, "MChestX",	"MChestY",	"MChestZ"
	, "LChestX", "LChestY", "LChestZ"
	, "LShoX", "LShoY", "LShoZ" 
	, "RPinkieX",	"RPinkieY",	"RPinkieZ" };
	int nMarkers = wSMarkers.size();
	// nx3
	int nPts = countRowsCSV(fileName, wristMarkers[0]);
	MatrixXd inPtsAlongRows = MatrixXd::Zero(nPts,nMarkers);
	std::string metaFileName;
	std::string reachCtr = getPartID(fileName) + "P"
		+fileName.at ((unsigned)(fileName.length()-5));
	// Calculating the length of string 
	metaFileName = 
		fileName.substr(0,fileName.find("_clean"))+".csv";	
	string headerLine = 
		getHeaderRowCSV(metaFileName);

	int ctr = 0;
	while ((ctr < nMarkers)) {
		int startMarkerDelim = 
			headerLine.find(wSMarkers[ctr]);
		if (startMarkerDelim!=std::string::npos) { 
			fillMatCSV(ctr, metaFileName
				, wSMarkers[ctr], wSMarkers[ctr+1]
				, wSMarkers[ctr+2], inPtsAlongRows);
		}
		ctr +=3;
	}

	MatrixXd inPtsRPi = inPtsAlongRows.rightCols(3).transpose()
		, inPtsRTh = inPtsAlongRows.block(0,3,nPts,3).transpose()
		, inPtsRWr = inPtsAlongRows.leftCols(3).transpose()
		, inPtsRLA = inPtsAlongRows.block(0,6,nPts,3).transpose()
		, inPtsREl = inPtsAlongRows.block(0,9,nPts,3).transpose()
		, inPtsRUA = inPtsAlongRows.block(0,12,nPts,3).transpose()
		, inPtsRSh = inPtsAlongRows.block(0,15,nPts,3).transpose()
		,	inPtsRCh = inPtsAlongRows.block(0,18,nPts,3).transpose()
		,	inPtsMCh = inPtsAlongRows.block(0,21,nPts,3).transpose()
		,	inPtsLCh = inPtsAlongRows.block(0,24,nPts,3).transpose()
		,	inPtsLSh = inPtsAlongRows.block(0,27,nPts,3).transpose();

	/*** END - Reading input CSV Motion Files ***/
	/************************************************/
	/**		 A.1 Finding Nominal Plane Fit 		**/ 
	//(a) ShoToUa
	Plane bfpUa(cartToHom(inPtsRUA));
	Vector3d bfpUaV1 = bfpUa.getPlaneV1N()
		, bfpUaV2 = bfpUa.getPlaneV2N()
		, bfpUaV3 = bfpUa.getPlaneNormalVectN();
	/** ELLIPSE FITTING **/
	//(a) ShoToUa
	Ellipse3D e3dUA(cartToHom(inPtsRUA));
	MatrixXd thetasUA = e3dUA.findThetas();
	MatrixXd outPtsRUA = homToCart(
		e3dUA.getPointAtThetasH(fixThetas(thetasUA)));
	//(b) shoToEl
	Ellipse3D e3dEl(cartToHom(inPtsREl));
	MatrixXd thetasEl = e3dEl.findThetas();
	MatrixXd outPtsREl = homToCart(
		e3dEl.getPointAtThetasH(fixThetas(thetasEl)));
	//(c) shoToLa
	Ellipse3D e3dLA(cartToHom(inPtsRLA));
	MatrixXd thetasLA = e3dLA.findThetas();
	MatrixXd outPtsRLA = homToCart(
		e3dLA.getPointAtThetasH(fixThetas(thetasLA)));
	//(d) shoToWr
	Ellipse3D e3dWr(cartToHom(inPtsRWr));
	MatrixXd thetasWr = e3dWr.findThetas();
	MatrixXd outPtsRWr = homToCart(
		e3dWr.getPointAtThetasH(fixThetas(thetasWr)));
	//(e) shoToTh
	Ellipse3D e3dTh(cartToHom(inPtsRTh));
	MatrixXd thetasTh = e3dTh.findThetas();
	MatrixXd outPtsRTh = homToCart(
		e3dTh.getPointAtThetasH(fixThetas(thetasTh)));
	//(e) shoToPi
	Ellipse3D e3dPi(cartToHom(inPtsRPi));
	MatrixXd thetasPi = e3dPi.findThetas();
	MatrixXd outPtsRPi = homToCart(
		e3dPi.getPointAtThetasH(fixThetas(thetasPi)));


	/** For print in Mathematica **/
	/// inPts ///
	//(a) for human
	printEigenMathematica( inPtsRWr.transpose(), cout
		, "inPtsRWr"+ reachCtr);	
	printEigenMathematica( inPtsRTh.transpose(), cout
		, "inPtsRTh"+ reachCtr);	
	printEigenMathematica( inPtsRPi.transpose(), cout
		, "inPtsRPi"+ reachCtr);	
	printEigenMathematica( inPtsRLA.transpose(), cout
		, "inPtsRLA"+ reachCtr);	
	printEigenMathematica( inPtsREl.transpose(), cout
		, "inPtsREl"+ reachCtr);	
	printEigenMathematica( inPtsRUA.transpose(), cout
		, "inPtsRUA"+ reachCtr);	
	printEigenMathematica( inPtsRSh.transpose(), cout
		, "inPtsRSh"+ reachCtr);	
	printEigenMathematica( inPtsRCh.transpose(), cout
		, "inPtsRCh"+ reachCtr);	
	printEigenMathematica( inPtsMCh.transpose(), cout
		, "inPtsMCh"+ reachCtr);	
	printEigenMathematica( inPtsLCh.transpose(), cout
		, "inPtsLCh"+ reachCtr);	
	// (1) outPts
	printEigenMathematica( outPtsRWr.transpose(), cout
		, "outPtsRWr"+ reachCtr);	
	printEigenMathematica( outPtsRTh.transpose(), cout
		, "outPtsRTh"+ reachCtr);	
	printEigenMathematica( outPtsRPi.transpose(), cout
		, "outPtsRPi"+ reachCtr);	
	printEigenMathematica( outPtsRLA.transpose(), cout
		, "outPtsRLA"+ reachCtr);	
	printEigenMathematica( outPtsREl.transpose(), cout
		, "outPtsREl"+ reachCtr);	
	printEigenMathematica( outPtsRUA.transpose(), cout
		, "outPtsRUA"+ reachCtr);	

	return 0;
}

