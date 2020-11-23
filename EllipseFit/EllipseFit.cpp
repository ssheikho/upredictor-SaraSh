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

#include "DifferentialIKErrorTerm.h"
#include "DikProblem.h"
#include "DikSolver.h"

#include <math.h> 
#include <Eigen/Core>
#include <Eigen/Dense>
//#include <Eigen/SVD>

//#include <cmath>
#include <iostream>
class DikProblem;

DikProblem *dIkProbInChest;

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
	std::string metaFileName;
//	std::string reachCtr = getPartID(fileName) + "P"
//		+fileName.at ((unsigned)(fileName.length()-5));

	int phase = 1;
//whichPhase(fileName);		
	// Calculating the length of string 
	metaFileName = 
		fileName.substr(0,fileName.find("_clean"))+".csv";	
	string headerLine = 
		getHeaderRowCSV(metaFileName);


	int nPts = countRowsCSV(fileName, wristMarkers[0]);
	MatrixXd inPtsAlongRows = MatrixXd::Zero(nPts,nMarkers);
	int nPtsTotal = countRowsCSV(
		metaFileName, wSMarkers[0]);
	MatrixXd inPtsAlongRowsTotal = 
		MatrixXd::Zero(nPtsTotal,nMarkers);

	int ctr = 0;
	while ((ctr < nMarkers)) {
		int startMarkerDelim = 
			headerLine.find(wSMarkers[ctr]);
		if (startMarkerDelim!=std::string::npos) { 
			fillMatCSV(ctr, metaFileName
				, wSMarkers[ctr], wSMarkers[ctr+1]
				, wSMarkers[ctr+2], inPtsAlongRowsTotal);
		}
		ctr +=3;
	}


	if (phase ==1)	
		inPtsAlongRows = inPtsAlongRowsTotal.topRows(nPts);
	
	else	{
		inPtsAlongRows = inPtsAlongRowsTotal.bottomRows(nPts);
		cout <<"here"<<endl;
	}
	// positions in meter
	MatrixXd inPts = 0.01 * inPtsAlongRows.transpose();
	MatrixXd inPtsRCh = inPts.block(Markers::RCh,0,3,nPts)
		,	inPtsMCh = inPts.block(Markers::MCh,0,3,nPts)
		,	inPtsLCh = inPts.block(Markers::LCh,0,3,nPts)
		,	inPtsRPi = inPts.block(Markers::RPi,0,3,nPts)
		,	inPtsRTh = inPts.block(Markers::RTh,0,3,nPts)
		,	inPtsRWr = inPts.block(Markers::RWr,0,3,nPts)
		,	inPtsRLA = inPts.block(Markers::RLA,0,3,nPts)
		,	inPtsREl = inPts.block(Markers::REl,0,3,nPts)
		,	inPtsRUA = inPts.block(Markers::RUA,0,3,nPts)
		,	inPtsRSh = inPts.block(Markers::RSh,0,3,nPts);

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
	
	MatrixXd efPts = MatrixXd::Zero(nMarkers,nPts);
	efPts = inPts;
	efPts.block(Markers::RUA,0,3,nPts) = outPtsRUA;
	efPts.block(Markers::REl,0,3,nPts) = outPtsREl;
	efPts.block(Markers::RLA,0,3,nPts) = outPtsRLA;
	efPts.block(Markers::RWr,0,3,nPts) = outPtsRWr;
	efPts.block(Markers::RTh,0,3,nPts) = outPtsRTh;
	efPts.block(Markers::RPi,0,3,nPts) = outPtsRPi;

	/** For print in Mathematica **/
	/// inPts ///
	//(a) for human
/*
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
*/

		/************************************************/
		/**		 A.1 Finding Nominal Plane Fit 		**/ 

	int TotreachCount = 1;
	int reachCtr = 1;		
	// Best Fit Plane pars ONLY FOR WRIST//

	// start and end of reach
	MatrixXd inStartPtsAlongRows = 
 		 MatrixXd::Zero(1,nMarkers);
	MatrixXd inTargetPtsAlongRows =
		 MatrixXd::Zero(1,nMarkers);


	MatrixXd v1sAlongRows = MatrixXd::Zero(1,3)
		, v2sAlongRows = MatrixXd::Zero(1,3)
		, v3sAlongRows = MatrixXd::Zero(1,3);

	
			// collect start and target positions of the reach //
			inStartPtsAlongRows =
				 inPtsAlongRows.row(0);
			inTargetPtsAlongRows =
				 inPtsAlongRows.bottomRows(1);

			// Best Fit Plane model //
			Plane bfp(cartToHom(inPts.topRows(3)));
			v1sAlongRows = 
				bfp.getPlaneV1N().transpose();
			v2sAlongRows = 
				bfp.getPlaneV2N().transpose();
			v3sAlongRows =
				bfp.getPlaneNormalVectN().transpose();
			// Angular model plane - Yaw/Dip/Roll in PCA //
			MatrixXd planeAngleToZAlongRows
				= MatrixXd::Zero(TotreachCount,1)
				, pYawZAlongRows = MatrixXd::Zero(TotreachCount,1)
				, pDipXAlongRows = MatrixXd::Zero(TotreachCount,1)
				, pTiltYAlongRows = MatrixXd::Zero(TotreachCount,1);

			planeAngleToZAlongRows(reachCtr-1,0) = bfp.getPlaneNtoZ()
			,	pYawZAlongRows(reachCtr-1,0) = bfp.getPlaneYawZ()
			,	pDipXAlongRows(reachCtr-1,0) = bfp.getPlaneDipX()
			,	pTiltYAlongRows(reachCtr-1,0) = bfp.getPlaneTiltY();

cout <<"here"<<endl;
	/*** Bayesian linear regression Model  ***/ 				

	///**** 		nominal Plane Parameters 		***///

	/**		(A) TEST planeNtoZ /Dip w.r.t Y-axis		
			  const term + start position (x,y,z) + PC1 normalized
			  + travelled distance + target position (x,y,z) 		**/

	/* input matrix X for the linear model XW=y */

	int inXdim =  1 + 3 * nMarkers + nMarkers / 3;
	//m.n
	MatrixXd inXmatplAngModel = MatrixXd::Zero(inXdim,reachCtr);

	// const term
	// 1 constraint
	inXmatplAngModel.row(0) = MatrixXd::Ones(1,reachCtr);

	// start position (x,y,z) //
	// nMarkers constraints
	inXmatplAngModel.block(1,0,nMarkers,reachCtr) = 
		inStartPtsAlongRows.topRows(reachCtr).transpose();

	// PC1 normalized//
	// nMarkers constraints
	MatrixXd reachDist = MatrixXd::Zero(reachCtr,nMarkers);
	reachDist = (inTargetPtsAlongRows.topRows(reachCtr) 
		- inStartPtsAlongRows.topRows(reachCtr));

	for (int i = 0; i < nMarkers; i += 3) 
		inXmatplAngModel.block(nMarkers +1 + i,0,3,reachCtr) 
			= reachDist.block(0,i,reachCtr,3).rowwise()
				.normalized().transpose();

	//  shortest travelled distance from start to goal //
	// nMarkers/3 constraints
	for (int i = 0, row = 0; i < nMarkers; i += 3, row++) 
		inXmatplAngModel.row(2 * nMarkers +1 + row) = 
			reachDist.block(0,i,reachCtr,3)
				.rowwise().norm().transpose();

	// target position (x,y,z) //
	// nMarkers constraints
	inXmatplAngModel.bottomRows(nMarkers) = 
		inTargetPtsAlongRows.transpose();
	// Y input model // 
	//d.n
	MatrixXd inYmatplAngModel = MatrixXd::Zero(3,1);
	inYmatplAngModel = v3sAlongRows.transpose();
cout <<"here"<<endl;
	/** Build Bayesian linear regression model of PlaneNtoZ	**/
	//Likelihood
	vector<double> ALiPlNtoZ = findALi(inXmatplAngModel
			, inYmatplAngModel);
	cout << "ALiPlNtoZ: " << endl;
	for (int i = 0; i < ALiPlNtoZ.size(); i++) 
		cout << ALiPlNtoZ[i] << endl;
	/** Evidence **/
	vector<double> AEvPlNtoZ = findAEv(inXmatplAngModel
			, inYmatplAngModel);
	cout << "AEvPlNtoZ: " << endl;
	for (int i = 0; i < AEvPlNtoZ.size(); i++) 
		cout << AEvPlNtoZ[i] << endl;

	MatrixXd planeV3hatEv = 
		predictY(inXmatplAngModel, AEvPlNtoZ);

	Vector3d zAxis(3,1);
	zAxis << 0.0, 0.0, 1.0;
	MatrixXd planeAngleToZhatEv = MatrixXd::Zero(reachCtr,1);
	for (int i = 0; i < reachCtr; i++) {
		Vector3d planeV3 = planeV3hatEv.col(i);
		double dotProduct = planeV3.dot(zAxis);
		Vector3d crossProduct = planeV3.cross(zAxis);
		planeAngleToZhatEv(i,0) = 
			atan2(crossProduct.norm(), dotProduct);
	}

	MatrixXd trainEplaneNToZEv = 
		( planeAngleToZAlongRows.topRows(reachCtr) 
		- planeAngleToZhatEv).cwiseAbs();
	cout << "trainEplaneNToZEv mean = " 
			 << trainEplaneNToZEv.mean() << endl;




	/** output vector y(n.1) predictions **/
	MatrixXd planeV3hatLi = 
		predictY(inXmatplAngModel, ALiPlNtoZ);

	/** training error **/
	MatrixXd planeAngleToZhatLi = MatrixXd::Zero(reachCtr,1);
	for (int i = 0; i < reachCtr; i++) {
		Vector3d planeV3 = planeV3hatLi.col(i);
		double dotProduct = planeV3.dot(zAxis);
		Vector3d crossProduct = planeV3.cross(zAxis);
		planeAngleToZhatLi(i,0) 
			= atan2(crossProduct.norm(), dotProduct);
	}

	MatrixXd trainEplaneNToZ = 
		(planeAngleToZAlongRows.topRows(reachCtr) 
		- planeAngleToZhatLi).cwiseAbs();
	cout << "trainEplaneNToZ mean = " 
			 << trainEplaneNToZ.mean() << endl;

	/* MATHEMATICA OUTPUT: 	*/
	// Plane angles
	//printEigenMathematica( EfAlphaAlongRows.topRows(reachCtr)
	//	, cout, "EfAlphaAlongRows");	
	printEigenMathematica( pDipXAlongRows.topRows(reachCtr)
		, cout, "pDipXAlongRows");	
	printEigenMathematica( pTiltYAlongRows.topRows(reachCtr)
		, cout, "pTiltYAlongRows");	
	printEigenMathematica( pYawZAlongRows
		, cout, "pYawZAlongRows");	
	printEigenMathematica( planeAngleToZAlongRows.topRows(reachCtr)
		,	cout,	"planeAngleToZAlongRows");
	//LR 	model 
	printEigenMathematica( planeAngleToZhatLi.transpose()
		, cout, "planeAngleToZhatLi");	
	printEigenMathematica( trainEplaneNToZ
		, cout, "trainEplaneNToZ");	
/* LR 	Ev model */
	printEigenMathematica( planeAngleToZhatEv.transpose()
		, cout, "planeAngleToZhatEv");	
	printEigenMathematica( trainEplaneNToZEv
		, cout, "trainEplaneNToZEv");	



	/*** END-Bayesian linear regression Model  ***/ 				
	printEigenMathematica( inPts.transpose(), cout
		, "inPts");	
	printEigenMathematica( efPts.transpose(), cout
		, "efPts");	

	return 0;
}

