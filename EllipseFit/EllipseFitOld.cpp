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
	srand(time(NULL));

/*************************************/
/* H1.A. Ellipse Fitting Algorith */

/************************************/
/*** Itterate through all N reach 
			samples per participant ***/

	/**********************************************/
	/*** Reading input CSV Motion Files ***/
	int colNo = 0;
	std::vector<std::string> wristMarkers =
 		{"x", "y", "z"};

	std::vector<int> nCounterPerN;

	int TotreachCount = 0;
	int nRowsJoint = 0;

	std::vector<std::string> partIdsV;
	std::vector<char> reachCatsV;
	std::vector<std::string> metaFileNames;
	

	for (int i=1; i<argc; i++){

		// select a motion segment; i.e. P1 or P2 
		string fileName = argv[i];
		size_t whichStudy = fileName.find("study_4");		
		size_t whichPhase = fileName.find("phase1");

		if (whichPhase!= string::npos && whichStudy != string::npos) {
			reachCatsV.push_back(
				fileName.at ((unsigned)(fileName.length()-5)));
			partIdsV.push_back(getPartID(fileName));

			// Calculating the length of string 
			metaFileNames.push_back (
				fileName.substr(0,fileName.length()-17)+".csv" );
			//Add marker names to csv header
			//setHeaderCSV(metaFileNames.back());

			// nx3
			int nRowsTemp = countRowsCSV
				(argv[i], wristMarkers[colNo]);
			nCounterPerN.push_back(nRowsTemp);
			
			nRowsJoint += nRowsTemp;
			TotreachCount ++;
		}
	}

	/** Copy input marker data into Eigen Matrix **/
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

	// start and end of reach
	MatrixXd inStartPtsAlongRows = 
 		 MatrixXd::Zero(TotreachCount,nMarkers);
	MatrixXd inTargetPtsAlongRows =
		 MatrixXd::Zero(TotreachCount,nMarkers);

	// Best Fit Plane pars ONLY FOR WRIST//
	MatrixXd v1sAlongRows = MatrixXd::Zero(TotreachCount,3)
		, v2sAlongRows = MatrixXd::Zero(TotreachCount,3)
		, v3sAlongRows = MatrixXd::Zero(TotreachCount,3);

	MatrixXd planeAngleToZAlongRows = 
		MatrixXd::Zero(TotreachCount,1)
		, pYawZAlongRows = MatrixXd::Zero(TotreachCount,1)
		, pDipXAlongRows = MatrixXd::Zero(TotreachCount,1)
		, pTiltYAlongRows = MatrixXd::Zero(TotreachCount,1);

	// ELLIPSE FITTING
	MatrixXd e3dCAlongRows = MatrixXd::Zero(TotreachCount,3);
	//2d
	MatrixXd EfAlphaAlongRows = MatrixXd::Zero(TotreachCount,1)
	, EfCxAlongRows = MatrixXd::Zero(TotreachCount,1)
	, EfCyAlongRows = MatrixXd::Zero(TotreachCount,1)
	, EfAaxisAlongRows = MatrixXd::Zero(TotreachCount,1)
	, EfBaxisAlongRows = MatrixXd::Zero(TotreachCount,1)
	, EfAxesRatioAlongRows = MatrixXd::Zero(TotreachCount,1);


	int reachCtr = 0, EfpassCtr = 0;
	for (std::vector<int>::const_iterator i 
		= nCounterPerN.begin(); i != nCounterPerN.end(); ++i){

		int nRows = *i;
		MatrixXd inPtsAlongRows = MatrixXd::Zero(nRows,nMarkers);
		MatrixXd inPtsAlongCols = MatrixXd::Zero(nMarkers,nRows);

		string headerLine = 
			getHeaderRowCSV(metaFileNames[reachCtr]);
		int ctr = 0;

		while ((ctr < nMarkers)) {

			int startMarkerDelim = headerLine.find(wSMarkers[ctr]);
			if (startMarkerDelim!=std::string::npos) 
				fillMatCSV(ctr, metaFileNames[reachCtr]
					, wSMarkers[ctr], wSMarkers[ctr+1]
					, wSMarkers[ctr+2], inPtsAlongRows);

			ctr +=3;
		}
		inPtsAlongCols = inPtsAlongRows.transpose();
	/*** END - Reading input CSV Motion Files ***/
	/************************************************/
		
		/** ELLIPSE FITTING **/
		//Ellipse3D e3Dw(cartToHom(inPtsAlongCols.topRows(3)));
		double wDis = (inPtsAlongRows.block(0,0,1,3)
			- inPtsAlongRows.block(nRows-1,0,1,3)).norm()/2.0;
		//double EfCost = e3Dw.getTotalXYCost();
		//if (EfCost > 0 && (e3Dw.getEllipse().getA() > 0.75*wDis) 
		//	&& (e3Dw.getEllipse().getA() < 1.5*wDis)) {

		/************************************************/
		/**		 A.1 Finding Nominal Plane Fit 		**/ 
		
			// collect start and target positions of the reach //
			inStartPtsAlongRows.row(reachCtr) =
				 inPtsAlongRows.row(0);
			inTargetPtsAlongRows.row(reachCtr) =
				 inPtsAlongRows.bottomRows(1);

			// Best Fit Plane model //
			Plane bfp(cartToHom(inPtsAlongCols.topRows(3)));
			v1sAlongRows.row(reachCtr) = bfp.getPlaneV1N().transpose();
			v2sAlongRows.row(reachCtr) = bfp.getPlaneV2N().transpose();
			v3sAlongRows.row(reachCtr) =
				bfp.getPlaneNormalVectN().transpose();

			// Angular model plane - Yaw/Dip/Roll in PCA //
			planeAngleToZAlongRows(reachCtr,0) = bfp.getPlaneNtoZ();
			pYawZAlongRows(reachCtr,0) = bfp.getPlaneYawZ();
			pDipXAlongRows(reachCtr,0) = bfp.getPlaneDipX();
			pTiltYAlongRows(reachCtr,0) = bfp.getPlaneTiltY();

		/** ELLIPSE FITTING 
			e3dCAlongRows.row(reachCtr) = 
				homToCart(e3Dw.getCenter()).transpose();
			// 2d
			Ellipse ellipse = e3Dw.getEllipse();
			EfAlphaAlongRows(reachCtr,0) = ellipse.getRT2D().getTheta();
			EfCxAlongRows(reachCtr,0) = ellipse.getRT2D().getCX();
			EfCyAlongRows(reachCtr,0) = ellipse.getRT2D().getCY();
			EfAaxisAlongRows(reachCtr,0) = ellipse.getA();
			EfBaxisAlongRows(reachCtr,0) = ellipse.getB();
			EfAxesRatioAlongRows(reachCtr,0) = 
				ellipse.getA()/ellipse.getB();

			//output pts
			MatrixXd thetas = e3Dw.findThetas();
			MatrixXd outPtsWAlongRows = homToCart(
				e3Dw.getPointAtThetasH(fixThetas(thetas))).transpose();
**/
	/** For Plotting in Mathematica 
			// simulated pts on ellipse
			size_t nPts = 200;
			MatrixXd simEwAlongRows = homToCart(
				e3Dw.ellipticalInterpolator(nPts)).transpose();
			
			cout << "(* Part-" + std::to_string(reachCtr) 
					 << " *) " << endl;

			cout << "pNToZ = " << planeAngleToZAlongRows(reachCtr,0) 
					 << ", " << "pYawZ = " << pYawZAlongRows(reachCtr,0)
 					 << ", " << "pTiltY = " << pTiltYAlongRows(reachCtr,0) 
					 << endl;
					
			MatrixXd plDirVecs = MatrixXd::Zero(2,3);
			plDirVecs.row(0) = v1sAlongRows.row(reachCtr);
			plDirVecs.row(1) = v2sAlongRows.row(reachCtr);
			printEigenMathematica( bfp.getPlaneOrigin().transpose()
				, cout, "PlaneO"+ std::to_string(reachCtr));		
			printEigenMathematica( plDirVecs , cout
				, "PlaneDir"+ std::to_string(reachCtr));		
			printEigenMathematica( bfp.getPlaneNormalVectN()
				, cout , "PlaneN"+ std::to_string(reachCtr));		
**/
/*		
			

			cout << "Ellipse" + std::to_string(reachCtr) << ": " 
					 << "Alpha= " << EfAlphaAlongRows(reachCtr,0) << ", "
					 << "Cx= " << EfCxAlongRows(reachCtr,0) << ", "
					 << "Cy= " << EfCyAlongRows(reachCtr,0) << ", "
					 << "A= " << EfAaxisAlongRows(reachCtr,0) << ", "
					 << "B= " << EfBaxisAlongRows(reachCtr,0) << endl;

			cout << "Ellipse3d" + std::to_string(reachCtr) 
					 << ": " << endl;
*/
/*
			printEigenMathematica( e3dCAlongRows.row(reachCtr)
				, cout, "center"+ std::to_string(reachCtr));		
			printEigenMathematica( e3Dw.getA()
				, cout, "A"+ std::to_string(reachCtr));	
			printEigenMathematica( e3Dw.getB()
				, cout, "B"+ std::to_string(reachCtr));	


			printEigenMathematica( inStartPtsAlongRows.block
				(reachCtr,0,1,3), cout
				, "inStartPtW"+ std::to_string(reachCtr));	
			printEigenMathematica( inTargetPtsAlongRows.block
				(reachCtr,0,1,3), cout
				, "inTargetPtsW"+ std::to_string(reachCtr));

			printEigenMathematica(inPtsAlongCols.topRows(3).transpose()
				, cout, "inPtsW"+ std::to_string(reachCtr));		
			printEigenMathematica(simEwAlongRows
				, cout, "simEW"+ std::to_string(reachCtr));		
			printEigenMathematica(outPtsWAlongRows
				, cout, "outPtsW"+ std::to_string(reachCtr));	
*/
	/** END - For Plotting in Mathematica **/
			reachCtr ++;
	//	}	
	}


	/*** Bayesian linear regression Model  ***/ 				

	cout << "reachCount: " << reachCtr << endl;
	cout << "TotreachCount: " << TotreachCount << endl;

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
		inXmatplAngModel.block(nMarkers +1 + i,0,3,reachCtr) = 
			reachDist.block(0,i,reachCtr,3).rowwise().normalized().transpose();

	//  shortest travelled distance from start to goal //
	// nMarkers/3 constraints
	for (int i = 0, row = 0; i < nMarkers; i += 3, row++) 
		inXmatplAngModel.row(2 * nMarkers +1 + row) = 
			reachDist.block(0,i,reachCtr,3).rowwise().norm().transpose();

	// target position (x,y,z) //
	// nMarkers constraints
	inXmatplAngModel.bottomRows(nMarkers) = 
		inTargetPtsAlongRows.topRows(reachCtr).transpose();
	// Y input model // 
	//d.n
	MatrixXd inYmatplAngModel = MatrixXd::Zero(3,reachCtr);
	inYmatplAngModel = v3sAlongRows.topRows(reachCtr).transpose();

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

	return 0;
}

