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

		if (whichPhase!= string::npos && whichStudy != string::npos) 			{
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
	{	"RThumbX",	"RThumbY",	"RThumbZ" 
	, "RLArmX",	"RLArmY",	"RLArmZ"		
	, "RElbX",	"RElbY",	"RElbZ"				
	, "RUArmX",	"RUArmY",	"RUArmZ"		
	, "RShoX",	"RShoY",	"RShoZ"
	, "RChestX",	"RChestY",	"RChestZ"	
	/*, "MChestX",	"MChestY",	"MChestZ"
	, "LChestX", "LChestY", "LChestZ"*/ };
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
		Ellipse3D e3Dw(cartToHom(inPtsAlongCols.topRows(3)));
		MatrixXd thetas = fixThetas(e3Dw.findThetas());
		double wDis = (inPtsAlongRows.block(0,0,1,3)
			- inPtsAlongRows.block(nRows-1,0,1,3)).norm()/2.0;
		double EfCost = e3Dw.getTotalXYCost();
		if (EfCost > 0 && (e3Dw.getEllipse().getA() > 0.75*wDis) 
		&& (e3Dw.getEllipse().getA() < 1.5*wDis)) {

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

		/** ELLIPSE FITTING **/
			Ellipse ellipse = e3Dw.getEllipse();
			EfAlphaAlongRows(reachCtr,0) = ellipse.getRT2D().getTheta();
			EfCxAlongRows(reachCtr,0) = ellipse.getRT2D().getCX();
			EfCyAlongRows(reachCtr,0) = ellipse.getRT2D().getCY();
			EfAaxisAlongRows(reachCtr,0) = ellipse.getA();
			EfBaxisAlongRows(reachCtr,0) = ellipse.getB();
			EfAxesRatioAlongRows(reachCtr,0) = 
				ellipse.getA()/ellipse.getB();

	/** For Plotting in Mathematica **/
			MatrixXd e3dCenter = e3Dw.getCenter();
			MatrixXd e3dA = e3Dw.getA();
			MatrixXd e3dB = e3Dw.getB();
			// simulated pts on ellipse
			size_t nPts = 200;
			MatrixXd simEwAlongRows = homToCart(
				e3Dw.ellipticalInterpolator(nPts)).transpose();
			MatrixXd outPtsWAlongRows = homToCart(
				e3Dw.getPointAtThetasH(thetas)).transpose();
			

			cout << "Plane" + std::to_string(reachCtr) << "= " 
					 << "pNToZ= " << planeAngleToZAlongRows(reachCtr,0) << ", "
					 << "pYawZ= " << pYawZAlongRows(reachCtr,0) << endl;

			cout << "Ellipse" + std::to_string(reachCtr) << ": " 
					 << "Alpha= " << EfAlphaAlongRows(reachCtr,0) << ", "
					 << "Cx= " << EfCxAlongRows(reachCtr,0) << ", "
					 << "Cy= " << EfCyAlongRows(reachCtr,0) << ", "
					 << "A= " << EfAaxisAlongRows(reachCtr,0) << ", "
					 << "B= " << EfBaxisAlongRows(reachCtr,0) << endl;

			cout << "Ellipse3d" + std::to_string(reachCtr) << ": " << endl;
			printEigenMathematica(  homToCart(e3dCenter).transpose()
				, cout, "center"+ std::to_string(reachCtr));		
			printEigenMathematica( e3dA
				, cout, "A"+ std::to_string(reachCtr));	
			printEigenMathematica( e3dB
				, cout, "B"+ std::to_string(reachCtr));	



			printEigenMathematica( inStartPtsAlongRows.block
				(reachCtr,0,1,3), cout
				, "inStartPtW"+ std::to_string(reachCtr));	
			printEigenMathematica( inTargetPtsAlongRows.block
				(reachCtr,0,1,3), cout
				, "inTargetPtsW"+ std::to_string(reachCtr));


	printEigenMathematica(inPtsAlongRows.block(0,3,nRows,3)
		, cout, "inPtsUA"+ std::to_string(reachCtr));	
	printEigenMathematica(inPtsAlongRows.leftCols(3)
		, cout, "inPtsW"+ std::to_string(reachCtr));	
	
	printEigenMathematica(simEwAlongRows
		, cout, "simEW"+ std::to_string(reachCtr));		
	printEigenMathematica(outPtsWAlongRows
		, cout, "outPtsW"+ std::to_string(reachCtr));	
	
	/** END - For Plotting in Mathematica **/
			reachCtr ++;
		}	
	}

	return 0;
}

