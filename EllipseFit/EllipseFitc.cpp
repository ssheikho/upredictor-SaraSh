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
/*	
	string fileName = argv[1];
	string  partId = getPartID(fileName);
	char phase = fileName.at ((unsigned)(fileName.length()-5));

	// nx3
	int nRows = countRowsCSV (argv[1], wristMarkers[colNo]);

	nCounterPerN.push_back(nRows);
*/
	std::vector<int> nCounterPerN;

	int reachCount = 0;
	int nRowsJoint = 0;

	std::vector<std::string> partIdsV;
	std::vector<char> reachCatsV;
	string fileName = "";


	MatrixXd inWSPtsAlongRows = MatrixXd::Zero(100,6);		
		fillMatCSVWS(argv[1], 0, 0
			, 100, inWSPtsAlongRows);
cout << inWSPtsAlongRows<< endl;

	for (int i=1; i<argc; i++){

		fileName = argv[i];
		reachCatsV.push_back(
			fileName.at ((unsigned)(fileName.length()-5)));

		// select a motion segment; i.e. P1 or P2 
		if (reachCatsV[i-1] == '1'){
			partIdsV.push_back(getPartID(fileName));

			// nx3
			int nRowsTemp = countRowsCSV
				(argv[i], wristMarkers[colNo]);
			nCounterPerN.push_back(nRowsTemp);

			nRowsJoint += nRowsTemp;
			reachCount ++;
		}
	}

	/*** END - Reading input CSV Motion Files ***/
	/************************************************/

	/************************************************/
	/**		 A.1 Finding Nominal Plane Fit 		**/ 

	int reachCountTot = reachCount;
	MatrixXd inWPtsAlongRowsJoint(nRowsJoint,3);		
	MatrixXd inWPtsPCAalongRowsJoint(nRowsJoint,3);		

	MatrixXd eigValsRowWise = 
		MatrixXd::Zero(nCounterPerN.size(),3);

	MatrixXd v1sAlongRows = 
		MatrixXd::Zero(nCounterPerN.size(),3);
	MatrixXd v2sAlongRows = 
		MatrixXd::Zero(nCounterPerN.size(),3);
	MatrixXd v3sAlongRows = 
		MatrixXd::Zero(nCounterPerN.size(),3);

	
	MatrixXd pYawZAlongRows = 
		MatrixXd::Zero(nCounterPerN.size(),1);
	MatrixXd pDipXAlongRows = 
		MatrixXd::Zero(nCounterPerN.size(),1);
	MatrixXd pTiltYAlongRows = 
		MatrixXd::Zero(nCounterPerN.size(),1);
	
	//FOR TEST
	MatrixXd planeAngleToZAlongRows = 
		MatrixXd::Zero(nCounterPerN.size(),1);

	MatrixXd HypErrShellAxesRowWise = 
		MatrixXd::Zero(nCounterPerN.size(),3);

	MatrixXd fixedLEcenterOffsetsAlongRows = 
		MatrixXd::Zero(nCounterPerN.size(),1);
	MatrixXd fixedLEllipsoidAxesRowWise = 
		MatrixXd::Zero(nCounterPerN.size(),3);

	std::vector<std::pair<double, double>>
		 maxMinAngularErrsV;


	MatrixXd inWstartPtsAlongRows = 
 		 MatrixXd::Zero(nCounterPerN.size(),3);
	MatrixXd inWtargetPtsAlongRows =
		 MatrixXd::Zero(nCounterPerN.size(),3);

	// ELLIPSE FITTING
	MatrixXd EfAlphaAlongRows =
 		 MatrixXd::Zero(nCounterPerN.size(),1);

	int nRowsPrev = 0;
	int reachCounter = 1;

	Vector3d zAxis (0.0, 0.0, 1.0) ;  

	
	/***	 The nominal Plane in PCA 	***/
	for (std::vector<int>::const_iterator i 
		= nCounterPerN.begin(); i != nCounterPerN.end(); ++i){
		int nRows = *i;

		MatrixXd inWPtsAlongRows = MatrixXd::Zero(nRows,3);
		fillMatCSV(colNo, argv[reachCounter]
			, wristMarkers[0], wristMarkers[1]
			, wristMarkers[2], inWPtsAlongRows);
		MatrixXd inWPtsAlongCols = inWPtsAlongRows.transpose();

		inWstartPtsAlongRows.row(reachCounter-1) =
			 inWPtsAlongRows.row(0);
			
		inWtargetPtsAlongRows.row(reachCounter-1) =
			 inWPtsAlongRows.bottomRows(1);
			

		/*** Best Fit Plane model ***/
		Plane bfp(cartToHom(inWPtsAlongCols));
		v1sAlongRows.row(reachCounter-1) = 
			bfp.getPlaneV1N().transpose();
		v2sAlongRows.row(reachCounter-1) = 
			bfp.getPlaneV2N().transpose();
		v3sAlongRows.row(reachCounter-1) =
			bfp.getPlaneNormalVectN().transpose();

		/***Eigen System Solver***/
	 	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> 
			eigenSys = bfp.getPlaneEigenSys();

		Eigen::MatrixXd EigVectsMatRowStacked =
			MatrixXd::Zero(3,3);
		EigVectsMatRowStacked =
			getOrthEigVectsMatRowStacked (eigenSys);	
	
		Eigen::MatrixXd diagEigValsMat = MatrixXd::Zero(3,3);
		diagEigValsMat = getDiagEigMat(eigenSys);
	
		eigValsRowWise(reachCounter-1,0) = diagEigValsMat(0,0);
		eigValsRowWise(reachCounter-1,1) = diagEigValsMat(1,1);
		eigValsRowWise(reachCounter-1,2) = diagEigValsMat(2,2);

		/**** Angular model plane - 
						Yaw/Dip/Roll in PCA ****/
		planeAngleToZAlongRows(reachCounter-1,0) = 				
			bfp.getPlaneNtoZ();
/*
			bfp.getPlaneNtoZ()> M_PI/2.0 ?  
			(bfp.getPlaneNtoZ()-M_PI) : bfp.getPlaneNtoZ();
*/	
		pYawZAlongRows(reachCounter-1,0) 
			= bfp.getPlaneYawZ();
		pDipXAlongRows(reachCounter-1,0)
			= bfp.getPlaneDipX();
		pTiltYAlongRows(reachCounter-1,0) 
			= bfp.getPlaneTiltY();

		/***	 Allign points to PCA 	***/
		MatrixXd inWPtsPCAaligned = 
			MatrixXd::Zero(nRows,3);
		inWPtsPCAaligned	= allignPtsToPCA
		 (inWPtsAlongRows, EigVectsMatRowStacked);

		/*** Error space - Hyperbloid Axes ***/    
		Eigen::MatrixXd noiseVarEigVects = 
			buildPCAGausNoiseVarEigVects(nRows, eigenSys);
		Eigen::MatrixXd eigVectsErrs = 
			buildEigVectsErrs(nRows, noiseVarEigVects);

		HypErrShellAxesRowWise.row(reachCounter-1) = 
			 getHypErrShellAxes(eigVectsErrs
			 , diagEigValsMat).transpose();

		/*** the hyperbolic error bounds ***/
		MatrixXd HypErrBounds = MatrixXd::Zero(nRows,2);
		HypErrBounds = getHypErrBounds(HypErrShellAxesRowWise.row
			(reachCounter-1).transpose(), inWPtsPCAaligned);

		/*** Fixed-Length Normal Vector A.2 Error space - 
				 ellipsoid major axes with a center offset
				 √2h3 from the origin along the 3 axis. ***/    
		fixedLEcenterOffsetsAlongRows(reachCounter-1,0) = 
			sqrt(2.0) *	HypErrShellAxesRowWise(reachCounter-1,2);

		fixedLEllipsoidAxesRowWise.row(reachCounter-1) = 
			 getFixedLEllipsoidAxes(HypErrShellAxesRowWise.row
			 (reachCounter-1).transpose()).transpose();


		/*** Max & Min Angular errors ***/
		/**** provides the angular width of the error 
					distribution aligned with the major axes of
					the dataset. ****/
		maxMinAngularErrsV.push_back(getMaxMinAngularErrs
			 (HypErrShellAxesRowWise.row(reachCounter-1)
			 .transpose()));

		/** ELLIPSE FITTING
		
		MatrixXd thetas = fixThetas(e3D.findThetas());
		double cX = ellipse.getRT2D().getCX();
		double cY = ellipse.getRT2D().getCY();
		double a = ellipse.getA();
		double b = ellipse.getB();
		double EfAxisRatio = b/a;
		 **/
		Ellipse3D e3D(cartToHom(inWPtsAlongCols));
		Ellipse ellipse = e3D.getEllipse();
		EfAlphaAlongRows(reachCounter-1,0) 
			= ellipse.getRT2D().getTheta();

	  /**** JOINT - combine all reach samples to one
									big data matrix, each centered arount
									its own mean ****/
/*



	printEigenMathematica(
			HypErrShellAxesRowWise.row(reachCounter-1), cout,
			 "HypErrShellAxes"+ std::to_string(reachCounter));	


		inWPtsAlongRowsJoint.block
				(nRowsPrev, 0, nRows, 3) = inWPtsAlongRows;
*/


		nRowsPrev += nRows;
		reachCounter ++;
	}


	/*** Bayesian linear regression for Plane PC1	, 
				i.e. in-plnae rotation of the ellipse ***/

	/* input matrix X for the linear model XW=y */
	MatrixXd inXmatPyawZ = MatrixXd::Zero(9,nCounterPerN.size());

	// start position (x,y,z)
	inXmatPyawZ.topRows(3) = inWstartPtsAlongRows
		.rowwise().normalized().transpose();
	inXmatPyawZ.row(3) = inWstartPtsAlongRows
		.rowwise().norm().transpose();

	//PC1
	inXmatPyawZ.block(4,0,3,nCounterPerN.size()) = 
		(inWtargetPtsAlongRows - inWstartPtsAlongRows)
		.rowwise().normalized().transpose();

	// Path shortest distance from start to goal
	inXmatPyawZ.row(7) = 
		(inWtargetPtsAlongRows - inWstartPtsAlongRows)
		.rowwise().norm().transpose();

	// target position (x,y,z)
	inXmatPyawZ.row(8) =inWtargetPtsAlongRows
		.rowwise().norm().transpose();
	//inXmatPyawZ.bottomRows(3) = inWtargetPtsAlongRows.transpose();
	//	.rowwise().normalized().transpose();

	/* output vector pYawZAlongRows(n.1) predictions */
	vector<double> WLiPyawZ = findALi
		(inXmatPyawZ, pYawZAlongRows.transpose());

	MatrixXd pYawZhatLi = 
		predictY(inXmatPyawZ, WLiPyawZ);

	/* training error */
	MatrixXd trainEpYawZLi = (pYawZAlongRows 
		- pYawZhatLi.transpose()).cwiseAbs();
	
	cout << trainEpYawZLi.mean() << endl;


	/*** Bayesian linear regression for PlaneNtoZ	***/

	/* input matrix X for the linear model XW=y */
	MatrixXd inXmat = MatrixXd::Zero(6,nCounterPerN.size());
	// start position (x,y,z)
	inXmat.topRows(3) = inWstartPtsAlongRows.transpose();
	//inXmat.row(3) = 
	//	inWstartPtsAlongRows.rowwise().norm().transpose();
	//PC1
	inXmat.block(3,0,3,nCounterPerN.size()) = 
		(inWtargetPtsAlongRows - inWstartPtsAlongRows).transpose();
	//	.rowwise().normalized().transpose();
	// Path shortest distance from start to goal
	//inXmat.row(6) = 
	//	(inWtargetPtsAlongRows - inWstartPtsAlongRows)
	//	.rowwise().norm().transpose();

	// target position (x,y,z)
	//inXmat.row(8) =
	//	inWtargetPtsAlongRows.rowwise().norm().transpose();
	//inXmat.bottomRows(3) = inWtargetPtsAlongRows.transpose();

	/*** compute angles (dot prod) between N-Joint 
			 plane, and N  of the nominal plane of each
			 individual reach sample	 

	for (int i=0; i<v3sAlongRows.rows(); i++)
		inXmat(7,i) = (v3sAlongRows.row(i)
			.transpose().dot(zAxis));
	***/

	/* output vector y(n.1) predictions */
	vector<double> ALi = 
			findALi(inXmat, planeAngleToZAlongRows.transpose());

	MatrixXd planeAngleToZhatLi = 
		predictY(inXmat, ALi);

	/* training error */
	MatrixXd trainEplaneNToZLi = (planeAngleToZAlongRows 
		- planeAngleToZhatLi.transpose()).cwiseAbs();
	
	cout << trainEplaneNToZLi.mean() << endl;

	/* MATHEMATICA OUTPUT: 
	

	printEigenMathematica( planeAngleToZhatLi.transpose()
		, cout, "planeAngleToZhatLi");	
	printEigenMathematica( planeAngleToZAlongRows
		,	cout,	"planeAngleToZAlongRows");
	printEigenMathematica( trainEplaneNToZLi
		, cout, "trainEplaneNToZLi");	

		
	printEigenMathematica( EfAlphaAlongRows
		, cout, "EfAlphaAlongRows");	
	printEigenMathematica( EfAlphaAlongRows
		, cout, "pYawZAlongRows");	
	printEigenMathematica( EfAlphaAlongRows
		, cout, "pTiltYAlongRows");	

	*/ 

	/************************************/
	/**	 joint fit planar fit to all
					 N reach samples per participant 	**/

/*
 	/// Cluster the start of motion, and analyse consistensy 
			//of plane structure for all reach samples per cluster /// 
	
	//fit k planes to each cluster
	int K = 2;
	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> KCluster = 
		KmeansCluster	(inWstartPtsAlongRows, K);
	MatrixXd yClustersV = MatrixXd::Zero(inWstartPtsAlongRows.rows(),1);
	yClustersV = KCluster.first;
	MatrixXd nPerCluster = KCluster.second;
	
	//sor inPts by cluster order
	MatrixXd inWPtsAlongRowsJTracker =
		 MatrixXd::Zero(nRowsJoint,3);

	MatrixXd inWstartPtsAlongRowsJTracker =
		 MatrixXd::Zero(inWstartPtsAlongRows.rows(),3);
	MatrixXd inWtargetPtsAlongRowsJTracker =
		 MatrixXd::Zero(inWtargetPtsAlongRows.rows(),3);

	MatrixXd reachCountTracker = MatrixXd::Zero(reachCountTot,1);
	for (int i=0; i<reachCountTot; i++)
		reachCountTracker(i,0) = i;

	sortMatAscWkey(reachCountTracker,yClustersV);
	sortMatAscColWise(inWstartPtsAlongRows);
	int nRowsJF = 0; 
	for (int i = 0; i<reachCountTracker.rows(); i++){
		int	nRowsJTracker = 0;
		for (int counter=0; counter<reachCountTracker(i); counter++)
			nRowsJTracker +=  nCounterPerN[counter];

		inWPtsAlongRowsJTracker.block(nRowsJF,0, nCounterPerN
			[reachCountTracker(i)], 3) = inWPtsAlongRowsJoint.block
			(nRowsJTracker, 0, nCounterPerN[reachCountTracker(i)], 3);
		
		
		inWstartPtsAlongRowsJTracker.row(i) =
			inWPtsAlongRowsJTracker.row(nRowsJF);
		inWtargetPtsAlongRowsJTracker.row(i) =
			inWPtsAlongRowsJTracker.row(nRowsJF + nCounterPerN
				[reachCountTracker(i)]-1);

		nRowsJF += nCounterPerN[reachCountTracker(i)];

	}

	int nRowsJCluster = 0, Ntot = 0;
	for (int k = 0; k<nPerCluster.rows(); k++){
		int	nRowsCBTracker = 0;

		for (int counter=Ntot; counter 
				< Ntot + nPerCluster(k,0); counter++)
			nRowsCBTracker +=  nCounterPerN[reachCountTracker(counter)];
		

		MatrixXd inWPtsAlongRowsCJ = MatrixXd::Zero(nRowsCBTracker,3);
		inWPtsAlongRowsCJ = inWPtsAlongRowsJTracker.block
			(nRowsJCluster, 0, nRowsCBTracker, 3);

		MatrixXd inWstartPtsAlongRowsCJ =
			inWstartPtsAlongRowsJTracker.block(Ntot,0,nPerCluster(k,0),3);
		MatrixXd inWtargetPtsAlongRowsCJ = 
			inWtargetPtsAlongRowsJTracker.block(Ntot,0,nPerCluster(k,0),3);


		//// JOINT - Eigen System Solver ////
		MatrixXd covMatJoint =
			buildCovMatPCA(inWPtsAlongRowsCJ);
		std::pair<Eigen::MatrixXd, Eigen::MatrixXd> eigenSysJoint 
			= solveEigensystem(covMatJoint); 		

		Eigen::MatrixXd EigVectsMatRowStackedJoint =
			MatrixXd::Zero(3,3);
		EigVectsMatRowStackedJoint =
			getOrthEigVectsMatRowStacked (eigenSysJoint);	
	
		Eigen::MatrixXd diagEigValsMatJoint = MatrixXd::Zero(3,3);
		diagEigValsMatJoint = getDiagEigMat(eigenSysJoint);

		////	JOINT - @Q The nominal Plane in PCA ////
		Vector3d v1Joint = 
			EigVectsMatRowStackedJoint.row(0).transpose();
		Vector3d v2Joint = 
			EigVectsMatRowStackedJoint.row(1).transpose();
		Vector3d v3Joint = 
			EigVectsMatRowStackedJoint.row(2).transpose();

		//// JOINT - Angular model plane - 
						Yaw/Dip/Roll in PCA /////
		double planeYawJointZ = atan2(v3Joint(0,0)
			, v3Joint(1,0)) - M_PI/2.0;

		double planeDipJointX = 
			acos(v3Joint(2,0)/v3Joint.norm());

		Vector3d crossProduct = v3Joint.cross(zAxis);
		double planeTiltJointY = acos(crossProduct.dot(v2Joint));
*/

		/*** JOINT -  Projection to PCA-aligned bases ***/
			/**** to keep track and discard unparallel motions ****/
/*	
		// nx3
		MatrixXd inWPtsPCAaligned = 
			MatrixXd::Zero(nRowsCBTracker,3);
		inWPtsPCAaligned	= allignPtsToPCA
		 (inWPtsAlongRowsCJ , EigVectsMatRowStackedJoint);

		//// JOINT - Error space - Hyperbloid Axes ////    
		Eigen::MatrixXd HypErrShellAxesJoint = 
			MatrixXd::Zero(3,1);

		Eigen::MatrixXd noiseVarEigVectsJoint = 
			buildPCAGausNoiseVarEigVects(nRowsCBTracker, eigenSysJoint);
		Eigen::MatrixXd eigVectsErrsJoint = 
			buildEigVectsErrs(nRowsCBTracker, noiseVarEigVectsJoint);

		//3.1
		HypErrShellAxesJoint = getHypErrShellAxes
			 (eigVectsErrsJoint, diagEigValsMatJoint);
*/

		/*** JOINT - Fixed-Length Normal Vector A.2
			Error space - ellipsoid major axes with a center
		  offset √2h3 from the origin along the 3 axis. ***/    

/*
		double fixedLEcenterOffsetJoint = 
			sqrt(2.0) *	HypErrShellAxesJoint(2,0);
		Eigen::MatrixXd fixedLEllipsoidAxesJoint =
			 MatrixXd::Zero(3,1);
		fixedLEllipsoidAxesJoint =	
			getFixedLEllipsoidAxes(HypErrShellAxesJoint);
*/

		/*** JOINT - the hyperbolic error bounds ***/
/*
		MatrixXd HypErrBounds = MatrixXd::Zero(nRowsCBTracker,2);
		HypErrBounds = getHypErrBounds
			(HypErrShellAxesJoint, inWPtsPCAaligned);

*/
		/*** JOINT - Max & Min Angular errors ***/
		/**** provides the angular width of the error 
					distribution aligned with the major axes of
					the dataset. ****/
/*
		std::pair<double, double> maxMinAngularErrsJoint = 
			getMaxMinAngularErrs(HypErrShellAxesJoint);

		nRowsJCluster += nRowsCBTracker;
		Ntot += nPerCluster(k,0);

		printEigenMathematica (inWPtsAlongRowsCJ, cout
			,"inWPtsAlongRowsCJ"+ std::to_string(k));		
		printEigenMathematica (inWPtsPCAaligned, cout
			,"inWPtsPCAaligned"+ std::to_string(k));		
		printEigenMathematica (HypErrBounds, cout
			, "HypErrBounds"+ std::to_string(k));		
		printEigenMathematica (inWstartPtsAlongRowsCJ, cout
			,"inWstartPtsAlongRows"+ std::to_string(k));		
		printEigenMathematica (inWtargetPtsAlongRowsCJ, cout
			,"inWtargetPtsAlongRows"+ std::to_string(k));		

		printEigenMathematica
	 		 (v1Joint, cout,"v1Joint"+ std::to_string(k));		
		printEigenMathematica
	 		 (v2Joint, cout,"v2Joint"+ std::to_string(k));		
		printEigenMathematica
	 		 (v3Joint, cout,"v3Joint"+ std::to_string(k));

	}
*/



	/**	 test and elliminate reach samples 
					with high planar fit residual 	**/

/*
		MatrixXd v3DotV3JalongRows =
		 MatrixXd::Zero(partIdsV.size(),1);
*/

	/*** compute angles (dot prod) between N-Joint 
			 plane, and N  of the nominal plane of each
			 individual reach sample	 ***/
/*
	for (int i=0; i<partIdsV.size(); i++){

		Vector3d v3 = v3sAlongRows.block(i,0,1,3).transpose();
		v3DotV3JalongRows(i,0) = abs(v3Joint.dot(v3));


		if (v3DotV3JalongRows(i,0) >= 0.8){
			inWPtsAlongRowsJTracker.block(nRowsJF,0
				 , nCounterPerN[i], 3) = inWPtsAlongRowsJoint.block
				 (nRowsTracker, 0, nCounterPerN[i], 3);

			reachCountTrackerFV.push_back(i);
			nRowsJF += nCounterPerN[i];
		
		}
		nRowsTracker += nCounterPerN[i];

	}
*/

	/*** JOINT F-  Projection to PCA-aligned bases ***/
/*
	MatrixXd inWPtsPCAalignedJF = 
		 MatrixXd::Zero(nRowsJF,3);
	inWPtsPCAalignedJF	= allignPtsToPCA
		 (inWPtsAlongRowsJF , EigVectsMatRowStackedJoint);
*/
	/*** JOINT F- the hyperbolic error bounds ***/
/*
	MatrixXd HypErrBoundsJF = MatrixXd::Zero(nRowsJF,2);
	HypErrBoundsJF = getHypErrBounds
		(HypErrShellAxesJoint, inWPtsPCAalignedJF);
*/

	/**		 A.1 END Finding Nominal Plane Fit 		**/ 
	/************************************************/


	/**********************************/
	/** 			Terminal Output 			**/ 
   
   /**  joint fit planar fit to all
                     N reach samples per participant    **/
/*
   cout << *partIdsV.begin() << ","
        << *reachCatsV.begin() << ","
        << reachCountTot << "/"
        << reachCountTrackerFV.size()<< ", , "
 
        << diagEigValsMatJoint(0,0) << ", "
        << diagEigValsMatJoint(1,1) << ", "
        << diagEigValsMatJoint(2,2) << ", "
         
        // 3 Eigen Vectors v1,v2,v3 as regression pars //
        << v1Joint(0,0) << ", "
        << v1Joint(1,0) << ", "
		    << v1Joint(2,0) << ", "
		    << v2Joint(0,0) << ", "
		    << v2Joint(1,0) << ", "
		    << v2Joint(2,0) << ", "
		    << v3Joint(0,0) << ", "
		    << v3Joint(1,0) << ", "
		    << v3Joint(2,0) << ", , "
		      
 
        // Plane angular structure //
        << planeYawJointZ << ", "
        << planeDipJointX << ", "
        << planeTiltJointY << ", "
 
 
        /// The Cartesian error space of fitted 
             planar orientation measurements ///
              
        << maxMinAngularErrsJoint.first << ", "
        << maxMinAngularErrsJoint.second << ", "
 
        // Maximum Ellipsoid (error) shell h1,h2,h3-  // 
        << HypErrShellAxesJoint(0,0) << ", "
        << HypErrShellAxesJoint(1,0) << ", "
        << HypErrShellAxesJoint(2,0) << ", "
 
        //Fixed Length N-vector
        << fixedLEcenterOffsetJoint << ", "
        << fixedLEllipsoidAxesJoint(0,0) << ", "
        << fixedLEllipsoidAxesJoint(1,0) << ", "
        << fixedLEllipsoidAxesJoint(2,0) << endl; 
*/

/*  
    for (std::vector<int>::const_iterator i 
        = reachCountTrackerFV.begin(); i 
				!= reachCountTrackerFV.end(); ++i){

	for (int i=0; i<partIdsV.size(); i++){
*/
        /*** compute angles (dot prod) between N-Joint 
                 plane, and N  of the nominal plane of each
                 individual reach sample     ***/
/*
        Vector3d v3 = v3sAlongRows.block(i,0,1,3).transpose();
        //double v3DotV3Joint = abs(v3Joint.dot(v3));

        cout << partIdsV[i] << "," 
        //cout << i+1 << ","
             << reachCatsV[i] << ","
             << nCounterPerN[i] << ","
						 //<< yClustersV(i,0) << ","

             << eigValsRowWise(i,0) << ", "
             << eigValsRowWise(i,1) << ", "
             << eigValsRowWise(i,2) << ", "
 
        // 3 Eigen Vectors v1,v2,v3 as regression pars //
             << v1sAlongRows(i,0) << ", "
             << v1sAlongRows(i,1) << ", "
             << v1sAlongRows(i,2) << ", "
             << v2sAlongRows(i,0) << ", "
             << v2sAlongRows(i,1) << ", "
             << v2sAlongRows(i,2) << ", "
             << v3sAlongRows(i,0) << ", "
             << v3sAlongRows(i,1) << ", "
             << v3sAlongRows(i,2) << ", "
            // << v3DotV3Joint << ", " 

        // Plane angular structure //
		         << pYawZAlongRows(i,0) << ", "
		         << pDipXAlongRows(i,0) << ", "
		         << pTiltYAlongRows(i,0) << ", "

		         << planeAngleToZAlongRows(i,0) << ", "



        /// The Cartesian error space of fitted 
           //  planar orientation measurements ///
              
             << maxMinAngularErrsV[i].first << ", "
             << maxMinAngularErrsV[i].second << ", "

            // Maximum Ellipsoid (error) shell h1,h2,h3- //
              
             << HypErrShellAxesRowWise(i,0) << ", "
             << HypErrShellAxesRowWise(i,1) << ", "
             << HypErrShellAxesRowWise(i,2) << ", "

            //Fixed Length N-vector
             << fixedLEcenterOffsetsAlongRows(i,0) << ", "
             << fixedLEllipsoidAxesRowWise(i,0) << ", "
             << fixedLEllipsoidAxesRowWise(i,1) << ", "
             << fixedLEllipsoidAxesRowWise(i,2) << endl; 
    }
*/

		/*** 
				MATHEMATICA OUTPUT: 
		***/ 	
	/*
		printEigenMathematica
	 		 (inWPtsAlongRowsJF, cout,"inWPtsAlongRowsJF");		
		printEigenMathematica
	 		 (inWPtsAlongRowsJoint, cout,"inWPtsAlongRowsJ");		

		printEigenMathematica
	 		 (v1Joint, cout,"v1Joint");		
		printEigenMathematica
	 		 (v2Joint, cout,"v2Joint");		
		printEigenMathematica
	 		 (v3Joint, cout,"v3Joint");
			
		printEigenMathematica
	 		 (planeAngleToZAlongRows, cout,"planeAngleToZAlongRows");		

		printEigenMathematica
	 		 (pYawZAlongRows, cout,
			 "pYawZAlongRows");	
		printEigenMathematica
	 		 (pDipXAlongRows, cout,
			 "pDipXAlongRows");	
		printEigenMathematica
	 		 (pTiltYAlongRows, cout,
			 "pTiltYAlongRows");	

		printEigenMathematica
	 		 (inWPtsPCAaligned, cout,"inWPtsPCAaligned");		
		printEigenMathematica
			(HypErrBounds, cout, "HypErrBounds");	

		cout << "fixedLEcenterOffsetJoint: " 
				 << fixedLEcenterOffsetJoint << endl;	
	 	
*/

	/*** END MATHEMATICA OUTPUT: 
	***/ 

	return 0;
}


/**** END ****/

	
/****** My method - 
		Plane Normal Angles to XY *****/
	/*
	Vector3d PlaneNormalVectN = 
		e3D.getEIF().getPlane().getPlaneNormalVectN();
	double dotProduct = PlaneNormalVectN.dot(zAxis);
	Vector3d crossProductt = 
		PlaneNormalVectN.cross(zAxis);

	double planeAngleToXY = 
		atan2(crossProductt.norm(), dotProduct);
	*/
//************************* Bayes Fitting *************************/
	
/*

	int maxPolyOrder = 6;

	
	MatrixXd x_1Vec = (inWPtsAlongCols.block(0,0,1,n)).transpose();	
	MatrixXd x_2Vec = (inWPtsAlongCols.block(1,0,1,n)).transpose();	
	MatrixXd x_3Vec = (inWPtsAlongCols.block(2,0,1,n)).transpose();	
	MatrixXd xVec = MatrixXd::Zero(n,1);	
	MatrixXd yVec = MatrixXd::Zero(n,1);

	
	for(size_t i = 0; i < n; i++) { 
		xVec(i,0) = x_1Vec(i,0)/x_3Vec(i,0);
		yVec(i,0) = x_2Vec(i,0)/x_3Vec(i,0);	
	}

	MatrixXd evidenceMat(maxPolyOrder+1,1);		
	MatrixXd AMetaLi = MatrixXd::Zero(maxPolyOrder+1,maxPolyOrder+1);	
	MatrixXd AMetaEv = MatrixXd::Zero(maxPolyOrder+1,maxPolyOrder+1);
	MatrixXd likelihoodMat(maxPolyOrder+1,1);
	Eigen::MatrixXd yVecHatMetaLi	= 
			MatrixXd::Zero(maxPolyOrder+1,n);
	
	//2d
	MatrixXd evidenceMat2d(maxPolyOrder+1,1);		
	MatrixXd AMetaLi2d = MatrixXd::Zero(maxPolyOrder+1,maxPolyOrder+1);	
	MatrixXd AMetaEv2d = MatrixXd::Zero(maxPolyOrder+1,maxPolyOrder+1);
	MatrixXd likelihoodMat2d(maxPolyOrder+1,1);
	Eigen::MatrixXd yVecHatMetaLi2d	= 
			MatrixXd::Zero(maxPolyOrder+1,n);

	cout << argv[1] << ": " << endl;
		for(size_t k = 0; k < maxPolyOrder; k++) {
			Eigen::MatrixXd Xmat = MatrixXd::Zero(k+1,n);
			Xmat = buildXMatIn(xVec, k);			
			Eigen::MatrixXd Xmat2d = MatrixXd::Zero(k+1,n);
			Xmat2d = buildXMatIn(xVec2d, k);

			vector<double> AvecLi = findALi(Xmat, yVec.transpose());	
			vector<double> AvecLi2d = findALi(Xmat2d, yVec2d.transpose());	

			for(size_t j = 0; j < AvecLi.size(); j++) {
				AMetaLi(k,j) = AvecLi[j];
				AMetaLi2d(k,j) = AvecLi2d[j];
			}
		
			Eigen::MatrixXd yVecHatLi = 
				predictY(Xmat, AMetaLi.block(k,0,1,k+1));
			Eigen::MatrixXd yVecHatLi2d = 
				predictY(Xmat2d, AMetaLi2d.block(k,0,1,k+1));
			
			yVecHatMetaLi2d.block(k,0,1,n) = yVecHatLi2d;

		
			evidenceMat(k,0) = 
				computePy_xvalpha(Xmat, yVec.transpose());
			likelihoodMat(k,0) = 
				computePy_xav(Xmat, yVec.transpose());
			evidenceMat2d(k,0) = 
				computePy_xvalpha(Xmat, yVec2d.transpose());
			likelihoodMat2d(k,0) = 
				computePy_xav(Xmat2d, yVec2d.transpose());


	MatrixXd outPtsAlongColsBays = MatrixXd::Zero(3,n);
*/
	/*
		for(size_t i = 0; i < n; i++) { 
			//outPtsAlongColsBays(0,i) = xVec(i,0) * x_3Vec(i,0);
			yVecHatMetaLi(k,i) = yVecHatLi(0,i) * x_3Vec(i,0);	
			//outPtsAlongColsBays(2,i) = x_3Vec(i,0);	
		}

		
		//cout<<outPtsAlongColsBays.block(0,0,1,n).transpose()<<endl;
		//cout<<outPtsAlongColsBays.block(1,0,1,n).transpose()<<endl;
		//cout<<outPtsAlongColsBays.block(2,0,1,n).transpose()<<endl;

	//printEigenMathematicaSci(evidenceMat, cout, "evidenceMat");	
	//printEigenMathematicaSci(outPtsAlongColsBays, cout, 		"outPtsAlongColsBays");

	}

	//printEigenMathematicaSci(yVecHatMetaLi, cout, "yVecHatMetaLi");
	//printEigenMathematicaSci(xVec, cout, "xVec");
	//cout << evidenceMat(maxPolyOrder-1,0) << endl;	

		
	//cout << "evidence for " << argv[1] <<":\n" << evidenceHom << endl;
	//cout << "likelihood for " << argv[1] <<":\n" << likelihoodMat(k,0) << endl;
*/
	
	/*********** End - Bayes Fitting ***********/



	/*********** Bayes Fitting ***********/
/*
		int n = inPts2dRotatedAlongCols.cols();
		MatrixXd xVec = (inPts2dRotatedAlongCols.block(0,0,1,n)).transpose();
		MatrixXd yVec = (inPts2dRotatedAlongCols.block(1,0,1,n)).transpose();	
		int maxPolyOrder = 5;
		MatrixXd evidenceMat(maxPolyOrder+1,1);		
		MatrixXd AMetaLi = MatrixXd::Zero(maxPolyOrder+1,maxPolyOrder+1);	
		MatrixXd AMetaEv = MatrixXd::Zero(maxPolyOrder+1,maxPolyOrder+1);
		MatrixXd likelihoodMat(maxPolyOrder+1,1);
		Eigen::MatrixXd yVecHatMetaLi	= 
			MatrixXd::Zero(maxPolyOrder+1,n);
		Eigen::MatrixXd yVecHatMetaEv	= 
			MatrixXd::Zero(maxPolyOrder+1,n);
size_t k = 4;
		//for(size_t k = 0; k < maxPolyOrder; k++) {
			Eigen::MatrixXd Xmat = MatrixXd::Zero(k+1,n);
			Xmat = buildXMatIn(xVec, k);
			
			vector<double> AvecLi = findALi(Xmat, yVec.transpose());	
		
			for(size_t j = 0; j < AvecLi.size(); j++) 
				AMetaLi(k,j) = AvecLi[j];

		
		Eigen::MatrixXd yVecHatLi = 
				predictY(Xmat, AMetaLi.block(k,0,1,k+1));
		yVecHatMetaLi.block(k,0,1,n) = yVecHatLi;

		
		evidenceMat(k,0) = 
				computePy_xvalpha(Xmat, yVec.transpose());
		likelihoodMat(k,0) = 
				computePy_xav(Xmat, yVec.transpose());
	//}

*/
	/*********** End - Bayes Fitting ***********/

	
	/*********** Mathematica Plotting fncs ***********/
	/*
	printEigenMathematica(inWPtsAlongCols, cout, "inPts");
	printEigenMathematica(e3D.getEIF().getPtsXYAlongCols()
		, cout, "inPtsRotated");
	printEigenMathematica(ptsOnEllipsePseudo2D
		, cout, "outPts2d");
	printEigenMathematica(ptsOnEllipsePseudo3D
		, cout, "outPts3d");
	printEigenMathematica(fixedSignThetas
		, cout, "Thetas");
	printEigenMathematica(speed
		, cout, "speed");
	printEigenMathematica(acceleration
		, cout, "acceleration");
	printEigenMathematica(jerk
		, cout, "jerk");	*/

	/********* End - Mathematica Plotting fncs *********/

