#include "csv.h"
#include "ParseCSV.h"

#include "BallCentroidsOnRamp.h"
#include <UBCUtil.h>

std::vector<std::string> BallCentroidsOnRamp::_rampMarkers = {
	"rampM1_X", "rampM1_Y", "rampM1_Z",
	"rampM2_X", "rampM2_Y", "rampM2_Z",
	"rampM3_X", "rampM3_Y", "rampM3_Z"
};

BallCentroidsOnRamp::BallCentroidsOnRamp(string inLine):
	_inLine(inLine)
	, _nCols(9)
	, _rampMarkerLoc(MatrixXd::Zero(3,3))
	, _rampBallCentroids(MatrixXd::Zero(3,5)) {

	setMarkerLocs();
	//setBallCentroids();
}

BallCentroidsOnRamp::~BallCentroidsOnRamp(){}


void BallCentroidsOnRamp::setMarkerLocs(){
	
	int colNo = 0;
	int nRows = countRowsCSV(_inLine, _rampMarkers[colNo]);

	MatrixXd outMatRamp(nRows,_nCols);
	MatrixXd rampMarkerLocRow(_nCols,1);	
	
	
	while (colNo < _nCols-2){
		fillMatCSV(colNo, _inLine, _rampMarkers[colNo], 
			_rampMarkers[colNo + 1], _rampMarkers[colNo + 2], outMatRamp);
		
		colNo = colNo + 3;
	}

	for (int col = 0; col < outMatRamp.cols(); col++){
		
		double sum = 0.0;
		double ctr = 0.0;
		double mean = 0.0;
		
<<<<<<< HEAD
		for (int row = 2; row < outMatRamp.rows(); row++){
			if(outMatRamp(row,col)!=0.0 
				/*&& abs(outMatRamp(row,col)-outMatRamp(row-1,col))<1.5*/){
=======
		for (int row = 1; row < outMatRamp.rows(); row++){
			if(outMatRamp(row,col)!=0.0){
>>>>>>> e070f2e51ee0c62744927929a7a1819317e943b7
				sum += outMatRamp(row,col); 
				ctr++;
			}
		}
		mean = sum/double(ctr);
		rampMarkerLocRow(col,0) = mean;
	}
	colNo = 0;
	int mult = 0;
	for (int col = 0; col < _rampMarkerLoc.cols(); col++){
		
		_rampMarkerLoc.block(0,col,3,1) = 
		rampMarkerLocRow.block(mult,0,3,1);

		//_rampMarkerLoc(0,col) = 
		//rampMarkerLocRow(0,(outMatRamp.cols()-mult)-3);

		//_rampMarkerLoc(1,col) = 
		//rampMarkerLocRow(0,(outMatRamp.cols()-mult)-2);
		
		//_rampMarkerLoc(2,col) = 
		//rampMarkerLocRow(0,(outMatRamp.cols()-mult)-1);

		mult += 3;
	}
	
}

void BallCentroidsOnRamp::setBallCentroids(){
	// ramp X direction/width 
	// total:22.25"/56.5cm space:3.25"/8.25cm devider-wall thickness: 1"/2.54cm 
	// move by 2.125"/5.4cm in X-direction

	// ramp Y direction/length
	// total:36"+2 3/4"/98.5cm space:3.25"/8.25cm devider-wall thickness: 1"/2.54cm 
	MatrixXd Offsets(3,1);	
	Offsets << -7.0,
	           -7.0,
							0.0;	

	for (int col = 0; col < _rampBallCentroids.cols(); col++){
		_rampBallCentroids.block(0,col,3,1) = 
			_rampMarkerLoc.block(0,col+1,3,1) + Offsets;
	}
}

// markers 1&2 are the lower markers 
MatrixXd BallCentroidsOnRamp::getMarkerLocs(){
	return _rampMarkerLoc;
	}

MatrixXd BallCentroidsOnRamp::getBallCentroids(){
	return _rampBallCentroids;
	}





