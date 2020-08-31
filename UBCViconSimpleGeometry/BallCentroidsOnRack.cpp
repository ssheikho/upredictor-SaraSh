#include "BallCentroidsOnRack.h"

std::vector<std::string> BallCentroidsOnRack::_rackMarkers = {
	"rackM1_X", "rackM1_Y", "rackM1_Z",
	"rackM2_X", "rackM2_Y", "rackM2_Z",
	"rackM3_X", "rackM3_Y", "rackM3_Z",
	"rackM4_X", "rackM4_Y", "rackM4_Z",
	"rackM5_X", "rackM5_Y", "rackM5_Z",
	"rackM6_X", "rackM6_Y", "rackM6_Z"
};

BallCentroidsOnRack::BallCentroidsOnRack(string inLine):
	_inLine(inLine)
	, _nCols(18)
	, _rackMarkerLoc(MatrixXd::Zero(3,6))
	, _ballCentroids(MatrixXd::Zero(3,4)) {

	setMarkerLocs();
	setBallCentroids();
}

BallCentroidsOnRack::~BallCentroidsOnRack(){}


void BallCentroidsOnRack::setMarkerLocs(){
	
	int colNo = 0;
	int nRows = countRowsCSV(_inLine, _rackMarkers[colNo]);

	MatrixXd outMatRack(nRows,_nCols);
	MatrixXd rackMarkerLocRow(1,_nCols);	
	
	
	while (colNo < _nCols-2){
		fillMatCSV(colNo, _inLine, _rackMarkers[colNo], 
			_rackMarkers[colNo + 1], 			_rackMarkers[colNo + 2], outMatRack);
		
		colNo = colNo + 3;
	}

	for (int col = 0; col < outMatRack.cols(); col++){
		
		double sum = 0.0;
		double ctr = 0.0;
		double mean = 0.0;
		
		for (int row = 1; row < outMatRack.rows(); row++){
<<<<<<< HEAD
			if(outMatRack(row,col)!=0.0
				/*&& abs(outMatRack(row,col)-outMatRack(row-1,col))<1.5*/){
=======
			if(outMatRack(row,col)!=0.0){
>>>>>>> e070f2e51ee0c62744927929a7a1819317e943b7
				sum += outMatRack(row,col); 
				ctr++;
			}
		}
		mean = sum/double(ctr);
		rackMarkerLocRow(0,col) = mean;
	}
	colNo = 0;
	int mult = 0;
	for (int col = 0; col < _rackMarkerLoc.cols(); col++){
		
		_rackMarkerLoc(0,col) = 
		rackMarkerLocRow(0,(outMatRack.cols()-mult)-3);

		_rackMarkerLoc(1,col) = 
		rackMarkerLocRow(0,(outMatRack.cols()-mult)-2);
		
		_rackMarkerLoc(2,col) = 
		rackMarkerLocRow(0,(outMatRack.cols()-mult)-1);

		mult += 3;
	}
	
}

void BallCentroidsOnRack::setBallCentroids(){
	MatrixXd Offsets(3,1);	
	Offsets << 6.0,
	           1.0,
                       0.0;	

	for (int col = 0; col < _ballCentroids.cols(); col++){
		_ballCentroids.block(0,col,3,1) = 
			_rackMarkerLoc.block(0,col+1,3,1) + Offsets;
	}
}

MatrixXd BallCentroidsOnRack::getMarkerLocs(){
	return _rackMarkerLoc;
	}

MatrixXd BallCentroidsOnRack::getBallCentroids(){
	return _ballCentroids;
	}





