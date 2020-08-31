#include "WristDisplacement.h"

std::vector<std::string> WristDisplacement::_wristMarkers = {
	"Rwrist_X",	"Rwrist_Y",	"Rwrist_Z"
};

WristDisplacement::WristDisplacement(string inLine):
	_inLine(inLine)
	, _nCols(3)
	, _wristDisplacementAcrossCols(MatrixXd::Zero(_nCols, 
			countRowsCSV(_inLine, _wristMarkers[0]))){

	setWristDisplacementAcrossCols();
}

WristDisplacement::~WristDisplacement(){}


void WristDisplacement::setWristDisplacementAcrossCols(){
	
	int colNo = 0;
	int nRows = countRowsCSV(_inLine, _wristMarkers[colNo]);
	MatrixXd outMatWristDis(nRows,_nCols);
	fillMatCSV(0, _inLine, _wristMarkers[colNo], 
			_wristMarkers[colNo + 1], _wristMarkers[colNo + 2], outMatWristDis);
	
	_wristDisplacementAcrossCols = outMatWristDis.transpose();
	
}

MatrixXd WristDisplacement::getWristDisplacementAcrossCols(){
	return _wristDisplacementAcrossCols;
	}




