#ifndef WRIST_DISPLACEMENT_H
#define WRIST_DISPLACEMENT_H

#include <Eigen/Core>
#include <string>

#include "ParseCSV.h"
#include <UBCUtil.h>
using namespace Eigen;
using namespace std;

class WristDisplacement {
public:
	WristDisplacement(string inLine);
	~WristDisplacement();

	MatrixXd getWristDisplacementAcrossCols();
	void setWristDisplacementAcrossCols();

	static std::vector<std::string> _wristMarkers;	

protected:
	
	string _inLine;
	int _nCols;	
	MatrixXd _wristDisplacementAcrossCols;
};

#endif
