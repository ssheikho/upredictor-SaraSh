
#include "csv.h"
#include "ParseCSV.h"

#include "BallCentroidsOnRack.h"
#include "BallCentroidsOnRamp.h"
#include <UBCUtil.h>

#include <Eigen/Dense>

#include <iostream>

using namespace Eigen;
using namespace std;

int main(int argc, char **argv) {
<<<<<<< HEAD
	srand(time(NULL));
	string inLine = argv[1];
	BallCentroidsOnRack bcRack(inLine);
	MatrixXd ballsCentroidsbcRack = bcRack.getBallCentroids();

	printEigenMathematica(bcRack.getMarkerLocs(), cout, "MarkerLocs");
	//BallCentroidsOnRamp bcRamp(inLine);

	//cout << argv[1] <<endl; 
	//cout << bcRack.getMarkerLocs()<< endl;

	//cout <<"markers Centroids on Ramp: " <<endl; 
	//cout << bcRamp.getMarkerLocs() << endl;
=======

	string inLine = "2016_10_28_empty_Setup2.csv";
	BallCentroidsOnRack bcRack(inLine);
	
	MatrixXd ballsCentroidsbcRack = bcRack.getBallCentroids();
	cout <<"balls Centroids on Rack: " <<endl; 
	cout << bcRack.getBallCentroids() << endl;

	BallCentroidsOnRamp bcRamp(inLine);
	cout <<"markers Centroids on Ramp: " <<endl; 
	cout << bcRamp.getMarkerLocs() << endl;
>>>>>>> e070f2e51ee0c62744927929a7a1819317e943b7
	
	return 0;
}

