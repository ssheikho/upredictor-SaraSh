
#include "csv.h"
#include <fstream>
#include <iostream>
#include "WristDisplacement.h"

#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <map>

using namespace std;

int main(int argc, char **argv) {
	string inLine = "2016_10_26_participant1_1_90_407.csv";
 	WristDisplacement wd(inLine);
	//ViconData *vd = loadVD(f);

	MatrixXd wPosAcrossCols = wd.getWristDisplacementAcrossCols();
	//cout <<"Wrist Displacement Across Cols: " << endl; 
	//cout << wPosAcrossCols << endl;

	MatrixXd speed = simpleGradientAcrossCols(wPosAcrossCols);\


	printEigenMathematica(wPosAcrossCols, cout, "position");
	printEigenMathematica(speed, cout, "speed");

	return 0;
}
