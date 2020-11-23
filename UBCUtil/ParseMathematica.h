#ifndef PARSE_MATHEMATICA_H
#define PARSE_MATHEMATICA_H

#include <Eigen/Core>
#include <string>

using namespace Eigen;
using namespace std;

string getWholeFileMathematica(string fileName);
int countColsMathematica(string inLine);
int countRowsMathematica(string inMat);
void fillMat(int inRowNo, string inLine, MatrixXd &mat);
MatrixXd parseMathematica(string inMat);

#endif
