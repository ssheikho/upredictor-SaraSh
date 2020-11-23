#ifndef PARSE_CSV_H
#define PARSE_CSV_H

#include <Eigen/Core>
#include "csv.h"
#include <string>

using namespace Eigen;
using namespace std;

#define stringify( name ) # name

enum markerLocs {
	
	LChestX, LChestY, LChestZ,		/* 30 */
	LElbX, LElbY, LElbZ,					/* 39 */
	LLArmX,	LLArmY,	LLArmZ,				/* 42 */
	LPinkieX,	LPinkieY,	LPinkieZ,	/* 51 */
	LShoX,	LShoY,	LShoZ,				/* 33 */
	LThumbX,	LThumbY,	LThumbZ,	/* 48 */
	LUArmX,	LUArmY,	LUArmZ,				/* 36 */
	LWristX,	LWristY,	LWristZ,	/* 45 */	

	MChestX,	MChestY,	MChestZ,	/* 27 */

	RChestX,	RChestY,	RChestZ,	/* 24 */
	RElbX,	RElbY,	RElbZ,				/* 15 */
	RLArmX,	RLArmY,	RLArmZ,				/* 12 */
	RPinkieX,	RPinkieY,	RPinkieZ,	/* 3 */
	RShoX,	RShoY,	RShoZ,				/* 21 */	
	RThumbX,	RThumbY,	RThumbZ,	/* 6 */	
	RUArmX,	RUArmY,	RUArmZ,				/* 18 */	
	RWristX, RWristY, RWristZ 		/* 9 */

};
const int NumMarkers = 51;


// templated version of my_equal so it could work with both char and wchar_t
template<typename charT>
struct my_equal {
    my_equal( const std::locale& loc ) : loc_(loc) {}
    bool operator()(charT ch1, charT ch2) {
        return std::toupper(ch1, loc_) == std::toupper(ch2, loc_);
    }
private:
    const std::locale& loc_;
};

// find substring (case insensitive)
template<typename T>
int ci_find_substr( const T& str1, const T& str2, const std::locale& loc = std::locale() )
{
    typename T::const_iterator it = std::search( str1.begin(), str1.end(), 
        str2.begin(), str2.end(), my_equal<typename T::value_type>(loc) );
    if ( it != str1.end() ) return it - str1.begin();
    else return -1; // not found
}

string getWholeFile(string fileName);
string getHeaderRowCSV(string fileName);


string getCharfromString(string fileName);
string getPartID(string fileName);
int whichPhase(string fileName);

int countCols(string inLine);
int countRowsCSV(string inMat, string inCol);
int countRowsFColCSV(string inMat); 


void setHeaderCSV(string inLine);
int countEmptyHeadColsCSV (string inLine);
std::vector<std::string> getRecordedMarkers(string inLine);


void fillMatCSV(int inColNo, string inLine, string inX,
	 string inY, string inZ, MatrixXd &mat);
void fillWholeMatCSV(int inColNo, string inLine, MatrixXd &mat);



MatrixXd parseCSV(string inMat);

#endif
