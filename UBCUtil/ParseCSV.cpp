#include "ParseCSV.h"

#include <fstream>

#include <cstdlib>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

string getWholeFile(string fileName) {
	ifstream infile(fileName);
	string s, out = "";
	while (!infile.eof()) {	
		getline(infile, s);
		out += s;
	}
	infile.close();
	return out;
}

string getHeaderRowCSV(string fileName) {
	/* READ FILE */
	ifstream infile(fileName);
	string s, out = "";
	getline(infile, s);
	out += s;
	infile.close();
	return out;
}

string getPartID(string fileName) {

	/* position of "Part" in fileName; where ID numbering begins 
		in the current file naming structure */

	// box packing Study 1 	
	size_t startPos = fileName.find("participant");
	if (startPos!=std::string::npos) 
		fileName = fileName.substr(startPos + 11);
	
	// box packing Study 2
	else { startPos = fileName.find("Part");

		if (startPos!=std::string::npos) 
			fileName = fileName.substr(startPos + 4);
	}
	int len = fileName.find('_');
	string retS = fileName.substr(0,len);

	return retS;
}


void getStartFrame(string fileName, int& startF, int& endF) {

	size_t endPos = fileName.find("_clean");	
	int slashDelim = fileName.find("/");	

	while((slashDelim >= 0) && (slashDelim < endPos)) {
		fileName = fileName.substr(slashDelim+1,endPos );
		slashDelim = fileName.find("/");	
		endPos = fileName.find("_clean");	
	}

	size_t startPos = fileName.find("Part");	
	if (startPos!=std::string::npos) 
		fileName = fileName.substr(startPos + 4);

	int dashDelim = fileName.find('_');
	if(dashDelim !=std::string::npos) {
		fileName = fileName.substr(dashDelim + 1);
		dashDelim = fileName.find('_');

		fileName = fileName.substr(dashDelim + 1);
	}

	int len = fileName.find('_');
	string startFrame = fileName.substr(0,len);		
	startF = stoi(startFrame);
	//end frame
	fileName = fileName.substr(len + 1);
	len = fileName.find('_');
	string endFrame = fileName.substr(0,len);
	endF = stoi(endFrame);

}

int whichPhase(string fileName) {


		string retS = fileName.substr(
			((unsigned)(fileName.length()-5)), 1);
		int phaseNum = stoi(retS);
		return phaseNum;
}

int countCols(string inLine) {
	inLine = inLine.substr(1);
	int commaDelim = inLine.find(',');
	int endBracketDelim = inLine.find('}');
	int ctr = 0;
	while((commaDelim >= 0) && (commaDelim < endBracketDelim)) {
		//cout << inLine.substr(0, commaDelim) << endl;		
		inLine = inLine.substr(commaDelim + 1);
		commaDelim = inLine.find(',');
		endBracketDelim = inLine.find('}');
		ctr++;
	}
	ctr++;
	//cout << inLine.substr(0, endBracketDelim) << endl;
	return ctr;
}

int countRowsCSV(string inMat, string inCol) {

	io::CSVReader<1> in(inMat);

	in.read_header(io::ignore_extra_column, inCol);
	in.next_line();
	int ctr = 0;
	while(in.read_row(inCol)){
		ctr++;

	}
	return ctr - 1;

}

int countRowsFColCSV(string inMat) {

	io::CSVReader<3> in(inMat);
	//in.read_header(io::ignore_extra_column, inCol);
	int ctr = 0 ;
	double x, y, z;
	while(in.read_row(x, y, z)) ctr++;
	return ctr;
}

void fillMatCSV(int inColNo, string inLine, string inX, 
	string inY, string inZ, MatrixXd &mat) {

	io::CSVReader<3> in(inLine);
	in.read_header(io::ignore_extra_column, inX, inY, inZ);
	double inX_, inY_, inZ_;
	in.next_line();
	int rowNo = 0;
	while(in.read_row(inX_, inY_, inZ_) && rowNo < mat.rows()){	
		// inX_, inY_, inZ_: contain the value from the file
		mat(rowNo, inColNo) = inX_;
		mat(rowNo, inColNo+1) = inY_;
		mat(rowNo, inColNo+2) = inZ_;

		rowNo++;
	}

}

void setHeaderCSV(string inLine) {

	std::vector<std::string> wristMarkers = {
	"LChestX", "LChestY", "LChestZ",		/* 30 */
	"LElbX", "LElbY", "LElbZ",					/* 39 */
	"LLArmX",	"LLArmY",	"LLArmZ",				/* 42 */
	"LPinkX",	"LPinkY",	"LPinkZ",	/* 51 */
	"LShoX",	"LShoY",	"LShoZ",				/* 33 */
	"LThumbX",	"LThumbY",	"LThumbZ",	/* 48 */
	"LUArmX",	"LUArmY",	"LUArmZ",				/* 36 */
	"LWristX",	"LWristY",	"LWristZ",	/* 45 */	

	"MChestX",	"MChestY",	"MChestZ",	/* 27 */

	"RChestX",	"RChestY",	"RChestZ",	/* 24 */
	"RElbX",	"RElbY",	"RElbZ",				/* 15 */
	"RLArmX",	"RLArmY",	"RLArmZ",				/* 12 */
	"RPinkX",	"RPinkY",	"RPinkZ",	/* 3 */
	"RShoX",	"RShoY",	"RShoZ",				/* 21 */	
	"RThumbX",	"RThumbY",	"RThumbZ",	/* 6 */	
	"RUArmX",	"RUArmY",	"RUArmZ",				/* 18 */	
	"RWristX", "RWristY", "RWristZ" 		/* 9 */
	};
	int ctr = countEmptyHeadColsCSV(inLine);

	std::vector<std::string> foundMarkers = 
		getRecordedMarkers(inLine);

	// write FILE //
	fstream file;
	file.open(inLine);
	
	for (int i = 0; i < ctr; i++) 
		file << " ,";

	for (int i = 0; i < foundMarkers.size(); i++) {
		file << foundMarkers[i] << "X" << " ,"
				 << foundMarkers[i] << "Y" << " ," 
				 << foundMarkers[i] << "Z" ;		
		if (i < foundMarkers.size()-1) file << " ,";
	}
	file << "\n";

	file.close();

}


int countEmptyHeadColsCSV (string inLine) {


	/* READ FILE */
	string headerLine = getHeaderRowCSV(inLine);
	int commaDelim = headerLine.find(',');
	int endRowDelim = headerLine.find("LChest");

	//NO HEADER cols
	int ctr = 0;
	while((commaDelim >= 0) && (commaDelim < endRowDelim)) {

		headerLine = headerLine.substr(commaDelim + 1);
		commaDelim = headerLine.find(',');
		//start of the marker pos recordings
		endRowDelim = headerLine.find("LChest");
		ctr++;
	}
	return ctr;
}

std::vector<std::string> getRecordedMarkers(string inLine) {

	std::vector<std::string> retVecMarkers;
	/* READ FILE */
	string headerLine = getHeaderRowCSV(inLine);

	// Left side markers
	int startRowDelim = headerLine.find("LChest");
	int endRowDelim = headerLine.find("MChest");
	while ((startRowDelim >=0) && 
		(startRowDelim < endRowDelim)) {
		headerLine = headerLine.substr(startRowDelim);

		int markerLen = headerLine.find(',');		
		string markerName = headerLine.substr(0,markerLen);
		retVecMarkers.push_back(markerName);

		startRowDelim =	headerLine.substr(markerLen).find("L") 
			+ markerLen;
		endRowDelim =  headerLine.find("MChest");
	}
	// right side markers
	retVecMarkers.push_back("MChest");
	startRowDelim =	headerLine.find("R");

	while ((startRowDelim >=0) && (endRowDelim >=0)
		&& (retVecMarkers.size() < NumMarkers/3)) {
		headerLine = headerLine.substr(startRowDelim);

		int markerLen = headerLine.find(',');		
		string markerName = headerLine.substr(0,markerLen);
		retVecMarkers.push_back(markerName);
		endRowDelim = headerLine.substr(markerLen).find("R");
		startRowDelim =	endRowDelim	+ markerLen;
	}
	return retVecMarkers;
}



