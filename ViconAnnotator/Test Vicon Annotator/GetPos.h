#pragma once
#ifndef GETPOS_H
#define GETPOS_H

#include <deque>
#include <fstream>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <map>
#include <list>

using namespace std;

class GetPos {
	typedef deque <deque <double> > ::iterator record_iterator;
	typedef deque        <double>   ::iterator field_iterator;
	typedef map<int, string> markers;


public:
	struct Trajectory;
	struct Marker;
	int numFrames;
	const Trajectory pos(const string& filename, const int& Frame);
	map<int, string> column(const string& filename);
	vector<Marker> MarkerList;
	float* g_vertex_buffer_data;
	void parsedata(const string& filename);
	void parsedata2(const string& filename, int numFrames);
	const vector<double> get(const string& filename, const int& Frame);

};

struct GetPos::Trajectory { // struct for coordinate x,y,z of the marker
	string x;
	string y;
	string z;
};

struct GetPos::Marker { // struct for coordinate x,y,z of the marker
	deque<string> x, y, z, TS;
	string MarkerName;
};



#endif
