#include "AmpGen/Binning.h"

#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/Types.h"

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include <iostream>

using namespace AmpGen;

//test function
std::string AmpGen::binning(){
	return "message to print ('')> ('')>";
}


// read the binning scheme from a .txt file
// put the informaiton into three vectors: m+ coord, m- coord and the bin
void AmpGen::readBinning(const std::string& binningFile, std::vector<real_t>& s01, std::vector<real_t>& s02, std::vector<int>& bins){

	std::ifstream in_File{binningFile};

	if (!in_File.good()){
		ERROR("Unable to open binning file");
		return;
	}

	real_t first{0}, second{0}, third{0}; // variables to put each number into to add to the vectors

	while( in_File >> first >> second >> third ){
		s01.push_back(first);
		s02.push_back(second);
		bins.push_back(int(third));
	}
	in_File.close();
	return;
}


// find the bin of an event
// taken from Jake's function
// by finding the nearest dalitz coordinate in the list from the text file
int AmpGen::nearestBinIndex(Event& event, const std::vector<real_t>& s01_list, const std::vector<real_t>& s02_list, const std::vector<int>& binList){
	// distance = distance ^2 on the dalitz plane between event and each point in the binning scheme lists

	real_t s01 = event.s(0, 1); // aka m-
	real_t s02 = event.s(0, 2); // aka m+

	int minIndex{0};
	real_t minDistance{100};

	// #pragma omp parallel for // doesn't actually seem to speed things up
	for (int i=0;i<s01_list.size();i++){
		real_t distance = std::pow(s01_list[i] - s01, 2) + std::pow(s02_list[i] - s02, 2);

		if (distance < minDistance) {
			minIndex = i;
			minDistance = distance;
		}
	}

	int bin{binList[minIndex]};
	if (s02 < s01){
		bin = -bin;
	}
		
	return bin;
	}
