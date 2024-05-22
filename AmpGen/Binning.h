#ifndef AMPGEN_BINNING_H
#define AMPGEN_BINNING_H

#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/Types.h"

#include <string>
#include <vector>

namespace AmpGen
{
	std::string binning();


	void readBinning(const std::string& binning_file, std::vector<real_t>& s01, std::vector<real_t>& s02, std::vector<int>& bins);

	int nearestBinIndex(Event& event, const std::vector<real_t>& s01_list, const std::vector<real_t>& s02_list, const std::vector<int>& binList);
}

#endif