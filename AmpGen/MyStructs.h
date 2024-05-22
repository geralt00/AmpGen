#ifndef AMPGEN_MYSTRUCTS_H
#define AMPGEN_MYSTRUCTS_H

#include "AmpGen/NamedParameter.h"

#include <string>
#include <vector>

namespace AmpGen
{	
	// STRUCT FOR GENERATION - so can loop through same code twice to make B- then B+ events
	struct BtypeInfo
	{
		size_t nTypes{2};
		std::vector<int> signs{-1, 1}; 
		std::vector<std::string> prefixes{"Bminus_", "Bplus_"};	
	};

	// STRUCT FOR GENERATION - as above but for generating BES typed data
	struct BESeventInfo
	{
		size_t nTags{5};
		std::vector<std::string> tagNames{"KK", "KsPi0", "KminusPiplus", "KplusPiminus", "Kspipi"};
	};

	// STRUCT FOR BINNED/FITTING - so don't have to use -1/1 to specify which B and bind we are doing
	struct signs
	{
		int BminusSign{-1};
		int BplusSign{1};

		int negBinSign{-1};
		int posBinSign{1};

		//since it's CP odd in tag side, with psi(3770) = -1, we need to flip the sign of the CP even component
		int CPoddSign{1};
		int CPevenSign{-1};

	};

	// STRUCT FOR AMPLITUDE INFO - so can create it once then store without having to use AmplitudeEvaluator or .GetValNoCache etc
	struct amplitudeInfo
	{
		std::vector<real_t> A;
		std::vector<real_t> Abar;
		std::vector<real_t> deltaD;
		std::vector<int> bins;
	};
	struct bkgfraction
	{
	std::vector<double> frac;
	std::vector<std::string> tag;
	};
	struct normeff
	{
	std::vector<double> frac;
	std::vector<std::string> tag;
	};
	
	// struct dalitzPair{
	// 	real_t plus{0};
	// 	real_t minus{0};
	// };

	template <typename T>
	struct dalitzPair{
		T plus;
		T minus;
	};


}

#endif