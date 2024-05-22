#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/Binning.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/Formalism.h"
#include "AmpGen/Generator.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PhaseCorrection.h"
#include "AmpGen/ProfileClock.h"

#if ENABLE_AVX
  using EventList_type = AmpGen::EventListSIMD;
#else
  using EventList_type = AmpGen::EventList; 
#endif

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>


#include "TRandom3.h"


using namespace AmpGen;
using RealVec_t = std::vector<real_t>;

template <typename T> void createZeros(std::vector<T>& list, const size_t& nItems){
  for (size_t i{0}; i < nItems;i++){
    list.push_back(0);
  }
  return;
}


int main(int argc, char * argv[]){

	//******* PARSE THE ARGS FROM COMMAND LINE AND THE OPTIONS FILE*********
	OptionsParser::setArgs( argc, argv );

	EventType signalType{NamedParameter<std::string>("EventType", "", "Signal Type to generate, D0 -> Ks0pi-pi+")}; // notes -> MINUS COMES FIRST

	const size_t      nThreads = NamedParameter<size_t>("nCores", 8, "Number of threads to use" );
	#ifdef _OPENMP
	omp_set_num_threads( nThreads );
	INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
	omp_set_dynamic( 0 );
	#endif

	const std::string binningFile = NamedParameter<std::string>("Binning", "KsPiPi_equal.txt", "text file with list of binning scheme points");
	
	const size_t      blocksize = NamedParameter<size_t>("Blocksize", 100000, "Number of events to loop through at a time" );

	const std::string modelFile = NamedParameter<std::string>("Model", "Kspipi.opt", "text file with plain D -> f model");

	const std::string outputFilename = NamedParameter<std::string>("Output", "output", "text file to print the output values - NO .txt needed");
	
	const size_t      nBins = NamedParameter<size_t>("nBins", 8, "Number of bins in scheme, default 8" );

	const size_t      nInt = NamedParameter<size_t>("nInt", 1000000, "Total number of events to generate" );

	const size_t      seed = NamedParameter<size_t>("Seed", 0, "First seed to use in generation, dfault 0" );



//*************** READ IN THE BINNING SCHEME **************
	INFO("Reading in binning scheme");
	std::vector<real_t> mMinus;
	std::vector<real_t> mPlus;
	std::vector<int> bins;

	readBinning(binningFile, mMinus, mPlus, bins); //pushes back the lists from text file into these vectors


// *********** SET UP THE MPS  /  A ***********
	MinuitParameterSet MPS_Kspipi;
	MPS_Kspipi.loadFromFile(modelFile);
	AddCPConjugate(MPS_Kspipi);

	CoherentSum A(signalType, MPS_Kspipi);
	CoherentSum Abar(signalType.conj(true), MPS_Kspipi);
	A.prepare(); Abar.prepare();

// ********** SET UP THE PHASE CORRECTION **********
	MinuitParameterSet MPS_PC;
	MPS_PC.loadFromStream();

	PhaseCorrection phaseCorrection{MPS_PC};
	phaseCorrection.compilePolynomialExpressions(signalType);


// loop over nInt / blocksize times through generating and summing blocksize events
	size_t nLoops = nInt / blocksize;
	for(size_t n{0}; n < nLoops; n++){
	ProfileClock TIME_totalLoop; TIME_totalLoop.start();

	// ************** GENERATE IN EVENTS ***************
		INFO("Generating " << blocksize << " events, " << n+1 << "th time of " << nLoops);

		TRandom3 rndm(seed+n); 
		EventList_type events{ Generator<>(signalType, &rndm).generate(nInt) };

		auto evalA = A.amplitudeEvaluator( &events ); 
		auto evalAbar = Abar.amplitudeEvaluator( &events ); 


	//****** PREPARE LISTS OF EACH VALUE WE WANT *******
		// I know it looks like many variables, but I think is easier having names than lists
		// So have pos/neg bins, and 4 variables summing over
		RealVec_t sumA_posBins, sumAbar_posBins, sumAAbarReal_posBins, sumAAbarImag_posBins;
		createZeros(sumA_posBins, nBins+1);	createZeros(sumAbar_posBins, nBins+1);	createZeros(sumAAbarReal_posBins, nBins+1);	createZeros(sumAAbarImag_posBins, nBins+1);
		RealVec_t sumA_negBins, sumAbar_negBins, sumAAbarReal_negBins, sumAAbarImag_negBins;
		createZeros(sumA_negBins, nBins+1);	createZeros(sumAbar_negBins, nBins+1);	createZeros(sumAAbarReal_negBins, nBins+1);	createZeros(sumAAbarImag_negBins, nBins+1);

		

	//******* PREPARE LAMBDAS TO FILL THE THINGS ******
		auto fillSums = [&evalA, &evalAbar](Event& event, int& binNumber, real_t& correction, RealVec_t& FSums, RealVec_t& FbarSums, RealVec_t& cSums, RealVec_t& sSums)
		{
			complex_t complexA = evalA(event);
			complex_t complexAbar = evalAbar(event);


			real_t modA{ std::abs(complexA) };
			real_t modAbar{ std::abs(complexAbar) };
				//convert CP convention
	double temp = std::arg(complexA * std::conj(complexAbar))-M_PI;
	while (temp < -M_PI)
	{
		temp += 2.0 * M_PI;
	}
	while (temp > M_PI)
	{
		temp -= 2.0 * M_PI;
	}
			real_t totalPhase { temp + correction};

			FSums[binNumber] += std::pow(modA, 2);
			FbarSums[binNumber] += std::pow(modAbar, 2);
			cSums[binNumber] += modA * modAbar * cos(totalPhase);
			sSums[binNumber] += modA * modAbar * sin(totalPhase);
			return;
		};



	// ******** FILL THE SUMS ********
		ProfileClock CLOCK_summing;
		CLOCK_summing.start();
		// #pragma omp parallel for // think this causes problems with the pointers into the tree
		INFO("looping through events");
		for (auto event : events){
			int bin = nearestBinIndex(event, mMinus, mPlus, bins);
			real_t deltaC{ 0 };
			if (NamedParameter<size_t>("PhaseCorrection::Order", 0, "Order of polynomial of phase correction") != 0){
				deltaC = phaseCorrection.eval(event);
			}else{
				deltaC = phaseCorrection.evalBias(event);
			}

			if (bin < 0){
				int unsignedBin{ abs(bin) };
				fillSums(event, unsignedBin, deltaC, sumA_negBins, sumAbar_negBins, sumAAbarReal_negBins, sumAAbarImag_negBins);
			}else{
				fillSums(event, bin, deltaC, sumA_posBins, sumAbar_posBins, sumAAbarReal_posBins, sumAAbarImag_posBins);
			}
		}
		
		CLOCK_summing.stop();
		INFO("Took " << CLOCK_summing << "ms to sum over the events");


	// ********* WRITE THE RESULTS TO FILE **********

		std::ofstream outFile{ outputFilename+std::to_string(n)+".txt" };
		if ( !outFile.good()) {
			ERROR("Unable to open output file to print result");
			return 1;
		}

		INFO("Writing results to options file " << outputFilename);
		// default for options output is Fs by counting number of events and ci/si from summing amplitudes only
		for (size_t i{1}; i < nBins+1 ; i++){ 
				outFile 	<< "sumA_ " << i << "     Fix     " << sumA_posBins[i] << "     0\n"
							<< "sumAbar_ " << i << "  Fix     " << sumAbar_posBins[i] << "     0\n"
							<< "sumA_- " << i << "     Fix     " << sumA_negBins[i] << "     0\n"
							<< "sumAbar_- " << i << "  Fix     " << sumAbar_negBins[i] << "     0\n"
							<< "sumCos_ " << i << "     Fix     " << sumAAbarReal_posBins[i]<< "     0\n"
							<< "sumSin_ " << i << "     Fix     " << sumAAbarImag_posBins[i] << "     0\n"
							<< "sumCos_- " << i << "    Fix     " << sumAAbarReal_negBins[i] << "     0\n"
							<< "sumSin_- " << i << "    Fix     " << sumAAbarImag_negBins[i] << "     0\n";
		}
		outFile.close();

		TIME_totalLoop.stop();
		INFO("Took " << TIME_totalLoop << "[ms] to do " << n << "th iteration");

	}

	if (nInt != nLoops*blocksize){
		WARNING("Summed over " << nLoops*blocksize << " events rather than " << nInt);
	}

	return 0;
}