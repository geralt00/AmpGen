#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/Binning.h"
#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/EventType.h"
#include "AmpGen/Formalism.h"
#include "AmpGen/Generator.h"
#include "AmpGen/MyStructs.h"
#include "AmpGen/PhaseCorrection_KLpipi.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PhaseCorrection.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/QMI.h"

#if ENABLE_AVX
	using EventList_type = AmpGen::EventListSIMD;
#else
	using EventList_type = AmpGen::EventList; 
#endif

#include <fstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TRandom3.h"

using namespace AmpGen;


int main(int argc, char * argv[]){

	OptionsParser::setArgs( argc, argv );

	const std::string binningFile = NamedParameter<std::string>("Binning", "KsPiPi_equal.txt", "text file with list of binning scheme points");

	const std::string modelFile = NamedParameter<std::string>("Model", "Kspipi.opt", "Options file with the amplitude model in it");

	const size_t      nInt = NamedParameter<size_t>("nInt", 1e7, "Number of events to calculate normalisation - should be large");

	const size_t      seed = NamedParameter<size_t>("Seed", 0, "Random seed for generation");
	TRandom3 rndm(seed); 

	const std::string tupleFilename = NamedParameter<std::string>("Tuple", "tuple.root", "name of root file to write data to");

	const std::string textFilename = NamedParameter<std::string>("Amplitudes", "amplitudes.txt", "name of text file to write amplitudes to");

	EventType signalType{NamedParameter<std::string>("EventType", "", "Signal Type to generate, i.e. D0 K0S0 pi- pi+")}; // notes -> MINUS COMES FIRST


	const size_t      nThreads = NamedParameter<size_t>("nCores", 8, "Number of threads to use" );
	#ifdef _OPENMP
	omp_set_num_threads( nThreads );
	INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
	omp_set_dynamic( 0 );
	#endif

	// make events
	ProfileClock TIME_MCgeneration; TIME_MCgeneration.start();
	EventList_type events =  Generator<>(signalType, &rndm).generate(nInt);
	TIME_MCgeneration.stop();
	INFO("Took " << TIME_MCgeneration << "ms to generate " << nInt << " flat phase space events");
	
	MinuitParameterSet MPS_Kspipi;
	MPS_Kspipi.loadFromFile(modelFile);
	AddCPConjugate(MPS_Kspipi);


// SET UP KLpipi model(2018)
	std::vector<Expression> Phi = QMI::dalitz(signalType);
    std::vector<CompiledExpression<real_t(const real_t*, const real_t*)> > cPhi;
    for (int i=0;i<Phi.size();i++){
            CompiledExpression<real_t(const real_t*, const real_t*)> cPhi_i(Phi[i], "Phi_" + std::to_string(i) , signalType.getEventFormat(), &MPS_Kspipi);
            cPhi_i.prepare(); cPhi_i.compile();
            cPhi.push_back(cPhi_i);
    }
	D0ToKLpipi2018 tKLpipi;tKLpipi.init();

	amplitudeInfo amplitudes = fillAmplitudeInfo_KLpipi(events, tKLpipi, cPhi);


	// Read in binning scheme
	ProfileClock TIME_binningReading;
	TIME_binningReading.start();
	std::vector<real_t> mMinus;
	std::vector<real_t> mPlus;
	std::vector<int> bins;

	const size_t nBins = NamedParameter<size_t>("nBins", 8, "Number of bins in binning file, assumed 8" );

	readBinning(binningFile, mMinus, mPlus, bins);
	TIME_binningReading.stop();
	INFO("Took " << TIME_binningReading << "ms to read in the binning scheme");
	

	// bin events -> save to struct
	ProfileClock TIME_doBinning;
	TIME_doBinning.start();
	std::vector<int> binList; // each will have a bin number for each event
	int eventBin{0};

	for (auto& event:events){
		binList.push_back(nearestBinIndex(event, mMinus, mPlus, bins));
	}
	
	TIME_doBinning.stop();
	INFO("Took " << TIME_doBinning << "ms to bin the events");
	amplitudes.bins = binList;


	// save events to root file
	INFO("Writing Dalitz Event List to tree");
	TFile* outputFile = TFile::Open(tupleFilename.c_str(), "RECREATE"); 
	outputFile->cd();
	events.tree("DalitzEventList")->Write();
	outputFile->Close();



	// save struct informaiton to text file
	std::ofstream outFile{textFilename};
	for (size_t i{0}; i < nInt; i++){
		outFile << amplitudes.A[i] << "  " << amplitudes.Abar[i] << "  " << amplitudes.deltaD[i] << "  " << amplitudes.bins[i] << std::endl;
	}

	outFile.close();


	return 0;
}