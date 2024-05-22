#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/Formalism.h"
#include "AmpGen/Generator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PhaseCorrection.h"
#include "AmpGen/MyStructs.h"
#include "AmpGen/WriteToTree.h"
#include "extern/DtoKpipiAmplitude.h"
#include "extern/Belle2010Amplitude.h"
#include "extern/D0ToKLpipi2018.h"

#include "extern/TDalitz.h"
#include "AmpGen/QMI.h"


#include "AmpGen/Particle.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/Generator.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/AddCPConjugate.h"


#ifdef _OPENMP
	#include <omp.h>
#endif
#if ENABLE_AVX
	using EventList_type = AmpGen::EventListSIMD;
#else
	using EventList_type = AmpGen::EventList; 
#endif


#include <string>
#include <vector>

#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"



using namespace AmpGen;


// TRandom 3 unaccepted input to the setRandom function so trying a separate function to mimic the syntax it gets passed around with in AmpGen.cpp. Not sure why this works and just calling it doesn't. 
void setTRandom(TRandom* rndm, auto& generator){
	generator.setRandom(rndm);
	return;
}


int main(int argc , char* argv[] ){

// ******* ARGPARSING *******
	OptionsParser::setArgs( argc, argv );

	const size_t      blockSize = NamedParameter<size_t>("blockSize", 1e5, "Blocksize for generation");

	EventType signalType{NamedParameter<std::string>("EventType", "", "Signal Type to generate")};

	const std::string modelFile = NamedParameter<std::string>("Model", "Kspipi.opt", "output filename, needs .root");
	
	const size_t      nBins = NamedParameter<size_t>("nBins", 100, "Number of bins for projection histograms");

	const size_t      nEvents = NamedParameter<size_t>("nEvents", 15000, "number of events to generate, multiplies by the BR of each tag");

	const std::string outputFileName = NamedParameter<std::string>("Output", "outputFile.root", "output filename, needs .root");

	const size_t      seed = NamedParameter<size_t>("Seed", 0, "Random seed for generation");
	TRandom3 rndm; rndm.SetSeed( seed );

	const size_t      nThreads = NamedParameter<size_t>("nCores", 8, "Number of threads to use" );
	#ifdef _OPENMP
		omp_set_num_threads( nThreads );
		INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
		omp_set_dynamic( 0 );
	#endif


//***********************************************************
//***************** READ IN THINGS AND SET UP ***************
//***********************************************************

// ***** SET UP THE PHASE CORRECTION *****
	MinuitParameterSet MPS_phaseCorrection; 
	MPS_phaseCorrection.loadFromStream();	// base mps with the phase correction info

	PhaseCorrection deltaCorrection{MPS_phaseCorrection}; 
	deltaCorrection.compilePolynomialExpressions(signalType);

// ***** SET UP SIGNAL (Kspipi) AMPLITUDE *****
	MinuitParameterSet MPS_Kspipi;
	MPS_Kspipi.loadFromFile(modelFile);
	AddCPConjugate(MPS_Kspipi);

	CoherentSum A(signalType, MPS_Kspipi);
	CoherentSum Abar(signalType.conj(true), MPS_Kspipi);
	A.prepare(); Abar.prepare();



//************************************************************************************
// ************** GO THROUGH EACH OF THE 3 GROUPED TAG TYPES AND GENERATE ************
//************************************************************************************
//Shenghui taged BELLE2010
    std::string kaon_type = "KS";

    std::cout << "Running KS-amplitude calculator\n";
    DtoKpipiAmplitude* amp;
    amp = new Belle2010Amplitude(kaon_type, "/publicfs/ucas/user/zengshh/LHCb/B2DK_D2K0S0PiPi/analysis/ampgen/extern/Belle2010/bin/cp_mult.txt");

	std::vector<Expression> Phi = QMI::dalitz(signalType);
    std::vector<CompiledExpression<real_t(const real_t*, const real_t*)> > cPhi;
    for (int i=0;i<Phi.size();i++){
            CompiledExpression<real_t(const real_t*, const real_t*)> cPhi_i(Phi[i], "Phi_" + std::to_string(i) , signalType.getEventFormat(), &MPS_Kspipi);
            cPhi_i.prepare(); cPhi_i.compile();
            cPhi.push_back(cPhi_i);
    }

	signs signs{};
	Generator<PhaseSpace> generator(signalType);
	setTRandom(&rndm, generator); // needed instead of straight using .setRandom - see explanation at function
	generator.setBlockSize(blockSize);
	generator.setNormFlag(true);

	real_t thisEventFraction;
	size_t thisNEvents;

//**************  CP:
	int CPsign;
	auto totalA_CP = [&amp, &cPhi, &deltaCorrection, &CPsign](Event event){
		return totalAmplitudeSquared_CP_test(CPsign, *amp, cPhi, deltaCorrection, event);
	};

	thisEventFraction = NamedParameter<real_t>("BEStag::CPeven:nEvents", 0.1, "fraction of nEvents for this tag");
	thisNEvents = thisEventFraction * nEvents;
	CPsign = signs.CPevenSign; // this should be a plus, even referring to the tag state
	EventList acceptedEvents_CPeven{signalType};
	generator.fillEventsUsingLambda(totalA_CP, acceptedEvents_CPeven, thisNEvents);

	thisEventFraction = NamedParameter<real_t>("BEStag::CPodd:nEvents", 0.1, "fraction of nEvents for this tag");
	thisNEvents = thisEventFraction * nEvents;
	CPsign = signs.CPoddSign; // should be minus, the tag is the odd state
	EventList acceptedEvents_CPodd{signalType};
	generator.fillEventsUsingLambda(totalA_CP, acceptedEvents_CPodd, thisNEvents);

//************** FLAVOUR:
	bool isFlavour; // true if for flavour, false if for flavourbar
	auto totalA_flavour = [&isFlavour, &amp, &cPhi](Event event){
		return totalAmplitudeSquared_flavour_test(isFlavour, *amp, cPhi, event);
	};

	thisEventFraction = NamedParameter<real_t>("BEStag::flavour:nEvents", 0.1, "fraction of nEvents for this tag");
	thisNEvents = thisEventFraction * nEvents;
	isFlavour = true;
	EventList acceptedEvents_flavour{signalType};
	generator.fillEventsUsingLambda(totalA_flavour, acceptedEvents_flavour, thisNEvents);

	thisEventFraction = NamedParameter<real_t>("BEStag::flavourBar:nEvents", 0.1, "fraction of nEvents for this tag");
	thisNEvents = thisEventFraction * nEvents;
	isFlavour = false;
	EventList acceptedEvents_flavourBar{signalType};
	generator.fillEventsUsingLambda(totalA_flavour, acceptedEvents_flavourBar, thisNEvents);

//************* Mixed (Kspipi tag):
	auto totalA_same = [&amp, &cPhi, &deltaCorrection](Event event_main, Event event_tag){
		return totalAmplitudeSquared_BES_test(*amp, cPhi, deltaCorrection, event_main, event_tag);
	};

	thisEventFraction = NamedParameter<real_t>("BEStag::Kspipi:nEvents", 0.1, "fraction of nEvents for this tag");
	thisNEvents = thisEventFraction * nEvents;
	EventList acceptedEvents_same_main{signalType};
	EventList acceptedEvents_same_tag{signalType};
	generator.fill2EventsUsingLambda(totalA_same, acceptedEvents_same_main, acceptedEvents_same_tag, thisNEvents);

//************* KlLpipi tag:

//	D0ToKLpipi2018 tKLpipi;tKLpipi.init();
//		auto totalA_mixed = [&amp, &tKLpipi, &cPhi, &deltaCorrection](Event event_main, Event event_tag){
//		return totalAmplitudeSquared_BES_KLpipi_test(*amp, tKLpipi, cPhi, deltaCorrection, event_main, event_tag);
//};

//	thisEventFraction = NamedParameter<real_t>("BEStag::KLpipi:nEvents", 0.1, "fraction of nEvents for this tag");
//	thisNEvents = thisEventFraction * nEvents;
//	EventList acceptedEvents_klpipi_main{signalType};
//	EventList acceptedEvents_klpipi_tag{signalType};
//	generator.fill2EventsUsingLambda(totalA_mixed, acceptedEvents_klpipi_main, acceptedEvents_klpipi_tag, thisNEvents);





//***********************************************************
// ************* WRITE IT OUT TO TREE ***********************
//***********************************************************
	TFile * outputFile = TFile::Open(outputFileName.c_str(), "RECREATE"); 
	outputFile->cd();
        auto my_dd = [&amp](Event& evt){
		dcomplex A = amp->get_amplitude(evt.s(0,2),evt.s(0,1));
		dcomplex Abar = amp->get_amplitude(evt.s(0,1), evt.s(0,2));
		
                return  (arg(A*conj(Abar))); // changed to make there ony one delta_D definition 19/9/23
        };
		auto my_A = [&amp](Event& evt){
		dcomplex A = amp->get_amplitude(evt.s(0,2),evt.s(0,1));
		dcomplex Abar = amp->get_amplitude(evt.s(0,1), evt.s(0,2));
		
                return  norm(A); // changed to make there ony one delta_D definition 19/9/23
    };
	auto writeAndPlot = [&nBins](EventList& acceptedEvents, std::string& label){
		INFO("Writing Dalitz Event List to tree for " << label << " events");
		acceptedEvents.tree((label + "__DalitzEventList").c_str())->Write();

		if (NamedParameter<bool>("DoPlots", true)){
			INFO("Making and writing plots for " << label << " too");
			writePlots(acceptedEvents, nBins, label);
		}
		return;
	};

	std::string labels[8] = {"CPeven", "CPodd", "flavour", "flavourBar", "same_signal", "same_tag"};
	writeAndPlot(acceptedEvents_CPeven, labels[0]);
	writeAndPlot(acceptedEvents_CPodd, labels[1]);
	writeAndPlot(acceptedEvents_flavour, labels[2]);
	writeAndPlot(acceptedEvents_flavourBar, labels[3]);
	writeAndPlot(acceptedEvents_same_main, labels[4]);
	writeAndPlot(acceptedEvents_same_tag, labels[5]);
	EventList_type flatEvents_CPeven = generator.generate(1000000);
	writeDalitzVariables(flatEvents_CPeven, "");
	QMI::writeValues(flatEvents_CPeven, my_dd, "dd");
	QMI::writeValues(flatEvents_CPeven, my_A, "amp_A");	




 outputFile->Close();


	return 0;
}
