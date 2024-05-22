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
#include "AmpGen/PhaseCorrection_KLpipi.h"
#include "AmpGen/MyStructs.h"
#include "AmpGen/WriteToTree.h"
#include "extern/D0ToKLpipi2018.h"
#include "extern/TDalitz.h"
#include "AmpGen/QMI.h"
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif
#if ENABLE_AVX
using EventList_type = AmpGen::EventListSIMD;
#else
using EventList_type = AmpGen::EventList;
#endif
#include <map>
#include <string>
#include <vector>

#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"

using namespace AmpGen;

int distribution(double mean) {
    std::random_device rd;
    std::mt19937 generator(rd());
    std::poisson_distribution<int> poisson(mean);
    return poisson(generator);
}

// TRandom 3 unaccepted input to the setRandom function so trying a separate function to mimic the syntax it gets passed around with in AmpGen.cpp. Not sure why this works and just calling it doesn't.
void setTRandom(TRandom *rndm, auto &generator)
{
	generator.setRandom(rndm);
	return;
}

int main(int argc, char *argv[])
{

	// ******* ARGPARSING *******
	OptionsParser::setArgs(argc, argv);

	const size_t blockSize = NamedParameter<size_t>("blockSize", 1e5, "Blocksize for generation");

	EventType signalType{NamedParameter<std::string>("EventType", "", "Signal Type to generate")};

	const std::string modelFile = NamedParameter<std::string>("Model", "Kspipi.opt", "output filename, needs .root");

	const size_t nBins = NamedParameter<size_t>("nBins", 100, "Number of bins for projection histograms");

	const size_t nEvents = NamedParameter<size_t>("nEvents", 15000, "number of events to generate, multiplies by the BR of each tag");

	const std::string outputFileName = NamedParameter<std::string>("Output", "outputFile.root", "output filename, needs .root");

	const size_t seed = NamedParameter<size_t>("Seed", 0, "Random seed for generation");

	const bool is_poisson = NamedParameter<bool>("is_poisson", false, "use poisson distribution for number of events");

	TRandom3 rndm;
	rndm.SetSeed(seed);

	const size_t nThreads = NamedParameter<size_t>("nCores", 8, "Number of threads to use");

	const bool plot_phase = NamedParameter<bool>("plot_phase", false, "plot phase of the amplitude");
#ifdef _OPENMP
	omp_set_num_threads(nThreads);
	INFO("Setting " << nThreads << " fixed threads for OpenMP");
	omp_set_dynamic(0);
#endif

	//***********************************************************
	//***************** READ IN THINGS AND SET UP ***************
	//***********************************************************

	// ***** SET UP THE PHASE CORRECTION *****
	MinuitParameterSet MPS_phaseCorrection;
	MPS_phaseCorrection.loadFromStream(); // base mps with the phase correction info

	PhaseCorrection deltaCorrection{MPS_phaseCorrection};
	deltaCorrection.compilePolynomialExpressions(signalType);
	PhaseCorrection_KLpipi deltaCorrection_KLpipi{MPS_phaseCorrection};
	deltaCorrection_KLpipi.compilePolynomialExpressions(signalType);
	// ***** SET UP SIGNAL (Kspipi) AMPLITUDE *****
	MinuitParameterSet MPS_Kspipi;
	MPS_Kspipi.loadFromFile(modelFile);
	AddCPConjugate(MPS_Kspipi);

	CoherentSum A(signalType, MPS_Kspipi);
	CoherentSum Abar(signalType.conj(true), MPS_Kspipi);
	A.prepare();
	Abar.prepare();

	// SET UP KLpipi model(2018)
	std::vector<Expression> Phi = QMI::dalitz(signalType);
	std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> cPhi;
	for (int i = 0; i < Phi.size(); i++)
	{
		CompiledExpression<real_t(const real_t *, const real_t *)> cPhi_i(Phi[i], "Phi_" + std::to_string(i), signalType.getEventFormat(), &MPS_Kspipi);
		cPhi_i.prepare();
		cPhi_i.compile();
		cPhi.push_back(cPhi_i);
	}
	D0ToKLpipi2018 tKLpipi;
	tKLpipi.init();
	//************************************************************************************
	// ************** GO THROUGH EACH OF THE 3 GROUPED TAG TYPES AND GENERATE ************
	//************************************************************************************
	signs signs{};
	Generator<PhaseSpace> generator(signalType);
	setTRandom(&rndm, generator); // needed instead of straight using .setRandom - see explanation at function
	generator.setBlockSize(blockSize);
	generator.setNormFlag(true);

	real_t thisEventFraction;
	size_t thisNEvents;

	//**************  CP:
	int CPsign;
	auto totalA_CP = [&A, &Abar, &deltaCorrection, &CPsign](Event event)
	{
		return totalAmplitudeSquared_CP(CPsign, A, Abar, deltaCorrection, event);
	};

	thisEventFraction = NamedParameter<real_t>("BEStag::CPeven:nEvents", 0.1, "fraction of nEvents for this tag");
	thisNEvents = thisEventFraction * nEvents;
	if (is_poisson)
	{
		thisNEvents = distribution(thisNEvents);
	}
	CPsign = signs.CPevenSign; // this should be a plus, even referring to the tag state
	EventList acceptedEvents_CPeven{signalType};
	generator.fillEventsUsingLambda(totalA_CP, acceptedEvents_CPeven, thisNEvents);

	thisEventFraction = NamedParameter<real_t>("BEStag::CPodd:nEvents", 0.1, "fraction of nEvents for this tag");
	thisNEvents = thisEventFraction * nEvents;
	if (is_poisson)
	{
		thisNEvents = distribution(thisNEvents);
	}
	CPsign = signs.CPoddSign; // should be minus, the tag is the odd state
	EventList acceptedEvents_CPodd{signalType};
	generator.fillEventsUsingLambda(totalA_CP, acceptedEvents_CPodd, thisNEvents);

	//************** FLAVOUR:
	bool isFlavour; // true if for flavour, false if for flavourbar
	auto totalA_flavour = [&isFlavour, &A, &Abar](Event event)
	{
		return totalAmplitudeSquared_flavour(isFlavour, A, Abar, event);
	};

	thisEventFraction = NamedParameter<real_t>("BEStag::flavour:nEvents", 0.1, "fraction of nEvents for this tag");
	thisNEvents = thisEventFraction * nEvents;
	if (is_poisson)
	{
		thisNEvents = distribution(thisNEvents);
	}
	isFlavour = true;
	EventList acceptedEvents_flavour{signalType};
	generator.fillEventsUsingLambda(totalA_flavour, acceptedEvents_flavour, thisNEvents);

	thisEventFraction = NamedParameter<real_t>("BEStag::flavourBar:nEvents", 0.1, "fraction of nEvents for this tag");
	thisNEvents = thisEventFraction * nEvents;
	if (is_poisson)
	{
		thisNEvents = distribution(thisNEvents);
	}
	isFlavour = false;
	EventList acceptedEvents_flavourBar{signalType};
	generator.fillEventsUsingLambda(totalA_flavour, acceptedEvents_flavourBar, thisNEvents);

	//************* SAME (Kspipi tag):
	auto totalA_same = [&A, &Abar, &deltaCorrection](Event event_main, Event event_tag)
	{
		return totalAmplitudeSquared_BES(A, Abar, deltaCorrection, event_main, event_tag);
	};

	thisEventFraction = NamedParameter<real_t>("BEStag::Kspipi:nEvents", 0.1, "fraction of nEvents for this tag");
	thisNEvents = thisEventFraction * nEvents;
	if (is_poisson)
	{
		thisNEvents = distribution(thisNEvents);
	}
	EventList acceptedEvents_same_main{signalType};
	EventList acceptedEvents_same_tag{signalType};
	generator.fill2EventsUsingLambda(totalA_same, acceptedEvents_same_main, acceptedEvents_same_tag, thisNEvents);

	/*
		auto totalA_Kspipi_KLpipi = [&A, &Abar, &tKLpipi, &cPhi, &deltaCorrection, &deltaCorrection_KLpipi](Event event_main, Event event_tag){
			return totalAmplitudeSquared_BES_KSpipi_KLpipi(A, Abar, tKLpipi, cPhi, deltaCorrection, deltaCorrection_KLpipi, event_main, event_tag);
		};

		thisEventFraction = NamedParameter<real_t>("BEStag::KLpipi:nEvents", 0.1, "fraction of nEvents for this tag");
		thisNEvents = thisEventFraction * nEvents;
		EventList acceptedEvents_Kspipi_KLpipi_main{signalType};
		EventList acceptedEvents_Kspipi_KLpipi_tag{signalType};
		generator.fill2EventsUsingLambda(totalA_Kspipi_KLpipi, acceptedEvents_Kspipi_KLpipi_main, acceptedEvents_Kspipi_KLpipi_tag, thisNEvents);
	*/

	//***********************************************************
	// ************* WRITE IT OUT TO TREE ***********************
	//***********************************************************

	TFile *outputFile = TFile::Open(outputFileName.c_str(), "RECREATE");
	outputFile->cd();

	auto writeAndPlot = [&nBins](EventList &acceptedEvents, std::string &label)
	{
		INFO("Writing Dalitz Event List to tree for " << label << " events");
		acceptedEvents.tree((label + "__DalitzEventList").c_str())->Write();

		if (NamedParameter<bool>("DoPlots", true))
		{
			INFO("Making and writing plots for " << label << " too");
			writePlots(acceptedEvents, nBins, label);
		}
		return;
	};

	std::string labels[8] = {"CPeven", "CPodd", "flavour", "flavourBar", "same_signal", "same_tag", "KS", "KL"};
	writeAndPlot(acceptedEvents_CPeven, labels[0]);
	writeAndPlot(acceptedEvents_CPodd, labels[1]);
	writeAndPlot(acceptedEvents_flavour, labels[2]);
	writeAndPlot(acceptedEvents_flavourBar, labels[3]);
	writeAndPlot(acceptedEvents_same_main, labels[4]);
	writeAndPlot(acceptedEvents_same_tag, labels[5]);
	/*
	writeAndPlot(acceptedEvents_Kspipi_KLpipi_main, labels[6]);
	writeAndPlot(acceptedEvents_Kspipi_KLpipi_tag, labels[7]);*/
	if (plot_phase)
	{
		auto my_dd = [&A, &Abar](Event &evt)
		{
			//return std::arg(A.getValNoCache(evt) * std::conj(Abar.getValNoCache(evt))) ;
			//return arg(A.getValNoCache(evt)*conj(-Abar.getValNoCache(evt)));//??? why is this negative
			return Deltadelta(evt, A, Abar);
		};
		auto amp_A = [&A, &Abar](Event &evt)
		{
			return norm(A.getValNoCache(evt));
		};
		auto my_bias = [&deltaCorrection](Event &event)
		{
			return deltaCorrection.evalBias(event);
		};
		EventList_type flatEvents = generator.generate(1000000);
		writeDalitzVariables(flatEvents, "");
		QMI::writeValues(flatEvents, my_bias, "dd_bias");
		QMI::writeValues(flatEvents, my_dd, "dd");
		QMI::writeValues(flatEvents, amp_A, "amp_A");
	}

	outputFile->Close();

	return 0;
}
