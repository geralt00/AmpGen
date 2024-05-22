// commented while testing phase correction class 6/2

#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/ArgumentPack.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventType.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Formalism.h"
#include "AmpGen/Generator.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MyStructs.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PhaseCorrection.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/Types.h"
#include "AmpGen/WriteToTree.h"

#if ENABLE_AVX
  #include "AmpGen/EventListSIMD.h"
  using EventList_type = AmpGen::EventListSIMD;
#else
  #include "AmpGen/EventList.h"
  using EventList_type = AmpGen::EventList; 
#endif

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TRandom3.h"


using namespace AmpGen;

int main(int argc, char * argv[]){


  //******* PARSE THE ARGS FROM COMMAND LINE AND THE OPTIONS FILE*********
	OptionsParser::setArgs( argc, argv );

	EventType signalType{NamedParameter<std::string>("EventType", "", "Signal Type to generate, i.e. D0 K0S0 pi- pi+")}; // notes -> (1) MINUS COMES FIRST

	const std::string inputData = NamedParameter<std::string>("Input", "data.root", "Root file containing B->D(->Kspipi)K data for both B+/-");

	const std::string logFile = NamedParameter<std::string>("Log", "Fit.log", "Name of the output log file for fit result");

	const std::string modelFile = NamedParameter<std::string>("Model", "Kspipi.opt", "Options file with the amplitude model in it");

	size_t		      nInt = NamedParameter<size_t>("nInt", 1e7, "Number of events to calculate normalisation - should be large");

	const std::string normalisationAmplitudes = NamedParameter<std::string>("NormalisationAmplitudes", "", "Text file with normalsiation amplitudes");

	const std::string normalisationTuple = NamedParameter<std::string>("NormalisationEvents", "", "Root file with normalsiation events");

	const std::string plotFile = NamedParameter<std::string>("Plots", "Result.root", "Name of the output root file to save fit result plots to");

	const size_t      seed = NamedParameter<size_t>("Seed", 0, "Random seed for generation");
	TRandom3 rndm(seed); 

	const size_t      nThreads = NamedParameter<size_t>("nCores", 8, "Number of threads to use" );
	#ifdef _OPENMP
	omp_set_num_threads( nThreads );
	INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
	omp_set_dynamic( 0 );
	#endif

	//******* LOAD THE DATA ******
	EventList_type events_Bminus((inputData + ":Bminus__DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"), WeightBranch_eff("eff_weight"));
	EventList_type events_Bplus((inputData + ":Bplus__DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"), WeightBranch_eff("eff_weight"));
	std::vector<EventList_type> Datalist_LHCb{events_Bminus, events_Bplus};
	// (gammagen syntax:)
	// EventList_type events_Bplus(inputData + ".root:Bp2DKp" , signalType);
	// EventList_type events_Bminus(inputData + ".root:Bm2DKm" , signalType);
	size_t nEvents { events_Bminus.size() + events_Bplus.size() }; // assume both files are the same size as they should be

	//******** LOAD THE MPS OF FITTING VARIABLES ********
	INFO("loading MPS");
	MinuitParameterSet MPS;
	MPS.loadFromStream();

	//******** CREATE POINTS IN FLAT PHASE SPACE FOR NORMALISATION ******
	// nInt points across space of defined decay and randomly distributed
	
	bool generatedEvents{ false };
	EventList_type eventsMC(signalType);

	if (normalisationTuple.empty()){
		INFO("No normalisation tuple found, generating " << nInt << " integration events");
		eventsMC =  Generator<>(signalType, &rndm).generate(nInt);
		generatedEvents = true;
	}else{
		ProfileClock TIME_readNormTree; TIME_readNormTree.start();
		ArgumentPack args;
		//eventsMC.loadFromFile( normalisationTuple+":Bplus_DalitzEventList" , args );
		EventList_type eventsMC((normalisationTuple + ":Bminus_DalitzEventList").c_str(), signalType, WeightBranch_eff("eff_weight"));
		TIME_readNormTree.stop();
		nInt = eventsMC.size();
	}


	// Prepare the phase correction
	PhaseCorrection phaseCorrection{MPS};
	phaseCorrection.compilePolynomialExpressions(signalType);

	//******** SORT OUT THE AMPLITUDE MODEL *******
	// load in amplitude model A = amplitude(D0 -> Ks0pipi) and Abar from a Dbar0
	INFO("about to load in kspipi model");
	MinuitParameterSet MPS_Kspipi;
	MPS_Kspipi.loadFromFile(modelFile);
	AddCPConjugate(MPS_Kspipi);

	CoherentSum A(signalType, MPS_Kspipi);
	CoherentSum Abar(signalType.conj(true), MPS_Kspipi);

	// connect the amplitude model and the points generated across phase space
	A.setMC(eventsMC);
	A.prepare();
	Abar.setMC(eventsMC);
	Abar.prepare();


	// none of the individual A will change so get the info into a usable format now for speed reasons.
	INFO("Collating amplitude information");
	amplitudeInfo amplitudesMC;
	std::ifstream inFile{normalisationAmplitudes};
	if (inFile.good() && !generatedEvents){ // case of the amplitudes and events are both fine
		INFO("reading in amplitude information from text file, if it doesn't match the tuple read in the fit may not be happy");
		double A_1, Abar_2, Dd_3, bin_4;
		nInt = 0; // need to now count the number of entries in the text file
		while ( inFile >> A_1 >> Abar_2 >> Dd_3 >> bin_4 ){
			amplitudesMC.A.push_back(A_1);
			amplitudesMC.Abar.push_back(Abar_2);
			amplitudesMC.deltaD.push_back(Dd_3);
			nInt++;
		}
	}
	else
	{ // other cases mean recalculating amplitude info, but for different reasons:
		// if regenerated events:
		if(generatedEvents){WARNING("Can't open file " << normalisationAmplitudes << " as have regenerated events. Calculating amplitude information now");}
		// if events were given but amplitudes weren't:
		else if (normalisationAmplitudes.empty()){WARNING("No amplitude file specified with --NormalisationAmplitudes. Calculating amplitude information now");}
		// if we read in the events fine but the amplitude file is bad
		else if (!inFile.good()){WARNING("Can't open file " << normalisationAmplitudes << ". Recalculating amplitude information now");}
		amplitudesMC = fillAmplitudeInfo(eventsMC, A, Abar);
	}

	amplitudeInfo amplitudes_Bminus = fillAmplitudeInfo(events_Bminus, A, Abar);
	amplitudeInfo amplitudes_Bplus = fillAmplitudeInfo(events_Bplus, A, Abar);
	// for normalisation - |A|^2 for later but when efficiency introduce the norm a shoud be weighted by the efficiency
	real_t normA{A.norm()};
	real_t normAbar{Abar.norm()};

	// struct for signs so don't have to use numbers
	signs signs{};
	//******** BUILD THE LOG LIKELIHOOD LAMBDAS **********
	//******get sWeight normalisation


	auto LL_total = [&events_Bminus, &events_Bplus, &eventsMC, 
					 &amplitudesMC, &amplitudes_Bminus, &amplitudes_Bplus, &nEvents,
					 &normA, &normAbar, &MPS,
					 &phaseCorrection, &signs](){

		phaseCorrection.updateCoeffs(MPS); // update the stored coefficients to latest in MPS
		// NORMALISATION:
		std::pair<real_t, real_t> normalisationCrossTerms = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC, eventsMC, phaseCorrection); // {cos term, sin term}
		real_t normalisation_Bplus = totalAmplitudeSquared_Integrated(signs.BplusSign, normA, normAbar, MPS, normalisationCrossTerms);
		real_t normalisation_Bminus = totalAmplitudeSquared_Integrated(signs.BminusSign, normA, normAbar, MPS, normalisationCrossTerms);
		real_t  ll_running{0};

		
		#pragma omp parallel for reduction (+:ll_running)
		// for (auto event : events_Bminus) // had to switch to below to make the omp run -> change back if you find out how?
		for (size_t i=0; i < events_Bminus.size(); i++)
		{
			// Event event{events_Bminus[i]}
			real_t correction{ phaseCorrection.eval(events_Bminus[i]) };
			real_t probability{totalAmplitudeSquared_XY(signs.BminusSign, amplitudes_Bminus.A[i], amplitudes_Bminus.Abar[i], amplitudes_Bminus.deltaD[i], correction, MPS) };
			ll_running += (events_Bminus[i].weight_bkg() * log(probability/normalisation_Bminus));			
			//ll_running += log(probability/normalisation_Bminus);
		}

		#pragma omp parallel for reduction (+:ll_running)
		// for (auto event : events_Bplus) // as above 
		for (size_t i=0; i < events_Bplus.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_Bplus[i]) };
			real_t probability{totalAmplitudeSquared_XY(signs.BplusSign, amplitudes_Bplus.A[i], amplitudes_Bplus.Abar[i], amplitudes_Bplus.deltaD[i], correction, MPS) };
			ll_running += (events_Bplus[i].weight_bkg() * log(probability/normalisation_Bplus));
			//ll_running += events_Bplus[i].weight_bkg() * log(probability/normalisation_Bplus);
		}
		return -2*ll_running;
	};

	INFO("built lambda");

	//******** CHECK HOW LONG LL TAKES TO CALCULATE ********
	if ( NamedParameter<bool>("DoTimeTestOnly", false, "Only do the time test for one LL?") ){
		ProfileClock  TIME_oneLL;
		TIME_oneLL.start();
		real_t testLL(LL_total());
		TIME_oneLL.stop();
		INFO("Start with "<<testLL<<" took "<<TIME_oneLL.t_duration<<"ms to calculate LL_total");
		return 0;
	}

	//********* DO THE FIT ***********
	ProfileClock TIME_fitting; TIME_fitting.start();
	Minimiser minimiser(LL_total, &MPS); // formerly mini

	if (NamedParameter<bool>("DoGradientTest", true)){
		minimiser.gradientTest();
	}
	
	minimiser.doFit();
	TIME_fitting.stop();
	// Now save it, note saving only into the logfile while getting things working
	FitResult* output = new FitResult(minimiser); // formerly fr
	output->writeToFile(logFile);
	INFO("Took " << TIME_fitting << "ms to fit to data");

	// *********** PLOT TIME ************
	if (NamedParameter<bool>("DoPlots", true)){
		INFO("Creating projections of fit and data");
		ProfileClock TIME_plots; TIME_plots.start();
		//writeUnbinnedFitPlot_LHCb(A, Abar, phaseCorrection, MPS, eventsMC, Datalist_LHCb);
		TIME_plots.stop();
		INFO("Took " << TIME_plots << "ms to make plots and write to file: " << plotFile);
	}
		
	return 0;
}