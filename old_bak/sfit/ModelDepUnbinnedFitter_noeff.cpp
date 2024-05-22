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

	const std::string inputData_DK = NamedParameter<std::string>("Input_DK", "lhcb_toy_b2dk_0.root", "Root file containing B->D(->Kspipi)K data for both B+/-");
	const std::string inputData_DPi = NamedParameter<std::string>("Input_DPi", "lhcb_toy_b2dpi_0.root", "Root file containing B->D(->Kspipi)pi data for both B+/-");

	const std::string logFile = NamedParameter<std::string>("Log", "Fit.log", "Name of the output log file for fit result");

	const std::string modelFile = NamedParameter<std::string>("Model", "Kspipi.opt", "Options file with the amplitude model in it");

	size_t		      nInt = NamedParameter<size_t>("nInt", 1e7, "Number of events to calculate normalisation - should be large");

	const std::string normalisationAmplitudes = NamedParameter<std::string>("NormalisationAmplitudes", "", "Text file with normalsiation amplitudes");

	const std::string normalisationTuple_DK = NamedParameter<std::string>("NormalisationEvents_DK", "", "Root file with normalsiation events");

	const std::string normalisationTuple_DPi = NamedParameter<std::string>("NormalisationEvents_DPi", "", "Root file with normalsiation events for KLpipi");

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
	EventList_type events_dk_Bminus((inputData_DK + ":Bminus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"), WeightBranch_eff("eff_weight"));
	EventList_type events_dk_Bplus((inputData_DK + ":Bplus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"), WeightBranch_eff("eff_weight"));
	EventList_type events_dpi_Bminus((inputData_DPi + ":Bminus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"), WeightBranch_eff("eff_weight"));
	EventList_type events_dpi_Bplus((inputData_DPi + ":Bplus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"), WeightBranch_eff("eff_weight"));


	std::vector<EventList_type> Datalist_LHCb{events_dk_Bminus, events_dk_Bplus, events_dpi_Bminus, events_dpi_Bplus};
	//******** LOAD THE MPS OF FITTING VARIABLES ********
	INFO("loading MPS");
	MinuitParameterSet MPS;
	MPS.loadFromStream();
	bool generatedEvents{ false };	
	EventList_type eventsMC(signalType);

	INFO("No normalisation tuple found, generating " << nInt << " integration events");
	eventsMC =  Generator<>(signalType, &rndm).generate(nInt);
	generatedEvents = true;


	//******** PREPARE THE PHASE CORRECTION ********
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


	amplitudeInfo amplitudes_dk_Bminus = fillAmplitudeInfo(events_dk_Bminus, A, Abar);
	amplitudeInfo amplitudes_dk_Bplus = fillAmplitudeInfo(events_dk_Bplus, A, Abar);
	amplitudeInfo amplitudes_dpi_Bminus = fillAmplitudeInfo(events_dpi_Bminus, A, Abar);
	amplitudeInfo amplitudes_dpi_Bplus = fillAmplitudeInfo(events_dpi_Bplus, A, Abar);

	// for normalisation - |A|^2 for later but when efficiency introduce the norm a shoud be weighted by the efficiency
	real_t normA{A.norm()};
	real_t normAbar{Abar.norm()};

	// for normalisation - |A|^2 for later but when efficiency introduce the norm a shoud be weighted by the efficiency
	real_t sum_sw[4] = {0, 0, 0, 0};
	int index_data=0;
	for(auto event : Datalist_LHCb){
		for (size_t i = 0; i <event.size() ; i++)
		{
			sum_sw[index_data] += event[i].weight_bkg();
		}
		index_data++;
	}

	INFO("sw_dk_m: "<<sum_sw[0]<<" sw_dk_p: "<<sum_sw[1]<<" sw_dpi_m: "<<sum_sw[2]<<" sw_dpi_p: "<<sum_sw[3]);
	// struct for signs so don't have to use numbers
	signs signs{};
	//******** BUILD THE LOG LIKELIHOOD LAMBDAS **********
	//******get sWeight normalisation
	auto dk_sWeight_norm = [&events_dk_Bminus, &events_dk_Bplus]()
	{
		double_t alpha{0};
		double_t  sw{0};
		double_t sw_sq{0};
		for (size_t i=0; i < events_dk_Bminus.size(); i++)
		{
			sw += events_dk_Bminus[i].weight_bkg();
			sw_sq += events_dk_Bminus[i].weight_bkg()*events_dk_Bminus[i].weight_bkg();
		}
		for (size_t i=0; i < events_dk_Bplus.size(); i++)
		{
			sw += events_dk_Bplus[i].weight_bkg();
			sw_sq += events_dk_Bplus[i].weight_bkg()*events_dk_Bplus[i].weight_bkg();
		}

		alpha = sw/sw_sq;	
		return alpha;
	};

	auto dpi_sWeight_norm = [&events_dpi_Bminus, &events_dpi_Bplus]()
	{
		double_t alpha{0};
		double_t  sw{0};
		double_t sw_sq{0};
		for (size_t i=0; i < events_dpi_Bminus.size(); i++)
		{
			sw += events_dpi_Bminus[i].weight_bkg();
			sw_sq += events_dpi_Bminus[i].weight_bkg()*events_dpi_Bminus[i].weight_bkg();
		}
		for (size_t i=0; i < events_dpi_Bplus.size(); i++)
		{
			sw += events_dpi_Bplus[i].weight_bkg();
			sw_sq += events_dpi_Bplus[i].weight_bkg()*events_dpi_Bplus[i].weight_bkg();
		}

		alpha = sw/sw_sq;	
		return alpha;
	};

	auto dk_alpha = dk_sWeight_norm();
	auto dpi_alpha = dpi_sWeight_norm();


	//B->DK 
	auto LL_DK = [&events_dk_Bminus, &events_dk_Bplus, &eventsMC, 
					 &amplitudesMC, &amplitudes_dk_Bminus, &amplitudes_dk_Bplus,
					 &normA, &normAbar, &MPS,
					 &phaseCorrection, &signs, & dk_alpha, &sum_sw](){

		phaseCorrection.updateCoeffs(MPS); // update the stored coefficients to latest in MPS
		// NORMALISATION:
		std::pair<real_t, real_t> normalisationCrossTerms_p = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC, eventsMC, phaseCorrection); // {cos term, sin term}
		std::pair<real_t, real_t> normalisationCrossTerms_m = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC, eventsMC, phaseCorrection); // {cos term, sin term}
		real_t normalisation_Bplus = totalAmplitudeSquared_Integrated(signs.BplusSign, normA, normAbar, MPS, normalisationCrossTerms_p);
		real_t normalisation_Bminus = totalAmplitudeSquared_Integrated(signs.BminusSign, normA, normAbar, MPS, normalisationCrossTerms_m);		//real_t normalisation_DPi_Bplus = totalAmplitudeSquared_DPi_Integrated(signs.BplusSign, normA[0], normAbar[0], MPS, normalisationCrossTerms_p);
		//real_t normalisation_DPi_Bminus = totalAmplitudeSquared_DPi_Integrated(signs.BminusSign, normA[1], normAbar[1], MPS, normalisationCrossTerms_m);
		real_t  ll_running{0};
		real_t ll_running_norm{0};
		#pragma omp parallel for reduction (+:ll_running)
		// for (auto event : events_Bminus) // had to switch to below to make the omp run -> change back if you find out how?
		for (size_t i=0; i < events_dk_Bminus.size(); i++)
		{
			// Event event{events_Bminus[i]}
			real_t correction{ phaseCorrection.eval(events_dk_Bminus[i]) };
			real_t probability{totalAmplitudeSquared_XY(signs.BminusSign, amplitudes_dk_Bminus.A[i], amplitudes_dk_Bminus.Abar[i], amplitudes_dk_Bminus.deltaD[i], correction, MPS) };
			//real_t probability_misid{totalAmplitudeSquared_DPi_XY(signs.BminusSign, amplitudes_dk_Bminus.A[i], amplitudes_dk_Bminus.Abar[i], amplitudes_dk_Bminus.deltaD[i], correction, MPS) };
			ll_running += (events_dk_Bminus[i].weight_bkg() * log((1.0*probability/normalisation_Bminus)));// + (0.0*probability_misid/normalisation_DPi_Bminus)));
			//ll_running += (events_dk_Bminus[i].weight_bkg() * log(1.0*probability));	

		}

		#pragma omp parallel for reduction (+:ll_running)
		// for (auto event : events_Bplus) // as above 
		for (size_t i=0; i < events_dk_Bplus.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_dk_Bplus[i]) };
			real_t probability{totalAmplitudeSquared_XY(signs.BplusSign, amplitudes_dk_Bplus.A[i], amplitudes_dk_Bplus.Abar[i], amplitudes_dk_Bplus.deltaD[i], correction, MPS) };
			//real_t probability_misid{totalAmplitudeSquared_DPi_XY(signs.BplusSign, amplitudes_dk_Bplus.A[i], amplitudes_dk_Bplus.Abar[i], amplitudes_dk_Bplus.deltaD[i], correction, MPS) };
			ll_running += (events_dk_Bplus[i].weight_bkg() * log((1.0*probability/normalisation_Bplus)));// + (0.0*probability_misid/normalisation_DPi_Bplus)));
		}

		return ll_running;
	};

	//B->DPi
	auto LL_DPi = [&events_dpi_Bminus, &events_dpi_Bplus, &eventsMC, 
					 &amplitudesMC, &amplitudes_dpi_Bminus, &amplitudes_dpi_Bplus,
					 &normA, &normAbar, &MPS,
					 &phaseCorrection, &signs, & dpi_alpha, & sum_sw](){

		phaseCorrection.updateCoeffs(MPS); // update the stored coefficients to latest in MPS
		// NORMALISATION:
		std::pair<real_t, real_t> normalisationCrossTerms_p = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC, eventsMC, phaseCorrection); // {cos term, sin term}
		std::pair<real_t, real_t> normalisationCrossTerms_m = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC, eventsMC, phaseCorrection); // {cos term, sin term}
		real_t normalisation_Bplus = totalAmplitudeSquared_DPi_Integrated(signs.BplusSign, normA, normAbar, MPS, normalisationCrossTerms_p);
		real_t normalisation_Bminus = totalAmplitudeSquared_DPi_Integrated(signs.BminusSign, normA, normAbar, MPS, normalisationCrossTerms_m);
		//real_t normalisation_DK_Bplus = totalAmplitudeSquared_Integrated(signs.BplusSign, normA[2], normAbar[2], MPS, normalisationCrossTerms_p);
		//real_t normalisation_DK_Bminus = totalAmplitudeSquared_Integrated(signs.BminusSign, normA[3], normAbar[3], MPS, normalisationCrossTerms_m);
		real_t  ll_running{0};
		#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_dpi_Bminus.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_dpi_Bminus[i]) };
			real_t probability{totalAmplitudeSquared_DPi_XY(signs.BminusSign, amplitudes_dpi_Bminus.A[i], amplitudes_dpi_Bminus.Abar[i], amplitudes_dpi_Bminus.deltaD[i], correction, MPS) };
			//real_t probability_misid{totalAmplitudeSquared_XY(signs.BminusSign, amplitudes_dpi_Bminus.A[i], amplitudes_dpi_Bminus.Abar[i], amplitudes_dpi_Bminus.deltaD[i], correction, MPS) };
			ll_running += (events_dpi_Bminus[i].weight_bkg() * log((1.0*probability/normalisation_Bminus)));//+(0.0*probability_misid/normalisation_DK_Bminus)));
			//ll_running += (events_dpi_Bminus[i].weight_bkg() * log(1.0*probability) + (1-events_dpi_Bminus[i].weight_bkg())*log(normalisation_Bminus));
		}

		#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_dpi_Bplus.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_dpi_Bplus[i]) };
			real_t probability{totalAmplitudeSquared_DPi_XY(signs.BplusSign, amplitudes_dpi_Bplus.A[i], amplitudes_dpi_Bplus.Abar[i], amplitudes_dpi_Bplus.deltaD[i], correction, MPS) };
			//real_t probability_misid{totalAmplitudeSquared_XY(signs.BplusSign, amplitudes_dpi_Bplus.A[i], amplitudes_dpi_Bplus.Abar[i], amplitudes_dpi_Bplus.deltaD[i], correction, MPS) };
			ll_running += (events_dpi_Bplus[i].weight_bkg() * log((1.0*probability/normalisation_Bplus)));//+(0.0*probability_misid/normalisation_DK_Bplus)));
			//ll_running += (events_dpi_Bplus[i].weight_bkg() * log(1.0*probability) + (1-events_dpi_Bplus[i].weight_bkg())*log(normalisation_Bplus));
		}

		return ll_running;
	};

	INFO("built lambda");

	auto LL_total = [&LL_DK, &LL_DPi, &dk_alpha, &dpi_alpha](){
		return (-2*dk_alpha*LL_DK()  -2*dpi_alpha*LL_DPi());
	};

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
//		writeUnbinnedFitPlot_LHCb(A, Abar, phaseCorrection, MPS, eventsMC, Datalist_LHCb);
		TIME_plots.stop();
		INFO("Took " << TIME_plots << "ms to make plots and write to file: " << plotFile);
	}
	
		
	return 0;
}









