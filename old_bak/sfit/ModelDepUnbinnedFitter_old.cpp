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
	size_t nInt_p{0}, nInt_m{0}, nInt_DK_p{0}, nInt_DK_m{0}, nInt_DPi_p{0}, nInt_DPi_m{0};
	//******** CREATE POINTS IN FLAT PHASE SPACE FOR NORMALISATION ******
	// nInt points across space of defined decay and randomly distributed
	bool generatedEvents{ false };	
	EventList_type eventsMC_DK_m((normalisationTuple_DK + ":Bminus_DalitzEventList").c_str(), signalType, WeightBranch_eff("eff_weight"));
	nInt_DK_m = eventsMC_DK_m.size();
	EventList_type eventsMC_DK_p((normalisationTuple_DK + ":Bplus_DalitzEventList").c_str(), signalType, WeightBranch_eff("eff_weight"));
	nInt_DK_p = eventsMC_DK_p.size();
	EventList_type eventsMC_DPi_m((normalisationTuple_DPi + ":Bminus_DalitzEventList").c_str(), signalType, WeightBranch_eff("eff_weight"));
	nInt_DPi_m = eventsMC_DPi_m.size();
	EventList_type eventsMC_DPi_p((normalisationTuple_DPi + ":Bplus_DalitzEventList").c_str(), signalType, WeightBranch_eff("eff_weight"));
	nInt_DPi_p = eventsMC_DPi_p.size();


	std::vector<size_t> nInts{ nInt_DK_p, nInt_DK_m, nInt_DPi_p, nInt_DPi_m };
	std::vector<EventList_type> eventsMC{eventsMC_DK_p, eventsMC_DK_m, eventsMC_DPi_p, eventsMC_DPi_m};
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
	A.setMC(eventsMC[0]);
	A.prepare();

	Abar.setMC(eventsMC[0]);
	Abar.prepare();

	// none of the individual A will change so get the info into a usable format now for speed reasons.
	INFO("Collating amplitude information");
	amplitudeInfo amplitudesMC_DK_p, amplitudesMC_DK_m, amplitudesMC_DPi_p, amplitudesMC_DPi_m;
	{ 
		amplitudesMC_DK_p = fillAmplitudeInfo(eventsMC[0], A, Abar);
		amplitudesMC_DK_m = fillAmplitudeInfo(eventsMC[1], A, Abar);
		amplitudesMC_DPi_p = fillAmplitudeInfo(eventsMC[2], A, Abar);
		amplitudesMC_DPi_m = fillAmplitudeInfo(eventsMC[3], A, Abar);
	}


	std::vector<amplitudeInfo> amplitudesMC{amplitudesMC_DK_p, amplitudesMC_DK_m, amplitudesMC_DPi_p, amplitudesMC_DPi_m};
	amplitudeInfo amplitudes_dk_Bminus = fillAmplitudeInfo(events_dk_Bminus, A, Abar);
	amplitudeInfo amplitudes_dk_Bplus = fillAmplitudeInfo(events_dk_Bplus, A, Abar);
	amplitudeInfo amplitudes_dpi_Bminus = fillAmplitudeInfo(events_dpi_Bminus, A, Abar);
	amplitudeInfo amplitudes_dpi_Bplus = fillAmplitudeInfo(events_dpi_Bplus, A, Abar);

	// for normalisation - |A|^2 for later but when efficiency introduce the norm a shoud be weighted by the efficiency
	real_t normA[4] = {0, 0, 0, 0};
	real_t normAbar[4] = {0, 0, 0, 0};
	int index_mc=0;
	for(auto event : eventsMC){
		real_t eff_weight_sum{0};
		for (size_t i = 0; i <nInts[index_mc] ; i++)
		{
			normA[index_mc] += amplitudesMC[index_mc].A[i]*amplitudesMC[index_mc].A[i]*eventsMC[index_mc][i].weight_eff();
			normAbar[index_mc] += amplitudesMC[index_mc].Abar[i]*amplitudesMC[index_mc].Abar[i]*eventsMC[index_mc][i].weight_eff();		
			eff_weight_sum += eventsMC[index_mc][i].weight_eff();
		}
		normA[index_mc] = normA[index_mc]/eff_weight_sum;
		normAbar[index_mc] = normAbar[index_mc]/eff_weight_sum;	
		index_mc++;
	}

	INFO("normA: "<<normA[0]<<" "<<normA[1]<<" "<<normA[2]<<" "<<normA[3]);

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
		std::pair<real_t, real_t> normalisationCrossTerms_p = totalAmplitudeSquared_Integrated_crossTerm(amplitudesMC[0], eventsMC[0], phaseCorrection); // {cos term, sin term}
		std::pair<real_t, real_t> normalisationCrossTerms_m = totalAmplitudeSquared_Integrated_crossTerm(amplitudesMC[1], eventsMC[1], phaseCorrection); // {cos term, sin term}

		real_t normalisation_Bplus = totalAmplitudeSquared_Integrated(signs.BplusSign, normA[0], normAbar[0], MPS, normalisationCrossTerms_p);
		real_t normalisation_Bminus = totalAmplitudeSquared_Integrated(signs.BminusSign, normA[1], normAbar[1], MPS, normalisationCrossTerms_m);
		//real_t normalisation_DPi_Bplus = totalAmplitudeSquared_DPi_Integrated(signs.BplusSign, normA[0], normAbar[0], MPS, normalisationCrossTerms_p);
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
			ll_running += (events_dk_Bminus[i].weight_bkg() *log(1.0*probability/normalisation_Bminus) );// + (0.0*probability_misid/normalisation_DPi_Bminus)));
			//ll_running += (events_dk_Bminus[i].weight_bkg() * log(1.0*probability));	

		}

		#pragma omp parallel for reduction (+:ll_running)
		// for (auto event : events_Bplus) // as above 
		for (size_t i=0; i < events_dk_Bplus.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_dk_Bplus[i]) };
			real_t probability{totalAmplitudeSquared_XY(signs.BplusSign, amplitudes_dk_Bplus.A[i], amplitudes_dk_Bplus.Abar[i], amplitudes_dk_Bplus.deltaD[i], correction, MPS) };
			//real_t probability_misid{totalAmplitudeSquared_DPi_XY(signs.BplusSign, amplitudes_dk_Bplus.A[i], amplitudes_dk_Bplus.Abar[i], amplitudes_dk_Bplus.deltaD[i], correction, MPS) };
			ll_running += (events_dk_Bplus[i].weight_bkg() *log(1.0*probability/normalisation_Bplus));// + (0.0*probability_misid/normalisation_DPi_Bplus)));
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
		std::pair<real_t, real_t> normalisationCrossTerms_p = totalAmplitudeSquared_Integrated_crossTerm(amplitudesMC[2], eventsMC[2], phaseCorrection); // {cos term, sin term}
		std::pair<real_t, real_t> normalisationCrossTerms_m = totalAmplitudeSquared_Integrated_crossTerm(amplitudesMC[3], eventsMC[3], phaseCorrection); // {cos term, sin term}

		real_t normalisation_Bplus = totalAmplitudeSquared_DPi_Integrated(signs.BplusSign, normA[2], normAbar[2], MPS, normalisationCrossTerms_p);
		real_t normalisation_Bminus = totalAmplitudeSquared_DPi_Integrated(signs.BminusSign, normA[3], normAbar[3], MPS, normalisationCrossTerms_m);
		//real_t normalisation_DK_Bplus = totalAmplitudeSquared_Integrated(signs.BplusSign, normA[2], normAbar[2], MPS, normalisationCrossTerms_p);
		//real_t normalisation_DK_Bminus = totalAmplitudeSquared_Integrated(signs.BminusSign, normA[3], normAbar[3], MPS, normalisationCrossTerms_m);
		real_t  ll_running{0};
		#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_dpi_Bminus.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_dpi_Bminus[i]) };
			real_t probability{totalAmplitudeSquared_DPi_XY(signs.BminusSign, amplitudes_dpi_Bminus.A[i], amplitudes_dpi_Bminus.Abar[i], amplitudes_dpi_Bminus.deltaD[i], correction, MPS) };
			//real_t probability_misid{totalAmplitudeSquared_XY(signs.BminusSign, amplitudes_dpi_Bminus.A[i], amplitudes_dpi_Bminus.Abar[i], amplitudes_dpi_Bminus.deltaD[i], correction, MPS) };
			ll_running += (events_dpi_Bminus[i].weight_bkg() *log(1.0*probability/normalisation_Bminus) );//+(0.0*probability_misid/normalisation_DK_Bminus)));
			//ll_running += (events_dpi_Bminus[i].weight_bkg() * log(1.0*probability) + (1-events_dpi_Bminus[i].weight_bkg())*log(normalisation_Bminus));
		}

		#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_dpi_Bplus.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_dpi_Bplus[i]) };
			real_t probability{totalAmplitudeSquared_DPi_XY(signs.BplusSign, amplitudes_dpi_Bplus.A[i], amplitudes_dpi_Bplus.Abar[i], amplitudes_dpi_Bplus.deltaD[i], correction, MPS) };
			//real_t probability_misid{totalAmplitudeSquared_XY(signs.BplusSign, amplitudes_dpi_Bplus.A[i], amplitudes_dpi_Bplus.Abar[i], amplitudes_dpi_Bplus.deltaD[i], correction, MPS) };
			ll_running += (events_dpi_Bplus[i].weight_bkg() *log(1.0*probability/normalisation_Bplus));//+(0.0*probability_misid/normalisation_DK_Bplus)));
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
		//writeUnbinnedFitPlot_LHCb(A, Abar, phaseCorrection, MPS, eventsMC, Datalist_LHCb);
		TIME_plots.stop();
		INFO("Took " << TIME_plots << "ms to make plots and write to file: " << plotFile);
	}
	
		
	return 0;
}








