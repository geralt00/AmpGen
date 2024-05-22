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
#include <stdexcept> // Include for std::runtime_error

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
double clip_log(double x, double epsilon = 1e-6) {
    if (x > epsilon) {
        // Direct computation if 'x' is above the threshold.
        return std::log(x);
    } else {
        // Apply Taylor approximation for 'x' near or below 'epsilon'.
		//INFO("log(x) is very small, using taylor approximation");
        double delta_x = x - epsilon;
        return std::log(epsilon) + delta_x / epsilon - std::pow(delta_x / epsilon, 2) / 2.0;

    }
}

double dom_weight (double x){
	if (x == 0) return 1.0;
	else return x;
}

using namespace AmpGen;

int main(int argc, char * argv[]){


  //******* PARSE THE ARGS FROM COMMAND LINE AND THE OPTIONS FILE*********
	OptionsParser::setArgs( argc, argv );

	EventType signalType{NamedParameter<std::string>("EventType", "", "Signal Type to generate, i.e. D0 K0S0 pi- pi+")}; // notes -> (1) MINUS COMES FIRST

	const std::string inputData_LL = NamedParameter<std::string>("Input_LL", "data.root", "Root file containing B->D(->Kspipi)K data for both B+/-");
	const std::string inputData_DD = NamedParameter<std::string>("Input_DD", "data.root", "Root file containing B->D(->Kspipi)K data for both B+/-");

	const std::string logFile = NamedParameter<std::string>("Log", "Fit.log", "Name of the output log file for fit result");

	const std::string modelFile = NamedParameter<std::string>("Model", "/shared/scratch/pc24403/cp_fit_test/sub/Kspipi.opt", "Options file with the amplitude model in it");

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
	EventList_type events_Bminus_LL((inputData_LL + ":Bminus_DalitzEventList").c_str() , signalType);
	EventList_type events_Bplus_LL((inputData_LL + ":Bplus_DalitzEventList").c_str() , signalType);

	EventList_type events_Bminus_DD((inputData_DD + ":Bminus_DalitzEventList").c_str() , signalType);
	EventList_type events_Bplus_DD((inputData_DD + ":Bplus_DalitzEventList").c_str() , signalType);	

	std::vector<EventList_type> Datalist_LHCb{events_Bminus_LL, events_Bminus_DD, events_Bplus_LL, events_Bplus_DD};
	size_t nEvents {events_Bminus_LL.size() + events_Bminus_DD.size() + events_Bplus_LL.size() + events_Bplus_DD.size()};

	//******** LOAD THE MPS OF FITTING VARIABLES ********
	INFO("loading MPS");
	MinuitParameterSet MPS;
	MPS.loadFromStream();

	//******** CREATE POINTS IN FLAT PHASE SPACE FOR NORMALISATION ******
	// nInt points across space of defined decay and randomly distributed
	bool generatedEvents{ false };	
	EventList_type eventsMC_DK_p_LL((normalisationTuple+ "_LL_p.root:DalitzEventList").c_str(), signalType);
	EventList_type eventsMC_DK_m_LL((normalisationTuple + "_LL_m.root:DalitzEventList").c_str(), signalType);
	EventList_type eventsMC_DK_p_DD((normalisationTuple+ "_DD_p.root:DalitzEventList").c_str(), signalType);
	EventList_type eventsMC_DK_m_DD((normalisationTuple + "_DD_m.root:DalitzEventList").c_str(), signalType);




	std::vector<EventList_type> eventsMC{eventsMC_DK_p_LL, eventsMC_DK_m_LL, eventsMC_DK_p_DD, eventsMC_DK_m_DD};


	// Prepare the phase correction
	PhaseCorrection phaseCorrection{MPS};
	phaseCorrection.compilePolynomialExpressions(signalType);

	//******** SORT OUT THE AMPLITUDE MODEL *******
	// load in amplitude model A = amplitude(D0 -> Ks0pipi) and Abar from a Dbar0
	INFO("about to load in kspipi model");
	MinuitParameterSet MPS_Kspipi;
	MPS_Kspipi.loadFromFile(modelFile);
	AddCPConjugate(MPS_Kspipi);


	//******** CREATE THE AMPLITUDE MODEL ********
	CoherentSum A_DK_DD_m(signalType, MPS_Kspipi);
	CoherentSum Abar_DK_DD_m(signalType.conj(true), MPS_Kspipi);
	CoherentSum A_DK_DD_p(signalType, MPS_Kspipi);
	CoherentSum Abar_DK_DD_p(signalType.conj(true), MPS_Kspipi);

	A_DK_DD_p.setMC(eventsMC[2]);
	A_DK_DD_p.prepare();
	Abar_DK_DD_p.setMC(eventsMC[2]);
	Abar_DK_DD_p.prepare();

	A_DK_DD_m.setMC(eventsMC[3]);
	A_DK_DD_m.prepare();
	Abar_DK_DD_m.setMC(eventsMC[3]);
	Abar_DK_DD_m.prepare();

	CoherentSum A_DK_LL_p(signalType, MPS_Kspipi);
	CoherentSum Abar_DK_LL_p(signalType.conj(true), MPS_Kspipi);
	CoherentSum A_DK_LL_m(signalType, MPS_Kspipi);
	CoherentSum Abar_DK_LL_m(signalType.conj(true), MPS_Kspipi);

	A_DK_LL_p.setMC(eventsMC[0]);
	A_DK_LL_p.prepare();
	Abar_DK_LL_p.setMC(eventsMC[0]);
	Abar_DK_LL_p.prepare();

	A_DK_LL_m.setMC(eventsMC[1]);
	A_DK_LL_m.prepare();
	Abar_DK_LL_m.setMC(eventsMC[1]);
	Abar_DK_LL_m.prepare();





	std::vector<CoherentSum> A{A_DK_LL_p, A_DK_LL_m, A_DK_DD_p,  A_DK_DD_m};
	std::vector<CoherentSum> Abar{Abar_DK_LL_p, Abar_DK_LL_m, Abar_DK_DD_p,  Abar_DK_DD_m};
	real_t normA[4] = {A_DK_LL_p.norm(), A_DK_LL_m.norm(), A_DK_DD_p.norm(), A_DK_DD_m.norm()};
	real_t normAbar[4] = {Abar_DK_LL_p.norm(), Abar_DK_LL_m.norm(), Abar_DK_DD_p.norm(), Abar_DK_DD_m.norm()};

	INFO("normA: "<<normA[0]<<" "<<normA[1]<<" "<<normA[2]<<" "<<normA[3]);

	// none of the individual A will change so get the info into a usable format now for speed reasons.
	INFO("Collating amplitude information");
	amplitudeInfo amplitudesMC_DK_LL_p, amplitudesMC_DK_LL_m, amplitudesMC_DK_DD_p, amplitudesMC_DK_DD_m;
	{ 
		amplitudesMC_DK_LL_p = fillAmplitudeInfo(eventsMC[0], A_DK_LL_p, Abar_DK_LL_p);
		amplitudesMC_DK_LL_m = fillAmplitudeInfo(eventsMC[1], A_DK_LL_m, Abar_DK_LL_m);
		amplitudesMC_DK_DD_p = fillAmplitudeInfo(eventsMC[2], A_DK_DD_p, Abar_DK_DD_p);
		amplitudesMC_DK_DD_m = fillAmplitudeInfo(eventsMC[3], A_DK_DD_m, Abar_DK_DD_m);

	}

	std::vector<amplitudeInfo> amplitudesMC{amplitudesMC_DK_LL_p, amplitudesMC_DK_LL_m, amplitudesMC_DK_DD_p, amplitudesMC_DK_DD_m};
	amplitudeInfo amplitudes_Bplus_LL = fillAmplitudeInfo(events_Bplus_LL, A[0], Abar[0]);
	amplitudeInfo amplitudes_Bminus_LL = fillAmplitudeInfo(events_Bminus_LL, A[1], Abar[1]);
	amplitudeInfo amplitudes_Bplus_DD = fillAmplitudeInfo(events_Bplus_DD, A[2], Abar[2]);
	amplitudeInfo amplitudes_Bminus_DD = fillAmplitudeInfo(events_Bminus_DD, A[3], Abar[3]);

	// for normalisation - |A|^2 for later but when efficiency introduce the norm a shoud be weighted by the efficiency
	// struct for signs so don't have to use numbers
	signs signs{};
	//******** BUILD THE LOG LIKELIHOOD LAMBDAS **********
	auto NLL_LL = [&events_Bminus_LL, &events_Bplus_LL, &eventsMC, 
					 &amplitudesMC, &amplitudes_Bminus_LL, &amplitudes_Bplus_LL, &nEvents,
					 &normA, &normAbar, &MPS,
					 &phaseCorrection, &signs]() -> double {

		phaseCorrection.updateCoeffs(MPS); // update the stored coefficients to latest in MPS
		// NORMALISATION:
		std::pair<real_t, real_t> normalisationCrossTerms_p = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC[0], eventsMC[0], phaseCorrection); // {cos term, sin term}
		std::pair<real_t, real_t> normalisationCrossTerms_m = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC[1], eventsMC[1], phaseCorrection); // {cos term, sin term}
		real_t normalisation_Bplus = totalAmplitudeSquared_Integrated(signs.BplusSign, normA[0], normAbar[0], MPS, normalisationCrossTerms_p);
		real_t normalisation_Bminus = totalAmplitudeSquared_Integrated(signs.BminusSign, normA[1], normAbar[1], MPS, normalisationCrossTerms_m);
		std::vector<real_t>  ll_running_LL{0, 0};
		std::vector<real_t> ll_mc{0, 0};
		ll_mc[0] = log(normalisation_Bminus);
		ll_mc[1] = log(normalisation_Bplus);
		
		//#pragma omp parallel for reduction (+:ll_running)
		// for (auto event : events_Bminus) // had to switch to below to make the omp run -> change back if you find out how?
		for (size_t i=0; i < events_Bminus_LL.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_Bminus_LL[i]) };
			real_t probability{totalAmplitudeSquared_XY(signs.BminusSign, amplitudes_Bminus_LL.A[i], amplitudes_Bminus_LL.Abar[i], amplitudes_Bminus_LL.deltaD[i], correction, MPS) };
			real_t amp_s2 = probability;
			real_t ln_data = log(amp_s2) - log(normalisation_Bminus);
			ll_running_LL[0] += ln_data;
		}

		//#pragma omp parallel for reduction (+:ll_running_LL)
		// for (auto event : events_Bplus) // as above 
		for (size_t i=0; i < events_Bplus_LL.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_Bplus_LL[i]) };
			real_t probability{totalAmplitudeSquared_XY(signs.BplusSign, amplitudes_Bplus_LL.A[i], amplitudes_Bplus_LL.Abar[i], amplitudes_Bplus_LL.deltaD[i], correction, MPS) };
			real_t amp_s2 = probability;
			real_t ln_data = log(amp_s2) - log(normalisation_Bplus);
			ll_running_LL[1] += ln_data;
		}


		return  -2*(ll_running_LL[0] + ll_running_LL[1]);
	};

	auto NLL_DD = [&events_Bminus_DD, &events_Bplus_DD, &eventsMC, 
					 &amplitudesMC, &amplitudes_Bminus_DD, &amplitudes_Bplus_DD, &nEvents,
					 &normA, &normAbar, &MPS,
					 &phaseCorrection, &signs]() -> double {

		phaseCorrection.updateCoeffs(MPS); // update the stored coefficients to latest in MPS
		// NORMALISATION:
		std::pair<real_t, real_t> normalisationCrossTerms_p = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC[2], eventsMC[2], phaseCorrection); // {cos term, sin term}
		std::pair<real_t, real_t> normalisationCrossTerms_m = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC[3], eventsMC[3], phaseCorrection); // {cos term, sin term}
		real_t normalisation_Bplus = totalAmplitudeSquared_Integrated(signs.BplusSign, normA[2], normAbar[2], MPS, normalisationCrossTerms_p);
		real_t normalisation_Bminus = totalAmplitudeSquared_Integrated(signs.BminusSign, normA[3], normAbar[3], MPS, normalisationCrossTerms_m);
		std::vector<real_t>  ll_running_DD{0, 0};
		std::vector<real_t> ll_mc{0, 0};
		
		//#pragma omp parallel for reduction (+:ll_running_DD)
		// for (auto event : events_Bminus) // had to switch to below to make the omp run -> change back if you find out how?
		for (size_t i=0; i < events_Bminus_DD.size(); i++)
		{
			// Event event{events_Bminus[i]}
			real_t correction{ phaseCorrection.eval(events_Bminus_DD[i]) };
			real_t probability{totalAmplitudeSquared_XY(signs.BminusSign, amplitudes_Bminus_DD.A[i], amplitudes_Bminus_DD.Abar[i], amplitudes_Bminus_DD.deltaD[i], correction, MPS) };
			real_t amp_s2 = probability;
			real_t ln_data = log(amp_s2)- log(normalisation_Bminus);

			ll_running_DD[0] += ln_data;


		}
		//#pragma omp parallel for reduction (+:ll_running_DD)
		// for (auto event : events_Bplus) // as above 
		for (size_t i=0; i < events_Bplus_DD.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_Bplus_DD[i]) };
			real_t probability{totalAmplitudeSquared_XY(signs.BplusSign, amplitudes_Bplus_DD.A[i], amplitudes_Bplus_DD.Abar[i], amplitudes_Bplus_DD.deltaD[i], correction, MPS) };
			real_t amp_s2 = probability;
			real_t ln_data = log(amp_s2) - log(normalisation_Bplus);
			ll_running_DD[1] += ln_data;

		}

		ll_mc[0] = log(normalisation_Bminus);
		ll_mc[1] = log(normalisation_Bplus);

		return  -2*(ll_running_DD[0] + ll_running_DD[1]);

	};

	INFO("built lambda");

	auto LL_total = [&NLL_LL, &NLL_DD](){
		return    NLL_DD() + NLL_LL();
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
	if (NamedParameter<bool>("DoGradientTest", false)){
		minimiser.gradientTest();
	}

	if (NamedParameter<bool>("DoScan", true)){

		INFO("Starting Scan");
		const double xmin = -1;
		const double xmax = 1;
		int step = 5000;
		const double stepsize = (xmax-xmin)/step;
		TGraph* rt = minimiser.scan(MPS["xPlus"], xmin, xmax, stepsize);
		TFile* f = new TFile("scan_rt.root", "RECREATE");
		rt->Write("rt");
		f->Close();
		TIME_fitting.stop();
		INFO("Took " << TIME_fitting << "ms to scan data");

	}
	else
	{
		minimiser.doFit();
		TIME_fitting.stop();
		// Now save it, note saving only into the logfile while getting things working
		FitResult* output = new FitResult(minimiser); // formerly fr
		output->writeToFile(logFile);
		INFO("Took " << TIME_fitting << "ms to fit to data");
	}
	// *********** PLOT TIME ************
	if (NamedParameter<bool>("DoPlots", true)){
		INFO("Creating projections of fit and data");
		ProfileClock TIME_plots; TIME_plots.start();
		writeUnbinnedFitPlot_LHCb(A, Abar, phaseCorrection, MPS, eventsMC, Datalist_LHCb, "LL");
		TIME_plots.stop();
		INFO("Took " << TIME_plots << "ms to make plots and write to file: " << plotFile);
	}
		
	return 0;
}