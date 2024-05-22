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
int main(int argc, char * argv[]){


  //******* PARSE THE ARGS FROM COMMAND LINE AND THE OPTIONS FILE*********
	OptionsParser::setArgs( argc, argv );

	EventType signalType{NamedParameter<std::string>("EventType", "", "Signal Type to generate, i.e. D0 K0S0 pi- pi+")}; // notes -> (1) MINUS COMES FIRST

	const std::string inputData_DK_LL = NamedParameter<std::string>("Input_DK_LL", "lhcb_toy_b2dk_0.root", "Root file containing B->D(->Kspipi)K data for both B+/-");
	const std::string inputData_DK_DD = NamedParameter<std::string>("Input_DK_DD", "lhcb_toy_b2dk_0.root", "Root file containing B->D(->Kspipi)K data for both B+/-");
	const std::string inputData_DPi_LL = NamedParameter<std::string>("Input_DPi_LL", "lhcb_toy_b2dpi_0.root", "Root file containing B->D(->Kspipi)pi data for both B+/-");
	const std::string inputData_DPi_DD = NamedParameter<std::string>("Input_DPi_DD", "lhcb_toy_b2dpi_0.root", "Root file containing B->D(->Kspipi)pi data for both B+/-");


	const std::string logFile = NamedParameter<std::string>("Log", "Fit.log", "Name of the output log file for fit result");

	const std::string modelFile = NamedParameter<std::string>("Model", "/software/pc24403/PCBPGGSZ/cpfit/sub/Kspipi.opt", "Options file with the amplitude model in it");

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
	EventList_type events_dk_LL_Bminus((inputData_DK_LL + ":Bminus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"));
	EventList_type events_dk_LL_Bplus((inputData_DK_LL + ":Bplus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"));
	EventList_type events_dk_DD_Bminus((inputData_DK_DD + ":Bminus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"));
	EventList_type events_dk_DD_Bplus((inputData_DK_DD + ":Bplus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"));

	EventList_type events_dpi_LL_Bminus((inputData_DPi_LL + ":Bminus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"));
	EventList_type events_dpi_LL_Bplus((inputData_DPi_LL + ":Bplus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"));
	EventList_type events_dpi_DD_Bminus((inputData_DPi_DD + ":Bminus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"));
	EventList_type events_dpi_DD_Bplus((inputData_DPi_DD + ":Bplus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw"));


	std::vector<EventList_type> Datalist_LHCb_LL{events_dk_LL_Bminus, events_dk_LL_Bplus, events_dpi_LL_Bminus, events_dpi_LL_Bplus};
	std::vector<EventList_type> Datalist_LHCb_DD{events_dk_DD_Bminus, events_dk_DD_Bplus, events_dpi_DD_Bminus, events_dpi_DD_Bplus};

	//******** LOAD THE MPS OF FITTING VARIABLES ********
	INFO("loading MPS");
	MinuitParameterSet MPS;
	MPS.loadFromStream();
	size_t nInt_p{0}, nInt_m{0}, nInt_DK_p{0}, nInt_DK_m{0}, nInt_DPi_p{0}, nInt_DPi_m{0};
	//******** CREATE POINTS IN FLAT PHASE SPACE FOR NORMALISATION ******
	// nInt points across space of defined decay and randomly distributed
	bool generatedEvents{ false };	
	EventList_type eventsMC_DK_LL_m((normalisationTuple_DK + "_LL_m.root:DalitzEventList").c_str(), signalType);
	EventList_type eventsMC_DK_LL_p((normalisationTuple_DK + "_LL_p.root:DalitzEventList").c_str(), signalType);
	EventList_type eventsMC_DK_DD_m((normalisationTuple_DK + "_DD_m.root:DalitzEventList").c_str(), signalType);
	EventList_type eventsMC_DK_DD_p((normalisationTuple_DK + "_DD_p.root:DalitzEventList").c_str(), signalType);

	EventList_type eventsMC_DPi_LL_m((normalisationTuple_DPi + "_LL_m.root:DalitzEventList").c_str(), signalType);
	EventList_type eventsMC_DPi_LL_p((normalisationTuple_DPi + "_LL_p.root:DalitzEventList").c_str(), signalType);
	EventList_type eventsMC_DPi_DD_m((normalisationTuple_DPi + "_DD_m.root:DalitzEventList").c_str(), signalType);
	EventList_type eventsMC_DPi_DD_p((normalisationTuple_DPi + "_DD_p.root:DalitzEventList").c_str(), signalType);


	std::vector<EventList_type> eventsMC_LL{eventsMC_DK_LL_p, eventsMC_DK_LL_m, eventsMC_DPi_LL_p, eventsMC_DPi_LL_m};
	std::vector<EventList_type> eventsMC_DD{eventsMC_DK_DD_p, eventsMC_DK_DD_m, eventsMC_DPi_DD_p, eventsMC_DPi_DD_m};

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

	//******** CREATE THE AMPLITUDE MODEL ********
	CoherentSum A_DK_LL_p(signalType, MPS_Kspipi);
	CoherentSum Abar_DK_LL_p(signalType.conj(true), MPS_Kspipi);
	CoherentSum A_DK_LL_m(signalType, MPS_Kspipi);
	CoherentSum Abar_DK_LL_m(signalType.conj(true), MPS_Kspipi);
	CoherentSum A_DPi_LL_p(signalType, MPS_Kspipi);
	CoherentSum Abar_DPi_LL_p(signalType.conj(true), MPS_Kspipi);
	CoherentSum A_DPi_LL_m(signalType, MPS_Kspipi);
	CoherentSum Abar_DPi_LL_m(signalType.conj(true), MPS_Kspipi);
	A_DK_LL_p.setMC(eventsMC_LL[0]);
	A_DK_LL_p.prepare();
	Abar_DK_LL_p.setMC(eventsMC_LL[0]);
	Abar_DK_LL_p.prepare();
	A_DK_LL_m.setMC(eventsMC_LL[1]);
	A_DK_LL_m.prepare();
	Abar_DK_LL_m.setMC(eventsMC_LL[1]);
	Abar_DK_LL_m.prepare();
	A_DPi_LL_p.setMC(eventsMC_LL[2]);
	A_DPi_LL_p.prepare();
	Abar_DPi_LL_p.setMC(eventsMC_LL[2]);
	Abar_DPi_LL_p.prepare();
	A_DPi_LL_m.setMC(eventsMC_LL[3]);
	A_DPi_LL_m.prepare();
	Abar_DPi_LL_m.setMC(eventsMC_LL[3]);
	Abar_DPi_LL_m.prepare();

	CoherentSum A_DK_DD_m(signalType, MPS_Kspipi);
	CoherentSum Abar_DK_DD_m(signalType.conj(true), MPS_Kspipi);
	CoherentSum A_DK_DD_p(signalType, MPS_Kspipi);
	CoherentSum Abar_DK_DD_p(signalType.conj(true), MPS_Kspipi);
	CoherentSum A_DPi_DD_m(signalType, MPS_Kspipi);
	CoherentSum Abar_DPi_DD_m(signalType.conj(true), MPS_Kspipi);
	CoherentSum A_DPi_DD_p(signalType, MPS_Kspipi);
	CoherentSum Abar_DPi_DD_p(signalType.conj(true), MPS_Kspipi);
	A_DK_DD_p.setMC(eventsMC_DD[0]);
	A_DK_DD_p.prepare();
	Abar_DK_DD_p.setMC(eventsMC_DD[0]);
	Abar_DK_DD_p.prepare();
	A_DK_DD_m.setMC(eventsMC_DD[1]);
	A_DK_DD_m.prepare();
	Abar_DK_DD_m.setMC(eventsMC_DD[1]);
	Abar_DK_DD_m.prepare();
	A_DPi_DD_p.setMC(eventsMC_DD[2]);
	A_DPi_DD_p.prepare();
	Abar_DPi_DD_p.setMC(eventsMC_DD[2]);
	Abar_DPi_DD_p.prepare();
	A_DPi_DD_m.setMC(eventsMC_DD[3]);
	A_DPi_DD_m.prepare();
	Abar_DPi_DD_m.setMC(eventsMC_DD[3]);
	Abar_DPi_DD_m.prepare();


	std::vector<CoherentSum> A_LL{A_DK_LL_p, A_DK_LL_m, A_DPi_LL_p, A_DPi_LL_m};
	std::vector<CoherentSum> A_DD{A_DK_DD_p, A_DK_DD_m, A_DPi_DD_p, A_DPi_DD_m};
	std::vector<CoherentSum> Abar_LL{Abar_DK_LL_p, Abar_DK_LL_m, Abar_DPi_LL_p, Abar_DPi_LL_m};
	std::vector<CoherentSum> Abar_DD{Abar_DK_DD_p, Abar_DK_DD_m, Abar_DPi_DD_p, Abar_DPi_DD_m};

	real_t normA_LL[4] = {A_DK_LL_p.norm(), A_DK_LL_m.norm(), A_DPi_LL_p.norm(), A_DPi_LL_m.norm()};
	real_t normAbar_LL[4] = {Abar_DK_LL_p.norm(), Abar_DK_LL_m.norm(), Abar_DPi_LL_p.norm(), Abar_DPi_LL_m.norm()};
	real_t normA_DD[4] = {A_DK_DD_p.norm(), A_DK_DD_m.norm(), A_DPi_DD_p.norm(), A_DPi_DD_m.norm()};
	real_t normAbar_DD[4] = {Abar_DK_DD_p.norm(), Abar_DK_DD_m.norm(), Abar_DPi_DD_p.norm(), Abar_DPi_DD_m.norm()};

	// none of the individual A will change so get the info into a usable format now for speed reasons.
	INFO("Collating amplitude information");
	amplitudeInfo amplitudesMC_DK_LL_p, amplitudesMC_DK_LL_m, amplitudesMC_DPi_LL_p, amplitudesMC_DPi_LL_m;
	amplitudeInfo amplitudesMC_DK_DD_p, amplitudesMC_DK_DD_m, amplitudesMC_DPi_DD_p, amplitudesMC_DPi_DD_m;
	{ 

		amplitudesMC_DK_LL_p = fillAmplitudeInfo(eventsMC_LL[0], A_DK_LL_p, Abar_DK_LL_p);
		amplitudesMC_DK_LL_m = fillAmplitudeInfo(eventsMC_LL[1], A_DK_LL_m, Abar_DK_LL_m);
		amplitudesMC_DPi_LL_p = fillAmplitudeInfo(eventsMC_LL[2], A_DPi_LL_p, Abar_DPi_LL_p);
		amplitudesMC_DPi_LL_m = fillAmplitudeInfo(eventsMC_LL[3], A_DPi_LL_m, Abar_DPi_LL_m);

		amplitudesMC_DK_DD_p = fillAmplitudeInfo(eventsMC_DD[0], A_DK_DD_p, Abar_DK_DD_p);
		amplitudesMC_DK_DD_m = fillAmplitudeInfo(eventsMC_DD[1], A_DK_DD_m, Abar_DK_DD_m);
		amplitudesMC_DPi_DD_p = fillAmplitudeInfo(eventsMC_DD[2], A_DPi_DD_p, Abar_DPi_DD_p);
		amplitudesMC_DPi_DD_m = fillAmplitudeInfo(eventsMC_DD[3], A_DPi_DD_m, Abar_DPi_DD_m);
	}


 	std::vector<amplitudeInfo>  amplitudesMC_LL{amplitudesMC_DK_LL_p, amplitudesMC_DK_LL_m, amplitudesMC_DPi_LL_p, amplitudesMC_DPi_LL_m};
	std::vector<amplitudeInfo>  amplitudesMC_DD{amplitudesMC_DK_DD_p, amplitudesMC_DK_DD_m, amplitudesMC_DPi_DD_p, amplitudesMC_DPi_DD_m};

	
	amplitudeInfo amplitudes_dk_LL_Bminus, amplitudes_dk_LL_Bplus, amplitudes_dpi_LL_Bminus, amplitudes_dpi_LL_Bplus;
	amplitudeInfo amplitudes_dk_DD_Bminus, amplitudes_dk_DD_Bplus, amplitudes_dpi_DD_Bminus, amplitudes_dpi_DD_Bplus;
	{
		amplitudes_dk_LL_Bminus = fillAmplitudeInfo(events_dk_LL_Bminus, A_DK_LL_m, Abar_DK_LL_m);
		amplitudes_dk_LL_Bplus = fillAmplitudeInfo(events_dk_LL_Bplus, A_DK_LL_p, Abar_DK_LL_p);
		amplitudes_dpi_LL_Bminus = fillAmplitudeInfo(events_dpi_LL_Bminus, A_DPi_LL_m, Abar_DPi_LL_m);
		amplitudes_dpi_LL_Bplus = fillAmplitudeInfo(events_dpi_LL_Bplus, A_DPi_LL_p, Abar_DPi_LL_p);

		amplitudes_dk_DD_Bminus = fillAmplitudeInfo(events_dk_DD_Bminus, A_DK_DD_m, Abar_DK_DD_m);
		amplitudes_dk_DD_Bplus = fillAmplitudeInfo(events_dk_DD_Bplus, A_DK_DD_p, Abar_DK_DD_p);
		amplitudes_dpi_DD_Bminus = fillAmplitudeInfo(events_dpi_DD_Bminus, A_DPi_DD_m, Abar_DPi_DD_m);
		amplitudes_dpi_DD_Bplus = fillAmplitudeInfo(events_dpi_DD_Bplus, A_DPi_DD_p, Abar_DPi_DD_p);
	}


	std::vector<amplitudeInfo> amplitudes_LL{amplitudes_dk_LL_Bplus, amplitudes_dk_LL_Bminus, amplitudes_dpi_LL_Bplus, amplitudes_dpi_LL_Bminus};
	std::vector<amplitudeInfo> amplitudes_DD{amplitudes_dk_DD_Bplus, amplitudes_dk_DD_Bminus, amplitudes_dpi_DD_Bplus, amplitudes_dpi_DD_Bminus};
	


	// for normalisation - |A|^2 for later but when efficiency introduce the norm a shoud be weighted by the efficiency
	/*
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
	*/
	// struct for signs so don't have to use numbers
	signs signs{};
	//******** BUILD THE LOG LIKELIHOOD LAMBDAS **********
	//******get sWeight normalisation

	// get the normalisation for the sWeights
	auto dk_sWeight_norm = [](EventList_type& events_dk_Bminus, EventList_type& events_dk_Bplus) {
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

	auto dpi_sWeight_norm = [](EventList_type& events_dpi_Bminus, EventList_type& events_dpi_Bplus) {
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


/*	
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
	*/
	auto dk_LL_alpha = dk_sWeight_norm(Datalist_LHCb_LL[0], Datalist_LHCb_LL[1]);
	auto dpi_LL_alpha = dpi_sWeight_norm(Datalist_LHCb_LL[2], Datalist_LHCb_LL[3]);
	auto dk_DD_alpha = dk_sWeight_norm(Datalist_LHCb_DD[0], Datalist_LHCb_DD[1]);
	auto dpi_DD_alpha = dpi_sWeight_norm(Datalist_LHCb_DD[2], Datalist_LHCb_DD[3]);

	INFO("dk_LL_alpha: "<<dk_LL_alpha<<" dpi_LL_alpha: "<<dpi_LL_alpha<<" dk_DD_alpha: "<<dk_DD_alpha<<" dpi_DD_alpha: "<<dpi_DD_alpha);

	std::vector<real_t> dk_alpha{dk_LL_alpha, dk_DD_alpha};
	std::vector<real_t> dpi_alpha{dpi_LL_alpha, dpi_DD_alpha};

	//Initialize the bad points finder
	bool flag_dk_LL_bp[200000]; bool flag_dk_LL_bm[200000];
	bool flag_dk_DD_bp[200000]; bool flag_dk_DD_bm[200000];
	bool flag_dpi_LL_bp[200000]; bool flag_dpi_LL_bm[200000];
	bool flag_dpi_DD_bp[200000]; bool flag_dpi_DD_bm[200000];


	//inintialise the bools to true
	for (size_t i = 0; i < 200000; i++)
	{
		flag_dk_LL_bp[i] = true;
		flag_dk_LL_bm[i] = true;
		flag_dk_DD_bp[i] = true;
		flag_dk_DD_bm[i] = true;
		flag_dpi_LL_bp[i] = true;
		flag_dpi_LL_bm[i] = true;
		flag_dpi_DD_bp[i] = true;
		flag_dpi_DD_bm[i] = true;
	}

	std::vector<bool*> flag_dk_LL{flag_dk_LL_bp, flag_dk_LL_bm};
	std::vector<bool*> flag_dk_DD{flag_dk_DD_bp, flag_dk_DD_bm};
	std::vector<bool*> flag_dpi_LL{flag_dpi_LL_bp, flag_dpi_LL_bm};
	std::vector<bool*> flag_dpi_DD{flag_dpi_DD_bp, flag_dpi_DD_bm};

	//B->DK LL
	auto NLL_DK_LL = [&events_dk_LL_Bminus, &events_dk_LL_Bplus, &eventsMC_LL,
						&amplitudes_LL, &amplitudesMC_LL,
						  &normA_LL, &normAbar_LL, &MPS,
					 &phaseCorrection, &signs, &flag_dk_LL, & dk_alpha]() -> double {

		phaseCorrection.updateCoeffs(MPS); // update the stored coefficients to latest in MPS
		// NORMALISATION:
		std::pair<real_t, real_t> normalisationCrossTerms_p = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC_LL[0], eventsMC_LL[0], phaseCorrection); // {cos term, sin term}
		std::pair<real_t, real_t> normalisationCrossTerms_m = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC_LL[1], eventsMC_LL[1], phaseCorrection); // {cos term, sin term}

		
		real_t normalisation_Bplus = totalAmplitudeSquared_Integrated(signs.BplusSign, normA_LL[0], normAbar_LL[0], MPS, normalisationCrossTerms_p);
		real_t normalisation_Bminus = totalAmplitudeSquared_Integrated(signs.BminusSign, normA_LL[1], normAbar_LL[1], MPS, normalisationCrossTerms_m);
		real_t normalisation_DPi_Bplus = totalAmplitudeSquared_DPi_Integrated(signs.BplusSign, normA_LL[2], normAbar_LL[2], MPS, normalisationCrossTerms_p);
		real_t normalisation_DPi_Bminus = totalAmplitudeSquared_DPi_Integrated(signs.BminusSign, normA_LL[3], normAbar_LL[3], MPS, normalisationCrossTerms_m);
		real_t frac = 0.22;
		real_t  ll_running{0};
		//#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_dk_LL_Bminus.size(); i++)
		{
			if(flag_dk_LL[1][i]==false) continue;
			real_t correction{ phaseCorrection.eval(events_dk_LL_Bminus[i]) };
			real_t probability{totalAmplitudeSquared_XY(signs.BminusSign, amplitudes_LL[1].A[i], amplitudes_LL[1].Abar[i], amplitudes_LL[1].deltaD[i], correction, MPS) };
			real_t probability_misid{totalAmplitudeSquared_DPi_XY(signs.BminusSign, amplitudes_LL[1].A[i], amplitudes_LL[1].Abar[i], amplitudes_LL[1].deltaD[i], correction, MPS) };
	
			if(probability<1e-6&probability!=0){
				flag_dk_LL[1][i] = false;
				INFO("Bad point "<<i<<" in B- DK LL gives probability "<<probability);
				return 0;
			}
			real_t ll_data = log((probability/normalisation_Bminus)*(1-frac) + (probability_misid/normalisation_DPi_Bminus)*frac);
			ll_running += events_dk_LL_Bminus[i].weight_bkg()*ll_data;
		}

		//#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_dk_LL_Bplus.size(); i++)
		{
			if(flag_dk_LL[0][i]==false) continue;
			real_t correction{ phaseCorrection.eval(events_dk_LL_Bplus[i]) };
			real_t probability{totalAmplitudeSquared_XY(signs.BplusSign, amplitudes_LL[0].A[i], amplitudes_LL[0].Abar[i], amplitudes_LL[0].deltaD[i], correction, MPS) };
			real_t probability_misid{totalAmplitudeSquared_DPi_XY(signs.BplusSign, amplitudes_LL[0].A[i], amplitudes_LL[0].Abar[i], amplitudes_LL[0].deltaD[i], correction, MPS) };
			if(probability<1e-6&probability!=0){
				flag_dk_LL[0][i] = false;
				INFO("Bad point "<<i<<" in B+ DK LL gives probability "<<probability);
				return 0;
			}
			real_t ll_data = log((probability/normalisation_Bplus)*(1-frac) + (probability_misid/normalisation_DPi_Bplus)*frac);
			ll_running += events_dk_LL_Bplus[i].weight_bkg()*ll_data;
		}

		return -2*dk_alpha[0]*ll_running;
	};

	//B->DK DD
	auto NLL_DK_DD = [&events_dk_DD_Bminus, &events_dk_DD_Bplus, &eventsMC_DD,
						&amplitudes_DD, &amplitudesMC_DD,
						  &normA_DD, &normAbar_DD, &MPS,
					 &phaseCorrection, &signs, &flag_dk_DD, & dk_alpha]() -> double {

		phaseCorrection.updateCoeffs(MPS); // update the stored coefficients to latest in MPS
		// NORMALISATION:
		std::pair<real_t, real_t> normalisationCrossTerms_p = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC_DD[0], eventsMC_DD[0], phaseCorrection); // {cos term, sin term}
		std::pair<real_t, real_t> normalisationCrossTerms_m = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC_DD[1], eventsMC_DD[1], phaseCorrection); // {cos term, sin term}
		
		real_t normalisation_Bplus = totalAmplitudeSquared_Integrated(signs.BplusSign, normA_DD[0], normAbar_DD[0], MPS, normalisationCrossTerms_p);
		real_t normalisation_Bminus = totalAmplitudeSquared_Integrated(signs.BminusSign, normA_DD[1], normAbar_DD[1], MPS, normalisationCrossTerms_m);
		real_t normalisation_DPi_Bplus = totalAmplitudeSquared_DPi_Integrated(signs.BplusSign,  normA_DD[2], normAbar_DD[2], MPS, normalisationCrossTerms_p);
		real_t normalisation_DPi_Bminus = totalAmplitudeSquared_DPi_Integrated(signs.BminusSign, normA_DD[3], normAbar_DD[3], MPS, normalisationCrossTerms_m);
		real_t frac = 0.21;
		real_t  ll_running{0};
		//#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_dk_DD_Bminus.size(); i++)
		{
			if(flag_dk_DD[1][i]==false) continue;
			real_t correction{ phaseCorrection.eval(events_dk_DD_Bminus[i]) };
			real_t probability{totalAmplitudeSquared_XY(signs.BminusSign, amplitudes_DD[1].A[i], amplitudes_DD[1].Abar[i], amplitudes_DD[1].deltaD[i], correction, MPS) };
			real_t probability_misid{totalAmplitudeSquared_DPi_XY(signs.BminusSign, amplitudes_DD[1].A[i], amplitudes_DD[1].Abar[i], amplitudes_DD[1].deltaD[i], correction, MPS) };
			if(probability<1e-6&probability!=0){
				flag_dk_DD[1][i] = false;
				INFO("Bad point "<<i<<" in B- DK DD gives probability "<<probability);
				return 0;
			}
			real_t ll_data = log((probability/normalisation_Bminus)*(1-frac) + (probability_misid/normalisation_DPi_Bminus)*frac);
			ll_running += events_dk_DD_Bminus[i].weight_bkg()*ll_data;
		}

		//#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_dk_DD_Bplus.size(); i++)
		{
			if(flag_dk_DD[0][i]==false) continue;
			real_t correction{ phaseCorrection.eval(events_dk_DD_Bplus[i]) };
			real_t probability{totalAmplitudeSquared_XY(signs.BplusSign, amplitudes_DD[0].A[i], amplitudes_DD[0].Abar[i], amplitudes_DD[0].deltaD[i], correction, MPS) };
			real_t probability_misid{totalAmplitudeSquared_DPi_XY(signs.BplusSign, amplitudes_DD[0].A[i], amplitudes_DD[0].Abar[i], amplitudes_DD[0].deltaD[i], correction, MPS) };
			if(probability<1e-6&probability!=0){
				flag_dk_DD[0][i] = false;
				INFO("Bad point "<<i<<" in B+ DK DD gives probability "<<probability);
				return 0;
			}
			real_t ll_data = log((probability/normalisation_Bplus)*(1-frac) + (probability_misid/normalisation_DPi_Bplus)*frac);
			ll_running += events_dk_DD_Bplus[i].weight_bkg()*ll_data;
		}

		return -2*dk_alpha[1]*ll_running;
	};


	//B->DPi LL
	auto NLL_DPi_LL = [&events_dpi_LL_Bminus, &events_dpi_LL_Bplus, &eventsMC_LL,
						&amplitudes_LL, &amplitudesMC_LL,
						  &normA_LL, &normAbar_LL, &MPS,
					 &phaseCorrection, &signs, &flag_dpi_LL, & dpi_alpha]() -> double {

		phaseCorrection.updateCoeffs(MPS); // update the stored coefficients to latest in MPS
		// NORMALISATION:
		std::pair<real_t, real_t> normalisationCrossTerms_p = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC_LL[2], eventsMC_LL[2], phaseCorrection); // {cos term, sin term}
		std::pair<real_t, real_t> normalisationCrossTerms_m = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC_LL[3], eventsMC_LL[3], phaseCorrection); // {cos term, sin term}

		
		real_t normalisation_Bplus = totalAmplitudeSquared_DPi_Integrated(signs.BplusSign, normA_LL[2], normAbar_LL[2], MPS, normalisationCrossTerms_p);
		real_t normalisation_Bminus = totalAmplitudeSquared_DPi_Integrated(signs.BminusSign, normA_LL[3], normAbar_LL[3], MPS, normalisationCrossTerms_m);
		//real_t normalisation_DPi_Bplus = totalAmplitudeSquared_DPi_Integrated(signs.BplusSign, normA[0], normAbar[0], MPS, normalisationCrossTerms_p);
		//real_t normalisation_DPi_Bminus = totalAmplitudeSquared_DPi_Integrated(signs.BminusSign, normA[1], normAbar[1], MPS, normalisationCrossTerms_m);
		real_t  ll_running{0};
		//#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_dpi_LL_Bminus.size(); i++)
		{
			if(flag_dpi_LL[1][i]==false) continue;
			real_t correction{ phaseCorrection.eval(events_dpi_LL_Bminus[i]) };
			real_t probability{totalAmplitudeSquared_DPi_XY(signs.BminusSign, amplitudes_LL[3].A[i], amplitudes_LL[3].Abar[i], amplitudes_LL[3].deltaD[i], correction, MPS) };
			//real_t probability_misid{totalAmplitudeSquared_DPi_XY(signs.BminusSign, amplitudes_dk_Bminus.A[i], amplitudes_dk_Bminus.Abar[i], amplitudes_dk_Bminus.deltaD[i], correction, MPS) };	
			if(probability<1e-6&probability!=0){
				flag_dpi_LL[1][i] = false;
				INFO("Bad point "<<i<<" in B- DPi LL gives probability "<<probability);
				return 0;

			}
			ll_running += events_dpi_LL_Bminus[i].weight_bkg()*(clip_log(probability) - (normalisation_Bminus));
		}

		//#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_dpi_LL_Bplus.size(); i++)
		{
			if(flag_dpi_LL[0][i]==false) continue;
			real_t correction{ phaseCorrection.eval(events_dpi_LL_Bplus[i]) };
			real_t probability{totalAmplitudeSquared_DPi_XY(signs.BplusSign, amplitudes_LL[2].A[i], amplitudes_LL[2].Abar[i], amplitudes_LL[2].deltaD[i], correction, MPS) };
			//real_t probability_misid{totalAmplitudeSquared_DPi_XY(signs.BplusSign, amplitudes_dk_Bplus.A[i], amplitudes_dk_Bplus.Abar[i], amplitudes_dk_Bplus.deltaD[i], correction, MPS) };
			if(probability<1e-6&probability!=0){
				flag_dpi_LL[0][i] = false;
				INFO("Bad point "<<i<<" in B+ DPi LL gives probability "<<probability);
				return 0;
			}
			ll_running += events_dpi_LL_Bplus[i].weight_bkg()*(clip_log(probability) - (normalisation_Bplus));
		}

		return -2*dpi_alpha[0]*ll_running;
	};

	//B->DPi DD
	auto NLL_DPi_DD = [&events_dpi_DD_Bminus, &events_dpi_DD_Bplus, &eventsMC_DD,
						&amplitudes_DD, &amplitudesMC_DD,
						  &normA_DD, &normAbar_DD, &MPS,
					 &phaseCorrection, &signs, &flag_dpi_DD, & dpi_alpha]() -> double {

		phaseCorrection.updateCoeffs(MPS); // update the stored coefficients to latest in MPS
		// NORMALISATION:
		std::pair<real_t, real_t> normalisationCrossTerms_p = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC_DD[2], eventsMC_DD[2], phaseCorrection); // {cos term, sin term}
		std::pair<real_t, real_t> normalisationCrossTerms_m = totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudesMC_DD[3], eventsMC_DD[3], phaseCorrection); // {cos term, sin term}

		
		real_t normalisation_Bplus = totalAmplitudeSquared_DPi_Integrated(signs.BplusSign, normA_DD[2], normAbar_DD[2], MPS, normalisationCrossTerms_p);
		real_t normalisation_Bminus = totalAmplitudeSquared_DPi_Integrated(signs.BminusSign, normA_DD[3], normAbar_DD[3], MPS, normalisationCrossTerms_m);
		//real_t normalisation_DPi_Bplus = totalAmplitudeSquared_DPi_Integrated(signs.BplusSign, normA[0], normAbar[0], MPS, normalisationCrossTerms_p);
		//real_t normalisation_DPi_Bminus = totalAmplitudeSquared_DPi_Integrated(signs.BminusSign, normA[1], normAbar[1], MPS, normalisationCrossTerms_m);
		real_t  ll_running{0};
		//#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_dpi_DD_Bminus.size(); i++)
		{
			if(flag_dpi_DD[1][i]==false) continue;
			real_t correction{ phaseCorrection.eval(events_dpi_DD_Bminus[i]) };
			real_t probability{totalAmplitudeSquared_DPi_XY(signs.BminusSign, amplitudes_DD[3].A[i], amplitudes_DD[3].Abar[i], amplitudes_DD[3].deltaD[i], correction, MPS) };
			//real_t probability_misid{totalAmplitudeSquared_DPi_XY(signs.BminusSign, amplitudes_dk_Bminus.A[i], amplitudes_dk_Bminus.Abar[i], amplitudes_dk_Bminus.deltaD[i], correction, MPS) };
			if(probability<1e-6&probability!=0){
				flag_dpi_DD[1][i] = false;
				INFO("Bad point "<<i<<" in B- DPi DD gives probability "<<probability);
				return 0;
			}
			ll_running += events_dpi_DD_Bminus[i].weight_bkg()*(clip_log(probability) - log(normalisation_Bminus));
		}

		//#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_dpi_DD_Bplus.size(); i++)
		{
			if(flag_dpi_DD[0][i]==false) continue;
			real_t correction{ phaseCorrection.eval(events_dpi_DD_Bplus[i]) };
			real_t probability{totalAmplitudeSquared_DPi_XY(signs.BplusSign, amplitudes_DD[2].A[i], amplitudes_DD[2].Abar[i], amplitudes_DD[2].deltaD[i], correction, MPS) };
			//real_t probability_misid{totalAmplitudeSquared_DPi_XY(signs.BplusSign, amplitudes_dk_Bplus.A[i], amplitudes_dk_Bplus.Abar[i], amplitudes_dk_Bplus.deltaD[i], correction, MPS) };
			if(probability<1e-6&probability!=0){
				flag_dpi_DD[0][i] = false;
				INFO("Bad point "<<i<<" in B+ DPi DD gives probability "<<probability);
				return 0;
			}
			ll_running += events_dpi_DD_Bplus[i].weight_bkg()*(clip_log(probability) - log(normalisation_Bplus));
		}

		return -2*dpi_alpha[1]*ll_running;
	};

	//******** BUILD THE TOTAL LL ********
	INFO("built lambda");

	auto LL_total = [&NLL_DK_LL, &NLL_DK_DD, &NLL_DPi_LL, &NLL_DPi_DD](){
		real_t ll_total{0};
		ll_total += NLL_DK_LL();
		ll_total += NLL_DK_DD();
		ll_total += NLL_DPi_LL();
		ll_total += NLL_DPi_DD();
		return ll_total;
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
	
	if (NamedParameter<bool>("DoScan", false))
	{
		INFO("Starting Scan");
		const double xmin = -1;
		const double xmax = 0;
		int step = 1000;
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

		Int_t allow_max = 20;
		Int_t num_bad_points_old = 0;
		for (Int_t j=0; j<allow_max; j++)
		{
			INFO("Starting fit for the "<<j<< "times");
			minimiser.doFit();
			Int_t num_bad_points = 0;
			INFO("The bad points are: ");
			for (size_t i=0; i < events_dk_LL_Bminus.size(); i++)
			{

				if(flag_dk_LL[1][i]==false)
				{
					INFO("B- DK LL: "<<i);
					num_bad_points++;
				}
			}
			for (size_t i=0; i < events_dk_LL_Bplus.size(); i++)
			{
				if(flag_dk_LL[0][i]==false)
				{
					INFO("B+ DK LL: "<<i);
					num_bad_points++;
				}
			}
			for (size_t i=0; i < events_dk_DD_Bminus.size(); i++)
			{
				if(flag_dk_DD[1][i]==false)
				{
					INFO("B- DK DD: "<<i);
					num_bad_points++;
				}
			}
			for (size_t i=0; i < events_dk_DD_Bplus.size(); i++)
			{
				if(flag_dk_DD[0][i]==false)
				{
					INFO("B+ DK DD: "<<i);
					num_bad_points++;
				}
			}
			for (size_t i=0; i < events_dpi_LL_Bminus.size(); i++)
			{
				if(flag_dpi_LL[1][i]==false)
				{
					INFO("B- DPi LL: "<<i);
					num_bad_points++;
				}
			}
			for (size_t i=0; i < events_dpi_LL_Bplus.size(); i++)
			{
				if(flag_dpi_LL[0][i]==false)
				{
					INFO("B+ DPi LL: "<<i);
					num_bad_points++;
				}
			}
			for (size_t i=0; i < events_dpi_DD_Bminus.size(); i++)
			{
				if(flag_dpi_DD[1][i]==false)
				{
					INFO("B- DPi DD: "<<i);
					num_bad_points++;
				}
			}
			for (size_t i=0; i < events_dpi_DD_Bplus.size(); i++)
			{
				if(flag_dpi_DD[0][i]==false)
				{
					INFO("B+ DPi DD: "<<i);
					num_bad_points++;
				}
			}
			if (num_bad_points>0)
			{
				INFO("There are "<<num_bad_points<<" bad points, rerunning fit");
			}
			else if(j != 0)
			{
				INFO("There are no more bad points, continuing with fit");

				break;
			}

			num_bad_points_old = num_bad_points;
			INFO("Final Fits with removal"<< num_bad_points_old<<" bad points");
			minimiser.doFit();
			TIME_fitting.stop();
		}

		// Now save it, note saving only into the logfile while getting things working
		FitResult* output = new FitResult(minimiser); // formerly fr
		output->writeToFile(logFile);
		INFO("Took " << TIME_fitting << "ms to fit to data");
	}


	// *********** PLOT TIME ************
	if (NamedParameter<bool>("DoPlots", true)){
		INFO("Creating projections of fit and data");
		ProfileClock TIME_plots; TIME_plots.start();
		writeUnbinnedFitPlot_LHCb(A_LL, Abar_LL, phaseCorrection, MPS, eventsMC_LL, Datalist_LHCb_LL, "LL");
		writeUnbinnedFitPlot_LHCb(A_DD, Abar_DD, phaseCorrection, MPS, eventsMC_DD, Datalist_LHCb_DD, "DD");
		TIME_plots.stop();
		INFO("Took " << TIME_plots << "ms to make plots and write to file: " << plotFile);
	}
	
		
	return 0;
}








