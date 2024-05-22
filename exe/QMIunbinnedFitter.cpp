#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/Binning.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventListSIMD.h"
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
#include "AmpGen/QMI.h"


#if ENABLE_AVX
	using EventList_type = AmpGen::EventListSIMD;
#else
	using EventList_type = AmpGen::EventList; 
#endif

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"



using namespace AmpGen;

// using CE=CompiledExpression<real_t(const real_t*, const real_t*)> ;
// need a way to make a vector of 0's
template <typename T> void createZeros(std::vector<T>& list, const size_t& nItems){
	for (size_t i{0}; i < nItems;i++){
		list.push_back(0);
	}
	return;
}

// workaround function so can keep all the event types together rather than having to load separately
void readinEventType(EventType& type, const char* string){
	EventType thisType{ NamedParameter<std::string>(string, " ") };
	type = thisType;
}


int main(int argc, char * argv[]){

	//******* PARSE THE ARGS FROM COMMAND LINE AND THE OPTIONS FILE*********
	OptionsParser::setArgs( argc, argv );
	
	EventType 		  signalType{NamedParameter<std::string>("EventType", "'D0 K0S0 pi- pi+'", "Signal Type to generate, should be for D0 -> Ks0pi-pi+")}; // notes -> (1) MINUS COMES FIRST 

	const std::string inputData_LHCb = NamedParameter<std::string>("InputLHCb", "data.root", "Root file containing B->D(->Kspipi)K data for both B+/-");

	const std::string inputData_BES = NamedParameter<std::string>("InputBES", "data.root", "Root file containing B->D(->Kspipi)K data for both B+/-");

	const std::string logFile = NamedParameter<std::string>("Log", "Fit.log", "Name of the output log file for fit result");

	const std::string modelFile = NamedParameter<std::string>("Model", "Kspipi.opt", "Options file with the amplitude model in it");

	size_t		 	  nInt = NamedParameter<size_t>("nInt", 1e7, "Number of normalisation events");

	const std::string normalisationAmplitudes = NamedParameter<std::string>("NormalisationAmplitudes", "", "Text file with normalsiation events amplitudes and bins");
	
	const std::string normalisationTuple = NamedParameter<std::string>("NormalisationEvents", "", "Root file with normalsiation events");

	const std::string plotFile = NamedParameter<std::string>("Plots", "Result.root", "Name of the output root file to save fit result plots to");

	const size_t      seed = NamedParameter<size_t>("Seed", 0, "Random seed for generation of normalisation events");


	//********Making structure for bkg fration*******
	auto bkg_fractions = NamedParameter<std::string>("Bkg_frac").getVector();
	auto eff_normalisation = NamedParameter<std::string>("Norm_eff").getVector();

	TRandom3 rndm(seed); 

	const size_t      nThreads = NamedParameter<size_t>("nCores", 8, "Number of threads to use" );
	#ifdef _OPENMP
		omp_set_num_threads( nThreads );
		INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
		omp_set_dynamic( 0 );
	#endif


//***********************************************************************************************************************
//*********************************************** READ IN / SET UP ******************************************************
//***********************************************************************************************************************

//******** CREATE THE BASE MPS ********
	MinuitParameterSet MPS; // this is the mps with the free parameters for fitting
	MPS.loadFromStream();

//****** SET UP PHASE CORRECTION
	PhaseCorrection phaseCorrection{MPS};
	phaseCorrection.compilePolynomialExpressions(signalType);


	std::vector<std::string> branches;
	branches.push_back("Bac_ID");

//******* READ IN EVENTS ********
    // 1/2: LHCb style
	EventList_type events_Bminus((inputData_LHCb + ":Bminus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw")); //WeightBranch_eff("weight_eff"),
	EventList_type events_Bplus((inputData_LHCb + ":Bplus_DalitzEventList").c_str() , signalType, WeightBranch_bkg("sig_sw")); // WeightBranch_eff("weight_eff"),
	const size_t nEvents_LHCb{events_Bminus.size()}; // assuming the same for each of B+/-
	INFO("LHCb  sWeight[0]= " << events_Bminus[0].weight_bkg());
    // 2/2: BES style:
    EventList_type events_CPodd((inputData_BES + ":CPodd__DalitzEventList").c_str(), signalType, WeightBranch_bkg("weight_bkg"), WeightBranch_eff("weight_eff"));
    EventList_type events_CPeven((inputData_BES + ":CPeven__DalitzEventList").c_str(), signalType, WeightBranch_bkg("weight_bkg"), WeightBranch_eff("weight_eff"));
	EventList_type events_same((inputData_BES + ":same_signal__DalitzEventList").c_str(), signalType, WeightBranch_bkg("weight_bkg"), WeightBranch_eff("weight_eff"));
	EventList_type events_tag((inputData_BES + ":same_tag__DalitzEventList").c_str(), signalType, WeightBranch_bkg("weight_bkg"), WeightBranch_eff("weight_eff"));// only events that we need the 'tagged' data for Kspipi
    EventList_type events_flavour((inputData_BES + ":flavour__DalitzEventList").c_str(), signalType, WeightBranch_bkg("weight_bkg"), WeightBranch_eff("weight_eff")); // not currently using
    EventList_type events_flavourBar((inputData_BES + ":flavourBar__DalitzEventList").c_str(), signalType, WeightBranch_bkg("weight_bkg"), WeightBranch_eff("weight_eff")); // not currently using
//    EventList_type events_same("/publicfs2/ucas/user/zengshh/LHCb/B2DK_D2K0S0PiPi/analysis/new_scripts/with_sim/besiii/data/same_weighted_signal.root:same_signal__DalitzEventList", signalType,WeightBranch_eff("weight_eff"), WeightBranch_bkg("weight_bkg"));
//    EventList_type events_tag("/publicfs2/ucas/user/zengshh/LHCb/B2DK_D2K0S0PiPi/analysis/new_scripts/with_sim/besiii/data/same_weighted_signal.root:same_tag__DalitzEventList", signalType,WeightBranch_eff("weight_eff"), WeightBranch_bkg("weight_bkg")); // only events that we need the 'tagged' data for Kspipi
	// 2/2 BES style: - TEMP: QCgen NAMES
	// EventList_type events_CPodd((inputData_BES + ":Signal_KK").c_str(), signalType);
    // EventList_type events_CPeven((inputData_BES + ":Signal_Kspi0").c_str(), signalType);
    // EventList_type events_flavour((inputData_BES + ":Signal_Kppim").c_str(), signalType);
    // EventList_type events_flavourBar((inputData_BES + ":Signal_Kmpip").c_str(), signalType);
    // EventList_type events_same((inputData_BES + ":Signal_Kspipi").c_str(), signalType);
    // EventList_type events_tag((inputData_BES + ":Tag_Kspipi").c_str(), signalType); // only events that we need the 'tagged' data for Kspipi

//---------GET BKGFRACTION-----------//
	bkgfraction bkg_frac;
	for (auto tag : bkg_fractions){
   	auto a = split(tag, ' ');
	double frac_temp =std::stod(a[1]);
	bkg_frac.frac.push_back(frac_temp);
	bkg_frac.tag.push_back(a[0]);
	}

	double *bkg_frac_val = get_bkg_fraction(bkg_frac);

	for(int i=0;i<7;i++){
		INFO("bkg_frac_val["<<i<<"] = "<<bkg_frac_val[i]);
	}

//--------GET EFFNORM----------------//
	normeff norm_eff;
	for (auto tag : eff_normalisation){
   	auto a = split(tag, ' ');
	double frac_temp =std::stod(a[1]);
	norm_eff.frac.push_back(frac_temp);
	norm_eff.tag.push_back(a[0]);
	}

	double *norm_eff_val = get_norm_eff(norm_eff);

	for(int i=0;i<7;i++){
		INFO("norm_eff_val["<<i<<"] = "<<norm_eff_val[i]);
	}
	

//********SAVE IN LIST **********

	std::vector<EventList_type> Datalist_LHCb{events_Bminus, events_Bplus};
	std::vector<EventList_type> Datalist_BESIII{events_flavour,events_flavourBar, events_CPodd, events_CPeven, events_same, events_tag};

//************ READ IN / GENERATE NORMALISATION EVENTS ************
	//EventList_type eventsMC(signalType);
	EventList_type eventsMC_CPodd(signalType);
	EventList_type eventsMC_CPeven(signalType);
	EventList_type eventsMC_same_signal(signalType);
	EventList_type eventsMC_same_tag(signalType);
	EventList_type eventsMC_flavour(signalType);
	EventList_type eventsMC_flavourBar(signalType);
	EventList_type eventsMC_Bplus(signalType);
	EventList_type eventsMC_Bminus(signalType);

	std::vector<EventList_type> eventsMC{eventsMC_CPodd, eventsMC_CPeven, eventsMC_same_signal, eventsMC_same_tag, eventsMC_flavour, eventsMC_flavourBar};
	int nInt_i[6] = {0,0,0,0,0,0};
	bool readEvents{ false };
	if (normalisationTuple.empty()){
		INFO("No normalisation tuple found, generating " << nInt << " integration events");
		for(auto eventsMC_i: eventsMC){
			eventsMC_i =  Generator<>(signalType, &rndm).generate(nInt);
		}
	}else{
		ArgumentPack args;//(WeightBranch_eff("weight_eff"));
		eventsMC[0].loadFromFile( normalisationTuple+":CPodd__DalitzEventList", args);
		eventsMC[1].loadFromFile( normalisationTuple+":CPeven__DalitzEventList", args);
		eventsMC[2].loadFromFile( normalisationTuple+":same_signal__DalitzEventList", args);
		eventsMC[3].loadFromFile( normalisationTuple+":same_tag__DalitzEventList", args);
		eventsMC[4].loadFromFile( normalisationTuple+":flavour__DalitzEventList", args);
		eventsMC[5].loadFromFile( normalisationTuple+":flavourBar__DalitzEventList", args);
		for(int i=0;i<6;i++){
			nInt_i[i] = eventsMC[i].size();
		}
		readEvents = true;
	}
	//temp for LHCb dataset
	eventsMC_Bplus =  Generator<>(signalType, &rndm).generate(nInt);

//********** READ IN MODEL / CREATE Kspipi AMPLITUDES ********
	// - model:
	MinuitParameterSet MPS_Kspipi;
	MPS_Kspipi.loadFromFile(modelFile);
	AddCPConjugate(MPS_Kspipi);
	
	// - Cohenrent sums A for the Kspipi amplitude
	CoherentSum A(signalType, MPS_Kspipi);
	CoherentSum Abar(signalType.conj(true), MPS_Kspipi);

	for(int i=0;i<6;i++){
		A.setMC(eventsMC[i]);  	A.prepare();
		Abar.setMC(eventsMC[i]);	Abar.prepare();
	}

	// |A|, |Abar|, Deltadelta for each of the normalisation events and input events B+/-
	// -> legacy way of doing things from jake because faster to get all at once at this point as these A's never change
    // LHCb section:
	amplitudeInfo amplitudes_Bminus = fillAmplitudeInfo(events_Bminus, A, Abar);
	amplitudeInfo amplitudes_Bplus = fillAmplitudeInfo(events_Bplus, A, Abar);

	amplitudeInfo amplitudes_flavour = fillAmplitudeInfo(events_flavour, A, Abar); // not currently using
	amplitudeInfo amplitudes_flavourBar = fillAmplitudeInfo(events_flavourBar, A, Abar); // as above 
	amplitudeInfo amplitudes_same = fillAmplitudeInfo(events_same, A, Abar);
	amplitudeInfo amplitudes_tag = fillAmplitudeInfo(events_tag, A, Abar);
	    // BES section:
	amplitudeInfo amplitudes_CPodd = fillAmplitudeInfo(events_CPodd, A, Abar);
	amplitudeInfo amplitudes_CPeven = fillAmplitudeInfo(events_CPeven, A, Abar);

//******* NORMALISATION EVENTS INFORMATION ***********

// **** READ IN / MAKE AMPLIDUDE INFO AND BIN LISTS
// if no file provided, or needed to generate the amplitude events, must recalculated the ampltidueInfo
	ProfileClock TIME_normalisationInformation; TIME_normalisationInformation.start();
	std::ifstream inFile{normalisationAmplitudes};
	amplitudeInfo amplitudesMC{};
	amplitudeInfo amplitudesMC_CPodd{}, amplitudesMC_CPeven{}, amplitudesMC_same_signal{}, amplitudesMC_same_tag{}, amplitudesMC_flavour{}, amplitudesMC_flavourBar{};
	std::vector<amplitudeInfo> eventsMC_amplitudes{amplitudesMC_CPodd, amplitudesMC_CPeven, amplitudesMC_same_signal, amplitudesMC_same_tag, amplitudesMC_flavour, amplitudesMC_flavourBar};
	if (inFile.good() && readEvents){
		// expect text file to be a list of number, one line with the value for one norm. event in order A, Abar, deltaD, binNumber
		double A_1, Abar_2, Dd_3, bin_4;
		while ( inFile >> A_1 >> Abar_2 >> Dd_3 >> bin_4 ){
			amplitudesMC.A.push_back(A_1);
			amplitudesMC.Abar.push_back(Abar_2);
			amplitudesMC.deltaD.push_back(Dd_3);
		}
	}else{
		if (!inFile.good() && readEvents){WARNING("Unable to open file " << normalisationAmplitudes << ", recalculating amplitudes and binning");}
		else if (inFile.good() && !readEvents){WARNING("Had to regenerate normalisation events, so not using amplitudes from " << normalisationAmplitudes);}		
		amplitudesMC = fillAmplitudeInfo(eventsMC_Bplus, A, Abar);
		for(int i=0;i<6;i++){
			eventsMC_amplitudes[i] = fillAmplitudeInfo(eventsMC[i], A, Abar);
		}
	}
	inFile.close();
	TIME_normalisationInformation.stop();
	INFO("Took " << TIME_normalisationInformation << "ms to get amplitude info for " << nInt << " normalisation events");

//***********************************************************************************************************************
//*********************************************** BUILD LIKELIHOODS *****************************************************
//***********************************************************************************************************************
    signs signs{};
	// variables to hold the various integration terms needed
	
	real_t normA{A.norm()};
	real_t normAbar{Abar.norm()};
	real_t integration_cosTerm, integration_sinTerm, integration_kspipiTerm;
	// - create the information we can outside of looping in the fit:
	real_t normA_flavour, normAbar_flavour, normA_flavourBar, normAbar_flavourBar, normA_CPodd, normAbar_CPodd, normA_CPeven, normAbar_CPeven, normA_same, normAbar_same, normA_tag, normAbar_tag;
	real_t integration_cosTerm_flavour, integration_sinTerm_flavour, integration_cosTerm_flavourBar, integration_sinTerm_flavourBar, integration_cosTerm_CPodd, integration_sinTerm_CPodd, integration_cosTerm_CPeven, integration_sinTerm_CPeven, integration_cosTerm_same, integration_sinTerm_same, integration_cosTerm_tag, integration_sinTerm_tag;
	real_t integration_kspipiTerm_flavour, integration_kspipiTerm_flavourBar, integration_kspipiTerm_CPodd, integration_kspipiTerm_CPeven, integration_kspipiTerm_same, integration_kspipiTerm_tag;

	std::vector<real_t> normA_i{normA_CPodd, normA_CPeven, normA_same, normA_tag, normA_flavour, normA_flavourBar};
	std::vector<real_t> normAbar_i{normAbar_CPodd, normAbar_CPeven, normAbar_same, normAbar_tag, normAbar_flavour, normAbar_flavourBar};
	std::vector<real_t> integration_cosTerm_i{integration_cosTerm_CPodd, integration_cosTerm_CPeven, integration_cosTerm_same, integration_cosTerm_tag, integration_cosTerm_flavour, integration_cosTerm_flavourBar};
	std::vector<real_t> integration_sinTerm_i{integration_sinTerm_CPodd, integration_sinTerm_CPeven, integration_sinTerm_same, integration_sinTerm_tag, integration_sinTerm_flavour, integration_sinTerm_flavourBar};
	std::vector<real_t> integration_kspipiTerm_i{integration_kspipiTerm_CPodd, integration_kspipiTerm_CPeven, integration_kspipiTerm_same, integration_kspipiTerm_tag, integration_kspipiTerm_flavour, integration_kspipiTerm_flavourBar};

	auto updateIntegrals = [&phaseCorrection, &eventsMC_amplitudes, &amplitudesMC, &eventsMC_Bplus, &eventsMC, &integration_cosTerm, &integration_sinTerm, &integration_kspipiTerm, &nInt, &normA_i, &normAbar_i, &integration_cosTerm_i, &integration_sinTerm_i, &integration_kspipiTerm_i, &nInt_i ]()
	{
		integration_cosTerm = 0; integration_sinTerm = 0; integration_kspipiTerm = 0;
		doUnbinnedIntegration(amplitudesMC, eventsMC_Bplus, phaseCorrection, integration_cosTerm, integration_sinTerm, integration_kspipiTerm);
		// divide by number of events
		integration_cosTerm = integration_cosTerm / nInt;
		integration_sinTerm = integration_sinTerm / nInt;
		integration_kspipiTerm = integration_kspipiTerm / nInt;
		for(int i=0;i<6;i++){
		doUnbinnedIntegration_temp(eventsMC_amplitudes[i], eventsMC[i], phaseCorrection, normA_i[i], normAbar_i[i], integration_cosTerm_i[i], integration_sinTerm_i[i], integration_kspipiTerm_i[i]);
			normA_i[i] = normA_i[i] / nInt_i[i];
			normAbar_i[i] = normAbar_i[i] / nInt_i[i];
			integration_cosTerm_i[i] = integration_cosTerm_i[i] / nInt_i[i];
			integration_sinTerm_i[i] =  integration_sinTerm_i[i] / nInt_i[i];
			integration_kspipiTerm_i[i] =  integration_kspipiTerm_i[i] / nInt_i[i];
		}
		return;

	};

	// not currently doing anything as without background considerations this pdf is always constant
     auto LL_flavour = [&events_flavour, &events_flavourBar, &amplitudes_flavour, &amplitudes_flavourBar, &normA_i, &normAbar_i, &bkg_frac_val, &norm_eff_val]()
     {
         real_t normalisation_flavour{ normAbar_i[4] }; // because the flavour = the tag state being D, and the overall amplitude is Abar for the signal being Dbar
         real_t normalisation_flavourBar{ normA_i[4] }; // as above: tag in Dbar state so signal in D state
         real_t  ll_running{0};
	 	#pragma omp parallel for reduction (+:ll_running)
	 	for (size_t i=0; i < events_flavour.size(); i++)
	 	{
	 		real_t probability{ totalAmplitudeSquared_flavour(amplitudes_flavour.Abar[i]) };
			real_t probability_bkg {bkg_frac_val[0]*events_flavour[i].weight_bkg()};
	 		ll_running += log((1-bkg_frac_val[0])*( probability / (normalisation_flavour)) + probability_bkg);
	 	}
	 	#pragma omp parallel for reduction (+:ll_running)
	 	for (size_t i=0; i < events_flavourBar.size(); i++)
	 	{
             real_t probability{ totalAmplitudeSquared_flavour(amplitudes_flavourBar.A[i]) };
			 real_t probability_bkg {bkg_frac_val[1]*events_flavourBar[i].weight_bkg()};
	 		ll_running += log((1-bkg_frac_val[1])*( probability / ( normalisation_flavourBar)) + probability_bkg);
	 	}
		     return ll_running;
     };
	
    auto LL_CP = [&events_CPeven, &events_CPodd, &eventsMC, &amplitudes_CPeven, &amplitudes_CPodd, 
				  &phaseCorrection, &signs, &normA_i, &normAbar_i, &integration_cosTerm_i, &bkg_frac_val, &norm_eff_val]()
    {
        real_t normalisation_odd{ (normA_i[0] + normAbar_i[0] + 2*integration_cosTerm_i[0]) };
        real_t normalisation_even{ (normA_i[1] + normAbar_i[1] - 2*integration_cosTerm_i[1]) };

        real_t  ll_even{0};
		
		#pragma omp parallel for reduction (+:ll_even)
		for (size_t i=0; i < events_CPeven.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_CPeven[i]) };
			real_t probability{ totalAmplitudeSquared_CP(signs.CPevenSign, amplitudes_CPeven.A[i], amplitudes_CPeven.Abar[i], amplitudes_CPeven.deltaD[i], correction) };
			real_t probability_bkg {bkg_frac_val[3]*events_CPeven[i].weight_bkg()};
			ll_even += log( (1-bkg_frac_val[3])*(probability / (normalisation_even)) + probability_bkg);
		}

		real_t ll_odd{0};
		
		#pragma omp parallel for reduction (+:ll_odd)
		for (size_t i=0; i < events_CPodd.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_CPodd[i]) };
			real_t probability{ totalAmplitudeSquared_CP(signs.CPoddSign, amplitudes_CPodd.A[i], amplitudes_CPodd.Abar[i], amplitudes_CPodd.deltaD[i], correction) };
			real_t probability_bkg {bkg_frac_val[2]*events_CPodd[i].weight_bkg()};
			ll_odd += log((1-bkg_frac_val[2]) *(probability / ( normalisation_odd))+ probability_bkg);
		}
        return ll_even + ll_odd;
    };
    auto LL_same = [&events_same, &events_tag, &eventsMC, &amplitudes_same, &amplitudes_tag, &phaseCorrection, &integration_kspipiTerm_i, &bkg_frac_val, &norm_eff_val]()
    {
        real_t normalisation{ integration_kspipiTerm_i[2] };
        real_t  ll_running{0};
		#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_same.size(); i++)
		{
            real_t deltaC{ phaseCorrection.eval(events_same[i]) };
            real_t deltaC_tag{ phaseCorrection.eval(events_tag[i]) };
			real_t probability{ totalAmplitudeSquared_BES(amplitudes_same.A[i], amplitudes_same.Abar[i], amplitudes_tag.A[i], amplitudes_tag.Abar[i], amplitudes_same.deltaD[i], amplitudes_tag.deltaD[i], deltaC, deltaC_tag) };
			real_t probability_bkg {bkg_frac_val[4]*events_same[i].weight_bkg()};
			ll_running += log( (1-bkg_frac_val[4]) *(probability / (normalisation)) + probability_bkg);
		}
        return ll_running;
    };


    auto LL_LHCb = [&events_Bminus, &events_Bplus, &eventsMC, 
					&amplitudesMC, &amplitudes_Bminus, &amplitudes_Bplus, &nEvents_LHCb,
					&normA, &normAbar, &integration_cosTerm, &integration_sinTerm, &MPS,
					&phaseCorrection, &signs, &bkg_frac_val](){
		
		std::pair<real_t, real_t> normalisationCrossTerms = {integration_cosTerm, integration_sinTerm}; // to be able to use the functions from SBF notation
		real_t normalisation_Bplus = totalAmplitudeSquared_Integrated(signs.BplusSign, normA, normAbar, MPS, normalisationCrossTerms) ;
		real_t normalisation_Bminus = totalAmplitudeSquared_Integrated(signs.BminusSign, normA, normAbar, MPS, normalisationCrossTerms) ;

		real_t  ll_running{0};
		#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_Bminus.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_Bminus[i]) };
			real_t probability{ totalAmplitudeSquared_XY(signs.BminusSign, amplitudes_Bminus.A[i], amplitudes_Bminus.Abar[i], amplitudes_Bminus.deltaD[i], correction, MPS)};
			ll_running += events_Bminus[i].weight_bkg()*log(probability / normalisation_Bminus);

		}

		#pragma omp parallel for reduction (+:ll_running)
		for (size_t i=0; i < events_Bplus.size(); i++)
		{
			real_t correction{ phaseCorrection.eval(events_Bplus[i]) };
			real_t probability{ totalAmplitudeSquared_XY(signs.BplusSign, amplitudes_Bplus.A[i], amplitudes_Bplus.Abar[i], amplitudes_Bplus.deltaD[i], correction, MPS)};
			ll_running += events_Bplus[i].weight_bkg()*log(probability / normalisation_Bplus);

		}
		return ll_running;
	};
	//******get sWeight normalisation
	auto sWeight_norm = [&events_Bminus, &events_Bplus](){
		real_t alpha{0};
		real_t  sw{0};
		real_t sw_sq{0};
		for (size_t i=0; i < events_Bminus.size(); i++)
		{
			sw += events_Bminus[i].weight_bkg();
			sw_sq += events_Bminus[i].weight_bkg()*events_Bminus[i].weight_bkg();
		}
		for (size_t i=0; i < events_Bplus.size(); i++)
		{
			sw += events_Bplus[i].weight_bkg();
			sw_sq += events_Bplus[i].weight_bkg()*events_Bplus[i].weight_bkg();
		}
		alpha = sw/(sw_sq);
		return alpha;
	};



     auto LL_total = [&LL_CP, &LL_flavour, &LL_same, &LL_LHCb, &phaseCorrection, &updateIntegrals, &sWeight_norm, &MPS]()
    {
        phaseCorrection.updateCoeffs(MPS);
		updateIntegrals();
        // return -2* (LL_CP() + LL_flavour() + LL_same() + LL_LHCb());
		return -2* (sWeight_norm()*LL_LHCb() + LL_same() + LL_CP() + LL_flavour());
        //return -2* (LL_CP() + LL_same() + LL_LHCb());

    };


//***********************************************************************************************************************
//*********************************************** DO FIT AND SAVE THINGS ************************************************
//***********************************************************************************************************************






//******** CHECK HOW LONG LL TAKES TO CALCULATE ********
	if ( NamedParameter<bool>("DoTimeTestOnly", false, "Only do the time test for one LL?")){
		ProfileClock  TIME_oneLL;
		TIME_oneLL.start();
		INFO(">> about to call LL_total");
		real_t testLL{LL_total()};
		TIME_oneLL.stop();
		INFO("Took "<<TIME_oneLL.t_duration << "ms to calculate one LogLikelihood lambda: " << testLL);
		return 0;
	}
	
//********* DO THE FIT ***********
	Minimiser minimiser(LL_total, &MPS); // formerly mini

	if ( NamedParameter<bool>("DoGradientTest", true, "do the gradient test?")){
		minimiser.gradientTest();
	}

	minimiser.doFit();

	//********* and minos errors
	if (NamedParameter<bool>("DoMINOS", false)){
		ProfileClock TIME_minos;
		TIME_minos.start();
		minimiser.minos(MPS["xPlus"]);
		minimiser.minos(MPS["xMinus"]);
		minimiser.minos(MPS["yPlus"]);
		minimiser.minos(MPS["yMinus"]);
		// minimiser.minos(MPS["rB"]);
		// minimiser.minos(MPS["dB"]);
		// minimiser.minos(MPS["gamma"]);
		TIME_minos.stop();
		INFO("Took " << TIME_minos.t_duration << "ms to calculate MINOS uncertainties");
	}

//************** SAVE THINGS ABOUT THE FIT TO FILE ***********
	// ********* Now save it
	FitResult* output = new FitResult(minimiser); // formerly fr
	INFO("writing result to logfile " << logFile);
	output->writeToFile(logFile);
	auto my_pc = [&phaseCorrection](Event& event){
				return phaseCorrection.eval(event) ;
    };
	
	//*********** plot results to root file - put this in WriteToTree?
	if (NamedParameter<bool>("DoPlots", true, "draw the plots to a root file?")){
		INFO("Creating projections of fit and data");
		ProfileClock TIME_plots; TIME_plots.start();
		//writeUnbinnedFitPlot_LHCb(A, Abar, phaseCorrection, MPS, eventsMC, Datalist_LHCb);
		writeUnbinnedFitPlot_BES(A, Abar, phaseCorrection, MPS, eventsMC, Datalist_BESIII);
		TIME_plots.stop();
		INFO("Took " << TIME_plots << "ms to make plots and write to file: " << plotFile);
	}


	std::cout << "\n bottom of code :D" << std::endl;

	return 0;
}