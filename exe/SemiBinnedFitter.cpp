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


int main(int argc, char * argv[]){

	//******* PARSE THE ARGS FROM COMMAND LINE AND THE OPTIONS FILE*********
	OptionsParser::setArgs( argc, argv );
	
	const std::string binningFile = NamedParameter<std::string>("Binning", "KsPiPi_equal.txt", "text file with list of binning scheme points");

	EventType 		  signalType{NamedParameter<std::string>("EventType", "'D0 K0S0 pi- pi+'", "Signal Type to generate, should be for D0 -> Ks0pi-pi+")}; // notes -> (1) MINUS COMES FIRST 

	const std::string inputData = NamedParameter<std::string>("Input", "data.root", "Root file containing B->D(->Kspipi)K data for both B+/-");

	const std::string logFile = NamedParameter<std::string>("Log", "Fit.log", "Name of the output log file for fit result");

	const std::string matrixFile = NamedParameter<std::string>("Matrix", "matrix.txt", "text file containing the inverse covariance matrix");

	const std::string modelFile = NamedParameter<std::string>("Model", "Kspipi.opt", "Options file with the amplitude model in it");

	const size_t 	  nBins = NamedParameter<size_t>("nBins", 8, "Number of bins in binning file, assumed 8" );

	size_t		 	  nInt = NamedParameter<size_t>("nInt", 1e7, "Number of normalisation events");

	const std::string normalisationAmplitudes = NamedParameter<std::string>("NormalisationAmplitudes", "", "Text file with normalsiation events amplitudes and bins");
	
	const std::string normalisationTuple = NamedParameter<std::string>("NormalisationEvents", "", "Root file with normalsiation events");

	const std::string plotFile = NamedParameter<std::string>("Plots", "Result.root", "Name of the output root file to save fit result plots to");

	const std::string cisiTrue = NamedParameter<std::string>("TrueHadParams", "", "text file with list for 'true' ci and si values");

	const std::string cisiResult = NamedParameter<std::string>("FitHadParams", "", "text file name for output fit ci and si values");

	const size_t      seed = NamedParameter<size_t>("Seed", 0, "Random seed for generation of normalisation events");
	TRandom3 rndm(seed); 

	const size_t      nThreads = NamedParameter<size_t>("nCores", 8, "Number of threads to use" );
	#ifdef _OPENMP
		omp_set_num_threads( nThreads );
		INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
		omp_set_dynamic( 0 );
	#endif


//***********************************************************************************************************************
//**************************************** READ IN ALL INFORMATION NEEDED ***********************************************
//******** CREATE THE MPS ********
	MinuitParameterSet MPS; // this is the mps with the free parameters for fitting
	MPS.loadFromStream();

//******** READ IN THE MODEL *****
	MinuitParameterSet MPS_Kspipi;
	MPS_Kspipi.loadFromFile(modelFile);
	AddCPConjugate(MPS_Kspipi);

//******** READ IN THE BINNING SCHEME *******
	ProfileClock TIME_binningReading;
	TIME_binningReading.start();
	std::vector<real_t> mMinus;
	std::vector<real_t> mPlus;
	std::vector<int> bins;


	readBinning(binningFile, mMinus, mPlus, bins);
	TIME_binningReading.stop();
	INFO("Took " << TIME_binningReading << "ms to read in the binning scheme");


	
//********* READ IN INV COV MATRIX ************
	real_t invCovMatrix[16][16]; //NOTE: doing this with 2*nBins threw an error, would need to change manually nBins not 8
	std::ifstream matrix{matrixFile};
	if (!matrix.good()){WARNING("Unable to read in matrix file " << matrixFile << ", constraint will end up being 0");}
	// empty variables to read the file in to, assuming form of matrix is {ci}{si} i E {1,...,8};
	real_t c1, c2, c3, c4, c5, c6, c7, c8, s1, s2, s3, s4, s5, s6, s7, s8;
	size_t nRow{0}; // to go through the rows of hte matrix
	while ( matrix >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7 >> c8 >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8 ){
		invCovMatrix[nRow][0] = c1;		invCovMatrix[nRow][8] = s1;
		invCovMatrix[nRow][1] = c2;		invCovMatrix[nRow][9] = s2;
		invCovMatrix[nRow][2] = c3;		invCovMatrix[nRow][10] = s3;
		invCovMatrix[nRow][3] = c4;		invCovMatrix[nRow][11] = s4;
		invCovMatrix[nRow][4] = c5;		invCovMatrix[nRow][12] = s5;
		invCovMatrix[nRow][5] = c6;		invCovMatrix[nRow][13] = s6;
		invCovMatrix[nRow][6] = c7;		invCovMatrix[nRow][14] = s7;
		invCovMatrix[nRow][7] = c8;		invCovMatrix[nRow][15] = s8;
		nRow++;
	}
	matrix.close();

//******** READ IN TRUE CI AND SI ************
	std::vector<real_t> trueCi;			std::vector<real_t> trueSi;
	createZeros(trueCi, nBins+1); 		createZeros(trueSi, nBins+1);
	MinuitParameterSet MPS_reference; 
	MPS_reference.loadFromFile(cisiTrue);
	for (size_t i{1}; i < nBins+1; i++){
		trueCi[i] = MPS_reference["c_" + std::to_string(i)]->mean();
		trueSi[i] = MPS_reference["s_" + std::to_string(i)]->mean();
	}


//******* READ IN EVENTS ********
	EventList_type events_Bminus((inputData + ":Bminus__DalitzEventList").c_str() , signalType); 
	EventList_type events_Bplus((inputData + ":Bplus__DalitzEventList").c_str() , signalType); 
	const size_t nEvents{events_Bminus.size()}; // assuming the same for each of B+/-


//********* OPTION: read in tuple made by gammagen:
	// EventList events_Bminus((inputData + ":Bm2DKm").c_str() , signalType); 
	// EventList events_Bplus((inputData + ":Bp2DKp").c_str() , signalType); 
	// const size_t nEvents{events_Bminus.size()}; // assuming the same for each of B+/-
//***

//****** SET UP PHASE CORRECTION
	PhaseCorrection phaseCorrection{MPS};
	phaseCorrection.compilePolynomialExpressions(signalType);

////*****************************************************************************************************************************
////**************************************** SET UP FOR UNBINNED/MD ASPECT OF FIT ***********************************************
//************ READ IN / GENERATE NORMALISATION EVENTS ************
	EventList_type eventsMC(signalType);
	bool readEvents{ false };

	if (normalisationTuple.empty()){
		INFO("No normalisation tuple found, generating " << nInt << " integration events");
		eventsMC =  Generator<>(signalType, &rndm).generate(nInt);
	}else{
		ArgumentPack args;
		eventsMC.loadFromFile( normalisationTuple+":DalitzEventList" , args );
		nInt = eventsMC.size();
		readEvents = true;
	}

//********** CREATE AMPLITUDES ********
	CoherentSum A(signalType, MPS_Kspipi);
	CoherentSum Abar(signalType.conj(true), MPS_Kspipi);

	A.setMC(eventsMC);
	A.prepare();
	Abar.setMC(eventsMC);
	Abar.prepare();

//********* SET UP AMPLITUDE INFORMATION VARIABLES FOR LATER ********
	// |A|^2 / nInt for normalistion:
	real_t normA{A.norm()};
	real_t normAbar{Abar.norm()};

	// |A|, |Abar|, Deltadelta for each of the normalisation events and input events B+/-
	// -> legacy way of doing things from jake because faster to get all at once at this point as these A's never change
	amplitudeInfo amplitudes_Bminus = fillAmplitudeInfo(events_Bminus, A, Abar);
	amplitudeInfo amplitudes_Bplus = fillAmplitudeInfo(events_Bplus, A, Abar);


////************************************************************************************************************************
////**************************************** SET UP FOR BINNED CI/SI ASPECT OF FIT *****************************************
//******* BIN THE EVENTS *******
	ProfileClock TIME_doBinning;
	TIME_doBinning.start();
	std::vector<int> binList_Bminus; // each will have a bin number for each event
	std::vector<int> binList_Bplus;
	int eventBin{0};

	for (auto& event:events_Bminus){
		binList_Bminus.push_back(nearestBinIndex(event, mMinus, mPlus, bins));
	}

	for (auto& event:events_Bplus){
		binList_Bplus.push_back(nearestBinIndex(event, mMinus, mPlus, bins));
	}
	TIME_doBinning.stop();
	INFO("Took " << TIME_doBinning << "ms to bin the events");

	amplitudes_Bminus.bins = binList_Bminus;
	amplitudes_Bplus.bins = binList_Bplus;

////************************************************************************************************************************
//*************************************************** NORMALISATION EVENTS INFORMATION *************************************

//****** READ IN / MAKE AMPLIDUDE INFO AND BIN LISTS
// if no file provided, or needed to generate the amplitude events, must recalculated the ampltidueInfo
	ProfileClock TIME_normalisationInformation; TIME_normalisationInformation.start();
	std::ifstream inFile{normalisationAmplitudes};
	amplitudeInfo amplitudesMC{};
	if (inFile.good() && readEvents){
		// expect text file to be a list of number, one line with the value for one norm. event in order A, Abar, deltaD, binNumber
		double A_1, Abar_2, Dd_3, bin_4;
		while ( inFile >> A_1 >> Abar_2 >> Dd_3 >> bin_4 ){
			amplitudesMC.A.push_back(A_1);
			amplitudesMC.Abar.push_back(Abar_2);
			amplitudesMC.deltaD.push_back(Dd_3);
			amplitudesMC.bins.push_back(bin_4);
		}
	}else{
		if (!inFile.good() && readEvents){WARNING("Unable to open file " << normalisationAmplitudes << ", recalculating amplitudes and binning");}
		else if (inFile.good() && !readEvents){WARNING("Had to regenerate normalisation events, so not using amplitudes from " << normalisationAmplitudes);}		
		amplitudesMC = fillAmplitudeInfo(eventsMC, A, Abar);
		for (auto& event:eventsMC){
			eventBin = nearestBinIndex(event, mMinus, mPlus, bins);
			amplitudesMC.bins.push_back(eventBin);
		}
	}
	inFile.close();
	TIME_normalisationInformation.stop();
	INFO("Took " << TIME_normalisationInformation << "ms to get amplitude info and bin " << nInt << " normalisation events");


//**************** CALCULATE Fi FOR USE LATER AS THEY ARE CONSTANT ******************
// not worth reading in with the other bits as doesn't take long (i think, unchecked so far)
	std::vector<real_t> Fi_posBins, Fbari_posBins, Fi_negBins, Fbari_negBins;
	createZeros(Fi_posBins, nBins+1);	createZeros(Fbari_posBins, nBins+1);
	createZeros(Fi_negBins, nBins+1);	createZeros(Fbari_negBins, nBins+1);

	ProfileClock CLOCK_summingFi;
	CLOCK_summingFi.start();
	for (size_t i=0; i < nInt; i++){
		int binNumber{ amplitudesMC.bins[i] };
		if (binNumber > 0){
			Fi_posBins[binNumber] += std::pow(amplitudesMC.A[i], 2);
			Fbari_posBins[binNumber] += std::pow(amplitudesMC.Abar[i], 2);
		}else{
			int unsignedBin = std::abs(binNumber);
			Fi_negBins[unsignedBin] += std::pow(amplitudesMC.A[i], 2);
			Fbari_negBins[unsignedBin] += std::pow(amplitudesMC.Abar[i], 2);
		}
	}
	CLOCK_summingFi.stop();
	INFO("Took " << CLOCK_summingFi << "ms to sum over the events and create Fi/Fbari");

////************************************************************************************************************************
////**************************************** MAKE FITTING LAMBDAS **********************************************************
	signs signs{}; // so don't need to hard write in -1 and +1 for B+/- and +/-ve bins

//********** GET INTEGRATION TERMS AND THINGS ************
	// pidgeon holes for 2 purposes: ci/si and cross terms in normalisation
	std::vector<real_t> currentCi, currentSi;
	createZeros(currentCi, nBins+1);	createZeros(currentSi, nBins+1);
	std::pair<real_t, real_t> crossTerms;

	auto updateIntegrals = [&currentCi, &currentSi, &crossTerms,
							&Fi_posBins, &Fi_negBins, &Fbari_posBins, &Fbari_negBins,
							&phaseCorrection, &amplitudesMC, &eventsMC,
							&nBins, &nInt]()
	{
		// create empty vectors
		std::vector<real_t> ci_posBins, ci_negBins, si_posBins, si_negBins;
		createZeros(ci_posBins, nBins+1);	createZeros(ci_negBins, nBins+1);
		createZeros(si_posBins, nBins+1);	createZeros(si_negBins, nBins+1);
		// and reset cross terms to zero:
		crossTerms.first = 0; crossTerms.second = 0;

		// do integration
		doBinnedIntegration(ci_posBins, si_posBins, ci_negBins, si_negBins, amplitudesMC, eventsMC, phaseCorrection);
		// normalise ci and si with Fi and put into pidgeon holes
		for (size_t i{1}; i < nBins+1; i++){
			currentCi[i] = 0.5*(ci_posBins[i] / std::pow(Fi_posBins[i]*Fbari_posBins[i] , 0.5 ) + ci_negBins[i] / std::pow( Fi_negBins[i]*Fbari_negBins[i] , 0.5 ));
			currentSi[i] = 0.5*(si_posBins[i] / std::pow(Fi_posBins[i]*Fbari_posBins[i] , 0.5 ) - si_negBins[i] / std::pow( Fi_negBins[i]*Fbari_negBins[i] , 0.5 ));
		// also while doing this sum pre-normalised ci and si for normalisation term
			crossTerms.first += (ci_posBins[i] + ci_negBins[i]) / nInt;
			crossTerms.second += (si_posBins[i] + si_negBins[i]) / nInt;
		}
		return;
	};


// // ******** GAUSSIAN CONSTRAINT **********
	// lambda to calculate current value for a given index on one side of (x-mu)T V (x-mu)
	auto gaussianConstraint = [&invCovMatrix, &trueCi, &trueSi, &currentCi, &currentSi, &nBins]()
	{
		// put each of x_i - x_i_true into a list
		std::vector<real_t> xSubtractMuVector;
		for (size_t i{1}; i < nBins+1; i++){
			xSubtractMuVector.push_back(currentCi[i] - trueCi[i]);
		}
		for (size_t i{1}; i < nBins+1; i++){
			xSubtractMuVector.push_back(currentSi[i] - trueSi[i]);
		} 

		// do the covariance matrix sum to get the exponent term	
		real_t ll{0};
		for (size_t i{0}; i < nBins*2; i++){
			for (size_t j{0}; j < nBins*2; j++){
				ll += xSubtractMuVector[i] * invCovMatrix[i][j] * xSubtractMuVector[j];
			}
		}

		return ll;
	};

	auto LL_total = [&nBins, &signs, &nEvents, &MPS,
					 &events_Bminus, &events_Bplus, &eventsMC,
					 &amplitudesMC, &amplitudes_Bminus, &amplitudes_Bplus, 
					 &normA, &normAbar, &currentCi, &currentSi, &crossTerms,
					 &gaussianConstraint, &updateIntegrals, 
					 &phaseCorrection]()
	{
		phaseCorrection.updateCoeffs(MPS); // update the stored coefficients to latest in MPS

		// ALL INTEGRATION HERE
		updateIntegrals();		
		
		// NORMALISATION NUMBERS:
		real_t normalisation_Bplus = totalAmplitudeSquared_Integrated(signs.BplusSign, normA, normAbar, MPS, crossTerms);
		real_t normalisation_Bminus = totalAmplitudeSquared_Integrated(signs.BminusSign, normA, normAbar, MPS, crossTerms);


		real_t  ll{0};
		#pragma omp parallel for reduction (+:ll)
		for (size_t i=0; i < nEvents; i++)
		{
			// Event event{events_Bminus[i]}
			real_t correction{ phaseCorrection.eval(events_Bminus[i]) };
			real_t probability{ totalAmplitudeSquared_XY(signs.BminusSign, amplitudes_Bminus.A[i], amplitudes_Bminus.Abar[i], amplitudes_Bminus.deltaD[i], correction, MPS) };
			ll += log(probability / normalisation_Bminus);
		}

		#pragma omp parallel for reduction (+:ll)
		for (size_t i=0; i < nEvents; i++)
		{
			real_t correction{ phaseCorrection.eval(events_Bplus[i]) };
			real_t probability{ totalAmplitudeSquared_XY(signs.BplusSign, amplitudes_Bplus.A[i], amplitudes_Bplus.Abar[i], amplitudes_Bplus.deltaD[i], correction, MPS) };
			ll += log(probability / normalisation_Bplus);
		}

		return -2*ll + gaussianConstraint();
	};



	
//******** CHECK HOW LONG LL TAKES TO CALCULATE ********
	if ( NamedParameter<bool>("DoTimeTestOnly", false, "Only do the time test for one LL?")){
		ProfileClock  TIME_oneLL;
		TIME_oneLL.start();
		real_t testLL{LL_total()}; // IGNORE UNUSED VARIABLE WARNING
		TIME_oneLL.stop();
		INFO("Took "<<TIME_oneLL.t_duration << "ms to calculate one LogLikelihood lambda: ");
		INFO("LL has initial value " << testLL << ", where " << gaussianConstraint() << " came from the constraint");
		return 0;
	}
	
//********* DO THE FIT ***********
	Minimiser minimiser(LL_total, &MPS); // formerly mini

	if ( NamedParameter<bool>("DoGradientTest", true, "do the gradient test?")){
		minimiser.gradientTest();
	}

	minimiser.doFit();

	INFO("final value of the gaussian constraint " << gaussianConstraint());

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


	//*********** plot results to root file - put this in WriteToTree?
	// NOTE: 4/8 added extra clause to only let the file draw if fit converged, needed for a specific pull study, uncommitted and don't need to keep in future
	if (output->status() == 0){
		if (NamedParameter<bool>("DoPlots", true, "draw the plots to a root file?")){
			INFO("Creating projections of fit and data");
			ProfileClock TIME_plots; TIME_plots.start();
			writeUnbinnedFitPlot(A, Abar, phaseCorrection, MPS, eventsMC, events_Bplus, events_Bminus);
			TIME_plots.stop();
			INFO("Took " << TIME_plots << "ms to make plots and write to file: " << plotFile);
		}
	}


	// save the ci and si values from the current fit - also put in the writeToTree header?
	if (NamedParameter<bool>("DoCiSiCalculation", true, "save the end ci and si values?")){
		std::ofstream outFile{cisiResult};
		for (int i{1}; i < nBins+1; i++){
			outFile << "c_" << i << " " << currentCi[i] << "\n"
					<< "s_" << i << " " << currentSi[i] << "\n";
		}
		outFile.close();
	}



	std::cout << "\n bottom of code :D" << std::endl;
	return 0;
}
