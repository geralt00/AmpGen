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
#include "AmpGen/NamedParameter.h"
// #include "AmpGen/PhaseCorrection.h" //NOTE - not used yet as contained within ci and si
#include "AmpGen/ProfileClock.h"
#include "AmpGen/MyStructs.h"
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

	EventType signalType{NamedParameter<std::string>("EventType", "", "Signal Type to generate, should be for D0 -> Ks0pi-pi+")}; // notes -> (1) MINUS COMES FIRST 

	const std::string fitType = NamedParameter<std::string>("Fitter", "poisson", "poisson or Chi2 minimiser to be used in fit");
	if (fitType != "poisson" && fitType != "chi2"){
		WARNING("Fit type '" << fitType << "' unknown, set to default poisson. Please specify poisson or chi2.");
	}

	const std::string inputData = NamedParameter<std::string>("Input", "data.root", "Root file containing B->D(->Kspipi)K data for both B+/-");

	const std::string logFile = NamedParameter<std::string>("Log", "Fit.log", "Name of the output log file for fit result");

	const std::string matrixFile = NamedParameter<std::string>("Matrix", "matrix.txt", "text file containing the inverse covariance matrix");

	const std::string plotFile = NamedParameter<std::string>("Plots", "Result.root", "Name of the output root file to save fit result plots to");

	const std::string resultFile = NamedParameter<std::string>("ResultFile", "", "text file for list of results of bin counts");

	const std::string cisiTrue = NamedParameter<std::string>("TrueHadParams", "", "text file with list for 'true' ci and si values");

	const size_t      nThreads = NamedParameter<size_t>("nCores", 8, "Number of threads to use" );
	#ifdef _OPENMP
		omp_set_num_threads( nThreads );
		INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
		omp_set_dynamic( 0 );
	#endif

//******** CREATE THE MPS********
	MinuitParameterSet MPS;
	MPS.loadFromStream();

// //******** READ IN THE BINNING SCHEME *******
	ProfileClock TIME_binningReading;
	TIME_binningReading.start();
	std::vector<real_t> mMinus;
	std::vector<real_t> mPlus;
	std::vector<int> bins;

	const size_t nBins = NamedParameter<size_t>("nBins", 8, "Number of bins in binning file, assumed 8" );

	readBinning(binningFile, mMinus, mPlus, bins);
	TIME_binningReading.stop();
	INFO("Took " << TIME_binningReading << "ms to read in the binning scheme");
	

// // ******* READ IN EVENTS ********
	EventList data_Bminus((inputData + ":Bminus__DalitzEventList").c_str() , signalType); 
	EventList data_Bplus((inputData + ":Bplus__DalitzEventList").c_str() , signalType); 
	const size_t nEvents{data_Bminus.size()}; // assuming the same for each of B+/-

// //********* TEMP: read in tuple made by gammagen:
// 	// EventList data_Bminus((inputData + ":Bm2DKm").c_str() , signalType); 
// 	// EventList data_Bplus((inputData + ":Bp2DKp").c_str() , signalType); 
// 	// const size_t nEvents{data_Bminus.size()}; // assuming the same for each of B+/-

// //******* BIN THE EVENTS *******
	ProfileClock TIME_doBinning;
	TIME_doBinning.start();
	std::vector<int> binList_Bminus; // each will have a bin number for each event
	std::vector<int> binList_Bplus;
	int eventBin{0};

	for (auto& event:data_Bminus){
		eventBin = nearestBinIndex(event, mMinus, mPlus, bins);
		binList_Bminus.push_back(eventBin);
	}

	for (auto& event:data_Bplus){
		eventBin = nearestBinIndex(event, mMinus, mPlus, bins);
		binList_Bplus.push_back(eventBin);
	}
	TIME_doBinning.stop();
	INFO("Took " << TIME_doBinning << "ms to bin the events");


// //******* COUNT NUMBER OF EVENTS IN EACH BIN *******
// 	// initialise some empty vectors to take no events -> nBins+1 entries, entry 0 should be 0
// 	INFO("Counting bin occupancy");
	std::vector<int> binCounts_posBins_Bminus, binCounts_negBins_Bminus;
	createZeros<int>(binCounts_posBins_Bminus, nBins+1);
	createZeros<int>(binCounts_negBins_Bminus, nBins+1);
	std::vector<int> binCounts_posBins_Bplus, binCounts_negBins_Bplus;
	createZeros<int>(binCounts_posBins_Bplus, nBins+1);
	createZeros<int>(binCounts_negBins_Bplus, nBins+1);
	int currentBin{0};

	// loop over each event's bin number
	for (size_t i{0}; i < nEvents; i++){
		// B- events
		currentBin = binList_Bminus[i];
		if (currentBin > 0){
			binCounts_posBins_Bminus[currentBin]++;
		}else{
			binCounts_negBins_Bminus[abs(currentBin)]++;
		}
		//B+ events
		currentBin = binList_Bplus[i];
		if (currentBin > 0){
			binCounts_posBins_Bplus[currentBin]++;
		}else{
			binCounts_negBins_Bplus[abs(currentBin)]++;
		}
	}


// ******** BUILD THE LOG LIKELIHOOD LAMBDA **********
	INFO("Building the LL");

	signs signs{}; // so don't need to hard write in -1 and +1 for B+/- and +/-ve bins

	std::vector<real_t> currentCi, currentSi, currentFi, currentFbari; // do we really need to do this with F too??
	createZeros(currentCi, nBins+1);	createZeros(currentSi, nBins+1);
	createZeros(currentFi, nBins+1);	createZeros(currentFbari, nBins+1);
	auto retrieveAllMPS = [ &MPS, &currentCi, &currentSi, &currentFi, &currentFbari , &nBins]()
	{
		for (size_t i{1}; i < nBins+1; i++){
			currentCi[i] = MPS["c_" + std::to_string(i)]->mean(); 
			currentSi[i] = MPS["s_" + std::to_string(i)]->mean();
			currentFi[i] = MPS["F_" + std::to_string(i)]->mean();
			currentFbari[i] = MPS["Fbar_" + std::to_string(i)]->mean();
		}
		return;
	};


// ******** gaussian constraint on ci and si *********

	// read in inverse covariance matrix
	real_t invCovMatrix[16][16]; //NOTE: doing this with 2*nBins threw an error, would need to change manually nBins not 8
	std::ifstream inFile{matrixFile};
	if (!inFile.good()){WARNING("unable to read in matrix file " << matrixFile << ", there will be no constraint in this fit");}
	// empty variables to read the file in to, assuming form of matrix is {ci}{si} i E {1,...,8};
	real_t c1, c2, c3, c4, c5, c6, c7, c8, s1, s2, s3, s4, s5, s6, s7, s8;
	size_t nRow{0}; // to go through the rows of hte matrix
	while ( inFile >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7 >> c8 >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8 ){
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


	// read in 'true' values for ci and si - should have a file in a named parameter containing them
	// note they may well be the same as the start values for the Free ci and si in the fit but this way they don't have to be
	std::vector<real_t> trueCi;			std::vector<real_t> trueSi;
	createZeros(trueCi, nBins+1); 		createZeros(trueSi, nBins+1);
	MinuitParameterSet MPS_reference; 
	MPS_reference.loadFromFile(cisiTrue);
	for (size_t i{1}; i < nBins+1; i++){
		trueCi[i] = MPS_reference["c_" + std::to_string(i)]->mean();
		trueSi[i] = MPS_reference["s_" + std::to_string(i)]->mean();
	}

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
//****

	// pidgeon holes for the decay widths so that can sum over them without losing
	std::vector<real_t> decayWidths_posBins_Bplus, decayWidths_posBins_Bminus, decayWidths_negBins_Bplus, decayWidths_negBins_Bminus;
	createZeros<real_t>(decayWidths_posBins_Bplus, nBins+1); createZeros<real_t>(decayWidths_posBins_Bminus, nBins+1);
	createZeros<real_t>(decayWidths_negBins_Bplus, nBins+1); createZeros<real_t>(decayWidths_negBins_Bminus, nBins+1);

	auto LL_total = [&nBins, &signs, &nEvents, &MPS, &fitType,
					 &retrieveAllMPS, &currentCi, &currentSi, &currentFi, &currentFbari,
					 &gaussianConstraint,
					 &binCounts_posBins_Bminus, &binCounts_negBins_Bminus, &binCounts_posBins_Bplus, &binCounts_negBins_Bplus,
					 &decayWidths_posBins_Bplus, &decayWidths_posBins_Bminus, &decayWidths_negBins_Bplus, &decayWidths_negBins_Bminus ]()
	{
		real_t ll{0};

		// Normalisation for the decay widths:
		real_t sumDecayWidths_Bplus = 0;
		real_t sumDecayWidths_Bminus = 0;


		// calculate partial width and sum them
		retrieveAllMPS();
		for (size_t i{1}; i < nBins+1; i++){ 
		
			decayWidths_negBins_Bminus[i] = decayWidth_Integrated_xy(signs.BminusSign, signs.negBinSign, currentFbari[i], currentFi[i], MPS, currentCi[i], currentSi[i]);
			decayWidths_negBins_Bplus[i] = decayWidth_Integrated_xy(signs.BplusSign, signs.negBinSign, currentFbari[i], currentFi[i], MPS, currentCi[i], currentSi[i]);

			decayWidths_posBins_Bminus[i] = decayWidth_Integrated_xy(signs.BminusSign, signs.posBinSign, currentFi[i], currentFbari[i], MPS, currentCi[i], currentSi[i]);
			decayWidths_posBins_Bplus[i] = decayWidth_Integrated_xy(signs.BplusSign, signs.posBinSign, currentFi[i], currentFbari[i], MPS, currentCi[i], currentSi[i]);

			sumDecayWidths_Bminus += decayWidths_negBins_Bminus[i] + decayWidths_posBins_Bminus[i]; 
			sumDecayWidths_Bplus += decayWidths_negBins_Bplus[i] + decayWidths_posBins_Bplus[i];
		}

		// normalise and calculate ll for each bin
		for (size_t i{1}; i < nBins+1; i++){ 
			real_t nExp_negBins_Bminus = nEvents * decayWidths_negBins_Bminus[i] / sumDecayWidths_Bminus;
			real_t nExp_posBins_Bminus = nEvents * decayWidths_posBins_Bminus[i] / sumDecayWidths_Bminus;
			real_t nExp_negBins_Bplus = nEvents * decayWidths_negBins_Bplus[i] / sumDecayWidths_Bplus;
			real_t nExp_posBins_Bplus = nEvents * decayWidths_posBins_Bplus[i] / sumDecayWidths_Bplus;
		
			// don't know how efficient this is to have the if here?
			if (fitType == "poisson"){
				ll -= 2*logPoisson(nExp_negBins_Bminus, binCounts_negBins_Bminus[i]);
				ll -= 2*logPoisson(nExp_posBins_Bminus, binCounts_posBins_Bminus[i]);
				ll -= 2*logPoisson(nExp_negBins_Bplus, binCounts_negBins_Bplus[i]);
				ll -= 2*logPoisson(nExp_posBins_Bplus, binCounts_posBins_Bplus[i]);
			}else{
				ll += Chi2(nExp_negBins_Bminus, binCounts_negBins_Bminus[i]);
				ll += Chi2(nExp_posBins_Bminus, binCounts_posBins_Bminus[i]);
				ll += Chi2(nExp_negBins_Bplus, binCounts_negBins_Bplus[i]);
				ll += Chi2(nExp_posBins_Bplus, binCounts_posBins_Bplus[i]);
			}
		}

		ll += gaussianConstraint();

		return ll;
	};


	
//******** CHECK HOW LONG LL TAKES TO CALCULATE ********
	if ( NamedParameter<bool>("DoTimeTestOnly", false, "Only do the time test for one LL?") ){
		ProfileClock  timeOneLL;

		timeOneLL.start();
		real_t testLL{LL_total()}; // IGNORE UNUSED VARIABLE WARNING
		timeOneLL.stop();
		INFO("Took "<<timeOneLL.t_duration << "ms to calculate one LogLikelihood lambda");
		INFO("LL has initial value " << testLL << ", where " << gaussianConstraint() << " came from the constraint");
		return 0;
	}
	
//********* DO THE FIT ***********
	Minimiser minimiser(LL_total, &MPS); // formerly mini

	if (NamedParameter<bool>("DoGradientTest", false)){
		minimiser.gradientTest();
	}

	minimiser.doFit();

	INFO("final value of the gaussian constraint " << gaussianConstraint());

	//********* and minos errors
	if ( NamedParameter<bool>("DoMINOS", false, "Calculate minos errors, default false") ){
		ProfileClock timeMINOS;
		timeMINOS.start();
		minimiser.minos(MPS["xPlus"]);
		minimiser.minos(MPS["xMinus"]);
		minimiser.minos(MPS["yPlus"]);
		minimiser.minos(MPS["yMinus"]);
		// minimiser.minos(MPS["rB"]);
		// minimiser.minos(MPS["dB"]);
		// minimiser.minos(MPS["gamma"]);
		for (size_t i{1}; i < nBins+1; i++){
			minimiser.minos(MPS["c_"+std::to_string(i)]);
			minimiser.minos(MPS["s_"+std::to_string(i)]);
		}
		timeMINOS.stop();
		INFO("Took " << timeMINOS.t_duration << "ms to calculate MINOS uncertainties");
	}

	// Now save it
	FitResult* output = new FitResult(minimiser); // formerly fr
	INFO("writing result to logfile " << logFile);
	output->writeToFile(logFile);


// *************** SAVE RAW RESULT ****************
	// print results to something easy to reimport and plot:
	if( !resultFile.empty()){
		writeFitResult(resultFile, MPS, nEvents, binCounts_posBins_Bminus, binCounts_posBins_Bplus, binCounts_negBins_Bminus, binCounts_negBins_Bplus);
	}
	
	// plot the result
	if (NamedParameter<bool>("DoPlots", true)){
		ProfileClock TIME_plots; TIME_plots.start();
		writeBinnedFitPlot(MPS, nEvents, binCounts_posBins_Bminus, binCounts_posBins_Bplus, binCounts_negBins_Bminus, binCounts_negBins_Bplus);
		TIME_plots.stop();
		INFO("Took " << TIME_plots << " to make and write plots to file " << plotFile);
	}

//*********

	std::cout << "\n bottom of code :D" << std::endl;
	return 0;
}
