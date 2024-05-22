#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/Binning.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/Expression.h"
#include "AmpGen/Formalism.h"
#include "AmpGen/Generator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/MyStructs.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PhaseCorrection.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/Projection.h"
#include "AmpGen/MyStructs.h"

#if ENABLE_AVX
  using EventList_type = AmpGen::EventListSIMD;
#else
  using EventList_type = AmpGen::EventList; 
#endif

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>


#include "TCanvas.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TRatioPlot.h"
#include "TStyle.h"
#include "TTree.h"


using namespace AmpGen;
using RealVec_t = std::vector<real_t>;

template <typename T> void createZeros(std::vector<T>& list, const size_t& nItems){
  for (size_t i{0}; i < nItems;i++){
    list.push_back(0);
  }
  return;
}

void readinEventType(EventType& type, const char* string){
	INFO("try for : " << NamedParameter<std::string>(string, " "));
	EventType thisType{ NamedParameter<std::string>(string, " ") };
	type = thisType;
}

int main(int argc, char * argv[]){

// ******* PARSE THE ARGS FROM COMMAND LINE AND THE OPTIONS FILE*********
	OptionsParser::setArgs( argc, argv );

	EventType signalType{NamedParameter<std::string>("EventType", "", "Signal Type to generate, default D0 -> Ks0pi-pi+")}; // notes -> (1) MINUS COMES FIRST

	const size_t      nThreads = NamedParameter<size_t>("nCores", 8, "Number of threads to use" );
	#ifdef _OPENMP
	omp_set_num_threads( nThreads );
	INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
	omp_set_dynamic( 0 );
	#endif

	const std::string inputData = NamedParameter<std::string>("Input", "", "Start of root file of B+- to D h+- events - NO .root NEEDED");

	const std::string binningFile = NamedParameter<std::string>("Binning", "../KsPiPi_equal.txt", "text file with list of binning scheme points");

	const std::string outputFilename = NamedParameter<std::string>("Output", "output.txt", "text file to print the output values - .txt needed");
	
	const std::string plotFilename = NamedParameter<std::string>("Plots", "", "file to print the plot - .pdf needed, if not specified, won't do plot");

	const size_t      nBins = NamedParameter<size_t>("nBins", 8, "Number of bins in scheme, default 8" );

		  size_t      nEvents = NamedParameter<size_t>("nEvents", 25000, "Number of events to generate/use for calculation, default 50000" );

	const size_t      seed = NamedParameter<size_t>("Seed", 0, "Random seed for generation");
	TRandom3 rndm(seed); 

	const size_t      nInt = NamedParameter<size_t>("nInt", 1e7, "Number of events to calculate normalisation - should be large");

	const std::string modelFile = NamedParameter<std::string>("Model", "Kspipi.opt", "Options file with the amplitude model in it");

//************ FOR ANY PLOTTING: CREATE THE CANVAS ***********
	// TCanvas* cv = new TCanvas("cv", "cv", 500, 500);
	// cv->SetLeftMargin(0.15);
	// cv->SetRightMargin(0.15);
	// cv->Divide(4, 2);
//************


//***************** MAKE INTEGRATION EVENTS AND INFORMAITON **********
	// make events
	TRandom3 random(seed);

	// ***** SET UP THE PHASE CORRECTION *****
	MinuitParameterSet MPS_phaseCorrection; 
	MPS_phaseCorrection.loadFromStream();	// base mps with the phase correction info

	PhaseCorrection deltaCorrection{MPS_phaseCorrection}; 
	deltaCorrection.compilePolynomialExpressions(signalType);
//	PhaseCorrection_KLpipi deltaCorrection_KLpipi{MPS_phaseCorrection}; 
//	deltaCorrection_KLpipi.compilePolynomialExpressions(signalType);

	// ***** SET UP SIGNAL (Kspipi) AMPLITUDE *****
	EventList_type data{ Generator<>(signalType, &random).generate(nInt) };
	// Event event{ data[2] };
	MinuitParameterSet MPS_Kspipi;
	MPS_Kspipi.loadFromFile(modelFile); 
	AddCPConjugate(MPS_Kspipi);

	// // get amlitude info
	CoherentSum A(signalType, MPS_Kspipi);
	CoherentSum Abar(signalType.conj(true), MPS_Kspipi);
	A.setEvents(data);		Abar.setEvents(data);
	A.prepare(); 			Abar.prepare();

	auto evalA = A.amplitudeEvaluator(&data);
	auto evalAbar = Abar.amplitudeEvaluator(&data);


	// std::ofstream outFile{outputFilename};
	// for ( auto event : data ){
	// 	outFile << std::arg(evalA(event)*std::conj(evalAbar(event))) - M_PI << ", ";
	// }
	// outFile.close();
	// amplitudeInfo amplitudesMC = fillAmplitudeInfo(data, A, Abar);
	// amplitudesMC.bins = {1, 5, -8, -4, 2, -3, 3, 8, -5, -1};
//*********************

//*********** FOR ANY PLOTTING: GRAPHICS ATTRIBUTES ********
	// int attEmptySquareMarker{25};
	// int attFillSquareMarker{21};
	// int attEmptyCircleMarker{24};
	// int attFillCircleMarker{20};
	// int attFillTriangleMarker{23};
	// int attEmptyTriangleMarker{32};

	// int colPurple{ 618 };
	// int colDarkGreen{ 419 };
	// int colGreen{ 409 };
	// int colBlue{ 867 };

//****************

//************ FOR ANY *GRAPH* PLOTTING: SET UP A FRAME ***********
	// // TH2* frame = new TH2D("frame", "",100, -8.5, 8.5, 100, 0, 0.19);
	// TH2* frame = new TH2D("frame", ";;;absolute bin number",100, -1.2, 1.2, 100, -1.2, 1.2);
	// // TH2* frame = new TH2D("frame", "",100, 0, 3.5, 100, 0, 3.5);
	// // // Prepare frame
	// // frame->SetMinimum(-1.2);
	// // frame->SetMaximum(1.2);

	// frame->SetTitle("c_{i} and s_{i} comparison");
	// frame->GetXaxis()->SetTitle("c_{i}");
	// frame->GetYaxis()->SetTitle("s_{i}");
	// // // // frame->SetTitle("F_{i} found from counting bin yields");
	// // // // frame->GetXaxis()->SetTitle("bin number, i");
	// // // // frame->GetYaxis()->SetTitle("F_{i}");
	// // // std::string title{ NamedParameter<std::string>("Title", "", "title of plot") };
	// // // frame->SetTitle(title.c_str());
	// // frame->GetXaxis()->SetTitle("m(K^{0}_{S} #pi^{-}) = m_{-}");
	// // frame->GetYaxis()->SetTitle("m(K^{0}_{S} #pi^{+}) = m_{+}");


	// frame->SetStats(0);

	// // // // ****** just draw:
	// // // // frame->Draw(" ");
	// // // ****** draw for colourbar settings (defntly not the proper way, but comes out as i want sooo):
	// frame->SetMinimum(0.5);
	// frame->SetMaximum(8.5);
	// double frameBins[nBins] = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
	// frame->SetContour(nBins, frameBins); // setting ticks to right numbers
	// frame->Fill(4, 5); // out of range point to make colourbar draw
	// frame->Draw("colz");
	// int colours[nBins+1] = {0, 2, 3, 4, 801, 6, 38, 419, 881};
	// gStyle->SetPalette(nBins+1, colours);

//************

//**************** READ IN A MPS FROM FILE *************
	// MinuitParameterSet MPS;
	// MPS.loadFromFile(NamedParameter<std::string>("MPSfile", "", "file to load for MPS"));
//**************************

//*************** READ IN MULTIPLE MPS FROM FILE *************
	// MinuitParameterSet MPS_opt1;
	// MPS_opt1.loadFromFile(NamedParameter<std::string>("parameterfile1", " "));
	// MinuitParameterSet MPS_opt2;
	// MPS_opt2.loadFromFile(NamedParameter<std::string>("parameterfile2", " "));
	// MinuitParameterSet MPS_opt3;
	// // MPS_opt3.loadFromFile(NamedParameter<std::string>("parameterfile3", " "));
	// std::string label1 = NamedParameter<std::string>("label1", " " );
	// std::string label2 = NamedParameter<std::string>("label2", " " );
	// // std::string label3 = NamedParameter<std::string>("label3", " " );
//**************************

//**************** READ IN MPS FROM STREAM *************
	// MinuitParameterSet MPS;
	// MPS.loadFromStream(); 
	// AddCPConjugate(MPS);
//**************************

//*************** PUT MPS VALUES FOR BIN YIELD INTO VARIABLES (nObs naming comvention) *******************
	// RealVec_t nObs_posBins_Bminus, nObs_negBins_Bminus, nObs_posBins_Bplus, nObs_negBins_Bplus;
	// createZeros<real_t>(nObs_posBins_Bminus, nBins+1); createZeros<real_t>(nObs_negBins_Bminus, nBins+1);
	// createZeros<real_t>(nObs_posBins_Bplus, nBins+1); createZeros<real_t>(nObs_negBins_Bplus, nBins+1);

	// RealVec_t nExp_posBins_Bminus, nExp_negBins_Bminus, nExp_posBins_Bplus, nExp_negBins_Bplus;
	// createZeros<real_t>(nExp_posBins_Bminus, nBins+1); createZeros<real_t>(nExp_negBins_Bminus, nBins+1);
	// createZeros<real_t>(nExp_posBins_Bplus, nBins+1); createZeros<real_t>(nExp_negBins_Bplus, nBins+1);

	// for (size_t i{1}; i < nBins+1; i++){
	// 	nObs_posBins_Bminus[i] = MPS["nObs_"+std::to_string(i) + "_Bminus"]->mean();
	// 	nObs_negBins_Bminus[i] = MPS["nObs_-"+std::to_string(i) + "_Bminus"]->mean(); 
	// 	nObs_posBins_Bplus[i] = MPS["nObs_"+std::to_string(i) + "_Bplus"]->mean();
	// 	nObs_negBins_Bplus[i] = MPS["nObs_-"+std::to_string(i) + "_Bplus"]->mean();

	// 	nExp_posBins_Bminus[i] = MPS["nExp_"+std::to_string(i) + "_Bminus"]->mean(); 
	// 	nExp_negBins_Bminus[i] = MPS["nExp_-"+std::to_string(i) + "_Bminus"]->mean(); 
	// 	nExp_posBins_Bplus[i] = MPS["nExp_"+std::to_string(i) + "_Bplus"]->mean();
	// 	nExp_negBins_Bplus[i] = MPS["nExp_-"+std::to_string(i) + "_Bplus"]->mean();
	// }
//**************************

//*************** AS ABOVE, BIN COUNT NAMING CONVENTION ***************
	// std::vector<size_t> binCounts_posBins_Bminus, binCounts_negBins_Bminus;
	// createZeros<size_t>(binCounts_posBins_Bminus, nBins+1);
	// createZeros<size_t>(binCounts_negBins_Bminus, nBins+1);
	// std::vector<size_t> binCounts_posBins_Bplus, binCounts_negBins_Bplus;
	// createZeros<size_t>(binCounts_posBins_Bplus, nBins+1);
	// createZeros<size_t>(binCounts_negBins_Bplus, nBins+1);

	// for (size_t i{1}; i < nBins+1; i++){
	// 	binCounts_posBins_Bminus[i] = MPS["Bminus_" + std::to_string(i)]->mean();
	// 	binCounts_negBins_Bminus[i] = MPS["Bminus_-" + std::to_string(i)]->mean(); 
	// 	binCounts_posBins_Bplus[i] = MPS["Bplus_" +  std::to_string(i)]->mean();
	// 	binCounts_negBins_Bplus[i] = MPS["Bplus_-" + std::to_string(i)]->mean();
	// }
//*******************

//*************** READ IN SOME VALUES FROM MULTIPLE MPS TO PLOT AND COMPARE **************
	// MinuitParameterSet MPS_calc;
	// MinuitParameterSet MPS_count;
	// // MPS_calc.loadFromFile("binnedFitParameters.opt");
	// MPS_calc.loadFromFile("FisSeparateDDbarEvents.opt");
	// MPS_count.loadFromFile("countedFisCheck.opt");

	// std::vector<real_t> F_calc_posbins, Fbar_calc_posbins, F_count_posbins, Fbar_count_posbins;
	// std::vector<real_t> F_calc_negbins, Fbar_calc_negbins, F_count_negbins, Fbar_count_negbins;
	// F_calc_posbins.push_back(0); Fbar_calc_posbins.push_back(0); F_count_posbins.push_back(0); Fbar_count_posbins.push_back(0);
	// F_calc_negbins.push_back(0);; Fbar_calc_negbins.push_back(0); F_count_negbins.push_back(0); Fbar_count_negbins.push_back(0);

	// for (int i{1}; i < nBins+1; i++){
	// 	F_calc_posbins.push_back(MPS_calc["F_"+std::to_string(i)]->mean());
	// 	Fbar_calc_posbins.push_back(MPS_calc["Fbar_"+std::to_string(i)]->mean());
	// 	F_calc_negbins.push_back(MPS_calc["F_-"+std::to_string(i)]->mean());
	// 	Fbar_calc_negbins.push_back(MPS_calc["Fbar_-"+std::to_string(i)]->mean());
	
	// 	F_count_posbins.push_back(MPS_count["F_"+std::to_string(i)]->mean());
	// 	Fbar_count_posbins.push_back(MPS_count["Fbar_"+std::to_string(i)]->mean());
	// 	F_count_negbins.push_back(MPS_count["F_-"+std::to_string(i)]->mean());
	// 	Fbar_count_negbins.push_back(MPS_count["Fbar_-"+std::to_string(i)]->mean());
	// }

	// INFO("just showing correct value loaded, " << MPS_calc["F_1"]->mean() << " should be 0.174263");

	// // SETUP
	// TCanvas* cv = new TCanvas("cv", "scatter plot", 700, 700);
	// TH1* frame = new TH1D("frame", "",100, -8.5, 8.5);

	// // Prepare frame
	// frame->SetMinimum(0);
	// frame->SetMaximum(0.2);
	// frame->SetTitle("F_{i} and #bar{F_{i}} values");
	// frame->GetXaxis()->SetTitle("bin number");
	// frame->GetYaxis()->SetTitle("F / #bar{F}");
	// frame->SetStats(0);
	// frame->Draw(" ");
	// TLegend* legend = new TLegend(0.1, 0.7, 0.4, 0.9);

	// // Draw each point of +/- bins x 2 MPS's
	// for (int i{1}; i<9; i++){
	// 	real_t x[2] ={-i, i};
	// 	INFO("F calc pos bins item " << i << ": " << F_calc_posbins[i]);
	// 	real_t y_calc[2] = {F_calc_negbins[i], F_calc_posbins[i]};
	// 	real_t y_calcbar[2] = {Fbar_calc_negbins[i], Fbar_calc_posbins[i]};
	// 	real_t y_count[2] = {F_count_negbins[i], F_count_posbins[i]};
	// 	real_t y_countbar[2] = {Fbar_count_negbins[i], Fbar_count_posbins[i]};

	// 	TGraph* g_calc = new TGraph(2, x, y_calc); 
	// 	g_calc->SetMarkerSize(3); g_calc->SetMarkerStyle(68); g_calc->SetMarkerColor(409); // light green
	// 	g_calc->Draw("SAME P");

	// 	TGraph* g_calcbar = new TGraph(2, x, y_calcbar); 
	// 	g_calcbar->SetMarkerSize(3); g_calcbar->SetMarkerStyle(70); g_calcbar->SetMarkerColor(419); // dark green
	// 	g_calcbar->Draw("SAME P");

	// 	TGraph* g_count = new TGraph(2, x, y_count); 
	// 	g_count->SetMarkerSize(3); g_count->SetMarkerStyle(68); g_count->SetMarkerColor(876); // light purple
	// 	g_count->Draw("SAME P");

	// 	TGraph* g_countbar = new TGraph(2, x, y_countbar); 
	// 	g_countbar->SetMarkerSize(3); g_countbar->SetMarkerStyle(70); g_countbar->SetMarkerColor(874); // dark purple
	// 	g_countbar->Draw("SAME P");

	// 	if (i==1){
	// 		// legend->AddEntry(g_calc, "F calculated with Amplitudes", "P");
	// 		// legend->AddEntry(g_calcbar, "#bar{F} calculated with Amplitudes", "P");
	// 		legend->AddEntry(g_calc, "F separate D#bar{D} events", "P");
	// 		legend->AddEntry(g_calcbar, "#bar{F} separate D#bar{D} events", "P");
	// 		legend->AddEntry(g_count, "F 1 set of event", "P");
	// 		legend->AddEntry(g_countbar, "#bar{F} 1 set of events", "P");
	// 	}

	// }
	// legend->Draw();
	// cv->SaveAs("CompareFs2.pdf");

	// return 0;
//************************************************************

//****************** CALCULATE EXPECTED nOBS **************
	// // create vectors to hold the nObs
	// RealVec_t binCounts_posBins_Bminus, binCounts_negBins_Bminus;
	// createZeros<real_t>(binCounts_posBins_Bminus, nBins+1);
	// createZeros<real_t>(binCounts_negBins_Bminus, nBins+1);
	// RealVec_t binCounts_posBins_Bplus, binCounts_negBins_Bplus;
	// createZeros<real_t>(binCounts_posBins_Bplus, nBins+1);
	// createZeros<real_t>(binCounts_negBins_Bplus, nBins+1);
	
	// std::vector<real_t> decayWidths_posBins_Bplus, decayWidths_posBins_Bminus, decayWidths_negBins_Bplus, decayWidths_negBins_Bminus;
	// createZeros<real_t>(decayWidths_posBins_Bplus, nBins+1); createZeros<real_t>(decayWidths_posBins_Bminus, nBins+1);
	// createZeros<real_t>(decayWidths_negBins_Bplus, nBins+1); createZeros<real_t>(decayWidths_negBins_Bminus, nBins+1);
	
	// real_t sumDecayWidths_Bplus = 0;
	// real_t sumDecayWidths_Bminus = 0;
	// signs signs{};
	
	// std::default_random_engine noiseGen;
	// noiseGen.seed(seed); // for adding some stat noise - NOTE without seed, running this script twice indep gives same random noise so same dist

	// for (size_t i{1}; i < nBins+1; i++){ 
	// 	real_t ci = MPS["c_" + std::to_string(i)]->mean(); 
	// 	real_t si = MPS["s_" + std::to_string(i)]->mean();
	// 	real_t Fi = MPS["F_" + std::to_string(i)]->mean();
	// 	real_t Fbari = MPS["Fbar_" + std::to_string(i)]->mean();

	// 	decayWidths_posBins_Bminus[i] = decayWidth_Integrated_xy(signs.BminusSign, signs.posBinSign, Fi, Fbari, MPS, ci, si);
	// 	decayWidths_posBins_Bplus[i] = decayWidth_Integrated_xy(signs.BplusSign, signs.posBinSign, Fi, Fbari, MPS, ci, si);
	// 	decayWidths_negBins_Bminus[i] = decayWidth_Integrated_xy(signs.BminusSign, signs.negBinSign, Fbari, Fi, MPS, ci, si);
	// 	decayWidths_negBins_Bplus[i] = decayWidth_Integrated_xy(signs.BplusSign, signs.negBinSign, Fbari, Fi, MPS, ci, si);


	// 	sumDecayWidths_Bminus += decayWidths_negBins_Bminus[i] + decayWidths_posBins_Bminus[i]; 
	// 	sumDecayWidths_Bplus += decayWidths_negBins_Bplus[i] + decayWidths_posBins_Bplus[i];
	// }

	// // add noise:
	// real_t sumNoisyBinCounts_Bminus{0}, sumNoisyBinCounts_Bplus{0};
	// for (size_t i{1}; i < nBins+1; i++){
	// 	binCounts_posBins_Bminus[i] =  nEvents * decayWidths_posBins_Bminus[i] / sumDecayWidths_Bminus;
	// 	binCounts_posBins_Bplus[i] = nEvents * decayWidths_posBins_Bplus[i] / sumDecayWidths_Bplus;
	// 	binCounts_negBins_Bminus[i] =  nEvents * decayWidths_negBins_Bminus[i] / sumDecayWidths_Bminus;
	// 	binCounts_negBins_Bplus[i] = nEvents * decayWidths_negBins_Bplus[i] / sumDecayWidths_Bplus;
		
	// 	// add poisson fluctuation
	// 	std::poisson_distribution<int> dist_pos_minus(binCounts_posBins_Bminus[i]);
	// 	std::poisson_distribution<int> dist_pos_plus(binCounts_posBins_Bplus[i]);
	// 	std::poisson_distribution<int> dist_neg_min(binCounts_negBins_Bminus[i]);
	// 	std::poisson_distribution<int> dist_neg_plus(binCounts_negBins_Bplus[i]);
		
	// 	binCounts_posBins_Bminus[i] = dist_pos_minus(noiseGen);
	// 	binCounts_posBins_Bplus[i] = dist_pos_plus(noiseGen);
	// 	binCounts_negBins_Bminus[i] = dist_neg_min(noiseGen);
	// 	binCounts_negBins_Bplus[i] = dist_neg_plus(noiseGen);

	// 	sumNoisyBinCounts_Bminus += binCounts_negBins_Bminus[i] + binCounts_posBins_Bminus[i];
	// 	sumNoisyBinCounts_Bplus += binCounts_negBins_Bplus[i] + binCounts_posBins_Bplus[i];
	// }

	// // final normalisation
	// for (size_t i{1}; i < nBins+1; i++){
	// 	binCounts_posBins_Bminus[i] =  nEvents * binCounts_posBins_Bminus[i] / sumNoisyBinCounts_Bminus;
	// 	binCounts_posBins_Bplus[i] = nEvents * binCounts_posBins_Bplus[i] / sumNoisyBinCounts_Bplus;
	// 	binCounts_negBins_Bminus[i] =  nEvents * binCounts_negBins_Bminus[i] / sumNoisyBinCounts_Bminus;
	// 	binCounts_negBins_Bplus[i] = nEvents * binCounts_negBins_Bplus[i] / sumNoisyBinCounts_Bplus;
	// }

//**********************
	
	
//*************** READ IN A BINNING SCHEME **************
	// INFO("Reading in binning scheme");
	// std::vector<real_t> mMinus;
	// std::vector<real_t> mPlus;
	// std::vector<int> bins;

	// readBinning(binningFile, mMinus, mPlus, bins); //pushes back the lists from text file into these vectors
//*****


// ************** MAKE FLAT EVENTS ***************
	// INFO("Generating events");
	// TRandom3 random(seed);
	// EventList_type data{ Generator<>(signalType, &random).generate(nInt) };
//******

// ************* READ EVENTS INTO EVENT LISTS ***************
	// INFO("Loading events");
	// EventList_type data_Bminus((inputData + ":Bminus__DalitzEventList").c_str() , signalType); 
	// EventList_type data_Bplus((inputData + ":Bplus__DalitzEventList").c_str() , signalType); 
	// EventList_type data(inputData.c_str() , signalType); 
	// nEvents = data.size();
//********

//***************** MAKE A COHERENT SUM INFO ********
	// CoherentSum A(signalType, MPS);
	// CoherentSum Abar(signalType.conj(true), MPS);
	// A.setEvents(data); Abar.setEvents(data);
	// A.prepare(); Abar.prepare();
	// auto evalA = A.amplitudeEvaluator( &data ); 
	// auto evalAbar = Abar.amplitudeEvaluator( &data ); 
//************

//************ CALC Ci and Si FROM THE FLAT DATA *************
	// RealVec_t F_posBins, Fbar_posBins, c_posBins, s_posBins;
	// createZeros(F_posBins, nBins+1);	createZeros(Fbar_posBins, nBins+1);	createZeros(c_posBins, nBins+1);	createZeros(s_posBins, nBins+1);
	// RealVec_t F_negBins, Fbar_negBins, c_negBins, s_negBins;
	// createZeros(F_negBins, nBins+1);	createZeros(Fbar_negBins, nBins+1);	createZeros(c_negBins, nBins+1);	createZeros(s_negBins, nBins+1);

	// // // formula 1 to try
	// auto fillSums = [&evalA, &evalAbar](Event& event, int& binNumber, RealVec_t& FSums, RealVec_t& FbarSums, RealVec_t& cSums, RealVec_t& sSums)
	// {
	// 	complex_t complexA = evalA(event);
	// 	complex_t complexAbar = evalAbar(event);


	// 	real_t modA{ std::abs(complexA) };
	// 	real_t modAbar{ std::abs(complexAbar) };
	// 	real_t Dd { std::arg(complexA * std::conj(complexAbar))};

	// 	FSums[binNumber] += std::pow(modA, 2);
	// 	FbarSums[binNumber] += std::pow(modAbar, 2);
	// 	cSums[binNumber] += modA * modAbar * cos(Dd);
	// 	sSums[binNumber] += modA * modAbar * sin(Dd);
	// 	return;
	// };
	// // // // formula 2 to cross check
	// // auto fillSums = [&evalA, &evalAbar](Event& event, int& binNumber, RealVec_t& FSums, RealVec_t& FbarSums, RealVec_t& cSums, RealVec_t& sSums)
	// // {
	// // 	complex_t A = evalA(event);
	// // 	complex_t Abar = evalAbar(event);

	// // 	FSums[binNumber] += A.real()*A.real() + A.imag()*A.imag();
	// // 	FbarSums[binNumber] += Abar.real()*Abar.real() + Abar.imag()*Abar.imag();
	// // 	cSums[binNumber] += A.real()*Abar.real() + A.imag()*Abar.imag();
	// // 	sSums[binNumber] += -A.real()*Abar.imag() + Abar.real()*A.imag();
	// // 	return;
	// // };

	// // loop through the flat events
	// INFO("looping through events");
	// for (auto event : data){
	// 	int bin = nearestBinIndex(event, mMinus, mPlus, bins);

	// 	if (bin < 0){
	// 		int unsignedBin{ abs(bin) };
	// 		fillSums(event, unsignedBin, F_negBins, Fbar_negBins, c_negBins, s_negBins);
	// 	}else{
	// 		fillSums(event, bin, F_posBins, Fbar_posBins, c_posBins, s_posBins);
	// 	}
	// }

	// // // normalise things
	// // INFO("normalising ci and si");
	// // real_t totalF{0}, totalFbar{0};
	// // for (size_t i{1}; i < nBins+1; i++){
	// // 	c_posBins[i] = c_posBins[i] / std::sqrt(F_posBins[i] * Fbar_posBins[i]);
	// // 	c_negBins[i] = c_negBins[i] / std::sqrt(F_negBins[i] * Fbar_negBins[i]);
	// // 	s_posBins[i] = s_posBins[i] / std::sqrt(F_posBins[i] * Fbar_posBins[i]);
	// // 	s_negBins[i] = s_negBins[i] / std::sqrt(F_negBins[i] * Fbar_negBins[i]);

	// // 	totalF += F_posBins[i] + F_negBins[i];
	// // 	totalFbar += Fbar_posBins[i] + Fbar_negBins[i];
	// // }


	// std::ofstream outFile{outputFilename};
	// INFO("Writing results to options file " << outputFilename);
	// // default for options output is Fs by counting number of events and ci/si from summing amplitudes only
	// for (size_t i{1}; i < nBins+1 ; i++){ 
	// 		outFile 	<< "sumA_ " << i << "     Fix     " << F_posBins[i]<< "     0\n"
	// 					<< "sumAbar_ "  << i << "  Fix     " << Fbar_posBins[i] << "     0\n"
	// 				 	<< "sumA_- " << i << "     Fix     " << F_negBins[i]  << "     0\n"
	// 					<< "sumAbar_- " << i << "  Fix     " << Fbar_negBins[i]  << "     0\n"
	// 					<< "sumCos_ " << i << "     Fix     " << c_posBins[i]<< "     0\n"
	// 					<< "sumSin_ " << i << "     Fix     " << s_posBins[i] << "     0\n"
	// 					<< "sumCos_- " << i << "    Fix     " << c_negBins[i] << "     0\n"
	// 					<< "sumSin_- " << i << "    Fix     " << s_negBins[i] << "     0\n";
	// }
	// outFile.close();

//***************



//************ MAKE PHASE CORRECTION *************
	// INFO("mps initialsation: ");
	// PhaseCorrection phaseCorrection{MPS};
	// phaseCorrection.compilePolynomialExpressions(signalType);
//*******************

// ******* PLOT SOMETHING: PHASE CORRECTION: coords ****** 
	// make compiled expressions for the coords:
	// dalitzPair<Expression> ex = phaseCorrection.transformedDalitzCoords(signalType);
	// CompiledExpression<real_t(const real_t*, const real_t*)> compiledZplus(ex.plus, "name", signalType.getEventFormat());
	// CompiledExpression<real_t(const real_t*, const real_t*)> compiledZminus(ex.minus, "name", signalType.getEventFormat());
	// compiledZplus.prepare();	compiledZminus.prepare();
	// compiledZplus.compile();	compiledZminus.compile();

	// // prepare the histogram:
	// TH2D* dalitz = new TH2D("dalitz", "#delta_{C}, phase correction;z''_{+};z'_{-};#delta_{C}/rad", 200, -1, 1, 200,-1, 1);
	// TH2D* dalitz = new TH2D("dalitz", ";;;#delta_{D}/rad", 200, -1, 1, 200,-1, 1);

	// for (auto event : data){
	// 	dalitz->Fill(compiledZplus(event.address()), compiledZminus(event.address()));
	// }


	// // draw the histo:
	// // dalitz->SetStats(false);
	// // dalitz->Draw("");


//**********

// ****** CHECK GAUSSIAN BIAS ********
	// INFO("checking coeffs");
	// phaseCorrection.printout();

	// INFO("checking a couplt of points:");
	// int count{1};
	// for ( auto event : data ){
	// 	// INFO("\n\n " << count);
	// 	INFO("********************************loop: " << count);
	// 	INFO("s02: " << event.s(0, 2) << ", s01: " <<  event.s(0, 1));
	// 	real_t erf = std::erf(10*(event.s(0, 2) - event.s(0, 1)));
	// 	real_t gauss1{0}, gauss2{0};
	// 	INFO("erf argument should be: " << (event.s(0, 2) - event.s(0, 1))*10 << ", leading to " << erf);
	// 	if (event.s(0, 2) > event.s(0, 1)){
	// 		gauss1 = std::pow(4*(event.s(0, 2)-1.0), 2) + std::pow(4*(event.s(0, 1)-1.25), 2);
	// 		gauss2 = std::pow(4*(event.s(0, 2)-2.5), 2) + std::pow(4*(event.s(0, 1)-1.25), 2);
	// 		INFO("Gauss arg bias 1 should be: " << std::pow(4*(event.s(0, 2)-1.0), 2) << " + " << std::pow(4*(event.s(0, 1)-1.25), 2) << " = " << gauss1);
	// 		INFO("Gauss arg bias 2 should be: " << std::pow(4*(event.s(0, 2)-2.5), 2) << " + " << std::pow(4*(event.s(0, 1)-1.25), 2) << " = " << gauss2);
	// 	}else{
	// 		gauss1 = std::pow(4*(event.s(0, 1)-1.0), 2) + std::pow(4*(event.s(0, 2)-1.25), 2);
	// 		gauss2 = std::pow(4*(event.s(0, 1)-2.5), 2) + std::pow(4*(event.s(0, 2)-1.25), 2);
	// 		INFO("Gauss arg bias 1 should be: " << std::pow(4*(event.s(0, 1)-1.0), 2) << " + " << std::pow(4*(event.s(0, 2)-1.25), 2) << " = " << gauss1);
	// 		INFO("Gauss arg bias 2 should be: " << std::pow(4*(event.s(0, 1)-2.5), 2) << " + " << std::pow(4*(event.s(0, 2)-1.25), 2) << " = " << gauss2);
	// 	}
	// 	INFO("so overall deltaC " << erf*std::exp(-gauss1) << " + " << erf*std::exp(-gauss2) << " = " << erf * (std::exp(-gauss1) - std::exp(-gauss2)));
	// 	phaseCorrection.evalBias(event);
	// 	count++;
	// }	



	// for (auto event : data){
	// 	INFO("*******LOOP " << count);
	// 	phaseCorrection.evalBias(event);
	// 	INFO("ran evalBias");
	// }

	// phaseCorrection.compilePolynomialExpressions(signalType);
	// INFO("empty initialsation: ");
	// PhaseCorrection phaseCorrection2;

//***************

//**************** PLOT SOMETHING: PHASE CORRECTION: DALITZ PLANE  *********
	size_t nHistoBins{ 500 };
	// TH2D* finalDalitz = new TH2D("dalitz", "#delta_{C}, phase correction;s_{-};s_{+};#delta_{C}/rad", nHistoBins, 0, 3.5, nHistoBins, 0, 3.5);
	TH2D* finalDalitz = new TH2D("dalitz", ";;;#delta_{D}/rad", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* finalBias = new TH2D("Bias", ";;;#delta_{D}^{bias}/rad", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);

	// TH2D* finalDalitz = new TH2D("dalitz", "#|delta_{D}|, absolute strong phase difference;s_{-};s_{+};#delta_{D}/rad", nHistoBins, 0, 3.5, nHistoBins, 0, 3.5);
	TH2D* fullCorrection = new TH2D("pc", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* fullPhaseDifference = new TH2D("pd", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);

	TH2D* scaleHisto = new TH2D("n", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);

	INFO("about to loop through data");
	for ( auto event : data ){
		// real_t PC{ phaseCorrection.eval(event) };
		real_t PD{ std::arg(evalA(event) * std::conj(evalAbar(event))) }; // actually delta_D
		real_t PC{ deltaCorrection.evalBias(event) };
		//if (PD < 0){PD += 2*M_PI;}
		//if (PC < 0){PC += 2*M_PI;}
		fullPhaseDifference->Fill(event.s(0, 1), event.s(0, 2), PD);
		fullCorrection->Fill(event.s(0, 1), event.s(0, 2), PC);
		scaleHisto->Fill(event.s(0, 1), event.s(0, 2));
	}
	INFO("filled histo");

	for ( size_t i{1}; i < nHistoBins; i++){
		for ( size_t j{1}; j < nHistoBins; j++){
			real_t scaler = scaleHisto->GetBinContent(i, j);
			real_t phaseDifference = fullPhaseDifference->GetBinContent(i, j);
			real_t correction = fullCorrection->GetBinContent(i, j);
			if (scaler == 0 && correction == 0){
				finalDalitz->SetBinContent(i, j, 0);
			}else{

				finalDalitz->SetBinContent(i, j, (phaseDifference+correction)/scaler);
				finalBias->SetBinContent(i, j, (correction)/scaler);
			}
		}

	}


	// // finalDalitz->GetXaxis()->SetTitle("s_{-}");
	// // finalDalitz->GetYaxis()->SetTitle("s_{+}");
	// // finalDalitz->SetStats(false);
	
	// // TH1D* xProjection = finalDalitz->ProjectionX();
	
	// // cv->cd(1);
	// // finalDalitz->Draw("COLZ1");
	// // cv->cd(2);
	// // fullCorrection->Draw("COLZ1");
	// // xProjection->Draw("");

	// //save to root file
	// TFile * outputFile = TFile::Open(plotFilename.c_str(), "RECREATE"); 
	// outputFile->cd();
	// finalDalitz->Write();
	// outputFile->Close();
//***********

// *********** PLOT SOMETHING: DALITZ PLANES OF B+/- TO COMPARE*********
	// TH2D* dalitz_Bminus = new TH2D("Bminus", ";;;fraction of Events", 200, 0.2, 3.1, 350, 0.2, 3.1);
	// TH2D* dalitz_Bplus = new TH2D("Bplus", ";;;fraction of Events", 200, 0.2, 3.1, 350, 0.2, 3.1);
	// // TH2D* dalitz_difference = new TH2D("difference", "B^{+} events - flipped B^{-} events", 200, 0, 3.5, 200, 0, 3.5);


	// for (int i{1}; i < 8; i++){
	// 	std::string filename{ inputData+ std::to_string(i)+".root"};
	// 	EventList_type data_Bminus((filename + ":Bminus__DalitzEventList").c_str() , signalType); 
	// 	EventList_type data_Bplus((filename + ":Bplus__DalitzEventList").c_str() , signalType); 

	// 	for ( auto event : data_Bminus ){
	// 		dalitz_Bminus->Fill(event.s(0, 1), event.s(0, 2));
	// 	}
	// 	for ( auto event : data_Bplus ){
	// 		dalitz_Bplus->Fill(event.s(0, 1), event.s(0, 2));
	// 	}


	// }


	// manually get difference between histograms as cant get Add function to work :(
	// for (int i{1}; i < 201; i++){
	// 	for (int j{1}; j < 201; j++){
	// 		int bin_Bminus = dalitz_Bminus->GetBinContent(j, i); // j i for flipping one plot
	// 		int bin_Bplus = dalitz_Bplus->GetBinContent(i, j);
	// 		// dalitz_difference->SetBinContent(i, j, bin_Bplus - bin_Bminus);
	// 	}
	// }


	// dalitz_Bminus->GetXaxis()->SetTitle("m(K^{0}_{S} #pi^{-}) = s_{01}");
	// dalitz_Bplus->GetXaxis()->SetTitle("m(K^{0}_{S} #pi^{-}) = s_{01}");
	// dalitz_Bminus->GetYaxis()->SetTitle("m(K^{0}_{S} #pi^{+}) = s_{02}");
	// dalitz_Bplus->GetYaxis()->SetTitle("m(K^{0}_{S} #pi^{+}) = s_{02}");

	// dalitz_difference->GetXaxis()->SetTitle("s_{01}(B^{-}) || s_{02}(B^{+})");
	// dalitz_difference->GetYaxis()->SetTitle("s_{02}(B^{-}) || s_{01}(B^{+})");
	// dalitz_difference->GetXaxis()->SetTitle("m(K^{0}_{S} #pi^{-}) = s_{01}");
	// dalitz_difference->GetYaxis()->SetTitle("m(K^{0}_{S} #pi^{+}) = s_{02}");
	

	// cv->cd(1);
	// dalitz_Bminus->SetStats(false);
	// dalitz_Bminus->Draw("COLZ");
	// cv->cd(2);
	// dalitz_Bplus->SetStats(false);
	// dalitz_Bplus->Draw("COLZ");
	// cv->cd(3);
	// dalitz_difference->SetStats(false);
	// dalitz_difference->Draw("colz1");
//*********

// ***************** BIN EVENTS IN FORM OF EVENT LIST *************
	// INFO("Binning events");
	// std::vector<int> binList_Bminus; 
	// std::vector<int> binList_Bplus; 
	// int eventBin{0};

	// for (auto& event:data_Bminus){
	// 	eventBin = nearestBinIndex(event, mMinus, mPlus, bins);
	// 	binList_Bminus.push_back(eventBin);
	// }
	// for (auto& event:data_Bplus){
	// 	eventBin = nearestBinIndex(event, mMinus, mPlus, bins);
	// 	binList_Bplus.push_back(eventBin);
	// }
//***********

// ******* ALT: BIN AND COUNT POPULATION IN ONE GO ***********
	// create vectors to hold the nObs
	// RealVec_t binCounts_posBins_Bminus, binCounts_negBins_Bminus;
	// createZeros<real_t>(binCounts_posBins_Bminus, nBins+1);
	// createZeros<real_t>(binCounts_negBins_Bminus, nBins+1);
	// RealVec_t binCounts_posBins_Bplus, binCounts_negBins_Bplus;
	// createZeros<real_t>(binCounts_posBins_Bplus, nBins+1);
	// createZeros<real_t>(binCounts_negBins_Bplus, nBins+1);

	// int currentBin{0};

	// ***** bin the events and count them ****
		// B- events
	// for (auto event : data_Bminus){
	// 	currentBin = nearestBinIndex(event, mMinus, mPlus, bins);
	// 	if (currentBin > 0){
	// 		binCounts_posBins_Bminus[currentBin]++;
	// 	}else{
	// 		binCounts_negBins_Bminus[abs(currentBin)]++;
	// 	}
	// }
		//B+ events
	// for (auto event : data_Bplus){
	// 	currentBin = nearestBinIndex(event, mMinus, mPlus, bins);
	// 	if (currentBin > 0){
	// 		binCounts_posBins_Bplus[currentBin]++;
	// 	}else{
	// 		binCounts_negBins_Bplus[abs(currentBin)]++;
	// 	}
	// }


//***************

//**************** COUNT BIN POPULATIONS - POS/NEG AND PLUS/MINUS SEPARATION ****************
	// // initialise some empty vectors to take no events -> nBins+1 entries, entry 0 should be 0
	// std::vector<size_t> binCounts_posBins_Bminus, binCounts_negBins_Bminus;
	// createZeros<size_t>(binCounts_posBins_Bminus, nBins+1);
	// createZeros<size_t>(binCounts_negBins_Bminus, nBins+1);
	// std::vector<size_t> binCounts_posBins_Bplus, binCounts_negBins_Bplus;
	// createZeros<size_t>(binCounts_posBins_Bplus, nBins+1);
	// createZeros<size_t>(binCounts_negBins_Bplus, nBins+1);
	// int currentBin{0};

	// // loop over each event's bin number
	// for (size_t i{0}; i < nEvents; i++){
	// 	// B- events
	// 	currentBin = binList_Bminus[i];
	// 	if (currentBin > 0){
	// 		binCounts_posBins_Bminus[currentBin]++;
	// 	}else{
	// 		binCounts_negBins_Bminus[abs(currentBin)]++;
	// 	}
	// 	//B+ events
	// 	currentBin = binList_Bplus[i];
	// 	if (currentBin > 0){
	// 		binCounts_posBins_Bplus[currentBin]++;
	// 	}else{
	// 		binCounts_negBins_Bplus[abs(currentBin)]++;
	// 	}
	// }
//********************

// //************ SAVE BIN FRACTIONS ************
	// std::ofstream outFile{outputFilename};
	// for ( size_t i{1}; i < nBins+1 ; i++){
	// 	outFile << "Bminus_" << i << "   Fix   " << binCounts_posBins_Bminus[i] << "  0\n"
	// 			<< "Bplus_" << i << "    Fix   " << binCounts_posBins_Bplus[i] << "  0\n\n"
	// 			<< "Bminus_-" << i << "  Fix   " << binCounts_negBins_Bminus[i] << "  0\n"
	// 			<< "Bplus_-" << i << "   Fix   " << binCounts_negBins_Bplus[i] << "  0\n\n\n";
	// };

	// outFile << "# counted bin fractions for each bin based off of: \n# "
	// 		<< NamedParameter<std::string>("Description" , " :( ", " explaination of xy input values used") << "\n\n\n";
//***************

//***************** CALCULATE SOMETHING: FI **************
	// std::vector<real_t> binCounts_posBins, binCounts_negBins;
	// createZeros<real_t>(binCounts_posBins, nBins+1);
	// createZeros<real_t>(binCounts_negBins, nBins+1);
	// real_t totalN{0};

	// for (int i{1}; i < 2; i++){

	// 	// EventList_type data((inputData + std::to_string(i) + ".root:DalitzEventList").c_str() , signalType); 
	// 	EventList_type data((inputData + ":DalitzEventList").c_str() , signalType); 
	// 	INFO("opening file: " << inputData );
	// 	totalN += data.size();
	// 	INFO("total nEvents is now: " << totalN);

	// 	for (auto event : data){
	// 	// for (int i{0}; i < 20; i++){
	// 		// Event event = data[i];
	// 		int currentBin = nearestBinIndex(event, mMinus, mPlus, bins);
	// 		if (currentBin > 0){
	// 			binCounts_posBins[currentBin]++;
	// 		}else{
	// 			binCounts_negBins[abs(currentBin)]++;
	// 		}
	// 	}
	// }

//********************

//************ SAVE SOMETHING: FI ********
	// std::ofstream outFile{outputFilename};
	// for ( size_t i{1}; i < nBins+1 ; i++){
	// 	outFile << "F_" << i << "   Fix   " << binCounts_posBins[i] / totalN << "  0\n"
	// 			<< "F_-" << i << "    Fix   " << binCounts_negBins[i] / totalN << "  0\n\n";
	// }

	// outFile.close();

//*********

//******** PLOT SOMETHING: ONLY A CERTAIN BIN ON DALITZ PLANE 	*******

	// EventList_type data((inputData+"_1.root:DalitzEventList").c_str() , signalType); 
	// INFO("opening file: " << inputData << " 1");


	// std::vector<int> binList;


	// for (auto event : data){
	// 	int currentBin = nearestBinIndex(event, mMinus, mPlus, bins);
	// 	binList.push_back(currentBin);
	// }
	
	// for (int b{1}; b < nBins+1; b++){
	// 	TH2D* histo_posBin = new TH2D("pos", ";s_{-};s_{+};", 200, 0, 3.5, 200, 0, 3.5);
	// 	TH2D* histo_negBin = new TH2D("neg", "", 200, 0, 3.5, 200, 0, 3.5);

	// 	for (int e{0}; e < nEvents; e++){
	// 		int thisBin = binList[e];
	// 		if (abs(thisBin) == b){
	// 			if (thisBin > 0){
	// 				histo_posBin->Fill(data[e].s(0, 1), data[e].s(0, 2));
	// 			}else{
	// 				histo_negBin->Fill(data[e].s(0, 1), data[e].s(0, 2));
	// 			}
	// 		}
	// 	}

	// 	cv->cd(b);
	// 	real_t max = std::max(histo_posBin->GetMaximum(), histo_negBin->GetMaximum());
	// 	histo_negBin->SetMaximum(max);
	// 	histo_posBin->SetMaximum(max);
	// 	histo_posBin->SetTitle(("Bin: "+std::to_string(b)).c_str());
	// 	histo_posBin->SetStats(false);
	// 	// histo_posBin->SetMarkerColor(colPurple);
	// 	// histo_negBin->SetMarkerColor(colBlue);
	// 	histo_posBin->Draw("colz");
	// 	histo_negBin->Draw("SAME colz");

	// 	INFO("histo for bin " << b << " has:  " << histo_posBin->GetEntries() << " entries");
	// 	INFO("histo for bin -" << b << " has: " << histo_negBin->GetEntries() << " entries");

	// }


//****************



// *************** PLOT SOMETHING: BINNING SCHEME **********
	// TH2* frame = new TH2D("frame", "",100, 0, 3.5, 100, 0, 3.5);
	// std::string title{ NamedParameter<std::string>("Title", "", "title of plot") };
	// frame->SetTitle(title.c_str());
	// frame->GetXaxis()->SetTitle("m(K^{0}_{S} #pi^{-}) = m_{-}");
	// frame->GetYaxis()->SetTitle("m(K^{0}_{S} #pi^{+}) = m_{+}");
	// frame->SetStats(0);

	// frame->SetMinimum(0.5);
	// frame->SetMaximum(8.5);
	// double frameBins[nBins] = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
	// frame->SetContour(nBins, frameBins); // setting ticks to right numbers
	// frame->Fill(4, 5); // out of range point to make colourbar draw
	// frame->Draw("colz");
	// int colours[nBins+1] = {0, 2, 3, 4, 801, 6, 38, 419, 881};
	// gStyle->SetPalette(nBins+1, colours);


	// INFO("binning and plotting")
	// for (auto event : data ){
	// 		int bin { nearestBinIndex(event, mMinus, mPlus, bins) };
	// 		double s01{event.s(0, 1)};
	// 		double s02{event.s(0, 2)};
	// 		double x[1] = {s01};
	// 		double y[1] = {s02};
	// 		TGraph* g = new TGraph(1, x, y); 
	// 		g->SetMarkerStyle(6); // g->SetMarkerSize(1);
	// 		g->SetMarkerColor(colours[abs(bin)]);
	// 		g->Draw("SAME P");
	// 	}

//*********

// ************* PLOT SOMTHING: A_CP  *****************
	// INFO("Plotting the graph");

	// note: currently asking for real_t inputs as last used from calculating decay widths,
	//		 would need to cahange and add a 1.0* if using int inputs from MPS or other
	// auto get_ACP = [](const real_t& Nplus, const real_t& Nminus){
	// 	return (Nplus - Nminus) / (Nplus + Nminus); 
	// };
	// auto get_err = [](const real_t& Nplus, const real_t& Nminus){
	// 	return std::sqrt(Nplus*Nminus / (Nplus + Nminus)); 
	// };

	// INFO("made lambdas");

	// TH1D* histo_ACP = new TH1D("henry", "A_{CP}", 17, -8.5, 8.5);

	// for (size_t i{1}; i < nBins+1; i++){
	// 	INFO("doing bin " << i << " with value " << get_ACP(binCounts_posBins_Bplus[i], binCounts_negBins_Bminus[i]));
	// 	histo_ACP->SetBinContent(9+i, get_ACP(binCounts_posBins_Bplus[i], binCounts_negBins_Bminus[i]));
	// 	histo_ACP->SetBinContent(9-i, get_ACP(binCounts_negBins_Bplus[i], binCounts_posBins_Bminus[i]));
	// 	INFO(" and error  " << get_err(binCounts_posBins_Bplus[i], binCounts_negBins_Bminus[i]) << "\n");
	// 	histo_ACP->SetBinError(9+i, get_err(binCounts_posBins_Bplus[i], binCounts_negBins_Bminus[i]));
	// 	histo_ACP->SetBinError(9-i, get_err(binCounts_negBins_Bplus[i], binCounts_posBins_Bminus[i]));
	// }

	// histo_ACP->SetStats(false);
	// histo_ACP->GetXaxis()->SetTitle("bin number");
	// histo_ACP->GetYaxis()->SetTitle("A_{CP}");
	// histo_ACP->SetLineColor(618);
	// histo_ACP->SetLineWidth(1.5);
	// histo_ACP->Draw("");

	// TLine* horizontal = new TLine(-8.5, 0, 8.5, 0);
	// horizontal->Draw("SAME");

//***************

// **************** PLOT SOMETHING: Ci / Si COMPARISON PLOT ***************
	// TEllipse* circle = new TEllipse(0, 0, 1, 1);
	// TLine* verticle = new TLine(0, -1.2, 0, 1.2);
	// TLine* horizontal = new TLine(-1.2, 0, 1.2, 0);
	// circle->SetLineStyle(9);
	// verticle->SetLineStyle(9);
	// horizontal->SetLineStyle(9);
	// circle->Draw("SAME");
	// verticle->Draw("SAME");
	// horizontal->Draw("SAME");

	// TLegend* legend = new TLegend(0.73 , 0.75, 0.85, 0.9);


	// for (int i{1}; i<9; i++){

	// 	// opt 1 pos bins
	// 	double x[1] = { MPS_opt1["c_"+std::to_string(i)]->mean() };
	// 	double y[1] = { MPS_opt1["s_"+std::to_string(i)]->mean() };
	// 	TGraph* g = new TGraph(1, x, y); 
	// 	g->SetMarkerStyle(attEmptyCircleMarker); //g->SetMarkerSize(4);  
	// 	g->SetMarkerColor(colours[i]); g->SetLineColorAlpha(colours[i], 1);
	// 	g->Draw("SAME P");

	// 	// opt 2 no errs pos bins
	// 	double x2[1] = {MPS_opt2["c_"+std::to_string(i)]->mean()};
	// 	double y2[1] = {MPS_opt2["s_"+std::to_string(i)]->mean()};
	// 	TGraph* g2 = new TGraph(1, x2, y2); 
	// 	g2->SetMarkerStyle(attEmptyTriangleMarker); //g2->SetMarkerSize(4); 
	// 	g2->SetMarkerColor(colours[i]); g2->SetLineColorAlpha(colours[i], 0);
	// 	g2->Draw("SAME P");

	// 	// opt 2 pos bins
	// 	// double x3[1] = {MPS_opt2["c_"+std::to_string(i)]->mean()};
	// 	// double y3[1] = {MPS_opt2["s_"+std::to_string(i)]->mean()};
	// 	// double ex[1] = {MPS_opt2["ERRc_"+std::to_string(i)]->mean()};
	// 	// double ey[1] = {MPS_opt2["ERRs_"+std::to_string(i)]->mean()};
	// 	// // TGraph* g3 = new TGraph(1, x3, y3);
	// 	// TGraph* g3 = new TGraphErrors(1, x3, y3, ex, ey); 
	// 	// g3->SetMarkerStyle(attFillCircleMarker); //g3->SetMarkerSize(4);
	// 	// g3->SetMarkerColor(colours[i]); g3->SetLineColorAlpha(colours[i], 1);
	// 	// g3->Draw("SAME P");

	// 	// opt 2 neg bins
	// 	// double x4[1] = {MPS_opt2["c_-"+std::to_string(i)]->mean()};
	// 	// double y4[1] = {MPS_opt2["s_-"+std::to_string(i)]->mean()};
	// 	// TGraph* g4 = new TGraph(1, x4, y4);
	// 	// // TGraph* g4 = new TGraphErrors(1, x4, y4, ex, ey); 
	// 	// g4->SetMarkerStyle(attEmptyCircleMarker); //g2->SetMarkerSize(4);
	// 	// g4->SetMarkerColor(colours[i]); g4->SetLineColorAlpha(colours[i], 0);
	// 	// g4->Draw("SAME P");

	// 	// ** compare a third set?
	// 		// opt 3 pos bins
	// 		// double x5[1] = {MPS_opt3["c_"+std::to_string(i)]->mean()};
	// 		// double y5[1] = {MPS_opt3["s_"+std::to_string(i)]->mean()};
	// 		// double ex5[1] = {MPS_opt3["ERRc_"+std::to_string(i)]->mean()};
	// 		// double ey5[1] = {MPS_opt3["ERRs_"+std::to_string(i)]->mean()};
	// 		// // TGraph* g5 = new TGraph(1, x5, y5);
	// 		// TGraph* g5 = new TGraphErrors(1, x5, y5, ex5, ey5); 
	// 		// g5->SetMarkerStyle(attFillSquareMarker); //g3->SetMarkerSize(4);
	// 		// g5->SetMarkerColor(colours[i]); g5->SetLineColorAlpha(colours[i], 1);
	// 		// g5->Draw("SAME P");
	// 	//**
	// }
	// INFO("adding entries to legend");
	// double x[1] = {1}; double y[1] = {1};
	// TGraph* g_forLegend = new TGraph(1, x, y);
	// TGraph* g3_forLegend = new TGraph(1, x, y);
	// // TGraph* g5_forLegend = new TGraph(1, x, y);
	// g_forLegend->SetMarkerStyle(attEmptyCircleMarker); g_forLegend->SetMarkerColor(1); g_forLegend->SetLineColorAlpha(1, 0);
	// g3_forLegend->SetMarkerStyle(attEmptyTriangleMarker); g3_forLegend->SetMarkerColor(1); g3_forLegend->SetLineColorAlpha(1, 0);
	// // g5_forLegend->SetMarkerStyle(attFillSquareMarker); g5_forLegend->SetMarkerColor(1); g5_forLegend->SetLineColorAlpha(1, 0);
	// legend->AddEntry(g_forLegend, (label1).c_str());
	// // legend->AddEntry(g2, (label1+" -ve bins").c_str());
	// legend->AddEntry(g3_forLegend, (label2).c_str());
	// // legend->AddEntry(g4, (label2+" -ve bins").c_str());
	// // legend->AddEntry(g5_forLegend, (label3).c_str());
	// // legend->AddEntry(g6, (label3+" -ve bins").c_str());
	
	// legend->Draw("SAME");
//*************

// ************** PLOT SOMETHING: Fi COMPARISON PLOT ***********

	// TLegend* legend = new TLegend(0.75 , 0.8, 0.9, 0.9);
	// TH1D* histo1 = new TH1D("harry", "Comparing F_{i}", 17, -8.5, 8.5);
	// TH1D* histo2 = new TH1D("harriet", "Comparing F_{i}", 17, -8.5, 8.5);
	// TH1D* histo3 = new TH1D("horatio", "Comparing F_{i}", 17, -8.5, 8.5);
	// TH1D* histo4 = new TH1D("hermentine", "Comparing F_{i}", 17, -8.5, 8.5);


	// for (int i{1}; i<nBins+1; i++){

	// 	histo1->SetBinContent(9-i, MPS_opt1["F_-"+std::to_string(i)]->mean());
	// 	histo1->SetBinContent(9+i, MPS_opt1["F_"+std::to_string(i)]->mean());

	// 	histo2->SetBinContent(9-i, MPS_opt2["F_-"+std::to_string(i)]->mean());
	// 	histo2->SetBinContent(9+i, MPS_opt2["F_"+std::to_string(i)]->mean());

	// 	// histo3->SetBinContent(9-i, MPS_opt1["Fbar_"+std::to_string(i)]->mean());
	// 	// histo3->SetBinContent(9+i, MPS_opt1["Fbar_-"+std::to_string(i)]->mean());

	// 	// histo4->SetBinContent(9-i, MPS_opt2["Fbar_"+std::to_string(i)]->mean());
	// 	// histo4->SetBinContent(9+i, MPS_opt2["Fbar_-"+std::to_string(i)]->mean());
	// }

	// histo1->SetLineColor(colPurple);
	// histo1->SetLineWidth(3);
	// histo2->SetLineColor(colBlue);
	// histo2->SetLineWidth(3); histo2->SetLineStyle(2);
	// histo3->SetLineColor(colBlue);
	// histo4->SetLineColor(colGreen);

	// histo1->SetStats(false);
	// histo1->GetXaxis()->SetTitle("bin Number, i");
	// histo1->GetYaxis()->SetTitle("F_{i}");
	// std::string title = NamedParameter<std::string>("plotTitle", "");
	// histo1->SetTitle(title.c_str());

	// histo1->Draw("");
	// histo2->Draw("SAME");
	// // histo3->Draw("SAME");
	// // histo4->Draw("SAME");

	// legend->AddEntry(histo1, (label1).c_str());
	// legend->AddEntry(histo2, (label2).c_str());
	// // legend->AddEntry(histo3, (label2).c_str());
	// // legend->AddEntry(histo4, (label2 + " #bar{F}_{-i}").c_str());

	// legend->Draw("SAME");

//*************

// ************** PLOT SOMETHING: Fi ERROR COMPARISON PLOT ***********

	// TLegend* legend = new TLegend(0.75 , 0.8, 0.9, 0.9);
	// TH1D* histo1 = new TH1D("harry", "title", 17, -8.5, 8.5);
	// TH1D* histo2 = new TH1D("harriet", "title", 17, -8.5, 8.5);


	// for (int i{1}; i<nBins+1; i++){

	// 	// histo1->SetBinContent(9-i, nEvents*MPS_opt1["Fbar_"+std::to_string(i)]->mean());
	// 	// histo2->SetBinContent(9-i, nEvents*MPS_opt1["Fbar_"+std::to_string(i)]->mean());

	// 	// histo1->SetBinContent(9+i, nEvents*MPS_opt1["F_"+std::to_string(i)]->mean());
	// 	// histo2->SetBinContent(9+i, nEvents*MPS_opt1["F_"+std::to_string(i)]->mean());

	// 	histo1->SetBinError(9-i, std::sqrt(nEvents * MPS_opt1["Fbar_"+std::to_string(i)]->mean()));
	// 	histo1->SetBinError(9+i, std::sqrt(nEvents * MPS_opt1["F_"+std::to_string(i)]->mean()));


	// 	histo2->SetBinError(9-i, nEvents*MPS_opt2["errFbar_"+std::to_string(i)]->mean());
	// 	histo2->SetBinError(9+i, nEvents*MPS_opt2["errF_"+std::to_string(i)]->mean());
	// }

	// histo1->SetLineColor(colPurple);
	// histo1->SetMarkerColor(colPurple);
	// histo2->SetLineColor(colBlue);
	// histo2->SetFillColor(colBlue);
	// histo2->SetMarkerColor(colBlue);

	// histo1->SetStats(false);
	// histo1->SetTitle(("Bin yield to compare uncertainty on F_{i}, nEvents = "+std::to_string(nEvents)).c_str());
	// histo1->GetXaxis()->SetTitle("bin Number, i");
	// histo1->GetYaxis()->SetTitle("uncertainty in bin yield");
	// histo2->SetStats(false);
	// histo2->SetTitle(("Bin yield to compare uncertainty on F_{i}, nEvents = "+std::to_string(nEvents)).c_str());
	// histo2->GetXaxis()->SetTitle("bin Number, i");
	// histo2->GetYaxis()->SetTitle("bin population");

	// histo1->Draw("E1");
	// histo2->Draw("SAME E2");
	// histo1->Draw("SAME E1");

	// legend->AddEntry(histo1, label1.c_str());
	// legend->AddEntry(histo2, label2.c_str());
	// legend->Draw("SAME");
//*********

//*************** PLOT SOMETHING: Fis x3 variations ************
	// // remember to get rid of the colourbar in the frame setup, and to check limits are right
	// TLegend* legend = new TLegend(0.1, 0.7, 0.4, 0.9);
	// for (int i{1}; i<9; i++){
	// 	real_t x[2] ={-i, i};
	// 	real_t y_sumA[2] = {MPS["Famps_-"+std::to_string(i)]->mean(), MPS["Famps_"+std::to_string(i)]->mean()};
	// 	real_t y_sumAbar[2] = {MPS["Fbaramps_"+std::to_string(i)]->mean(), MPS["Fbaramps_-"+std::to_string(i)]->mean()};
	// 	real_t y_count[2] = {MPS["Fbar_"+std::to_string(i)]->mean(), MPS["F_"+std::to_string(i)]->mean()};

	// 	TGraph* g_sumA = new TGraph(2, x, y_sumA); 
	// 	g_sumA->SetMarkerSize(3); g_sumA->SetMarkerStyle(2); g_sumA->SetMarkerColor(colDarkGreen); 
	// 	g_sumA->Draw("SAME P");

	// 	TGraph* g_sumAbar = new TGraph(2, x, y_sumAbar); 
	// 	g_sumAbar->SetMarkerSize(3); g_sumAbar->SetMarkerStyle(5); g_sumAbar->SetMarkerColor(colGreen);
	// 	g_sumAbar->Draw("SAME P");

	// 	TGraph* g_count = new TGraph(2, x, y_count); 
	// 	g_count->SetMarkerSize(1); g_count->SetMarkerStyle(24); g_count->SetMarkerColor(colPurple);
	// 	g_count->Draw("SAME P");


	// 	if (i==1){
	// 		g_sumA->Draw("SAME P");
	// 		legend->AddEntry(g_sumA, "F calculated with Amplitudes", "P");
	// 		legend->AddEntry(g_sumAbar, "#bar{F} calculated with Amplitudes", "P");
	// 		legend->AddEntry(g_count, "F by counting events", "P");
	// 	}
	// }
	// legend->Draw("SAME");
//*************

// ************* PLOT SOMETHING: FILL A HISTOGRAM OF BIN POPULATIONS (6/3/23) **********

	// TH1D* histo_posBins_Bminus = new TH1D("PBBM","Bminus", 8, 0.5, 8.5);
	// TH1D* histo_negBins_Bminus = new TH1D("NGBM","Bminus", 8, 0.5, 8.5);
	// TH1D* histo_posBins_Bplus = new TH1D("PBBP","Bplus", 8, 0.5, 8.5);
	// TH1D* histo_negBins_Bplus = new TH1D("NBBP","Bplus", 8, 0.5, 8.5);


	// for (size_t i{0}; i < nEvents; i++){
	// 	if (binList_Bminus[i] > 0){
	// 		histo_posBins_Bminus->Fill(binList_Bminus[i]);
	// 	}else{
	// 		histo_negBins_Bminus->Fill(abs(binList_Bminus[i]));
	// 	}
	// 	if (binList_Bplus[i] > 0){
	// 		histo_posBins_Bplus->Fill(binList_Bplus[i]);
	// 	}else{
	// 		histo_negBins_Bplus->Fill(abs(binList_Bplus[i]));
	// 	}
	// }



	// // - set up colours etc
	// histo_posBins_Bminus->SetLineColor(876);
	// histo_posBins_Bminus->SetLineWidth(4);
	// histo_posBins_Bplus->SetLineColor(822);
	// histo_posBins_Bplus->SetLineWidth(4);
	// histo_negBins_Bminus->SetLineColor(874);
	// histo_negBins_Bminus->SetLineWidth(4);
	// histo_negBins_Bminus->SetLineStyle(2);
	// histo_negBins_Bplus->SetLineColor(813);
	// histo_negBins_Bplus->SetLineWidth(4);
	// histo_negBins_Bplus->SetLineStyle(2);

	// // - legend
	// TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
	// legend->AddEntry(histo_posBins_Bminus, "B^{-}, +ve bins");
	// legend->AddEntry(histo_posBins_Bplus, "B^{+}, +ve bins");
	// legend->AddEntry(histo_negBins_Bminus, "B^{-}, -ve bins");
	// legend->AddEntry(histo_negBins_Bplus, "B^{+}, -ve bins");


	// // - pos bins B- serving as the frame, set up title etc here
	// histo_posBins_Bminus->GetXaxis()->SetTitle("abs(bin number)");
	// histo_posBins_Bminus->GetYaxis()->SetTitle("Number of events");
	// histo_posBins_Bminus->SetStats(false);

	// histo_posBins_Bminus->Draw();
	// histo_posBins_Bplus->Draw("SAME");
	// histo_negBins_Bminus->Draw("SAME");
	// histo_negBins_Bplus->Draw("SAME");
	// legend->Draw("SAME");


//**********


// ************* PLOT SOMETHING: FILL A HISTOGRAM OF *FIT RESULT* = BIN POPS (6/3/23) **********
	// // requires multiple pads:
	// TPad* pad_small = new TPad("paddy", "", 0.0, 0.0, 1.0, 0.2);
	// TPad* pad_main = new TPad("padricia", "", 0.0, 0.2, 1.0, 1.0);
	// pad_main->Divide(2, 1);		pad_main->Draw();
	// pad_small->Divide(2, 1);	pad_small->Draw();

	// TH1D* histo_Bminus_nObs = new TH1D("OBM","Bminus", 17, -8.5, 8.5);
	// TH1D* histo_Bminus_nExp = new TH1D("EBM","Bminus", 17, -8.5, 8.5);
	// TH1D* histo_Bplus_nObs = new TH1D("OBP","Bplus", 17, -8.5, 8.5);
	// TH1D* histo_Bplus_nExp = new TH1D("EBP","Bplus", 17, -8.5, 8.5);

	// TH1D* histo_Bminus_diff = new TH1D("DBM","", 17, -8.5, 8.5);
	// TH1D* histo_Bplus_diff = new TH1D("DBP","", 17, -8.5, 8.5);


	// for (int i{1}; i < nBins+1; i++){ // IGNORE SIGN COMPARISON WARNING
	// 		// fill nObs histograms for all signs of bins and B's
	// 		histo_Bminus_nObs->SetBinContent(9+i, nObs_posBins_Bminus[i]);
	// 		histo_Bminus_nObs->SetBinContent(9-i, nObs_negBins_Bminus[i]);

	// 		histo_Bplus_nObs->SetBinContent(9+i, nObs_posBins_Bplus[i]);
	// 		histo_Bplus_nObs->SetBinContent(9-i, nObs_negBins_Bplus[i]);

	// 		histo_Bminus_nExp->SetBinContent(9+i, nExp_posBins_Bminus[i]);
	// 		histo_Bminus_nExp->SetBinContent(9-i, nExp_negBins_Bminus[i]);

	// 		histo_Bplus_nExp->SetBinContent(9+i, nExp_posBins_Bplus[i]);
	// 		histo_Bplus_nExp->SetBinContent(9-i, nExp_negBins_Bplus[i]);

	// 		histo_Bplus_diff->SetBinContent(9+i, nObs_posBins_Bplus[i] / nExp_posBins_Bplus[i]);
	// 		histo_Bplus_diff->SetBinContent(9-i, nObs_negBins_Bplus[i] / nExp_negBins_Bplus[i]);

	// 		histo_Bminus_diff->SetBinContent(9+i, nObs_posBins_Bminus[i] / nExp_posBins_Bminus[i]);
	// 		histo_Bminus_diff->SetBinContent(9-i, nObs_negBins_Bminus[i] / nExp_negBins_Bminus[i]);
	// 	}

	// histo_Bminus_diff->SetBinContent(9, 1.);
	// histo_Bplus_diff->SetBinContent(9, 1.);

	// // - set up colours etc
	// histo_Bminus_nObs->SetLineColor(619);
	// histo_Bminus_nObs->SetFillColor(619);

	// histo_Bplus_nObs->SetLineColor(619);
	// histo_Bplus_nObs->SetFillColor(619);

	// histo_Bminus_nExp->SetLineColor(2);
	// histo_Bminus_nExp->SetLineWidth(2);
	// histo_Bminus_nExp->SetLineStyle(2);

	// histo_Bplus_nExp->SetLineColor(2);
	// histo_Bplus_nExp->SetLineWidth(2);
	// histo_Bplus_nExp->SetLineStyle(2);

	// histo_Bminus_diff->SetLineColor(2);
	// histo_Bplus_diff->SetLineColor(2);
	// histo_Bminus_diff->SetStats(false);
	// histo_Bplus_diff->SetStats(false);

	// histo_Bplus_diff->GetYaxis()->SetNdivisions(5, 0, 0);
	// histo_Bminus_diff->GetYaxis()->SetNdivisions(5, 0, 0);
	// histo_Bplus_diff->GetYaxis()->SetLabelSize(0.1);
	// histo_Bminus_diff->GetYaxis()->SetLabelSize(0.1);

	// // - legend
	// // TLegend* legend = new TLegend(0.425, 0.8, 0.575, 0.9); // rB = 1 angles 0
	// TLegend* legend = new TLegend(0.75, 0.8, 0.9, 0.9); // all others so far
	// legend->AddEntry(histo_Bminus_nObs, "DATA");
	// legend->AddEntry(histo_Bplus_nExp, "FIT");

	// // - a horizontal line for the difference
	// TLine* horizontal = new TLine(-8.5, 1, 8.5, 1);
	// horizontal->SetLineStyle(9);

	// // - pos bins B- serving as the frame, set up title etc here
	// histo_Bminus_nObs->GetXaxis()->SetTitle("abs(bin number)");
	// histo_Bminus_nObs->GetYaxis()->SetTitle("Number of events");
	// histo_Bminus_nObs->SetStats(false);
	// histo_Bplus_nObs->GetXaxis()->SetTitle("abs(bin number)");
	// histo_Bplus_nObs->GetYaxis()->SetTitle("Number of events");
	// histo_Bplus_nObs->SetStats(false);


	// pad_main->cd(1); // cv->cd(1);
	// histo_Bminus_nObs->Draw("");
	// histo_Bminus_nExp->Draw("SAME");
	// legend->Draw("SAME");
	
	// pad_small->cd(1);
	// histo_Bminus_diff->Draw("");
	// horizontal->Draw("SAME");
	
	// pad_main->cd(2); // cv->cd(2);
	// histo_Bplus_nObs->Draw("");
	// histo_Bplus_nExp->Draw("SAME");
	
	// pad_small->cd(2);
	// histo_Bplus_diff->Draw("");
	// horizontal->Draw("SAME");


//**********



//************ COUNT SOME FI'S FROM FLAT EVENTS *********
	// // ************* MAKE PIDGEON HOLES FOR THE COUNTS **********
	// RealVec_t F_posBins, Fbar_posBins;
	// createZeros(F_posBins, nBins+1);	createZeros(Fbar_posBins, nBins+1);	
	// RealVec_t F_negBins, Fbar_negBins;
	// createZeros(F_negBins, nBins+1);	createZeros(Fbar_negBins, nBins+1);	


	// // *************** OPTION TO LOOP OVER TWO EVENT LISTS **************
	// for ( auto flatData : flatDataVector) {
	// 	// //************** BIN THE EVENTS **************
		// INFO("Binning events");
		// std::vector<int> binList; // each will have a bin number for each event
		// int eventBin{0};

		// ProfileClock  CLOCK_binning;
		// CLOCK_binning.start();
		// for (auto& event:flatData){
		// 	eventBin = nearestBinIndex(event, mMinus, mPlus, bins);
		// 	binList.push_back(eventBin);
		// }


	// 	CLOCK_binning.stop();
	// 	INFO("took " << CLOCK_binning << "ms to bin " << nEvents << " events");

	// 	//****** PREPARE LISTS OF EACH VALUE WE WANT *******
		// std::vector<size_t> binCounts_posBins, binCounts_negBins;
		// createZeros<size_t>(binCounts_posBins, nBins+1);
		// createZeros<size_t>(binCounts_negBins, nBins+1);

		// INFO("counted bins to get Fs");
	// 	//******** COUNT THOSE BINS *********
		// for (size_t i{0}; i < nEvents; i++){
		// 	int currentBin = binList[i];
		// 	if (currentBin > 0){
		// 		binCounts_posBins[currentBin]++;
		// 	}else{
		// 		binCounts_negBins[abs(currentBin)]++;
		// 	}
	// 	}

	// 	//******** STORE FOR LATER ******
	// 	if (F_posBins[1] == 0){
	// 		for (int i{1}; i < nBins+1 ; i++){ //Will create a sign warning-> ignore: i needs to be an int
	// 			F_posBins[i] = 1.0 * binCounts_posBins[i] / nEvents;
	// 			F_negBins[i] = 1.0 * binCounts_negBins[i] / nEvents;
	// 		}
	// 	}else{
	// 		for (int i{1}; i < nBins+1 ; i++){ //Will create a sign warning-> ignore: i needs to be an int
	// 			Fbar_posBins[i] = 1.0 * binCounts_posBins[i] / nEvents;
	// 			Fbar_negBins[i] = 1.0 * binCounts_negBins[i] / nEvents;
	// 		}
	// 	}
	// }

	// // ********* WRITE THE RESULTS TO FILE **********
	// // NOTE: normalisation of F happens inline of output in this loop
	// std::ofstream outFile{ outputFilename };
	// if ( !outFile.good()) {
	// 	ERROR("Unable to open output file to print result");
	// 	return 1;
	// }


	// INFO("Writing results to options file " << outputFilename);
    // // note using definition Fbar_i = F_-i so that code in fitter won't need to change
	// for (int i{1}; i < nBins+1 ; i++){ //Will create a sign warning-> ignore: i needs to be an int

	// 	outFile << "F_" << i << "     Fix     " << F_posBins[i] << "     0\n"
	// 			<< "Fbar_" << i << "  Fix     " << Fbar_posBins[i] << "     0\n"
    //             << "F_" << -1*i << "    Fix     " << F_negBins[i] << "     0\n"
    //             << "Fbar_" << -1*i << " Fix     " << Fbar_negBins[i] << "     0\n\n";
	// }
	// outFile.close();
//*************

//*************** SAVE TO A ROOT FILE **************
	TFile * outputFile = TFile::Open(outputFilename.c_str(), "RECREATE"); 
	outputFile->cd();
	finalDalitz->Write();
	finalBias->Write();

	outputFile->Close();

//******************

//*************** FOR ANY PLOTTING: SAVE THE CANVAS ************
	// INFO("saving plot to file " << plotFilename);
	// cv->SaveAs((plotFilename).c_str());
//******************

	return 0;
}
