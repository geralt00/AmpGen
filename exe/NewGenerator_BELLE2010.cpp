#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/EventList.h"
#include "AmpGen/Formalism.h"
#include "AmpGen/Generator.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PhaseCorrection.h"
#include "AmpGen/MyStructs.h"
#include "AmpGen/WriteToTree.h"
#include "extern/DtoKpipiAmplitude.h"
#include "extern/Belle2010Amplitude.h"
#include "AmpGen/QMI.h"

#ifdef _OPENMP
	#include <omp.h>
#endif
#if ENABLE_AVX
	#include "AmpGen/EventListSIMD.h"
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
#include "TBranch.h"

using namespace AmpGen;

// TRandom 3 unaccepted input to the setRandom function so trying a separate function to mimic the syntax it gets passed around with in AmpGen.cpp. Not sure why this works and just calling it doesn't. 
void setTRandom(TRandom* rndm, auto& generator){
	generator.setRandom(rndm);
	return;
}


int main(int argc , char* argv[] ){

	// ******* ARGPARSING *******

	OptionsParser::setArgs( argc, argv );

	EventType eventType{NamedParameter<std::string>("EventType", "", "Signal Type to generate")};

	const size_t      blockSize = NamedParameter<size_t>("blockSize", 1e5, "Blocksize for generation");

	const size_t      nBins = NamedParameter<size_t>("nBins", 100, "Number of bins for projection histograms");

	const size_t      nEvents = NamedParameter<size_t>("nEvents", 15000, "number of events to generate, multiplies by the BR of each tag");
	
	const std::string outputFileName = NamedParameter<std::string>("Output", "outputFile.root", "output filename NO SUFFIX");

	const size_t      seed = NamedParameter<size_t>("Seed", 0, "Random seed for generation");
	TRandom3 rndm; rndm.SetSeed( seed );

	const size_t      nThreads = NamedParameter<size_t>("nCores", 8, "Number of threads to use" );
	#ifdef _OPENMP
		omp_set_num_threads( nThreads );
		INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
		omp_set_dynamic( 0 );
	#endif





	//******* CREATE AMPLITUDE PDF ********

	//** IMPORT MPS/AMPLITUDE MODEL **
	MinuitParameterSet MPS;
	MPS.loadFromStream(); // needs to have Kspipi model and event type and rB parameterised info
	AddCPConjugate(MPS);

	//** PREPARE AMPLITUDE ** 
	CoherentSum A(eventType, MPS);
	CoherentSum Abar(eventType.conj(true), MPS);
	A.prepare(); Abar.prepare();

	//** PREPARE AMPLITUDE OF KSPIPI_BEllE2010 **
    std::string kaon_type = "KS";

    std::cout << "Running KS-amplitude calculator\n";
    DtoKpipiAmplitude* amp;
    amp = new Belle2010Amplitude(kaon_type, "/publicfs/ucas/user/zengshh/LHCb/B2DK_D2K0S0PiPi/analysis/ampgen/extern/Belle2010/bin/cp_mult.txt");

	std::vector<Expression> Phi = QMI::dalitz(eventType);
    std::vector<CompiledExpression<real_t(const real_t*, const real_t*)> > cPhi;
    for (int i=0;i<Phi.size();i++){
            CompiledExpression<real_t(const real_t*, const real_t*)> cPhi_i(Phi[i], "Phi_" + std::to_string(i) , eventType.getEventFormat(), &MPS);
            cPhi_i.prepare(); cPhi_i.compile();
            cPhi.push_back(cPhi_i);
    }

	//** PREPARE THE PHASE CORRECTION **
	PhaseCorrection phaseCorrection{MPS}; 
    phaseCorrection.compilePolynomialExpressions(eventType);


	// ******** GENERATING EVENTS *********
	BtypeInfo Binfo{};

	TFile * outputFile = TFile::Open(outputFileName.c_str(), "RECREATE"); 
	outputFile->cd();

	// TEMP FOR GAMMAGEN MIMICKING:
	// std::string treeNames[2] = {"Bm2DKm", "Bp2DKp"};


    // loop through each tag and generate a root file of events
    for(unsigned i=0; i<Binfo.nTypes; i++)
	{

		int Bsign{Binfo.signs[i]};
		//** PREPARE PDF LAMBDA **

		auto AtotalCorrected = [&amp, &cPhi, &phaseCorrection, &MPS, &Bsign](Event event){
			real_t deltaC{0};
			if (phaseCorrection.doPolynomial()){
				deltaC = phaseCorrection.eval(event);
			}else if (phaseCorrection.doBias()){
				deltaC = phaseCorrection.evalBias(event);
			}
			return 0;
			//return totalAmplitudeSquared_rB_Belle2010(Bsign, event, *amp,cPhi, phaseCorrection, MPS);//Here is the problem
		};


		EventList acceptedEvents{eventType};
		
		Generator<PhaseSpace> generator(eventType);
		// generator.setRandom(rndm); // function from header - doesn't work :(
		setTRandom(&rndm, generator); // needed instead of straight using .setRandom - see explanation at function
		generator.setBlockSize(blockSize);
		generator.setNormFlag(true);

		generator.fillEventsUsingLambda(AtotalCorrected, acceptedEvents, nEvents);


		INFO("Writing Dalitz Event List to tree");
		acceptedEvents.tree((Binfo.prefixes[i]+"_DalitzEventList").c_str())->Write();
		// acceptedEvents.tree((treeNames[i]).c_str())->Write(); // gammaGen tree naming convention
		
		// ******* WRITING OUTPUT TO FILE ********
		if (NamedParameter<bool>("DoPlots", true)){
			INFO("Drawing plots to write to file");
			writePlots( acceptedEvents, nBins, Binfo.prefixes[i] );
			INFO( "Written Events");

			// do plots of phase correction? many plots at the moment, might delete later after used for diagnostics 26/9/23
			if (NamedParameter<bool>("DoPhaseCorrectionPlot", false)){
				INFO("Making a plot of the phase correction");
				EventList_type flatEvents = generator.generate(1000000);
				TH2D* PChisto = getPChisto_bias(phaseCorrection, flatEvents);
				PChisto->Write("deltaC_flat");
				INFO("Written plot to root file");
				TH2D* fullHisto = getFullPhaseHisto(phaseCorrection, A, Abar, flatEvents);
				fullHisto->Write("fullPhase_flat");
				TH2D* deltaDhisto = getDeltaDhisto( A, Abar, flatEvents);
				deltaDhisto->Write("deltaD_flat");
				
			}
			
		}

		if (NamedParameter<bool>("DoAmplitudes", false)){
			INFO("Writing the amplitudes to file");
			writeAmplitudes( A, Abar, acceptedEvents );
		}
		if (NamedParameter<bool>("DoPhasePlots", false)){
			INFO("Drawing the phase correction and strong phase difference plots");
			TH2D* histoPC = getPChisto_bias( phaseCorrection, acceptedEvents );
			TH2D* histoFull = getFullPhaseHisto( phaseCorrection, A, Abar, acceptedEvents );
			TH2D* histoDd = getDeltaDhisto( A, Abar, acceptedEvents );
			histoPC->Write();
			histoFull->Write();
			histoDd->Write();
		}





	}

	outputFile->Close();



	return 0;

}