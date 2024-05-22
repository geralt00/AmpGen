#include "AmpGen/WriteToTree.h"


#include "AmpGen/CoherentSum.h"
#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/Formalism.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MyStructs.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PhaseCorrection.h"
#include "AmpGen/PhaseCorrection_KLpipi.h"
#include "AmpGen/Types.h"


#if ENABLE_AVX
	using EventList_type = AmpGen::EventListSIMD;
#else
	using EventList_type = AmpGen::EventList; 
#endif

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"


#include <fstream>
#include <string>
#include <vector>


using namespace AmpGen;

template <typename T> void createZeros(std::vector<T>& list, const size_t& nItems){
	for (size_t i{0}; i < nItems;i++){
		list.push_back(0);
	}
	return;
}

// Must've been copied in from AmpGen, have edited QMI to look more like Tim's again
void AmpGen::writePlots(EventList_type& events, size_t nBins, std::string& prefix)
{		
	// make plots of each variable that's a 2 particle invariant mass
	auto plots = events.makeDefaultProjections(PlotOptions::Bins(nBins), PlotOptions::LineColor(kBlack), PlotOptions::Prefix(prefix)); 
	for ( auto& plot : plots ){
		plot->Write();}

	// amke a plot of each Dalitz option -> 3 cobinations of 3 inv masses
	if( NamedParameter<bool>("plots_2d",true) ){
		auto proj = events.eventType().defaultProjections(nBins);
		for( size_t i = 0 ; i < proj.size(); ++i ){
			for( size_t j = i+1 ; j < proj.size(); ++j ){ 
			events.makeProjection( Projection2D(proj[i], proj[j]), PlotOptions::LineColor(kBlack), PlotOptions::Prefix(prefix))->Write( ); 
			}
		}
	}

	return;
}

void AmpGen::writePlots(EventList& events, size_t nBins, std::string& prefix)
{	
	// make plots of each variable that's a 2 particle invariant mass
	auto plots = events.makeDefaultProjections(PlotOptions::Bins(nBins), PlotOptions::LineColor(kBlack), PlotOptions::Prefix(prefix)); 
	for ( auto& plot : plots ){
		plot->Write();}

	// amke a plot of each Dalitz option -> 3 cobinations of 3 inv masses
	if( NamedParameter<bool>("plots_2d",true) ){
		auto proj = events.eventType().defaultProjections(nBins);
		for( size_t i = 0 ; i < proj.size(); ++i ){
			for( size_t j = i+1 ; j < proj.size(); ++j ){ 
			events.makeProjection( Projection2D(proj[i], proj[j]), PlotOptions::LineColor(kBlack), PlotOptions::Prefix(prefix))->Write( ); 
			}
		}
	}

	return;
}

void AmpGen::writeAmplitudes(CoherentSum& A, CoherentSum& Abar, const EventList_type& acceptedEvents)
{
	TTree* tree = new TTree("AValues", "AValues");
	double AReal{0}, AImag{0}, AbarReal{0}, AbarImag{0};

	// A.prepare();
	// Abar.prepare();
	auto getA = A.amplitudeEvaluator( &acceptedEvents );
	auto getAbar = Abar.amplitudeEvaluator( &acceptedEvents );

	TBranch* bAReal = tree->Branch("A_real", &AReal);
	TBranch* bAImag = tree->Branch("A_imag", &AImag);
	TBranch* bAbarReal = tree->Branch("Abar_real", &AbarReal);
	TBranch* bAbarImag = tree->Branch("Abar_imag", &AbarImag);

	for (auto& leaf : acceptedEvents){
		AReal = getA(leaf).real();
		AImag = getA(leaf).imag();
		AbarReal = -1*getAbar(leaf).real();
		AbarImag = -1*getAbar(leaf).imag();
		tree->Fill();
	}
	tree->Write();
	return;
}


/*
void AmpGen::writeDalitzVariables(const EventList_type& events,const std::string& name="temp")
{
	// find inv masses
	auto s01 = [](Event& evt){ return evt.s(0, 1);};
	auto s02 = [](Event& evt){ return evt.s(0, 2);};
	auto s12 = [](Event& evt){ return evt.s(1, 2);};

	// lammbda to get the angles
	auto theta = [](Event& evt, size_t i, size_t j){ 
			double px_i = evt[4*i+0];
			double py_i = evt[4*i+1];
			double pz_i = evt[4*i+2];
			double px_j = evt[4*j+0];
			double py_j = evt[4*j+1];
			double pz_j = evt[4*j+2];
			double dot_ij = px_i * px_j + py_i * py_j + pz_i * pz_j;
			double p_i = std::pow(px_i * px_i + py_i * py_i + pz_i * pz_i, 0.5);
			double p_j = std::pow(px_j * px_j + py_j * py_j + pz_j * pz_j, 0.5);
			return std::acos(dot_ij/(p_i * p_j));
	} ;

	// new tree for the tfile
    TTree tree(("Dalitz_"+name).c_str(), ("Dalitz_"+ name).c_str());
	double x, y, z, theta12, theta13, theta23;
	auto bs01(tree.Branch("s01", &x));
	auto bs02(tree.Branch("s02", &y));
	auto bs12(tree.Branch("s12", &z));
	auto b_theta12(tree.Branch("theta_12", &theta12));
	auto b_theta13(tree.Branch("theta_13", &theta13));
	auto b_theta23(tree.Branch("theta_23", &theta23));

	// write each thing into the file
	for (auto &evt : events){
			x=s01(evt);
			y=s02(evt);
			z=s12(evt);
			theta12 = theta(evt, 1, 2);
			theta13 = theta(evt, 1, 3);
			theta23 = theta(evt, 2, 3);
			tree.Fill();
	}
	tree.Write();
	return;
}
*/
void AmpGen::writeDalitzVariables(const EventList_type& events,const std::string& name="temp")
{
	// find inv masses
	auto s01 = [](Event& evt){ return evt.s(0, 1);};
	auto s02 = [](Event& evt){ return evt.s(0, 2);};
	auto s12 = [](Event& evt){ return evt.s(1, 2);};

	// lammbda to get the angles
	auto theta = [](Event& evt, size_t i, size_t j){ 
			double px_i = evt[4*i+0];
			double py_i = evt[4*i+1];
			double pz_i = evt[4*i+2];
			double px_j = evt[4*j+0];
			double py_j = evt[4*j+1];
			double pz_j = evt[4*j+2];
			double dot_ij = px_i * px_j + py_i * py_j + pz_i * pz_j;
			double p_i = std::pow(px_i * px_i + py_i * py_i + pz_i * pz_i, 0.5);
			double p_j = std::pow(px_j * px_j + py_j * py_j + pz_j * pz_j, 0.5);
			return std::acos(dot_ij/(p_i * p_j));
	} ;

	// new tree for the tfile
    TTree tree(("Dalitz_"+name).c_str(), ("Dalitz_"+ name).c_str());
	double x, y, z, theta12, theta13, theta23;
	auto bs01(tree.Branch("s01", &x));
	auto bs02(tree.Branch("s02", &y));
	auto bs12(tree.Branch("s12", &z));
	auto b_theta12(tree.Branch("theta_12", &theta12));
	auto b_theta13(tree.Branch("theta_13", &theta13));
	auto b_theta23(tree.Branch("theta_23", &theta23));

	// write each thing into the file
	for (auto &evt : events){
			x=s01(evt);
			y=s02(evt);
			z=s12(evt);
			theta12 = theta(evt, 1, 2);
			theta13 = theta(evt, 1, 3);
			theta23 = theta(evt, 2, 3);
			tree.Fill();
	}
	tree.Write();
	return;
}

void AmpGen::writeFitResult(const std::string& filename, MinuitParameterSet& MPS, const size_t& nEvents,
							std::vector<int>& binCounts_posBins_Bminus, std::vector<int>& binCounts_posBins_Bplus , std::vector<int>& binCounts_negBins_Bminus , std::vector<int>& binCounts_negBins_Bplus)
{
	const size_t nBins = NamedParameter<size_t>("nBins", 8, "Number of bins in binning file, assumed 8" );
	signs signs{};

	std::ofstream outFile{filename};
	outFile << "nEvents    Fix    " << nEvents << "  0\n\n"
			<< "xPlus      Fix    " << MPS["xPlus"]->mean() << "  0\n" 
			<< "xMinus     Fix    " << MPS["xMinus"]->mean() << "  0\n" 
			<< "yPlus      Fix    " << MPS["yPlus"]->mean() << "  0\n" 
			<< "yMinus     Fix    " << MPS["yMinus"]->mean() << "  0\n\n\n\n";
		

	for (int i{1}; i < nBins+1; i++){
		outFile << "nObs_" << i << "_Bminus    Fix    " << binCounts_posBins_Bminus[i] << "  0\n"
				<< "nObs_" << i << "_Bplus     Fix    " << binCounts_posBins_Bplus[i] << "  0\n"
				<< "nObs_-" << i << "_Bminus   Fix    " << binCounts_negBins_Bminus[i] << "  0\n"
				<< "nObs_-" << i << "_Bplus    Fix    " << binCounts_negBins_Bplus[i] << "  0\n\n";

		real_t ci = MPS["c_" + std::to_string(i)]->mean(); 
		real_t si = MPS["s_" + std::to_string(i)]->mean();
		real_t Fi = MPS["F_" + std::to_string(i)]->mean();
		real_t Fbari = MPS["Fbar_" + std::to_string(i)]->mean();

		outFile << "nExp_" << i << "_Bminus    Fix    " << nEvents*decayWidth_Integrated_xy(signs.BminusSign, signs.posBinSign, Fi, Fbari, MPS, ci, si) << "  0\n"
				<< "nExp_" << i << "_Bplus     Fix    " << nEvents*decayWidth_Integrated_xy(signs.BplusSign, signs.posBinSign, Fi, Fbari, MPS, ci, si) << "  0\n";
		
		ci = MPS["c_-" + std::to_string(i)]->mean(); 
		si = MPS["s_-" + std::to_string(i)]->mean();
		Fi = MPS["F_-" + std::to_string(i)]->mean();
		Fbari = MPS["Fbar_-" + std::to_string(i)]->mean();

		outFile << "nExp_-" << i << "_Bminus   Fix    " << nEvents*decayWidth_Integrated_xy(signs.BminusSign, signs.negBinSign, Fbari, Fi, MPS, ci, si) << "  0\n"
				<< "nExp_-" << i << "_Bplus    Fix    " << nEvents*decayWidth_Integrated_xy(signs.BplusSign, signs.negBinSign, Fbari, Fi, MPS, ci, si) << "  0\n\n\n";

	}
	outFile.close();
	return;
}

void AmpGen::writeBinnedFitPlot(MinuitParameterSet& MPS, const size_t& nEvents,
						  std::vector<int>& binCounts_posBins_Bminus, std::vector<int>& binCounts_posBins_Bplus , std::vector<int>& binCounts_negBins_Bminus , std::vector<int>& binCounts_negBins_Bplus)
{
	const std::string plotFile = NamedParameter<std::string>("Plots", "Result.root", "Name of the output root file to save fit result plots to");
	const size_t nBins = NamedParameter<size_t>("nBins", 8, "Number of bins in binning file, assumed 8" );
	signs signs{};

	TH1D* histo_Bminus_nObs = new TH1D("obs_Bminus","Bminus", 17, -8.5, 8.5);
	TH1D* histo_Bminus_nExp = new TH1D("exp_Bminus","Bminus", 17, -8.5, 8.5);
	TH1D* histo_Bplus_nObs = new TH1D("obs_Bplus","Bplus", 17, -8.5, 8.5);
	TH1D* histo_Bplus_nExp = new TH1D("exp_Bplus","Bplus", 17, -8.5, 8.5);

	TH1D* histo_Bminus_diff = new TH1D("diff_Bminus","", 17, -8.5, 8.5);
	TH1D* histo_Bplus_diff = new TH1D("diff_Bplus","", 17, -8.5, 8.5);


	// pidgeon holes for the decay widths so that can sum over them without losing
	std::vector<real_t> decayWidths_posBins_Bplus, decayWidths_posBins_Bminus, decayWidths_negBins_Bplus, decayWidths_negBins_Bminus;
	createZeros<real_t>(decayWidths_posBins_Bplus, nBins+1); createZeros<real_t>(decayWidths_posBins_Bminus, nBins+1);
	createZeros<real_t>(decayWidths_negBins_Bplus, nBins+1); createZeros<real_t>(decayWidths_negBins_Bminus, nBins+1);
	real_t sumDecayWidths_Bplus = 0;
	real_t sumDecayWidths_Bminus = 0;
	for (size_t i{1}; i < nBins+1; i++){ 
		real_t ci = MPS["c_" + std::to_string(i)]->mean(); 
		real_t si = MPS["s_" + std::to_string(i)]->mean();
		real_t Fi = MPS["F_" + std::to_string(i)]->mean();
		real_t Fbari = MPS["Fbar_" + std::to_string(i)]->mean();
		decayWidths_posBins_Bminus[i] = decayWidth_Integrated_xy(signs.BminusSign, signs.posBinSign, Fi, Fbari, MPS, ci, si);
		decayWidths_posBins_Bplus[i] = decayWidth_Integrated_xy(signs.BplusSign, signs.posBinSign, Fi, Fbari, MPS, ci, si);

		decayWidths_negBins_Bminus[i] = decayWidth_Integrated_xy(signs.BminusSign, signs.negBinSign, Fbari, Fi, MPS, ci, si);
		decayWidths_negBins_Bplus[i] = decayWidth_Integrated_xy(signs.BplusSign, signs.negBinSign, Fbari, Fi, MPS, ci, si);


		sumDecayWidths_Bminus += decayWidths_negBins_Bminus[i] + decayWidths_posBins_Bminus[i]; 
		sumDecayWidths_Bplus += decayWidths_negBins_Bplus[i] + decayWidths_posBins_Bplus[i];
	}


	for (int i{1}; i < nBins+1; i++){ // IGNORE SIGN COMPARISON WARNING
		// fill nObs histograms for all signs of bins and B's
		histo_Bminus_nObs->SetBinContent(9+i, binCounts_posBins_Bminus[i]);
		histo_Bminus_nObs->SetBinContent(9-i, binCounts_negBins_Bminus[i]);

		histo_Bplus_nObs->SetBinContent(9+i, binCounts_posBins_Bplus[i]);
		histo_Bplus_nObs->SetBinContent(9-i, binCounts_negBins_Bplus[i]);
		
		real_t nExp_posBins_Bminus{ nEvents * decayWidths_posBins_Bminus[i] / sumDecayWidths_Bminus };
		real_t nExp_posBins_Bplus{ nEvents * decayWidths_posBins_Bplus[i] / sumDecayWidths_Bplus };
		real_t nExp_negBins_Bminus{ nEvents * decayWidths_negBins_Bminus[i] / sumDecayWidths_Bminus };
		real_t nExp_negBins_Bplus{ nEvents * decayWidths_negBins_Bplus[i] / sumDecayWidths_Bplus };

		histo_Bminus_nExp->SetBinContent(9+i, nExp_posBins_Bminus);
		histo_Bminus_nExp->SetBinContent(9-i, nExp_negBins_Bminus);

		histo_Bplus_nExp->SetBinContent(9+i, nExp_posBins_Bplus);
		histo_Bplus_nExp->SetBinContent(9-i, nExp_negBins_Bplus);

		histo_Bplus_diff->SetBinContent(9+i, (nExp_posBins_Bplus - binCounts_posBins_Bplus[i]) / std::sqrt(nExp_posBins_Bplus));
		histo_Bplus_diff->SetBinContent(9-i, (nExp_negBins_Bplus - binCounts_negBins_Bplus[i]) / std::sqrt(nExp_negBins_Bplus));

		histo_Bminus_diff->SetBinContent(9+i, (nExp_posBins_Bminus - binCounts_posBins_Bminus[i]) / std::sqrt(nExp_posBins_Bminus));
		histo_Bminus_diff->SetBinContent(9-i, (nExp_negBins_Bminus - binCounts_negBins_Bminus[i]) / std::sqrt(nExp_negBins_Bminus));
	}

	histo_Bminus_diff->SetBinContent(9, 0.);
	histo_Bplus_diff->SetBinContent(9, 0.);

	TFile* outputFile = TFile::Open(plotFile.c_str(), "RECREATE"); 
	outputFile->cd();
	histo_Bminus_nExp->Write();
	histo_Bminus_nObs->Write();
	histo_Bplus_nExp->Write();
	histo_Bplus_nObs->Write();
	histo_Bplus_diff->Write();
	histo_Bminus_diff->Write();
	outputFile->Close();


	return;
}

void AmpGen::writeAcpPlot(const std::string& filename, 
						  std::vector<int>& binCounts_posBins_Bminus, std::vector<int>& binCounts_posBins_Bplus , std::vector<int>& binCounts_negBins_Bminus , std::vector<int>& binCounts_negBins_Bplus)
{
	TCanvas* cv = new TCanvas("cv", "A_{CP}", 500, 500);

	const size_t nBins = NamedParameter<size_t>("nBins", 8, "Number of bins in binning file, assumed 8" );
	
	auto get_ACP = [](const int& Nplus, const int& Nminus){
		return 1.0*(Nplus - Nminus) / (Nplus + Nminus); 
	};
	auto get_err = [](const int& Nplus, const int& Nminus){
		return 2.0 * std::sqrt(1.0*Nplus*Nminus) / std::pow(Nplus + Nminus , 1.5); 
	};

	TH1D* histo_ACP = new TH1D("henry", "A_{CP}", 17, -8.5, 8.5);

	for (size_t i{1}; i < nBins+1; i++){
		histo_ACP->SetBinContent(9+i, get_ACP(binCounts_posBins_Bplus[i], binCounts_negBins_Bminus[i]));
		histo_ACP->SetBinContent(9-i, get_ACP(binCounts_negBins_Bplus[i], binCounts_posBins_Bminus[i]));

		histo_ACP->SetBinError(9+i, get_err(binCounts_posBins_Bplus[i], binCounts_negBins_Bminus[i]));
		histo_ACP->SetBinError(9-i, get_err(binCounts_negBins_Bplus[i], binCounts_posBins_Bminus[i]));
	}

	histo_ACP->SetStats(false);
	histo_ACP->GetXaxis()->SetTitle("bin number");
	histo_ACP->GetYaxis()->SetTitle("A_{CP}");
	histo_ACP->SetLineColor(618);
	histo_ACP->SetLineWidth(1.5);
	histo_ACP->Draw("");

	TLine* horizontal = new TLine(-8.5, 0, 8.5, 0);
	horizontal->Draw("SAME");

	cv->SaveAs((filename).c_str());

	return;
}
void AmpGen::writeUnbinnedFitPlot(CoherentSum& A, CoherentSum& Abar, PhaseCorrection& phaseCorrection, MinuitParameterSet& MPS, EventList_type& eventsMC,   EventList_type& events_Bplus, EventList_type& events_Bminus)
{
	signs signs{};
	EventType signalType{NamedParameter<std::string>("EventType", "", "Signal Type to generate, should be for D0 -> Ks0pi-pi+")}; 
	const std::string plotFile = NamedParameter<std::string>("Plots", "plot.root", "Name of the output root file for plots");
	
	// lambda with only event as input to plot pdf
	const std::function<double( const Event&)> pdfLambda_Bplus = [&signs, &A, &Abar, &phaseCorrection, &MPS]( const Event& event){
		return totalAmplitudeSquared_XY(signs.BplusSign, event, A, Abar, phaseCorrection, MPS);
	};
	const std::function<double( const Event&)> pdfLambda_Bminus = [&signs, &A, &Abar, &phaseCorrection, &MPS]( const Event& event){
		return totalAmplitudeSquared_XY(signs.BminusSign, event, A, Abar, phaseCorrection, MPS);
	};


	// projections of s01 = Kpi-, s01 = Kpi+
	auto projection_s01 = signalType.projection(100, {0, 1}); //100 = number of bins in the histos
	auto projection_s02 = signalType.projection(100, {0, 2});  // {0, x} = particles to take the inv mass of as the x axis

	// histograms of pdf/data, B+/-, pi+/- = 8 in total
	TH1D* histo_data_Bplus_piMinus = events_Bplus.makeProjection(projection_s01, PlotOptions::Prefix("Bplus_data"));
	TH1D* histo_pdf_Bplus_piMinus = eventsMC.makeProjection(projection_s01, PlotOptions::Prefix("Bplus_pdf"), WeightFunction(pdfLambda_Bplus));
	TH1D* histo_data_Bminus_piMinus = events_Bminus.makeProjection(projection_s01, PlotOptions::Prefix("Bminus_data"));
	TH1D* histo_pdf_Bminus_piMinus = eventsMC.makeProjection(projection_s01, PlotOptions::Prefix("Bminus_pdf"), WeightFunction(pdfLambda_Bminus));

	TH1D* histo_data_Bplus_piPlus = events_Bplus.makeProjection(projection_s02, PlotOptions::Prefix("Bplus_data"));
	TH1D* histo_pdf_Bplus_piPlus = eventsMC.makeProjection(projection_s02, PlotOptions::Prefix("Bplus_pdf"), WeightFunction(pdfLambda_Bplus));
	TH1D* histo_data_Bminus_piPlus = events_Bminus.makeProjection(projection_s02, PlotOptions::Prefix("Bminus_data"));
	TH1D* histo_pdf_Bminus_piPlus = eventsMC.makeProjection(projection_s02, PlotOptions::Prefix("Bminus_pdf"), WeightFunction(pdfLambda_Bminus));



	// save them all to a root file:
	TFile * outputFile = TFile::Open(plotFile.c_str(), "RECREATE"); 
	outputFile->cd();
	histo_data_Bplus_piMinus->Write();
	histo_pdf_Bplus_piMinus->Write();
	histo_data_Bminus_piMinus->Write();
	histo_pdf_Bminus_piMinus->Write();
	histo_data_Bplus_piPlus->Write();
	histo_pdf_Bplus_piPlus->Write();
	histo_data_Bminus_piPlus->Write();
	histo_pdf_Bminus_piPlus->Write();
	
	if (NamedParameter<bool>("DoPhaseCorrectionPlot", true, "save a histo of the final phase correction to the root file")){
		TH2D* histoPC = getPChisto(phaseCorrection, eventsMC);
		histoPC->Write();
	}
	
	
	outputFile->Close();

	return;
}
void AmpGen::writeUnbinnedFitPlot_LHCb(std::vector<CoherentSum>& A, std::vector<CoherentSum>& Abar, PhaseCorrection& phaseCorrection, MinuitParameterSet& MPS, std::vector<EventList_type>& eventsMC, std::vector<EventList_type>& Datalist_LHCb, std::string magtype)
{
	signs signs{};
	EventType signalType{NamedParameter<std::string>("EventType", "", "Signal Type to generate, should be for D0 -> Ks0pi-pi+")}; 
	const std::string plotFile = NamedParameter<std::string>("Plots", "plot.root", "Name of the output root file for plots");
	const std::function<double( const Event&)> pdfLambda_Bplus_dk = [&signs, &A, &Abar, &phaseCorrection, &MPS]( const Event& event){
		return totalAmplitudeSquared_XY(signs.BplusSign, event, A[0], Abar[0], phaseCorrection, MPS);
	};

	const std::function<double( const Event&)> pdfLambda_Bminus_dk = [&signs, &A, &Abar, &phaseCorrection, &MPS]( const Event& event){
		return totalAmplitudeSquared_XY(signs.BminusSign, event, A[1], Abar[1], phaseCorrection, MPS);
	};

	const std::function<double( const Event&)> pdfLambda_Bplus_dpi = [&signs, &A, &Abar, &phaseCorrection, &MPS]( const Event& event){
		return totalAmplitudeSquared_DPi_XY(signs.BplusSign, event, A[2], Abar[2], phaseCorrection, MPS);
	};
	
	const std::function<double( const Event&)> pdfLambda_Bminus_dpi = [&signs, &A, &Abar, &phaseCorrection, &MPS]( const Event& event){
		return totalAmplitudeSquared_DPi_XY(signs.BminusSign, event, A[3], Abar[3], phaseCorrection, MPS);
	};

	std::vector<function<double( const Event&)>> pdfLambdas{pdfLambda_Bplus_dk, pdfLambda_Bminus_dk, pdfLambda_Bplus_dpi, pdfLambda_Bminus_dpi};

	const std::function<double( const Event&)> sWeighted_data = []( const Event& event){
		return event.weight_bkg();
	};
	//
	//I'd like to make a dalitz variable rather than Histo
	
	// projections of s01 = Kpi-, s01 = Kpi+
	auto projection_s01 = signalType.projection(100, {0, 1}); //100 = number of bins in the histos
	auto projection_s02 = signalType.projection(100, {0, 2});  // {0, x} = particles to take the inv mass of as the x axis

	// histograms of pdf/data, B+/-, pi+/- = 8 in total
	TH1D* BDK_histo_data_Bplus_piMinus = Datalist_LHCb[1].makeProjection(projection_s01, PlotOptions::Prefix(magtype+"_BDK_plus_data"), WeightFunction(sWeighted_data));
	TH1D* BDK_histo_pdf_Bplus_piMinus = eventsMC[0].makeProjection(projection_s01, PlotOptions::Prefix(magtype+"_BDK_plus_pdf"), WeightFunction(pdfLambdas[0]));
	TH1D* BDK_histo_data_Bminus_piMinus = Datalist_LHCb[0].makeProjection(projection_s01, PlotOptions::Prefix(magtype+"_BDK_minus_data"), WeightFunction(sWeighted_data));
	TH1D* BDK_histo_pdf_Bminus_piMinus = eventsMC[1].makeProjection(projection_s01, PlotOptions::Prefix(magtype+"_BDK_minus_pdf"), WeightFunction(pdfLambdas[1]));

	TH1D* BDK_histo_data_Bplus_piPlus = Datalist_LHCb[1].makeProjection(projection_s02, PlotOptions::Prefix(magtype+"_BDK_plus_data"), WeightFunction(sWeighted_data));
	TH1D* BDK_histo_pdf_Bplus_piPlus = eventsMC[0].makeProjection(projection_s02, PlotOptions::Prefix(magtype+"_BDK_plus_pdf"), WeightFunction(pdfLambdas[0]));
	TH1D* BDK_histo_data_Bminus_piPlus = Datalist_LHCb[0].makeProjection(projection_s02, PlotOptions::Prefix(magtype+"_BDK_minus_data"), WeightFunction(sWeighted_data));
	TH1D* BDK_histo_pdf_Bminus_piPlus = eventsMC[1].makeProjection(projection_s02, PlotOptions::Prefix(magtype+"_BDK_minus_pdf"), WeightFunction(pdfLambdas[1]));
	
	// histograms of pdf/data, B+/-, pi+/- = 8 in total
	TH1D* BDPi_histo_data_Bplus_piMinus = Datalist_LHCb[3].makeProjection(projection_s01, PlotOptions::Prefix(magtype+"_BDPi_plus_data"), WeightFunction(sWeighted_data));
	TH1D* BDPi_histo_pdf_Bplus_piMinus = eventsMC[2].makeProjection(projection_s01, PlotOptions::Prefix(magtype+"_BDPi_plus_pdf"), WeightFunction(pdfLambdas[2]));
	TH1D* BDPi_histo_data_Bminus_piMinus = Datalist_LHCb[2].makeProjection(projection_s01, PlotOptions::Prefix(magtype+"_BDPi_minus_data"), WeightFunction(sWeighted_data));
	TH1D* BDPi_histo_pdf_Bminus_piMinus = eventsMC[3].makeProjection(projection_s01, PlotOptions::Prefix(magtype+"_BDPi_minus_pdf"), WeightFunction(pdfLambdas[3]));

	TH1D* BDPi_histo_data_Bplus_piPlus = Datalist_LHCb[3].makeProjection(projection_s02, PlotOptions::Prefix(magtype+"_BDPi_plus_data"), WeightFunction(sWeighted_data));
	TH1D* BDPi_histo_pdf_Bplus_piPlus = eventsMC[2].makeProjection(projection_s02, PlotOptions::Prefix(magtype+"_BDPi_plus_pdf"), WeightFunction(pdfLambdas[2]));
	TH1D* BDPi_histo_data_Bminus_piPlus = Datalist_LHCb[2].makeProjection(projection_s02, PlotOptions::Prefix(magtype+"_BDPi_minus_data"), WeightFunction(sWeighted_data));
	TH1D* BDPi_histo_pdf_Bminus_piPlus = eventsMC[3].makeProjection(projection_s02, PlotOptions::Prefix(magtype+"_BDPi_minus_pdf"), WeightFunction(pdfLambdas[3]));
	
	// save them all to a root file:
	TString create_type = "RECREATE";
	if (magtype == "DD") create_type = "UPDATE";
	TFile * outputFile = TFile::Open(plotFile.c_str(), create_type); 
	outputFile->cd();
	
	BDK_histo_data_Bplus_piMinus->Write();
	BDK_histo_pdf_Bplus_piMinus->Write();
	BDK_histo_data_Bminus_piMinus->Write();
	BDK_histo_pdf_Bminus_piMinus->Write();
	BDK_histo_data_Bplus_piPlus->Write();
	BDK_histo_pdf_Bplus_piPlus->Write();
	BDK_histo_data_Bminus_piPlus->Write();
	BDK_histo_pdf_Bminus_piPlus->Write();

	BDPi_histo_data_Bplus_piMinus->Write();
	BDPi_histo_pdf_Bplus_piMinus->Write();
	BDPi_histo_data_Bminus_piMinus->Write();
	BDPi_histo_pdf_Bminus_piMinus->Write();
	BDPi_histo_data_Bplus_piPlus->Write();
	BDPi_histo_pdf_Bplus_piPlus->Write();
	BDPi_histo_data_Bminus_piPlus->Write();
	BDPi_histo_pdf_Bminus_piPlus->Write();	


	if (NamedParameter<bool>("DoPhaseCorrectionPlot", true, "save a histo of the final phase correction to the root file")){
		TH2D* histoPC = getPChisto(phaseCorrection, eventsMC[0]);
		histoPC->Write();
	}
	
	
	
	outputFile->Close();

	return;
}

void AmpGen::writeUnbinnedFitPlot_BES(CoherentSum& A, CoherentSum& Abar, PhaseCorrection& phaseCorrection, MinuitParameterSet& MPS, std::vector<EventList_type>& eventsMC_list, std::vector<EventList_type>& Datalist_BESIII)
{
	signs signs{};
	EventType signalType{NamedParameter<std::string>("EventType", "", "Signal Type to generate, should be for D0 -> Ks0pi-pi+")}; 
	const std::string plotFile = NamedParameter<std::string>("Plots", "plot.root", "Name of the output root file for plots");
	const std::function<double( const Event&)> pdfLambda_CPodd = [&signs, &A, &Abar, &phaseCorrection]( const Event& event){

			return totalAmplitudeSquared_CP(signs.CPoddSign, A, Abar, phaseCorrection, event);

    };
	const std::function<double( const Event&)> pdfLambda_CPeven = [&signs, &A, &Abar, &phaseCorrection]( const Event& event){

			return totalAmplitudeSquared_CP(signs.CPevenSign, A, Abar, phaseCorrection, event);

    };
	auto pdfLambda_sig = [&signs, &A, &Abar, &phaseCorrection]( const Event& event_sig, const Event& event_tag){

			return totalAmplitudeSquared_BES(A, Abar, phaseCorrection, event_sig, event_tag);

    };
	const std::function<double( const Event&)> pdfLambda_flavour = [&signs, &A, &Abar]( const Event& event){


			return totalAmplitudeSquared_flavour(1, A, Abar, event);

	};
		const std::function<double( const Event&)> pdfLambda_flavourBar = [&signs, &A, &Abar]( const Event& event){


			return totalAmplitudeSquared_flavour(0, A, Abar, event);

	};


	//
	//I'd like to make a dalitz variable rather than Histo
	
	// projections of s01 = Kpi-, s01 = Kpi+
	auto projection_s01 = signalType.projection(100, {0, 1}); //100 = number of bins in the histos
	auto projection_s02 = signalType.projection(100, {0, 2});  // {0, x} = particles to take the inv mass of as the x axis
	// histograms of pdf/data, D2Kspipi, pi+/- = 4 in total
	TH1D* histo_data_flavour_piMinus = Datalist_BESIII[0].makeProjection(projection_s01, PlotOptions::Prefix("flavour_data"));
	TH1D* histo_pdf_flavour_piMinus = eventsMC_list[4].makeProjection(projection_s01, PlotOptions::Prefix("flavour_pdf"), WeightFunction(pdfLambda_flavour));
	TH1D* histo_data_flavour_piPlus = Datalist_BESIII[0].makeProjection(projection_s02, PlotOptions::Prefix("flavour_data"));
	TH1D* histo_pdf_flavour_piPlus = eventsMC_list[4].makeProjection(projection_s02, PlotOptions::Prefix("flavour_pdf"), WeightFunction(pdfLambda_flavour));

	TH1D* histo_data_flavlourBar_piMinus = Datalist_BESIII[1].makeProjection(projection_s01, PlotOptions::Prefix("flavourBar_data"));
	TH1D* histo_pdf_flavlourBar_piMinus = eventsMC_list[5].makeProjection(projection_s01, PlotOptions::Prefix("flavourBar_pdf"), WeightFunction(pdfLambda_flavourBar));
	TH1D* histo_data_flavlourBar_piPlus = Datalist_BESIII[1].makeProjection(projection_s02, PlotOptions::Prefix("flavourBar_data"));
	TH1D* histo_pdf_flavlourBar_piPlus = eventsMC_list[5].makeProjection(projection_s02, PlotOptions::Prefix("flavourBar_pdf"), WeightFunction(pdfLambda_flavourBar));

	TH1D* histo_data_CPodd_piMinus = Datalist_BESIII[2].makeProjection(projection_s01, PlotOptions::Prefix("CPodd_data"));
	TH1D* histo_pdf_CPodd_piMinus = eventsMC_list[0].makeProjection(projection_s01, PlotOptions::Prefix("CPodd_pdf"), WeightFunction(pdfLambda_CPodd));
	TH1D* histo_data_CPodd_piPlus = Datalist_BESIII[2].makeProjection(projection_s02, PlotOptions::Prefix("CPodd_data"));
	TH1D* histo_pdf_CPodd_piPlus = eventsMC_list[0].makeProjection(projection_s02, PlotOptions::Prefix("CPodd_pdf"), WeightFunction(pdfLambda_CPodd));
	TH1D* histo_data_CPeven_piMinus = Datalist_BESIII[3].makeProjection(projection_s01, PlotOptions::Prefix("CPeven_data"));
	TH1D* histo_pdf_CPeven_piMinus = eventsMC_list[1].makeProjection(projection_s01, PlotOptions::Prefix("CPeven_pdf"), WeightFunction(pdfLambda_CPeven));
	TH1D* histo_data_CPeven_piPlus = Datalist_BESIII[3].makeProjection(projection_s02, PlotOptions::Prefix("CPeven_data"));	
	TH1D* histo_pdf_CPeven_piPlus = eventsMC_list[1].makeProjection(projection_s02, PlotOptions::Prefix("CPeven_pdf"), WeightFunction(pdfLambda_CPeven));

	TH1D* histo_data_same_signal_piMinus = Datalist_BESIII[4].makeProjection(projection_s01, PlotOptions::Prefix("same_signal_data"));
	TH1D* histo_data_same_signal_piPlus = Datalist_BESIII[4].makeProjection(projection_s02, PlotOptions::Prefix("same_signal_data"));	
	TH1D* histo_data_same_tag_piMinus = Datalist_BESIII[5].makeProjection(projection_s01, PlotOptions::Prefix("same_tag_data"));
	TH1D* histo_data_same_tag_piPlus = Datalist_BESIII[5].makeProjection(projection_s02, PlotOptions::Prefix("same_tag_data"));	
	INFO("ALL the others are good");
	TH1D* histo_pdf_same_signal_piMinus = projection_s01.plot("same_signal_pdf");
	TH1D *histo_pdf_same_signal_piPlus = projection_s02.plot("same_signal_pdf");
	TH1D *histo_pdf_same_tag_piMinus = projection_s01.plot("same_tag_pdf");
	TH1D *histo_pdf_same_tag_piPlus = projection_s02.plot("same_tag_pdf");

	INFO("Events for same signal"<<eventsMC_list[2].size());
        for( unsigned i=0;i<eventsMC_list[3].size();++i ){
            Event evt_sig = eventsMC_list[2][i];
			Event evt_tag = eventsMC_list[3][i];
            //if( selection != nullptr && !selection(evt) ) continue;
			auto pos_sig_s01 = projection_s01(evt_sig);
			auto pos_sig_s02 = projection_s02(evt_sig);
			auto pos_tag_s01 = projection_s01(evt_tag);
			auto pos_tag_s02 = projection_s02(evt_tag);
            auto weight = pdfLambda_sig(eventsMC_list[2][i], eventsMC_list[3][i]);

            histo_pdf_same_signal_piMinus->Fill( pos_sig_s01, weight / evt_sig.genPdf() ); 
			histo_pdf_same_signal_piPlus->Fill( pos_sig_s02, weight / evt_sig.genPdf() );
			histo_pdf_same_tag_piMinus->Fill( pos_tag_s01,weight / evt_tag.genPdf() );
			histo_pdf_same_tag_piPlus->Fill( pos_tag_s02, weight / evt_tag.genPdf() );
        }

	INFO("ALL the dks are good");
	// save them all to a root file:
	TFile * outputFile = TFile::Open(plotFile.c_str(), "Update"); 
	outputFile->cd();
	//Since we utiliz all the Kspipi to obtain the phase correction, shall we do the chi2 fit with all the Kspipi events?

	histo_data_flavour_piMinus->Write();
	histo_pdf_flavour_piMinus->Write();
	histo_data_flavour_piPlus->Write();
	histo_pdf_flavour_piPlus->Write();
	histo_data_flavlourBar_piMinus->Write();
	histo_pdf_flavlourBar_piMinus->Write();
	histo_data_flavlourBar_piPlus->Write();
	histo_pdf_flavlourBar_piPlus->Write();
	histo_data_CPodd_piMinus->Write();
	histo_pdf_CPodd_piMinus->Write();
	histo_data_CPodd_piPlus->Write();
	histo_pdf_CPodd_piPlus->Write();
	histo_data_CPeven_piMinus->Write();
	histo_pdf_CPeven_piMinus->Write();
	histo_data_CPeven_piPlus->Write();
	histo_pdf_CPeven_piPlus->Write();

	histo_data_same_signal_piMinus->Write();
	histo_pdf_same_signal_piMinus->Write();
	histo_data_same_signal_piPlus->Write();
	histo_pdf_same_signal_piPlus->Write();
	histo_data_same_tag_piMinus->Write();
	histo_pdf_same_tag_piMinus->Write();
	histo_data_same_tag_piPlus->Write();
	histo_pdf_same_tag_piPlus->Write();
	outputFile->Close();

	return;
}

void AmpGen::writeUnbinnedFitPlot_KLpipi(CoherentSum& A, CoherentSum& Abar, PhaseCorrection& phaseCorrection, PhaseCorrection_KLpipi& PhaseCorrection_KLpipi, MinuitParameterSet& MPS, EventList_type& eventsMC, EventList_type& events_Bplus, EventList_type& events_Bminus)
{
	signs signs{};
	EventType signalType{NamedParameter<std::string>("EventType", "", "Signal Type to generate, should be for D0 -> Ks0pi-pi+")}; 
	const std::string plotFile = NamedParameter<std::string>("Plots", "plot.root", "Name of the output root file for plots");
	
	// lambda with only event as input to plot pdf
	const std::function<double( const Event&)> pdfLambda_Bplus = [&signs, &A, &Abar, &phaseCorrection, &MPS]( const Event& event){
		return totalAmplitudeSquared_XY(signs.BplusSign, event, A, Abar, phaseCorrection, MPS);
	};
	const std::function<double( const Event&)> pdfLambda_Bminus = [&signs, &A, &Abar, &phaseCorrection, &MPS]( const Event& event){
		return totalAmplitudeSquared_XY(signs.BminusSign, event, A, Abar, phaseCorrection, MPS);
	};


	// projections of s01 = Kpi-, s01 = Kpi+
	auto projection_s01 = signalType.projection(100, {0, 1}); //100 = number of bins in the histos
	auto projection_s02 = signalType.projection(100, {0, 2});  // {0, x} = particles to take the inv mass of as the x axis

	// histograms of pdf/data, B+/-, pi+/- = 8 in total
	TH1D* histo_data_Bplus_piMinus = events_Bplus.makeProjection(projection_s01, PlotOptions::Prefix("Bplus_data"));
	TH1D* histo_pdf_Bplus_piMinus = eventsMC.makeProjection(projection_s01, PlotOptions::Prefix("Bplus_pdf"), WeightFunction(pdfLambda_Bplus));
	TH1D* histo_data_Bminus_piMinus = events_Bminus.makeProjection(projection_s01, PlotOptions::Prefix("Bminus_data"));
	TH1D* histo_pdf_Bminus_piMinus = eventsMC.makeProjection(projection_s01, PlotOptions::Prefix("Bminus_pdf"), WeightFunction(pdfLambda_Bminus));

	TH1D* histo_data_Bplus_piPlus = events_Bplus.makeProjection(projection_s02, PlotOptions::Prefix("Bplus_data"));
	TH1D* histo_pdf_Bplus_piPlus = eventsMC.makeProjection(projection_s02, PlotOptions::Prefix("Bplus_pdf"), WeightFunction(pdfLambda_Bplus));
	TH1D* histo_data_Bminus_piPlus = events_Bminus.makeProjection(projection_s02, PlotOptions::Prefix("Bminus_data"));
	TH1D* histo_pdf_Bminus_piPlus = eventsMC.makeProjection(projection_s02, PlotOptions::Prefix("Bminus_pdf"), WeightFunction(pdfLambda_Bminus));



	// save them all to a root file:
	TFile * outputFile = TFile::Open(plotFile.c_str(), "RECREATE"); 
	outputFile->cd();
	histo_data_Bplus_piMinus->Write();
	histo_pdf_Bplus_piMinus->Write();
	histo_data_Bminus_piMinus->Write();
	histo_pdf_Bminus_piMinus->Write();
	histo_data_Bplus_piPlus->Write();
	histo_pdf_Bplus_piPlus->Write();
	histo_data_Bminus_piPlus->Write();
	histo_pdf_Bminus_piPlus->Write();
	
	if (NamedParameter<bool>("DoPhaseCorrectionPlot", true, "save a histo of the final phase correction to the root file")){
		TH2D* histoPC = getPChisto(phaseCorrection, eventsMC);
		TH2D* histoPC_KLpipi = getPChisto_KLpipi(PhaseCorrection_KLpipi, eventsMC);
		histoPC->Write();
		histoPC_KLpipi->Write();
	}
	
	
	outputFile->Close();

	return;
}

TH2D* AmpGen::getPChisto_bias(PhaseCorrection& PC, const EventList_type& events)
{
	size_t nHistoBins{ 500 };
	TH2D* finalDalitz = new TH2D("deltaC_bias", "#delta_{C}, phase correction;s_{-};s_{+};#delta_{C}/rad", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* fullCorrection = new TH2D("pc_bias", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* scaleHisto = new TH2D("n_bias", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);

	for ( auto event : events ){
		real_t deltaC{ PC.evalBias(event) };
		fullCorrection->Fill(event.s(0, 1), event.s(0, 2), deltaC);
		scaleHisto->Fill(event.s(0, 1), event.s(0, 2));
	}
	INFO("filled histo");

	for ( size_t i{1}; i < nHistoBins; i++){
		for ( size_t j{1}; j < nHistoBins; j++){
			real_t scaler = scaleHisto->GetBinContent(i, j);
			real_t correction = fullCorrection->GetBinContent(i, j);
			if (scaler == 0 && correction == 0){
				finalDalitz->SetBinContent(i, j, 0);
			}else{
				finalDalitz->SetBinContent(i, j, correction/scaler);
			}
		}

	}
	return finalDalitz;
}

TH2D* AmpGen::getPChisto_KLpipi_bias(PhaseCorrection_KLpipi& PC, const EventList_type& events)
{
	size_t nHistoBins{ 500 };
	TH2D* finalDalitz = new TH2D("deltaC_KLpipi", "#delta_{C}, phase correction;s_{-};s_{+};#delta_{C}/rad", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* fullCorrection = new TH2D("pc_KLpipi", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* scaleHisto = new TH2D("n_KLpipi", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);

	for ( auto event : events ){
		real_t deltaC{ PC.evalBias(event) };
		fullCorrection->Fill(event.s(0, 1), event.s(0, 2), deltaC);
		scaleHisto->Fill(event.s(0, 1), event.s(0, 2));
	}
	INFO("filled histo");

	for ( size_t i{1}; i < nHistoBins; i++){
		for ( size_t j{1}; j < nHistoBins; j++){
			real_t scaler = scaleHisto->GetBinContent(i, j);
			real_t correction = fullCorrection->GetBinContent(i, j);
			if (scaler == 0 && correction == 0){
				finalDalitz->SetBinContent(i, j, 0);
			}else{
				finalDalitz->SetBinContent(i, j, correction/scaler);
			}
		}

	}
	return finalDalitz;
}

TH2D* AmpGen::getPChisto(PhaseCorrection& PC, const EventList_type& events)
{
	size_t nHistoBins{ 500 };
	TH2D* finalDalitz = new TH2D("deltaC", "#delta_{C}, phase correction;s_{-};s_{+};#delta_{C}/rad", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* fullCorrection = new TH2D("pc", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* scaleHisto = new TH2D("n", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);

	for ( auto event : events ){
		real_t deltaC{ PC.eval(event) };
		fullCorrection->Fill(event.s(0, 1), event.s(0, 2), deltaC);
		scaleHisto->Fill(event.s(0, 1), event.s(0, 2));
	}
	INFO("filled histo");

	for ( size_t i{1}; i < nHistoBins; i++){
		for ( size_t j{1}; j < nHistoBins; j++){
			real_t scaler = scaleHisto->GetBinContent(i, j);
			real_t correction = fullCorrection->GetBinContent(i, j);
			if (scaler == 0 && correction == 0){
				finalDalitz->SetBinContent(i, j, 0);
			}else{
				finalDalitz->SetBinContent(i, j, correction/scaler);
			}
		}

	}
	return finalDalitz;
}

TH2D* AmpGen::getPChisto_KLpipi(PhaseCorrection_KLpipi& PC, const EventList_type& events)
{
	size_t nHistoBins{ 500 };
	TH2D* finalDalitz = new TH2D("deltaC_KLpipi", "#delta_{C}, phase correction;s_{-};s_{+};#delta_{C}/rad", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* fullCorrection = new TH2D("pc_KLpipi", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* scaleHisto = new TH2D("n_KLpipi", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);

	for ( auto event : events ){
		real_t deltaC{ PC.eval(event) };
		fullCorrection->Fill(event.s(0, 1), event.s(0, 2), deltaC);
		scaleHisto->Fill(event.s(0, 1), event.s(0, 2));
	}
	INFO("filled histo");

	for ( size_t i{1}; i < nHistoBins; i++){
		for ( size_t j{1}; j < nHistoBins; j++){
			real_t scaler = scaleHisto->GetBinContent(i, j);
			real_t correction = fullCorrection->GetBinContent(i, j);
			if (scaler == 0 && correction == 0){
				finalDalitz->SetBinContent(i, j, 0);
			}else{
				finalDalitz->SetBinContent(i, j, correction/scaler);
			}
		}

	}
	return finalDalitz;
}


TH2D* AmpGen::getFullPhaseHisto(PhaseCorrection& PC, CoherentSum& A, CoherentSum& Abar, const EventList_type& events)
{
	size_t nHistoBins{ 200 };
	TH2D* finalDalitz = new TH2D("totalDelta", "#delta_{C} + #delta_{D};s_{-};s_{+};phase/rad", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* fullPhase = new TH2D("pc", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* scaleHisto = new TH2D("n", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);

	for ( auto event : events ){
		real_t deltaC{ PC.eval(event) };
		real_t deltaD{ Deltadelta(event, A, Abar) };
		fullPhase->Fill(event.s(0, 1), event.s(0, 2), deltaC + deltaD);
		scaleHisto->Fill(event.s(0, 1), event.s(0, 2));
	}
	INFO("filled histo");

	for ( size_t i{1}; i < nHistoBins; i++){
		for ( size_t j{1}; j < nHistoBins; j++){
			real_t scaler = scaleHisto->GetBinContent(i, j);
			real_t correction = fullPhase->GetBinContent(i, j);
			if (scaler == 0 && correction == 0){
				finalDalitz->SetBinContent(i, j, 0);
			}else{
				finalDalitz->SetBinContent(i, j, correction/scaler);
			}
		}

	}
	return finalDalitz;
}


TH2D* AmpGen::getDeltaDhisto(CoherentSum& A, CoherentSum& Abar, const EventList_type& events)
{
	size_t nHistoBins{ 200 };
	TH2D* finalDalitz = new TH2D("deltaD", "#delta_{D};s_{-};s_{+};#delta_{D}/rad", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* fullPhase = new TH2D("pc", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);
	TH2D* scaleHisto = new TH2D("n", "", nHistoBins, 0.3, 3.2, nHistoBins, 0.3, 3.2);

	for ( auto event : events ){
		real_t deltaD{ Deltadelta(event, A, Abar) };
		fullPhase->Fill(event.s(0, 1), event.s(0, 2), deltaD);
		scaleHisto->Fill(event.s(0, 1), event.s(0, 2));
	}
	INFO("filled histo");

	for ( size_t i{1}; i < nHistoBins; i++){
		for ( size_t j{1}; j < nHistoBins; j++){
			real_t scaler = scaleHisto->GetBinContent(i, j);
			real_t correction = fullPhase->GetBinContent(i, j);
			if (scaler == 0 && correction == 0){
				finalDalitz->SetBinContent(i, j, 0);
			}else{
				finalDalitz->SetBinContent(i, j, correction/scaler);
			}
		}

	}
	return finalDalitz;
}
