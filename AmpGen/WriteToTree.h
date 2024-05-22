#ifndef AMPGEN_WRITETOTREE_H
#define AMPGEN_WRITETOTREE_H

#include "AmpGen/EventList.h"
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/PhaseCorrection.h"
#include "AmpGen/PhaseCorrection_KLpipi.h"
#include "AmpGen/Types.h"
#if ENABLE_AVX
	using EventList_type = AmpGen::EventListSIMD;
#else
	using EventList_type = AmpGen::EventList; 
#endif

#include <string>
#include <vector>

namespace AmpGen
{

	template <typename amp_t, typename events_t, typename TTree> void writeValues(events_t& list, amp_t& psi, TTree& t){
//	void writeValues(events_t& list, amp_t& psi, TTree& t){
		// lambda with only event as input to plot pdf
	std::vector<std::string> names{"Bplus", "Bminus", "CPodd", "CPeven", "same_sig", "same_tag"};


    //double x[6];
	for(int i=0;i<2;i++)
	{
		double x;
    	auto b(t.Branch(("weight_"+names[i]).c_str(), &x));
    	for (auto& evt : list){
    		x = psi[i](evt);
    	    t.Fill();
    	}
	}

        t.Write();

}
	void writePlots(EventList_type& events, size_t nBins, std::string& prefix);
	void writePlots(EventList& events, size_t nBins, std::string& prefix);

	void writeAmplitudes(CoherentSum& A, CoherentSum& Abar, const EventList_type& acceptedEvents);

	void writeDalitzVariables(const EventList_type& events,const std::string& name);
	void writeFitResult(const std::string& filename, MinuitParameterSet& MPS, const size_t& nEvents,
											std::vector<int>& binCounts_posBins_Bminus, std::vector<int>& binCounts_posBins_Bplus , std::vector<int>& binCounts_negBins_Bminus , std::vector<int>& binCounts_negBins_Bplus);

	void writeBinnedFitPlot(MinuitParameterSet& MPS, const size_t& nEvents,
					 					std::vector<int>& binCounts_posBins_Bminus, std::vector<int>& binCounts_posBins_Bplus , std::vector<int>& binCounts_negBins_Bminus , std::vector<int>& binCounts_negBins_Bplus);

	void writeAcpPlot(const std::string& filename, 
					  std::vector<int>& binCounts_posBins_Bminus, std::vector<int>& binCounts_posBins_Bplus , std::vector<int>& binCounts_negBins_Bminus , std::vector<int>& binCounts_negBins_Bplus);

	TH2D* getPChisto(PhaseCorrection& PC, const EventList_type& events);
	TH2D* getPChisto_KLpipi(PhaseCorrection_KLpipi& PC, const EventList_type& events);
	TH2D* getPChisto_bias(PhaseCorrection& PC, const EventList_type& events);
	TH2D* getPChisto_KLpipi_bias(PhaseCorrection_KLpipi& PC, const EventList_type& events);
	TH2D* getFullPhaseHisto(PhaseCorrection& PC, CoherentSum& A, CoherentSum& Abar, const EventList_type& events);
	TH2D* getDeltaDhisto(CoherentSum& A, CoherentSum& Abar, const EventList_type& events);
	void writeUnbinnedFitPlot_LHCb(std::vector<CoherentSum>& A, std::vector<CoherentSum>& Abar, PhaseCorrection& phaseCorrection, MinuitParameterSet& MPS, std::vector<EventList_type>& eventsMC,  std::vector<EventList_type>& Datalist_LHCb, std::string magtype);
	void writeUnbinnedFitPlot_BES(CoherentSum& A, CoherentSum& Abar, PhaseCorrection& phaseCorrection, MinuitParameterSet& MPS, std::vector<EventList_type>& eventsMC_list, std::vector<EventList_type>& Datalist_BESIII);


	void writeUnbinnedFitPlot(CoherentSum& A, CoherentSum& Abar, PhaseCorrection& phaseCorrection, MinuitParameterSet& MPS, EventList_type& eventsMC, EventList_type& events_Bplus, EventList_type& events_Bminus);
	void writeUnbinnedFitPlot_KLpipi(CoherentSum& A, CoherentSum& Abar, PhaseCorrection& phaseCorrection, PhaseCorrection_KLpipi& PhaseCorrection_KLpipi, MinuitParameterSet& MPS, EventList_type& eventsMC, EventList_type& events_Bplus, EventList_type& events_Bminus);

}




#endif