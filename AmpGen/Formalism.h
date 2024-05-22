#ifndef AMPGEN_FORMALISM_H
#define AMPGEN_FORMALISM_H

#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/PhaseCorrection.h"
#include "AmpGen/PhaseCorrection_KLpipi.h"
#include "AmpGen/MyStructs.h"
#include "AmpGen/Types.h"
#include "extern/DtoKpipiAmplitude.h"
#include "extern/Belle2010Amplitude.h"
#include "extern/D0ToKLpipi2018.h"
#include "AmpGen/D0ToKSpipi2018.h"



#if ENABLE_AVX
  using EventList_type = AmpGen::EventListSIMD;
#else
  using EventList_type = AmpGen::EventList; 
#endif

#include <map>
#include <string>
#include <vector>


namespace AmpGen
{
	using twoDVector_t = std::vector<std::vector<real_t>>;

	//Basic functions
	real_t logPoisson( real_t Exp, real_t Obs);
	real_t Chi2( real_t Exp,  real_t Obs);
	// used by fitters:
	amplitudeInfo fillAmplitudeInfo(EventList_type& events, CoherentSum& A, CoherentSum& Abar);
	amplitudeInfo fillAmplitudeInfo_KLpipi(EventList_type& events, D0ToKLpipi2018& tKLpipi, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi);

	// used to get ci and si:
	real_t Deltadelta(const Event& event,CoherentSum& A, CoherentSum& Abar);
	real_t Deltadelta_KL(complex<double>& A, complex<double>& Abar);
	// integrate over all events to get |A||Abar|cos/sin(deltaD) for MD fitter
	std::pair<real_t, real_t> totalAmplitudeSquared_Integrated_crossTerm(amplitudeInfo& amplitudes, EventList_type& events, PhaseCorrection& PC );
	std::pair<real_t, real_t> totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudeInfo& amplitudes, EventList_type& events, PhaseCorrection& PC );

	// integrate within bins to get ci and si for SBF
	void doBinnedIntegration(std::vector<real_t>& ci_posBins, std::vector<real_t>& si_posBins, std::vector<real_t>& ci_negBins, std::vector<real_t>& si_negBins, 
							 amplitudeInfo& amplitudes, EventList_type& events, PhaseCorrection& PC);
	void doUnbinnedIntegration(amplitudeInfo& amplitudes, EventList_type& events, PhaseCorrection& PC, real_t& cosTerm, real_t& sinTerm, real_t& complicatedTerm);
	void doUnbinnedIntegration_temp(amplitudeInfo& amplitudes, EventList_type& events, PhaseCorrection& PC, real_t &normA, real_t &normAbar, real_t& cosTerm, real_t& sinTerm, real_t& complicatedTerm);
	void doUnbinnedIntegration_KLpipi(amplitudeInfo& amplitudes, amplitudeInfo& amplitudes_klpipi, EventList_type& events,EventList_type& events_klpipi, PhaseCorrection& PC, PhaseCorrection_KLpipi& PC_KLipi, real_t &KLpipi_normA, real_t &KLpipi_normAbar, real_t& cosTerm, real_t& sinTerm, real_t& complicatedTerm);
	//void doUnbinnedIntegration_KLpipi(amplitudeInfo& amplitudes, amplitudeInfo& amplitudes_klpipi, EventList_type& events,EventList_type& events_klpipi,PhaseCorrection& PC, PhaseCorrection_KLpipi& PC_KLpipi, real_t& cosTerm, real_t& sinTerm, real_t& complicatedTerm);
	void getNorm(EventList_type& events, D0ToKLpipi2018& tKLpipi, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi, real_t& normA, real_t& normAbar);
	double* get_bkg_fraction(bkgfraction & bkg_frac);
	double* get_norm_eff(normeff & norm_eff);

	// used as the pdf for MD unbinned fit
	real_t totalAmplitudeSquared_XY(const int& BSign, real_t& modA, real_t& modAbar, real_t& Dd, real_t& correction, MinuitParameterSet& MPS);
	// wrapper to above but taking event as input for plotting
	real_t totalAmplitudeSquared_XY(const int& BSign, const Event& event, CoherentSum& A_coherentSum, CoherentSum& Abar_coherentSum,  PhaseCorrection& phaseCorrection, MinuitParameterSet& MPS);
	real_t totalAmplitudeSquared_DPi_XY(const int& BSign, const Event& event, CoherentSum& A_coherentSum, CoherentSum& Abar_coherentSum,  PhaseCorrection& phaseCorrection, MinuitParameterSet& MPS);
	// used as the pdf for MD unbinned fit
	real_t totalAmplitudeSquared_rB(const int& BSign, real_t& modA, real_t& modAbar, real_t& Dd, real_t& correction, MinuitParameterSet& MPS);
	// wrapper to above but taking event as input for plotting
	real_t totalAmplitudeSquared_rB(const int& BSign, const Event& event, CoherentSum& A_coherentSum, CoherentSum& Abar_coherentSum, PhaseCorrection& phaseCorrection, MinuitParameterSet& MPS);
	real_t totalAmplitudeSquared_rB_Belle2010(const int& BSign,CoherentSum& A, CoherentSum& Abar,  DtoKpipiAmplitude& amp, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi, real_t& deltaC, MinuitParameterSet& MPS,const Event& event);
	real_t totalAmplitudeSquared_DPi_XY(const int &BSign, real_t &modA, real_t &modAbar, real_t &Dd, real_t &correction, MinuitParameterSet &MPS);

	// used as the pdf for generation:
	real_t totalAmplitudeSquared_rB(const int& BSign, CoherentSum& A, CoherentSum& Abar, real_t& deltaC, MinuitParameterSet& MPS, Event& event);
	// used as the normalisation for MD unbinned fit
	real_t totalAmplitudeSquared_Integrated(const int& BSign,  real_t& normA, real_t& normAbar, MinuitParameterSet& MPS, std::pair<real_t, real_t> & crossTerm);

	// used as the normalisation for MD unbinned fit
	real_t totalAmplitudeSquared_DPi_Integrated(const int& BSign,  real_t& normA, real_t& normAbar, MinuitParameterSet& MPS, std::pair<real_t, real_t> & crossTerm);

	// as above but for rB MPS, NOTE - a bit convoluted to avoid putting the integration inside the normalisation for now, can change to not use x/y as an in between once integration inside anyways for the qmi/sb method
	real_t totalAmplitudeSquared_Integrated_rB(const int& BSign, real_t& normA, real_t& normAbar, MinuitParameterSet& MPS, std::pair<real_t, real_t>& crossTerm);

	// used in the MI binned fit as the pdf
	real_t decayWidth_Integrated_xy(const int& BSign, const int& binSign, real_t& Fi, real_t& Fbari, MinuitParameterSet& MPS, real_t& ci, real_t& si);
	real_t decayWidth_Integrated_rB(const int& BSign, const int& binSign, real_t& Fi, real_t& Fbari, MinuitParameterSet& MPS, real_t& ci, real_t& si);

	// bes style generation
	// CP:
	real_t totalAmplitudeSquared_CP(const int& CPsign, CoherentSum& A, CoherentSum& Abar, PhaseCorrection& PC, Event& event);
	real_t totalAmplitudeSquared_CP(const int& CPsign, CoherentSum& A, CoherentSum& Abar, PhaseCorrection& PC, const Event& event);
	real_t totalAmplitudeSquared_CP_test(const int& CPsign, DtoKpipiAmplitude& amp, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi, PhaseCorrection& PC, Event& event);

	// flavour:
	real_t totalAmplitudeSquared_flavour(const bool& isFlavour, CoherentSum& A, CoherentSum& Abar, Event& event);
	real_t totalAmplitudeSquared_flavour(const bool& isFlavour, CoherentSum& A, CoherentSum& Abar, const Event& event);
	real_t totalAmplitudeSquared_flavour_test(const bool& isFlavour, DtoKpipiAmplitude& amp, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi, Event& event);
	// BES Kspipi
	real_t totalAmplitudeSquared_BES_KSpipi(D0ToKSpipi2018& tKLpipi, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi, PhaseCorrection& PC, Event& event_sig, Event& event_tag);
	real_t totalAmplitudeSquared_flavour_KSpipi(const bool& isFlavour, D0ToKSpipi2018& tKLpipi, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi, Event& event);
	real_t totalAmplitudeSquared_CP_KSpipi(const int& CPsign, D0ToKSpipi2018& tKLpipi, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi, PhaseCorrection& PC, Event& event);

	// mixed:
	real_t totalAmplitudeSquared_BES(CoherentSum& A, CoherentSum& Abar, PhaseCorrection& PC, Event& event_main, Event& event_tag);
	real_t totalAmplitudeSquared_BES(CoherentSum& A, CoherentSum& Abar, PhaseCorrection& PC, const Event& event_main, const Event& event_tag);
	// Since there are no double KLpipi
	real_t totalAmplitudeSquared_BES_KLpipi(D0ToKSpipi2018& tKSpipi, D0ToKLpipi2018& tKLpipi, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi, PhaseCorrection& PC, PhaseCorrection_KLpipi& PC_KLipi, Event& event_sig, Event& event_tag);
	real_t totalAmplitudeSquared_flavour_KLpipi(const bool& isFlavour, D0ToKLpipi2018& tKLpipi, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi, Event& event);
	real_t totalAmplitudeSquared_CP_KLpipi(const int& CPsign, D0ToKLpipi2018& tKLpipi, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi, PhaseCorrection_KLpipi& PC_KLipi, Event& event);
	real_t totalAmplitudeSquared_BES_test(DtoKpipiAmplitude& amp, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi, PhaseCorrection& PC, Event& event_sig, Event& event_tag);
	real_t totalAmplitudeSquared_BES_KLpipi_test(DtoKpipiAmplitude& amp, D0ToKLpipi2018& tKLpipi, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi, PhaseCorrection& PC, Event& event_sig, Event& event_tag);
	real_t totalAmplitudeSquared_BES_KSpipi_KLpipi(CoherentSum& A, CoherentSum& Abar, D0ToKLpipi2018& tKLpipi, std::vector<CompiledExpression<real_t(const real_t*, const real_t*)>> &cPhi, PhaseCorrection& PC, PhaseCorrection_KLpipi& PC_KLipi, Event& event_sig, Event& event_tag);

	// old - hopefully can delete when above working
	real_t totalAmplitudeSquared_BES(bool& isKspipi, CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, PhaseCorrection& PC, Event& eventA, Event& eventB);
	
	// pdf for unbinned QMI fitting for each of the BES tags:
	// CP state the D->Kspipi is even or odd, sign = sign of the CP of this D, so an tag of Kspi0 for event and KK for odd
	real_t totalAmplitudeSquared_CP(const int& CPsign, real_t& modA, real_t& modAbar, real_t& Dd, real_t& correction);
	// flavour = the D is pure D^0 or Dbar^0 so a tag like K-pi+ from a Dbar or K+pi- from a D
	real_t totalAmplitudeSquared_flavour(real_t& modA);
	// same = both D's go to Kspipi so both sets of events need to be considered
	real_t totalAmplitudeSquared_BES(real_t& modA, real_t& modAbar, real_t& modA_tag, real_t& modAbar_tag, real_t& deltaD, real_t& deltaD_tag, real_t& deltaC, real_t& deltaC_tag);
	real_t totalAmplitudeSquared_BES_KLpipi(real_t& modA, real_t& modAbar, real_t& modA_tag, real_t& modAbar_tag, real_t& deltaD, real_t& deltaD_tag, real_t& deltaC, real_t& deltaC_tag);
	real_t totalAmplitudeSquared_BES(const Event &event, CoherentSum &A_coherentSum, CoherentSum &Abar_coherentSum, PhaseCorrection &phaseCorrection);

}


#endif
