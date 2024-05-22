#include "AmpGen/Formalism.h"

#include "AmpGen/Event.h"
#include "AmpGen/EventList.h"
#include "AmpGen/EventListSIMD.h"
#include "AmpGen/Expression.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MyStructs.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PhaseCorrection.h"
#include "AmpGen/PhaseCorrection_KLpipi.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/Types.h"
#include "extern/DtoKpipiAmplitude.h"
#include "extern/Belle2010Amplitude.h"
#include "extern/D0ToKLpipi2018.h"
#if ENABLE_AVX
using EventList_type = AmpGen::EventListSIMD;
#else
using EventList_type = AmpGen::EventList;
#endif

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <complex>

using namespace AmpGen;
using twoDVector_t = std::vector<std::vector<real_t>>;

real_t AmpGen::logPoisson(real_t Exp, real_t Obs)
{
	return -Exp + Obs*std::log(Exp) - std::lgamma(Obs+1);
}

// log of a chi^2 likelihood with expected number of events in bin E and observed O
real_t AmpGen::Chi2(real_t Exp, real_t Obs)
{
	return  std::pow(Exp-Obs, 2) / Exp;
}
amplitudeInfo AmpGen::fillAmplitudeInfo(EventList_type &events, CoherentSum &A, CoherentSum &Abar)
{
	ProfileClock TIME_amplitudes;
	TIME_amplitudes.start();

	amplitudeInfo amplitudes{};
	size_t nEvents{events.size()};
	std::vector<real_t> Dd(nEvents), absA(nEvents), absAbar(nEvents);

	auto evalA = A.amplitudeEvaluator(&events);
	auto evalAbar = Abar.amplitudeEvaluator(&events);

#pragma omp parallel for
	for (unsigned i = 0; i < nEvents; ++i)
	{
		complex_t thisA = evalA(events[i]);
		complex_t thisAbar = evalAbar(events[i]);
		absA[i] = std::abs(thisA);
		absAbar[i] = std::abs(thisAbar);
		// Dd[i] = std::arg( thisA * std::conj(thisAbar) );
		Dd[i] = Deltadelta(events[i], A, Abar); // changed to make there ony one delta_D definition 19/9/23
		//increasing the precision of the phase calculation
/*		if (i==1) {

			const Event event = events[i];
			event.print();
			INFO("Event " << i << "This A = " << thisA << "This Abar = " << thisAbar);
			INFO("Event " << i << " of " << nEvents << " has Dd = " << Dd[i] << " and |A| = " << absA[i] << " and |Abar| = " << absAbar[i]);
		}*/
	}
	amplitudes.A = absA;
	amplitudes.Abar = absAbar;
	amplitudes.deltaD = Dd;

	TIME_amplitudes.stop();
	INFO("Took " << TIME_amplitudes << "ms to collate amplitude information");

	return amplitudes;
}
double* AmpGen::get_bkg_fraction(bkgfraction & bkg_frac)
{
		double *bkg_frac_val = new double[7];
		INFO("bkg fraction size: " << bkg_frac.tag.size());
		for(size_t i{0}; i < bkg_frac.tag.size(); i++){
		if(bkg_frac.tag[i]=="kspipi"){
			bkg_frac_val[4] = bkg_frac.frac[i];
		}
		else if(bkg_frac.tag[i]=="CPeven")
		{
			bkg_frac_val[3] = bkg_frac.frac[i];
		}
		else if(bkg_frac.tag[i]=="CPodd")
		{
			bkg_frac_val[2] = bkg_frac.frac[i];
		}
		else if(bkg_frac.tag[i]=="flavourBar")
		{
			bkg_frac_val[1] = bkg_frac.frac[i];
		}
		else if(bkg_frac.tag[i]=="flavour")
		{
			bkg_frac_val[0] = bkg_frac.frac[i];
		}
		else if(bkg_frac.tag[i]=="Bplus")
		{
			bkg_frac_val[5] = bkg_frac.frac[i];
		}
		else if(bkg_frac.tag[i]=="Bminus")
		{
			bkg_frac_val[6] = bkg_frac.frac[i];
		}
		else
		{
			INFO("bkg fraction tag not found");
			break;
		}
	}
	return bkg_frac_val;

}
double* AmpGen::get_norm_eff(normeff & norm_eff)
{
		double *norm_eff_val = new double[7];
		for(size_t i{0}; i < norm_eff.tag.size(); i++){
		if(norm_eff.tag[i]=="kspipi"){
			norm_eff_val[4] = norm_eff.frac[i];
		}
		else if(norm_eff.tag[i]=="CPeven")
		{
			norm_eff_val[3] = norm_eff.frac[i];
		}
		else if(norm_eff.tag[i]=="CPodd")
		{
			norm_eff_val[2] = norm_eff.frac[i];
		}
		else if(norm_eff.tag[i]=="flavourBar")
		{
			norm_eff_val[1] = norm_eff.frac[i];
		}
		else if(norm_eff.tag[i]=="flavour")
		{
			norm_eff_val[0] = norm_eff.frac[i];
		}
		else if(norm_eff.tag[i]=="Bplus")
		{
			norm_eff_val[5] = norm_eff.frac[i];
		}
		else if(norm_eff.tag[i]=="Bminus")
		{
			norm_eff_val[6] = norm_eff.frac[i];
		}
		else
		{
			INFO("bkg fraction tag not found");
			break;
		}
	}
	return norm_eff_val;

}
void AmpGen::getNorm(EventList_type &events, D0ToKLpipi2018 &tKLpipi, std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> &cPhi, real_t &normA, real_t &normAbar)
{
	ProfileClock TIME_amplitudes_norm;
	TIME_amplitudes_norm.start();

	size_t nEvents{events.size()};
	complex<double> A;
	complex<double> Abar;
	INFO("Size: " << nEvents);
#pragma omp parallel for
	for (unsigned i = 0; i < nEvents; ++i)
	{
		vector<double> p1, p2, p3;
		for (int j = 0; j < 4; j++)
		{
			p1.push_back(cPhi[4 + j](events[i].address()));
			p2.push_back(cPhi[8 + j](events[i].address()));
			p3.push_back(cPhi[12 + j](events[i].address()));
		}
		// Sig Mention that all besiii model is with Ks pi+ pi- order
		complex<double> thisA = tKLpipi.Amp_PFT(p1, p3, p2);
		vector<double> reverse_p1;
		vector<double> reverse_p2;
		vector<double> reverse_p3;
		reverse_p1.clear();
		reverse_p2.clear();
		reverse_p3.clear();
		reverse_p1.push_back(-p1[0]);
		reverse_p1.push_back(-p1[1]);
		reverse_p1.push_back(-p1[2]);
		reverse_p1.push_back(p1[3]);
		reverse_p2.push_back(-p2[0]);
		reverse_p2.push_back(-p2[1]);
		reverse_p2.push_back(-p2[2]);
		reverse_p2.push_back(p2[3]);
		reverse_p3.push_back(-p3[0]);
		reverse_p3.push_back(-p3[1]);
		reverse_p3.push_back(-p3[2]);
		reverse_p3.push_back(p3[3]);
		complex<double> thisAbar = tKLpipi.Amp_PFT(reverse_p1, reverse_p2, reverse_p3);

		A += thisA;
		Abar += thisAbar;
	}
	normA = std::norm(A) / nEvents;
	normAbar = std::norm(Abar) / nEvents;

	TIME_amplitudes_norm.stop();
	INFO("Took " << TIME_amplitudes_norm << "ms to norm KLpipi amplitude information");

	return;
}
amplitudeInfo AmpGen::fillAmplitudeInfo_KLpipi(EventList_type &events, D0ToKLpipi2018 &tKLpipi, std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> &cPhi)
{
	ProfileClock TIME_amplitudes;
	TIME_amplitudes.start();

	amplitudeInfo amplitudes{};
	size_t nEvents{events.size()};
	std::vector<real_t> Dd(nEvents), absA(nEvents), absAbar(nEvents);

#pragma omp parallel for
	for (unsigned i = 0; i < nEvents; ++i)
	{
		vector<double> p1, p2, p3;
		for (int j = 0; j < 4; j++)
		{
			p1.push_back(cPhi[4 + j](events[i].address()));
			p2.push_back(cPhi[8 + j](events[i].address()));
			p3.push_back(cPhi[12 + j](events[i].address()));
		}
		// Sig Mention that all besiii model is with Ks pi+ pi- order
		complex<double> thisA = tKLpipi.Amp_PFT(p1, p2, p3);
		vector<double> reverse_p1;
		vector<double> reverse_p2;
		vector<double> reverse_p3;
		reverse_p1.clear();
		reverse_p2.clear();
		reverse_p3.clear();
		reverse_p1.push_back(-p1[0]);
		reverse_p1.push_back(-p1[1]);
		reverse_p1.push_back(-p1[2]);
		reverse_p1.push_back(p1[3]);
		reverse_p2.push_back(-p2[0]);
		reverse_p2.push_back(-p2[1]);
		reverse_p2.push_back(-p2[2]);
		reverse_p2.push_back(p2[3]);
		reverse_p3.push_back(-p3[0]);
		reverse_p3.push_back(-p3[1]);
		reverse_p3.push_back(-p3[2]);
		reverse_p3.push_back(p3[3]);
		complex<double> thisAbar = tKLpipi.Amp_PFT(reverse_p1, reverse_p3, reverse_p2);

		absA[i] = std::abs(thisA);
		absAbar[i] = std::abs(thisAbar);
		Dd[i] = Deltadelta_KL(thisA,thisAbar); // changed to make there only one delta_D definition 19/9/23
	}
	amplitudes.A = absA;
	amplitudes.Abar = absAbar;
	amplitudes.deltaD = Dd;

	TIME_amplitudes.stop();
	INFO("Took " << TIME_amplitudes << "ms to collate amplitude information");

	return amplitudes;
}

// return the phase difference between A and Abar -> Kspipi  for a single event in radians
real_t AmpGen::Deltadelta(const Event &event, CoherentSum &A, CoherentSum &Abar)
{

	//convert CP convention
	/*double temp = std::arg(A.getValNoCache(event) * std::conj(Abar.getValNoCache(event))) - M_PI;
	while (temp < -M_PI)
	{
		temp += 2.0 * M_PI;
	}
	while (temp > M_PI)
	{
		temp -= 2.0 * M_PI;
	}*/
	double temp = std::arg(A.getValNoCache(event) * std::conj(Abar.getValNoCache(event)));
	return temp;
}

// return the phase difference between A and Abar -> Kspipi  for a single event in radians
real_t AmpGen::Deltadelta_KL(complex<double> &A, complex<double> &Abar)
{

	//convert CP convention
	double temp = std::arg(A * std::conj(Abar)) - M_PI;
	
	while (temp < -M_PI)
	{
		temp += 2.0 * M_PI;
	}
	while (temp > M_PI)
	{
		temp -= 2.0 * M_PI;
	}
	return temp;
}

//***************************** INTEGRATIONS *****************
// INTEGRATE OVER EVERYONE
// return the cross term for normalisation of total amplitude squared -> this is the part that is integrated over/integration is done here
std::pair<real_t, real_t> AmpGen::totalAmplitudeSquared_Integrated_crossTerm(amplitudeInfo &amplitudes, EventList_type &events, PhaseCorrection &PC)
{

	real_t realAAbarStar{0}, imagAAbarStar{0};
	real_t phaseArgument{0};
	real_t weight_sum{0};
	// #pragma omp parallel for reduction(+:realAAbarStar, imagAAbarStar)
	// for (auto event : events){ // not sure how to use this syntax with the omp
	for (int i = 0; i < amplitudes.A.size(); i++)
	{
		phaseArgument = amplitudes.deltaD[i]; // + PC.eval(events[i]);
		realAAbarStar += events[i].weight_eff()* amplitudes.A[i] * amplitudes.Abar[i] * cos(phaseArgument);
		imagAAbarStar += events[i].weight_eff()* amplitudes.A[i] * amplitudes.Abar[i] * sin(phaseArgument);
		weight_sum += events[i].weight_eff();
	}
	return {realAAbarStar / weight_sum, imagAAbarStar / weight_sum};

}


std::pair<real_t, real_t> AmpGen::totalAmplitudeSquared_Integrated_crossTerm_noeff(amplitudeInfo &amplitudes, EventList_type &events, PhaseCorrection &PC)
{

	real_t realAAbarStar{0}, imagAAbarStar{0};
	real_t phaseArgument{0};
	real_t weight_sum{0};
	// #pragma omp parallel for reduction(+:realAAbarStar, imagAAbarStar)
	// for (auto event : events){ // not sure how to use this syntax with the omp
	for (int i = 0; i < amplitudes.A.size(); i++)
	{
		phaseArgument = amplitudes.deltaD[i]; // + PC.eval(events[i]);
		realAAbarStar += amplitudes.A[i] * amplitudes.Abar[i] * cos(phaseArgument);
		imagAAbarStar += amplitudes.A[i] * amplitudes.Abar[i] * sin(phaseArgument);
	}
	return {realAAbarStar / amplitudes.A.size(), imagAAbarStar / amplitudes.A.size()};

}

// INTEGRATE OVER EVERYONE AND BINS
// take a vector for ci and si and populate it as well as the overall normalisation terms
void AmpGen::doBinnedIntegration(std::vector<real_t> &ci_posBins, std::vector<real_t> &si_posBins, std::vector<real_t> &ci_negBins, std::vector<real_t> &si_negBins,
								 amplitudeInfo &amplitudes, EventList_type &events, PhaseCorrection &PC)
{
	for (int i{0}; i < events.size(); i++)
	{
		int currentBin{amplitudes.bins[i]};

		real_t totalPhase = amplitudes.deltaD[i] + PC.eval(events[i]);
		real_t totalA = amplitudes.A[i] * amplitudes.Abar[i];

		if (currentBin > 0)
		{
			ci_posBins[currentBin] += totalA * cos(totalPhase);
			si_posBins[currentBin] += totalA * sin(totalPhase);
		}
		else
		{
			int unsignedBin{std::abs(currentBin)};
			ci_negBins[unsignedBin] += totalA * cos(totalPhase);
			si_negBins[unsignedBin] += totalA * sin(totalPhase);
		}
	}

	return;
}

void AmpGen::doUnbinnedIntegration(amplitudeInfo &amplitudes, EventList_type &events, PhaseCorrection &PC, real_t &cosTerm, real_t &sinTerm, real_t &complicatedTerm)
{
	int nEvents = events.size();
	for (int i{0}; i < nEvents; i++)
	{
		int tagged_i{(i + nEvents / 2) % nEvents};

		real_t totalPhase = amplitudes.deltaD[i] + PC.eval(events[i]);
		real_t totalA = amplitudes.A[i] * amplitudes.Abar[i];

		real_t tagged_totalPhase = amplitudes.deltaD[tagged_i] + PC.eval(events[tagged_i]);
		real_t tagged_totalA = amplitudes.A[tagged_i] * amplitudes.Abar[tagged_i];

		cosTerm += totalA * cos(totalPhase);
		sinTerm += totalA * sin(totalPhase);

		complicatedTerm += std::pow(amplitudes.A[i] * amplitudes.Abar[tagged_i], 2) + std::pow(amplitudes.Abar[i] * amplitudes.A[tagged_i], 2) - 2 * totalA * tagged_totalA * cos(totalPhase - tagged_totalPhase);
	}

	return;
}

void AmpGen::doUnbinnedIntegration_temp(amplitudeInfo &amplitudes, EventList_type &events, PhaseCorrection &PC, real_t &normA, real_t &normAbar,real_t &cosTerm, real_t &sinTerm, real_t &complicatedTerm)
{
	int nEvents = events.size();
	for (int i{0}; i < nEvents; i++)
	{
		int tagged_i{(i + nEvents / 2) % nEvents};

		real_t totalPhase = amplitudes.deltaD[i] + PC.eval(events[i]);
		real_t totalA = amplitudes.A[i] * amplitudes.Abar[i];

		real_t tagged_totalPhase = amplitudes.deltaD[tagged_i] + PC.eval(events[tagged_i]);
		real_t tagged_totalA = amplitudes.A[tagged_i] * amplitudes.Abar[tagged_i];

		normA +=  std::pow(amplitudes.A[i], 2);
		normAbar +=  std::pow(amplitudes.Abar[i], 2);
		cosTerm +=  totalA * cos(totalPhase);
		sinTerm +=  totalA * sin(totalPhase);

		complicatedTerm += (std::pow(amplitudes.A[i] * amplitudes.Abar[tagged_i], 2) + std::pow(amplitudes.Abar[i] * amplitudes.A[tagged_i], 2) - 2 * totalA * tagged_totalA * cos(totalPhase - tagged_totalPhase));
	}

	return;
}

void AmpGen::doUnbinnedIntegration_KLpipi(amplitudeInfo &amplitudes, amplitudeInfo &amplitudes_klpipi, EventList_type &events, EventList_type &events_klpipi, PhaseCorrection &PC, PhaseCorrection_KLpipi &PC_KLipi, real_t &KLpipi_normA, real_t &KLpipi_normAbar, real_t &cosTerm, real_t &sinTerm, real_t &complicatedTerm)
{
	int nEvents = events.size();
	for (int i{0}; i < nEvents; i++)
	{
		int tagged_i{(i + nEvents / 2) % nEvents};

		real_t totalPhase = amplitudes_klpipi.deltaD[i] + PC_KLipi.eval(events_klpipi[i]);// do not need adding pi since the amplitude is already in the correct convention
		real_t totalA = amplitudes_klpipi.A[i] * amplitudes_klpipi.Abar[i];

		real_t tagged_totalPhase = amplitudes.deltaD[i] + PC.eval(events[i]);
		real_t tagged_totalA = amplitudes.A[i] * amplitudes.Abar[i];

		KLpipi_normA += std::pow(amplitudes_klpipi.A[i], 2);
		KLpipi_normAbar += std::pow(amplitudes_klpipi.Abar[i], 2);
		cosTerm += totalA * cos(totalPhase);
		sinTerm += totalA * sin(totalPhase);

		complicatedTerm += std::pow(amplitudes.A[i] * amplitudes_klpipi.Abar[i], 2) + std::pow(amplitudes.Abar[i] * amplitudes_klpipi.A[i], 2) - 2 * totalA * tagged_totalA * cos(totalPhase - tagged_totalPhase);
	}

	return;
}

// ***********************************************************************************************************************************
// ************************** DIFFERENT VARIATIONS ON THE OVERALL PDF FOR FITTING AND GENERATING *************************************
// ***********************************************************************************************************************************

// |Atot|^2 for one event, xy formalism
// used in MD unbinned fitting
real_t AmpGen::totalAmplitudeSquared_XY(const int &BSign, real_t &modA, real_t &modAbar, real_t &Dd, real_t &correction, MinuitParameterSet &MPS)
{
	real_t phase{Dd + correction};

	if (BSign == 1)
	{
		real_t xPlus{MPS["xPlus"]->mean()};
		real_t yPlus{MPS["yPlus"]->mean()};
		real_t rB2{std::pow(xPlus, 2) + std::pow(yPlus, 2)};

		return std::pow(modA, 2) * rB2 + std::pow(modAbar, 2) + 2 * modA * modAbar * (xPlus * cos(phase) - yPlus * sin(phase));
	}
	else
	{
		real_t xMinus{MPS["xMinus"]->mean()};
		real_t yMinus{MPS["yMinus"]->mean()};
		real_t rB2{std::pow(xMinus, 2) + std::pow(yMinus, 2)};

		return std::pow(modA, 2) + std::pow(modAbar, 2) * rB2 + 2 * modA * modAbar * (xMinus * cos(phase) + yMinus * sin(phase));
	}
}

real_t AmpGen::totalAmplitudeSquared_DPi_XY(const int &BSign, real_t &modA, real_t &modAbar, real_t &Dd, real_t &correction, MinuitParameterSet &MPS)
{
	real_t phase{Dd + correction};
	real_t xXi{MPS["xXi"]->mean()};
	real_t yXi{MPS["yXi"]->mean()};

	if (BSign == 1)
	{
		real_t xPlus{MPS["xPlus"]->mean()};
		real_t yPlus{MPS["yPlus"]->mean()};
		real_t xPlus_DPi = xPlus * xXi - yPlus * yXi;
		real_t yPlus_DPi = yPlus * xXi + xPlus * yXi;
		
		real_t rB2{std::pow(xPlus_DPi, 2) + std::pow(yPlus_DPi, 2)};

		return std::pow(modA, 2) * rB2 + std::pow(modAbar, 2) + 2 * modA * modAbar * (xPlus_DPi * cos(phase) - yPlus_DPi * sin(phase));
	}
	else
	{
		real_t xMinus{MPS["xMinus"]->mean()};
		real_t yMinus{MPS["yMinus"]->mean()};
		real_t xMinus_DPi = xMinus * xXi - yMinus * yXi;
		real_t yMinus_DPi = yMinus * xXi + xMinus * yXi;

		real_t rB2{std::pow(xMinus_DPi, 2) + std::pow(yMinus_DPi, 2)};

		return std::pow(modA, 2) + std::pow(modAbar, 2) * rB2 + 2 * modA * modAbar * ( xMinus_DPi * cos(phase) + yMinus_DPi * sin(phase));
	}
}

// as above but taking event as argument for plotting reasons
real_t AmpGen::totalAmplitudeSquared_XY(const int &BSign, const Event &event, CoherentSum &A_coherentSum, CoherentSum &Abar_coherentSum, PhaseCorrection &phaseCorrection, MinuitParameterSet &MPS)
{ // set up variables
	complex_t A{A_coherentSum.getValNoCache(event)};
	complex_t Abar{Abar_coherentSum.getValNoCache(event)};

	real_t modA{std::abs(A)}, modAbar{std::abs(Abar)};
	real_t Dd{Deltadelta(event, A_coherentSum, Abar_coherentSum)}; // changed to make there ony one delta_D definition 19/9/23
	real_t correction{phaseCorrection.eval(event)};

	return totalAmplitudeSquared_XY(BSign, modA, modAbar, Dd, correction, MPS);
}
real_t AmpGen::totalAmplitudeSquared_DPi_XY(const int &BSign, const Event &event, CoherentSum &A_coherentSum, CoherentSum &Abar_coherentSum, PhaseCorrection &phaseCorrection, MinuitParameterSet &MPS)
{ // set up variables
	complex_t A{A_coherentSum.getValNoCache(event)};
	complex_t Abar{Abar_coherentSum.getValNoCache(event)};

	real_t modA{std::abs(A)}, modAbar{std::abs(Abar)};
	real_t Dd{Deltadelta(event, A_coherentSum, Abar_coherentSum)}; // changed to make there ony one delta_D definition 19/9/23
	real_t correction{phaseCorrection.eval(event)};

	return totalAmplitudeSquared_DPi_XY(BSign, modA, modAbar, Dd, correction, MPS);
}
// as 2 above but fitting to rB dB and gamma
real_t AmpGen::totalAmplitudeSquared_rB(const int &BSign, real_t &modA, real_t &modAbar, real_t &Dd, real_t &correction, MinuitParameterSet &MPS)
{
	real_t phase{Dd + correction};
	real_t rB{MPS["rB"]->mean()};
	real_t dB{MPS["dB"]->mean() * M_PI / 180};
	real_t gamma{MPS["gamma"]->mean() * M_PI / 180};

	if (BSign == 1)
	{
		return std::pow(modA * rB, 2) + std::pow(modAbar, 2) + 2 * modA * modAbar * rB * cos(phase + (dB + gamma));
	}
	else
	{
		return std::pow(modA, 2) + std::pow(modAbar * rB, 2) + 2 * modA * modAbar * rB * cos(phase - (dB - gamma));
	}
}
// as above but taking event as argument for plotting reasons
real_t AmpGen::totalAmplitudeSquared_rB(const int &BSign, const Event &event, CoherentSum &A_coherentSum, CoherentSum &Abar_coherentSum, PhaseCorrection &phaseCorrection, MinuitParameterSet &MPS)
{ // set up variables
	complex_t A{A_coherentSum.getValNoCache(event)};
	complex_t Abar{Abar_coherentSum.getValNoCache(event)};

	real_t modA{std::abs(A)}, modAbar{std::abs(Abar)};
	real_t Dd{Deltadelta(event, A_coherentSum, Abar_coherentSum)}; // changed to make there ony one delta_D definition 19/9/23
	real_t correction{phaseCorrection.eval(event)};

	return totalAmplitudeSquared_rB(BSign, modA, modAbar, Dd, correction, MPS);
}

//Temp still need debug
real_t AmpGen::totalAmplitudeSquared_rB_Belle2010(const int &BSign, CoherentSum &A, CoherentSum &Abar, DtoKpipiAmplitude &amp, std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> &cPhi, real_t &deltaC, MinuitParameterSet &MPS, const Event &event)
{
	dcomplex CA = amp.get_amplitude(cPhi[1](event.address()), cPhi[0](event.address()));
	dcomplex CAbar = amp.get_amplitude(cPhi[0](event.address()), cPhi[1](event.address()));

	real_t modA{std::abs(A.getValNoCache(event))};
	real_t modAbar{std::abs(Abar.getValNoCache(event))};

	real_t rB{MPS["rB"]->mean()};
	real_t deltaB{MPS["dB"]->mean() * M_PI / 180};
	real_t gamma{MPS["gamma"]->mean() * M_PI / 180};

	real_t phaseDiff{-arg(CAbar * conj(CA))}; // will be in radians already
											  //	real_t phaseDiff{Deltadelta(event, A, Abar)};  // will be in radians already

	if (BSign == 1)
	{
		return std::pow(modA * rB, 2) + std::pow(modAbar, 2) + 2 * rB * modA * modAbar * cos(phaseDiff + deltaC + (deltaB + gamma));
		// 2*rB*std::real(A.getVal()* std::conjugate(Abar.getVal())) = above with deltaB=gamma=zero <- from jonas in meeting
	}
	else
	{
		return std::pow(modA, 2) + std::pow(modAbar * rB, 2) + 2 * rB * modA * modAbar * cos(phaseDiff + deltaC - (deltaB - gamma));
	}
}

// |Atot|^2 for one event, rB formalism -> as above but different input types
// used in generating
real_t AmpGen::totalAmplitudeSquared_rB(const int &BSign, CoherentSum &A, CoherentSum &Abar, real_t &deltaC, MinuitParameterSet &MPS, Event &event)
{
	real_t modA{std::abs(A.getValNoCache(event))};
	real_t modAbar{std::abs(Abar.getValNoCache(event))};

	real_t rB{MPS["rB"]->mean()};
	real_t deltaB{MPS["dB"]->mean() * M_PI / 180};
	real_t gamma{MPS["gamma"]->mean() * M_PI / 180};

	real_t phaseDiff{Deltadelta(event, A, Abar)}; // will be in radians already

	if (BSign == 1)
	{
		return std::pow(modA * rB, 2) + std::pow(modAbar, 2) + 2 * rB * modA * modAbar * cos(phaseDiff + deltaC + (deltaB + gamma));
		// 2*rB*std::real(A.getVal()* std::conjugate(Abar.getVal())) = above with deltaB=gamma=zero <- from jonas in meeting
	}
	else
	{
		return std::pow(modA, 2) + std::pow(modAbar * rB, 2) + 2 * rB * modA * modAbar * cos(phaseDiff + deltaC - (deltaB - gamma));
	}
}

// |Atot|^2 integrated for all events
// used in MD fitting normalisation - xy formalism
real_t AmpGen::totalAmplitudeSquared_DPi_Integrated(const int &BSign, real_t &normA, real_t &normAbar, MinuitParameterSet &MPS, std::pair<real_t, real_t> &crossTerm)
{
	real_t xXi{MPS["xXi"]->mean()};
	real_t yXi{MPS["yXi"]->mean()};

	if (BSign == 1)
	{
		real_t xPlus{MPS["xPlus"]->mean()};
		real_t yPlus{MPS["yPlus"]->mean()};
		real_t xPlus_DPi = xPlus * xXi - yPlus * yXi;
		real_t yPlus_DPi = yPlus * xXi + xPlus * yXi;

		real_t rB2{std::pow(xPlus_DPi, 2) + std::pow(yPlus_DPi, 2)};
		return normA * rB2 + normAbar + 2 * (xPlus_DPi * crossTerm.first - yPlus_DPi * crossTerm.second);
	}
	else
	{
		real_t xMinus{MPS["xMinus"]->mean()};
		real_t yMinus{MPS["yMinus"]->mean()};
		real_t xMinus_DPi = xMinus * xXi - yMinus * yXi;
		real_t yMinus_DPi = yMinus * xXi + xMinus * yXi;

		real_t rB2{std::pow(xMinus_DPi, 2) + std::pow(yMinus_DPi, 2)};

		return normA + normAbar * rB2 + 2 * (xMinus_DPi * crossTerm.first + yMinus_DPi * crossTerm.second);
	}
}
real_t AmpGen::totalAmplitudeSquared_Integrated(const int &BSign, real_t &normA, real_t &normAbar, MinuitParameterSet &MPS, std::pair<real_t, real_t> &crossTerm)
{
	
	if (BSign == 1)
	{
		real_t xPlus{MPS["xPlus"]->mean()};
		real_t yPlus{MPS["yPlus"]->mean()};
		real_t rB2{std::pow(xPlus, 2) + std::pow(yPlus, 2)};
		return normA * rB2 + normAbar + 2 * ((xPlus) * crossTerm.first - (yPlus) * crossTerm.second);
	}
	else
	{
		real_t xMinus{MPS["xMinus"]->mean()};
		real_t yMinus{MPS["yMinus"]->mean()};
		real_t rB2{std::pow(xMinus, 2) + std::pow(yMinus, 2)};

		return normA + normAbar * rB2 + 2 * ((xMinus) * crossTerm.first + (yMinus) * crossTerm.second);
	}
}
// |Atot|^2 integrated for all events
// used in MD fitting normalisation - rB formalism, bad way to do it to avoid putting integral in LL 14/6 -> CHANGE WHEN INTEGRAL INSIDE ANYWAYS
real_t AmpGen::totalAmplitudeSquared_Integrated_rB(const int &BSign, real_t &normA, real_t &normAbar, MinuitParameterSet &MPS, std::pair<real_t, real_t> &crossTerm)
{
	real_t rB{MPS["rB"]->mean()};
	real_t dB{MPS["dB"]->mean() * M_PI / 180};
	real_t gamma{MPS["gamma"]->mean() * M_PI / 180};

	if (BSign == 1)
	{
		real_t xPlus = rB * cos(dB + gamma);
		real_t yPlus = rB * sin(dB + gamma);

		return normA * std::pow(rB, 2) + normAbar + 2 * (xPlus * crossTerm.first - yPlus * crossTerm.second);
	}
	else
	{
		real_t xMinus = rB * cos(dB - gamma);
		real_t yMinus = rB * sin(dB - gamma);
	}
}

// |Atot|^2 integrated over a BIN - xy formalism
// used in MI binned fitting
real_t AmpGen::decayWidth_Integrated_xy(const int &BSign, const int &binSign, real_t &Fi, real_t &Fbari, MinuitParameterSet &MPS, real_t &ci, real_t &si)
{
	if (BSign == 1)
	{
		real_t xPlus{MPS["xPlus"]->mean()};
		real_t yPlus{MPS["yPlus"]->mean()};
		real_t rB2{std::pow(xPlus, 2) + std::pow(yPlus, 2)};

		return Fi * rB2 + Fbari + 2 * std::sqrt(Fi * Fbari) * (xPlus * ci - yPlus * si * binSign);
	}
	else
	{
		real_t xMinus{MPS["xMinus"]->mean()};
		real_t yMinus{MPS["yMinus"]->mean()};
		real_t rB2{std::pow(xMinus, 2) + std::pow(yMinus, 2)};

		return Fi + Fbari * rB2 + 2 * std::sqrt(Fi * Fbari) * (xMinus * ci + yMinus * si * binSign);
	}
}

// |Atot|^2 integrated over a BIN - rB formalism
// used in MI binned fitting
real_t AmpGen::decayWidth_Integrated_rB(const int &BSign, const int &binSign, real_t &Fi, real_t &Fbari, MinuitParameterSet &MPS, real_t &ci, real_t &si)
{
	real_t rB{MPS["rB"]->mean()};
	real_t dB{MPS["dB"]->mean() * M_PI / 180};
	real_t gamma{MPS["gamma"]->mean() * M_PI / 180};

	if (BSign == 1)
	{
		real_t xPlus{rB * cos(dB + gamma)};
		real_t yPlus{rB * sin(dB + gamma)};

		return Fi * std::pow(rB, 2) + Fbari + 2 * std::sqrt(Fi * Fbari) * (xPlus * ci - yPlus * si * binSign);
	}
	else
	{
		real_t xMinus{rB * cos(dB - gamma)};
		real_t yMinus{rB * sin(dB - gamma)};

		return Fi + Fbari * std::pow(rB, 2) + 2 * std::sqrt(Fi * Fbari) * (xMinus * ci + yMinus * si * binSign);
	}
}

// ***********************************************************************************************************************************
// ************************** BES STYLE FORMULATION THINGS FOR FITTING AND GENERATING ************************************************
// ***********************************************************************************************************************************

//************ GENERATING DATA TYPES:
// CP state
real_t AmpGen::totalAmplitudeSquared_CP(const int &CPsign, CoherentSum &A, CoherentSum &Abar, PhaseCorrection &PC, Event &event)
{
	real_t modA{std::abs(A.getValNoCache(event))};
	real_t modAbar{std::abs(Abar.getValNoCache(event))};
	real_t deltaD{Deltadelta(event, A, Abar)};
	real_t deltaC{PC.evalBias(event)}; // note will need to change so can implement a bias to ues evalBais also

	return std::pow(modA, 2) + std::pow(modAbar, 2) + 2 * CPsign * modA * modAbar * cos(deltaD + deltaC);
}

// flavour state - isFlavour = true if the SIGNAL is a D, false if SIGNAL is Dbar
real_t AmpGen::totalAmplitudeSquared_flavour(const bool &isFlavour, CoherentSum &A, CoherentSum &Abar, Event &event)
{
	if (isFlavour)
	{
		return std::norm(Abar.getValNoCache(event));
	}
	else
	{
		return std::norm(A.getValNoCache(event));
	}
}

// most general version, mainly for the Kspipi tag but could be attached anywhere
real_t AmpGen::totalAmplitudeSquared_BES(CoherentSum &A, CoherentSum &Abar, PhaseCorrection &PC, Event &event_main, Event &event_tag)
{
	// |A| for each of the amplitude and amplitude bar for the D that goes to Kspipi (A) and the D to the tag (B)
	real_t modA_main{std::abs(A.getValNoCache(event_main))};
	real_t modAbar_main{std::abs(Abar.getValNoCache(event_main))};
	real_t modA_tag{std::abs(A.getValNoCache(event_tag))};
	real_t modAbar_tag{std::abs(Abar.getValNoCache(event_tag))};

	// phase things
	real_t deltaC_main{PC.evalBias(event_main)};
	real_t deltaC_tag{PC.evalBias(event_tag)};
	real_t deltaD_main{Deltadelta(event_main, A, Abar)};
	real_t deltaD_tag{Deltadelta(event_tag, A, Abar)};

	return std::pow(modA_main * modAbar_tag, 2) + std::pow(modAbar_main * modA_tag, 2) - 2 * modA_main * modAbar_main * modA_tag * modAbar_tag * cos(deltaD_main + deltaC_main - (deltaD_tag + deltaC_tag));
}
//************ GENERATING DATA TYPES (KLpipi):
// CP state
real_t AmpGen::totalAmplitudeSquared_CP_KLpipi(const int &CPsign, D0ToKLpipi2018 &tKLpipi, std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> &cPhi, PhaseCorrection_KLpipi &PC_KLipi, Event &event)
{
	vector<double> p1, p2, p3;
	for (int i = 0; i < 4; i++)
	{
		p1.push_back(cPhi[4 + i](event.address()));
		p2.push_back(cPhi[8 + i](event.address()));
		p3.push_back(cPhi[12 + i](event.address()));
	}
	complex<double> A = tKLpipi.Amp_PFT(p1, p2, p3);
	vector<double> reverse_p1;
	vector<double> reverse_p2;
	vector<double> reverse_p3;
	reverse_p1.clear();
	reverse_p2.clear();
	reverse_p3.clear();
	reverse_p1.push_back(-p1[0]);
	reverse_p1.push_back(-p1[1]);
	reverse_p1.push_back(-p1[2]);
	reverse_p1.push_back(p1[3]);
	reverse_p2.push_back(-p2[0]);
	reverse_p2.push_back(-p2[1]);
	reverse_p2.push_back(-p2[2]);
	reverse_p2.push_back(p2[3]);
	reverse_p3.push_back(-p3[0]);
	reverse_p3.push_back(-p3[1]);
	reverse_p3.push_back(-p3[2]);
	reverse_p3.push_back(p3[3]);
	complex<double> Abar = (tKLpipi.Amp_PFT(reverse_p1, reverse_p3, reverse_p2));

	real_t modA{std::abs(A)};
	real_t modAbar{std::abs(Abar)};
	real_t deltaD{Deltadelta_KL(A,Abar)};
	real_t deltaC{PC_KLipi.evalBias(event)}; // note will need to change so can implement a bias to ues evalBais also

	return std::pow(modA, 2) + std::pow(modAbar, 2) + 2 * CPsign * modA * modAbar * cos(deltaD + deltaC);
}

// flavour state - isFlavour = true if the SIGNAL is a D, false if SIGNAL is Dbar
real_t AmpGen::totalAmplitudeSquared_flavour_KLpipi(const bool &isFlavour, D0ToKLpipi2018 &tKLpipi, std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> &cPhi, Event &event)
{
	vector<double> p1, p2, p3;
	for (int i = 0; i < 4; i++)
	{
		p1.push_back(cPhi[4 + i](event.address()));
		p2.push_back(cPhi[8 + i](event.address()));
		p3.push_back(cPhi[12 + i](event.address()));
	}
	complex<double> A = tKLpipi.Amp_PFT(p1, p2, p3);
	vector<double> reverse_p1;
	vector<double> reverse_p2;
	vector<double> reverse_p3;
	reverse_p1.clear();
	reverse_p2.clear();
	reverse_p3.clear();
	reverse_p1.push_back(-p1[0]);
	reverse_p1.push_back(-p1[1]);
	reverse_p1.push_back(-p1[2]);
	reverse_p1.push_back(p1[3]);
	reverse_p2.push_back(-p2[0]);
	reverse_p2.push_back(-p2[1]);
	reverse_p2.push_back(-p2[2]);
	reverse_p2.push_back(p2[3]);
	reverse_p3.push_back(-p3[0]);
	reverse_p3.push_back(-p3[1]);
	reverse_p3.push_back(-p3[2]);
	reverse_p3.push_back(p3[3]);
	complex<double> Abar = conj(tKLpipi.Amp_PFT(reverse_p1, reverse_p3, reverse_p2));

	if (isFlavour)
	{
		return std::norm(A);
	}
	else
	{
		return std::norm(Abar);
	}
}

// most general version, mainly for the KLpipi tag but could be attached anywhere
real_t AmpGen::totalAmplitudeSquared_BES_KSpipi(D0ToKSpipi2018 &tKLpipi, std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> &cPhi, PhaseCorrection &PC, Event &event_sig, Event &event_tag)
{

	vector<double> p1, p2, p3;
	vector<double> p4, p5, p6;
	for (int i = 0; i < 4; i++)
	{
		p1.push_back(cPhi[4 + i](event_sig.address()));
		p2.push_back(cPhi[8 + i](event_sig.address()));
		p3.push_back(cPhi[12 + i](event_sig.address()));
		p4.push_back(cPhi[4 + i](event_tag.address()));
		p5.push_back(cPhi[8 + i](event_tag.address()));
		p6.push_back(cPhi[12 + i](event_tag.address()));
	}
	// Sig Mention that all besiii model is with Ks pi+ pi- order
	complex<double> A_sig = tKLpipi.Amp_PFT(p1, p2, p3);
	vector<double> reverse_p1;
	vector<double> reverse_p2;
	vector<double> reverse_p3;
	reverse_p1.clear();
	reverse_p2.clear();
	reverse_p3.clear();
	reverse_p1.push_back(-p1[0]);
	reverse_p1.push_back(-p1[1]);
	reverse_p1.push_back(-p1[2]);
	reverse_p1.push_back(p1[3]);
	reverse_p2.push_back(-p2[0]);
	reverse_p2.push_back(-p2[1]);
	reverse_p2.push_back(-p2[2]);
	reverse_p2.push_back(p2[3]);
	reverse_p3.push_back(-p3[0]);
	reverse_p3.push_back(-p3[1]);
	reverse_p3.push_back(-p3[2]);
	reverse_p3.push_back(p3[3]);
	complex<double> Abar_sig = tKLpipi.Amp_PFT(reverse_p1, reverse_p3, reverse_p2);

	// Tag
	complex<double> A_tag = tKLpipi.Amp_PFT(p4, p5, p6);
	vector<double> reverse_p4;
	vector<double> reverse_p5;
	vector<double> reverse_p6;
	reverse_p4.clear();
	reverse_p5.clear();
	reverse_p6.clear();
	reverse_p4.push_back(-p4[0]);
	reverse_p4.push_back(-p4[1]);
	reverse_p4.push_back(-p4[2]);
	reverse_p4.push_back(p4[3]);
	reverse_p5.push_back(-p5[0]);
	reverse_p5.push_back(-p5[1]);
	reverse_p5.push_back(-p5[2]);
	reverse_p5.push_back(p5[3]);
	reverse_p6.push_back(-p6[0]);
	reverse_p6.push_back(-p6[1]);
	reverse_p6.push_back(-p6[2]);
	reverse_p6.push_back(p6[3]);
	complex<double> Abar_tag = tKLpipi.Amp_PFT(reverse_p4, reverse_p6, reverse_p5);

	// |A| for each of the amplitude and amplitude bar for the D that goes to Kspipi (A) and the D to the tag (B)
	real_t modA_sig{std::abs(A_sig)};
	real_t modAbar_sig{std::abs(Abar_sig)};
	real_t modA_tag{std::abs(A_tag)};
	real_t modAbar_tag{std::abs(Abar_tag)};

	// phase things
	real_t deltaC_sig{PC.evalBias(event_sig)};
	real_t deltaC_tag{PC.evalBias(event_tag)};
	real_t deltaD_sig{std::arg(A_sig * std::conj(Abar_sig))};
	real_t deltaD_tag{std::arg(A_tag * std::conj(Abar_tag))};

	return std::pow(modA_sig * modAbar_tag, 2) + std::pow(modAbar_sig * modA_tag, 2) - 2 * modA_sig * modAbar_sig * modA_tag * modAbar_tag * cos(deltaD_sig + deltaC_sig - (deltaD_tag + deltaC_tag));
}

//************ GENERATING DATA TYPES (KLpipi):
// CP state
real_t AmpGen::totalAmplitudeSquared_CP_KSpipi(const int &CPsign, D0ToKSpipi2018 &tKLpipi, std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> &cPhi, PhaseCorrection &PC, Event &event)
{
	vector<double> p1, p2, p3;
	for (int i = 0; i < 4; i++)
	{
		p1.push_back(cPhi[4 + i](event.address()));
		p2.push_back(cPhi[8 + i](event.address()));
		p3.push_back(cPhi[12 + i](event.address()));
	}
	complex<double> A = tKLpipi.Amp_PFT(p1, p2, p3);
	// Need to be covert since the exchange of model needed
	vector<double> reverse_p1;
	vector<double> reverse_p2;
	vector<double> reverse_p3;
	reverse_p1.clear();
	reverse_p2.clear();
	reverse_p3.clear();
	reverse_p1.push_back(-p1[0]);
	reverse_p1.push_back(-p1[1]);
	reverse_p1.push_back(-p1[2]);
	reverse_p1.push_back(p1[3]);
	reverse_p2.push_back(-p2[0]);
	reverse_p2.push_back(-p2[1]);
	reverse_p2.push_back(-p2[2]);
	reverse_p2.push_back(p2[3]);
	reverse_p3.push_back(-p3[0]);
	reverse_p3.push_back(-p3[1]);
	reverse_p3.push_back(-p3[2]);
	reverse_p3.push_back(p3[3]);
	complex<double> Abar = tKLpipi.Amp_PFT(reverse_p1, reverse_p3, reverse_p2);

	real_t modA{std::abs(A)};
	real_t modAbar{std::abs(Abar)};
	real_t deltaD{std::arg(A * std::conj(Abar))};
	real_t deltaC{PC.evalBias(event)}; // note will need to change so can implement a bias to ues evalBais also

	return std::pow(modA, 2) + std::pow(modAbar, 2) + 2 * CPsign * modA * modAbar * cos(deltaD + deltaC);
}

// flavour state - isFlavour = true if the SIGNAL is a D, false if SIGNAL is Dbar
real_t AmpGen::totalAmplitudeSquared_flavour_KSpipi(const bool &isFlavour, D0ToKSpipi2018 &tKLpipi, std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> &cPhi, Event &event)
{
	vector<double> p1, p2, p3;
	for (int i = 0; i < 4; i++)
	{
		p1.push_back(cPhi[4 + i](event.address()));
		p2.push_back(cPhi[8 + i](event.address()));
		p3.push_back(cPhi[12 + i](event.address()));
	}
	complex<double> A = tKLpipi.Amp_PFT(p1, p2, p3);
	vector<double> reverse_p1;
	vector<double> reverse_p2;
	vector<double> reverse_p3;
	reverse_p1.clear();
	reverse_p2.clear();
	reverse_p3.clear();
	reverse_p1.push_back(-p1[0]);
	reverse_p1.push_back(-p1[1]);
	reverse_p1.push_back(-p1[2]);
	reverse_p1.push_back(p1[3]);
	reverse_p2.push_back(-p2[0]);
	reverse_p2.push_back(-p2[1]);
	reverse_p2.push_back(-p2[2]);
	reverse_p2.push_back(p2[3]);
	reverse_p3.push_back(-p3[0]);
	reverse_p3.push_back(-p3[1]);
	reverse_p3.push_back(-p3[2]);
	reverse_p3.push_back(p3[3]);
	complex<double> Abar = tKLpipi.Amp_PFT(reverse_p1, reverse_p3, reverse_p2);

	if (isFlavour)
	{
		return std::norm(Abar);
	}
	else
	{
		return std::norm(A);
	}
}

// most general version, mainly for the KLpipi tag but could be attached anywhere
real_t AmpGen::totalAmplitudeSquared_BES_KLpipi(D0ToKSpipi2018 & tKSpipi,  D0ToKLpipi2018 &tKLpipi, std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> &cPhi, PhaseCorrection &PC, PhaseCorrection_KLpipi &PC_KLipi, Event &event_sig, Event &event_tag)
{

	vector<double> p1, p2, p3;
	vector<double> p4, p5, p6;
	for (int i = 0; i < 4; i++)
	{
		p1.push_back(cPhi[4 + i](event_sig.address()));
		p2.push_back(cPhi[8 + i](event_sig.address()));
		p3.push_back(cPhi[12 + i](event_sig.address()));
		
		p4.push_back(cPhi[4 + i](event_tag.address()));
		p5.push_back(cPhi[8 + i](event_tag.address()));
		p6.push_back(cPhi[12 + i](event_tag.address()));
	}
	// Sig
	vector<double> reverse_p1;
	vector<double> reverse_p2;
	vector<double> reverse_p3;
	reverse_p1.clear();
	reverse_p2.clear();
	reverse_p3.clear();
	reverse_p1.push_back(-p1[0]);
	reverse_p1.push_back(-p1[1]);
	reverse_p1.push_back(-p1[2]);
	reverse_p1.push_back(p1[3]);
	reverse_p2.push_back(-p2[0]);
	reverse_p2.push_back(-p2[1]);
	reverse_p2.push_back(-p2[2]);
	reverse_p2.push_back(p2[3]);
	reverse_p3.push_back(-p3[0]);
	reverse_p3.push_back(-p3[1]);
	reverse_p3.push_back(-p3[2]);
	reverse_p3.push_back(p3[3]);
	complex<double> A_sig = tKLpipi.Amp_PFT(p1, p2, p3);

	complex<double> Abar_sig = tKLpipi.Amp_PFT(reverse_p1, reverse_p3, reverse_p2);


	// Tag
	complex<double> A_tag = tKSpipi.Amp_PFT(p4, p5, p6);
	vector<double> reverse_p4;
	vector<double> reverse_p5;
	vector<double> reverse_p6;
	reverse_p4.clear();
	reverse_p5.clear();
	reverse_p6.clear();
	reverse_p4.push_back(-p4[0]);
	reverse_p4.push_back(-p4[1]);
	reverse_p4.push_back(-p4[2]);
	reverse_p4.push_back(p4[3]);
	reverse_p5.push_back(-p5[0]);
	reverse_p5.push_back(-p5[1]);
	reverse_p5.push_back(-p5[2]);
	reverse_p5.push_back(p5[3]);
	reverse_p6.push_back(-p6[0]);
	reverse_p6.push_back(-p6[1]);
	reverse_p6.push_back(-p6[2]);
	reverse_p6.push_back(p6[3]);
	complex<double> Abar_tag = tKSpipi.Amp_PFT(reverse_p4, reverse_p6, reverse_p5);

	// |A| for each of the amplitude and amplitude bar for the D that goes to Kspipi (A) and the D to the tag (B)
	real_t modA_sig{std::abs(A_sig)};
	real_t modAbar_sig{std::abs(Abar_sig)};
	real_t modA_tag{std::abs(A_tag)};
	real_t modAbar_tag{std::abs(Abar_tag)};

	// phase things
	real_t deltaC_sig{PC.evalBias(event_sig)};
	real_t deltaC_tag{PC.evalBias(event_tag)};
	real_t deltaD_sig{std::arg(A_sig * conj(Abar_sig))};
	real_t deltaD_tag{std::arg(A_tag * conj(Abar_tag))};

	return std::pow(modA_sig * modAbar_tag, 2) + std::pow(modAbar_sig * modA_tag, 2) - 2 * modA_sig * modAbar_sig * modA_tag * modAbar_tag * cos(deltaD_sig + deltaC_sig - (deltaD_tag + deltaC_tag));
}
// Shenghui added for using BELLE2010 model
real_t AmpGen::totalAmplitudeSquared_CP_test(const int &CPsign, DtoKpipiAmplitude &amp, std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> &cPhi, PhaseCorrection &PC, Event &event)
{
	dcomplex A = amp.get_amplitude(cPhi[1](event.address()), cPhi[0](event.address()));
	dcomplex Abar = amp.get_amplitude(cPhi[0](event.address()), cPhi[1](event.address()));

	real_t modA{std::abs(A)};
	real_t modAbar{std::abs(Abar)};
	real_t deltaD{std::arg(A * std::conj(Abar)) + M_PI};
	real_t deltaC{PC.evalBias(event)}; // note will need to change so can implement a bias to ues evalBais also

	return std::pow(modA, 2) + std::pow(modAbar, 2) + 2 * CPsign * modA * modAbar * cos(deltaD + deltaC);
}
// most general version, mainly for the Kspipi tag but could be attached anywhere
real_t AmpGen::totalAmplitudeSquared_BES_test(DtoKpipiAmplitude &amp, std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> &cPhi, PhaseCorrection &PC, Event &event_sig, Event &event_tag)
{
	dcomplex A_sig = amp.get_amplitude(cPhi[1](event_sig.address()), cPhi[0](event_sig.address()));
	dcomplex Abar_sig = amp.get_amplitude(cPhi[0](event_sig.address()), cPhi[1](event_sig.address()));
	dcomplex A_tag = amp.get_amplitude(cPhi[1](event_tag.address()), cPhi[0](event_tag.address()));
	dcomplex Abar_tag = amp.get_amplitude(cPhi[0](event_tag.address()), cPhi[1](event_tag.address()));

	// |A| for each of the amplitude and amplitude bar for the D that goes to Kspipi (A) and the D to the tag (B)
	real_t modA_sig{std::abs(A_sig)};
	real_t modAbar_sig{std::abs(Abar_sig)};
	real_t modA_tag{std::abs(A_tag)};
	real_t modAbar_tag{std::abs(Abar_tag)};

	// phase things
	real_t deltaC_sig{PC.eval(event_sig)};
	real_t deltaC_tag{PC.eval(event_tag)};
	real_t deltaD_sig{std::arg(A_sig * std::conj(Abar_sig)) + M_PI};
	real_t deltaD_tag{std::arg(A_tag * std::conj(Abar_tag)) + M_PI};

	return std::pow(modA_sig * modAbar_tag, 2) + std::pow(modAbar_sig * modA_tag, 2) - 2 * modA_sig * modAbar_sig * modA_tag * modAbar_tag * cos(deltaD_sig + deltaC_sig - (deltaD_tag + deltaC_tag));
}

// most general version, mainly for the KLpipi tag but could be attached anywhere
real_t AmpGen::totalAmplitudeSquared_BES_KSpipi_KLpipi(CoherentSum &A, CoherentSum &Abar, D0ToKLpipi2018 &tKLpipi, std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> &cPhi, PhaseCorrection &PC, PhaseCorrection_KLpipi &PC_KLipi, Event &event_sig, Event &event_tag)
{
	// |A| for each of the amplitude and amplitude bar for the D that goes to Kspipi (A) and the D to the tag (B)
	real_t modA_tag{std::abs(A.getValNoCache(event_tag))};
	real_t modAbar_tag{std::abs(Abar.getValNoCache(event_tag))};

	vector<double> p4, p5, p6;
	for (int i = 0; i < 4; i++)
	{
		p4.push_back(cPhi[4 + i](event_sig.address()));
		p5.push_back(cPhi[8 + i](event_sig.address()));
		p6.push_back(cPhi[12 + i](event_sig.address()));
	}

	// Tag
	complex<double> A_sig = tKLpipi.Amp_PFT(p4, p5, p6);
	vector<double> reverse_p4;
	vector<double> reverse_p5;
	vector<double> reverse_p6;
	reverse_p4.clear();
	reverse_p5.clear();
	reverse_p6.clear();
	reverse_p4.push_back(-p4[0]);
	reverse_p4.push_back(-p4[1]);
	reverse_p4.push_back(-p4[2]);
	reverse_p4.push_back(p4[3]);
	reverse_p5.push_back(-p5[0]);
	reverse_p5.push_back(-p5[1]);
	reverse_p5.push_back(-p5[2]);
	reverse_p5.push_back(p5[3]);
	reverse_p6.push_back(-p6[0]);
	reverse_p6.push_back(-p6[1]);
	reverse_p6.push_back(-p6[2]);
	reverse_p6.push_back(p6[3]);
	complex<double> Abar_sig = tKLpipi.Amp_PFT(reverse_p4, reverse_p6, reverse_p5);

	// |A| for each of the amplitude and amplitude bar for the D that goes to Kspipi (A) and the D to the tag (B)
	real_t modA_sig{std::abs(A_sig)};
	real_t modAbar_sig{std::abs(Abar_sig)};

	// phase things
	// phase things
	real_t deltaC_sig{PC_KLipi.evalBias(event_tag)};
	real_t deltaD_sig{Deltadelta_KL(A_sig,Abar_sig)};
	real_t deltaC_tag{PC.evalBias(event_sig)};
	real_t deltaD_tag{Deltadelta(event_tag, A, Abar)};

	// Not set a M_PI here since the phase convention is different from the BELLE2010 model

	return std::pow(modA_sig * modAbar_tag, 2) + std::pow(modAbar_sig * modA_tag, 2) - 2 * modA_sig * modAbar_sig * modA_tag * modAbar_tag * cos(deltaD_sig + deltaC_sig - (deltaD_tag + deltaC_tag));
}

// flavour state - isFlavour = true if the SIGNAL is a D, false if SIGNAL is Dbar
real_t AmpGen::totalAmplitudeSquared_flavour_test(const bool &isFlavour, DtoKpipiAmplitude &amp, std::vector<CompiledExpression<real_t(const real_t *, const real_t *)>> &cPhi, Event &event)
{
	dcomplex A = amp.get_amplitude(cPhi[1](event.address()), cPhi[0](event.address()));
	dcomplex Abar = amp.get_amplitude(cPhi[0](event.address()), cPhi[1](event.address()));

	if (isFlavour)
	{
		return std::norm(Abar);
	}
	else
	{
		return std::norm(A);
	}
}

// //TEMP: no longer used, can delete?
// real_t AmpGen::totalAmplitudeSquared_BES(bool& isKspipi, CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, PhaseCorrection& PC, Event& eventA, Event& eventB)
// {
// 	// |A| for each of the amplitude and amplitude bar for the D that goes to Kspipi (A) and the D to the tag (B)
// 	real_t modA{std::abs(A.getValNoCache(eventA))};
// 	real_t modAbar{std::abs(Abar.getValNoCache(eventA))};
// 	real_t modB{std::abs(B.getValNoCache(eventB))};
// 	real_t modBbar{std::abs(Bbar.getValNoCache(eventB))};

// 	// phase things
// 	real_t deltaC{PC.eval(eventA)};
// 	real_t deltaD_A{Deltadelta(eventA, A, Abar)};
// 	real_t deltaD_B{Deltadelta(eventB, B, Bbar)};

// 	if (isKspipi){
// 		// if the tag is also Kspipi, that one also needs a phase correction in the model
// 		real_t deltaC_B{PC.eval(eventB)};
// 		return  std::pow(modA*modBbar, 2) + std::pow(modAbar*modB, 2) - 2*modA*modAbar*modB*modBbar*cos(deltaD_A + deltaC - deltaD_B - deltaC_B);
// 	}else{
// 		return std::pow(modA*modBbar, 2) + std::pow(modAbar*modB, 2) - 2*modA*modAbar*modB*modBbar*cos(deltaD_A + deltaC - deltaD_B);
// 	}
// }

//********** FITTING DATA TYPES:

// CP state the D->Kspipi is even or odd, sign = sign of the CP of this D, so an tag of Kspi0 for even and KK for odd
real_t AmpGen::totalAmplitudeSquared_CP(const int &CPsign, real_t &modA, real_t &modAbar, real_t &Dd, real_t &correction)
{
	return std::pow(modA, 2) + std::pow(modAbar, 2) + 2 * CPsign * modA * modAbar * cos(Dd + correction);
}

// flavour = the D is pure D^0 or Dbar^0 so a tag like K-pi+ from a Dbar or K+pi- from a D
real_t AmpGen::totalAmplitudeSquared_flavour(real_t &modA)
{
	return std::pow(modA, 2);
}

// same = both D's go to Kspipi so both sets of events need to be considered
real_t AmpGen::totalAmplitudeSquared_BES(real_t &modA, real_t &modAbar, real_t &modA_tag, real_t &modAbar_tag, real_t &deltaD, real_t &deltaD_tag, real_t &deltaC, real_t &deltaC_tag)
{
	return std::pow(modA * modAbar_tag, 2) + std::pow(modAbar * modA_tag, 2) - 2 * modA * modAbar * modA_tag * modAbar_tag * cos(deltaD + deltaC - (deltaD_tag + deltaC_tag)); // should we split up the cosine like in the LHCb formula for statistical reasons?
																																											   // return std::pow(modA*modAbar_tag, 2) + std::pow(modAbar*modA_tag, 2) - 2*modA*modAbar*modA_tag*modAbar_tag*(cos(deltaD + deltaC)*cos(deltaD_tag + deltaC_tag) + sin(deltaD + deltaC)*sin(deltaD_tag + deltaC_tag));
}

real_t AmpGen::totalAmplitudeSquared_BES_KLpipi(real_t &modA, real_t &modAbar, real_t &modA_tag, real_t &modAbar_tag, real_t &deltaD, real_t &deltaD_tag, real_t &deltaC, real_t &deltaC_tag)
{
	return std::pow(modA * modAbar_tag, 2) + std::pow(modAbar * modA_tag, 2) - 2 * modA * modAbar * modA_tag * modAbar_tag * cos(deltaD + deltaC - (deltaD_tag + deltaC_tag)); // should we split up the cosine like in the LHCb formula for statistical reasons?
																																											   // return std::pow(modA*modAbar_tag, 2) + std::pow(modAbar*modA_tag, 2) - 2*modA*modAbar*modA_tag*modAbar_tag*(cos(deltaD + deltaC)*cos(deltaD_tag + deltaC_tag) + sin(deltaD + deltaC)*sin(deltaD_tag + deltaC_tag));
}
// as above but taking event as argument for plotting reasons
// For plot
real_t AmpGen::totalAmplitudeSquared_CP(const int &CPsign, CoherentSum &A, CoherentSum &Abar, PhaseCorrection &PC, const Event &event)
{
	real_t modA{std::abs(A.getValNoCache(event))};
	real_t modAbar{std::abs(Abar.getValNoCache(event))};
	real_t deltaD{Deltadelta(event, A, Abar)};
	real_t deltaC{PC.eval(event)}; // note will need to change so can implement a bias to ues evalBais also

	return std::pow(modA, 2) + std::pow(modAbar, 2) + 2 * CPsign * modA * modAbar * cos(deltaD + deltaC);
}

// flavour state - isFlavour = true if the SIGNAL is a D, false if SIGNAL is Dbar
real_t AmpGen::totalAmplitudeSquared_flavour(const bool &isFlavour, CoherentSum &A, CoherentSum &Abar, const Event &event)
{
	if (isFlavour)
	{
		return std::norm(Abar.getValNoCache(event));
	}
	else
	{
		return std::norm(A.getValNoCache(event));
	}
}
real_t AmpGen::totalAmplitudeSquared_BES(CoherentSum &A, CoherentSum &Abar, PhaseCorrection &PC, const Event &event_main, const Event &event_tag)
{
	// |A| for each of the amplitude and amplitude bar for the D that goes to Kspipi (A) and the D to the tag (B)
	real_t modA_main{std::abs(A.getValNoCache(event_main))};
	real_t modAbar_main{std::abs(Abar.getValNoCache(event_main))};
	real_t modA_tag{std::abs(A.getValNoCache(event_tag))};
	real_t modAbar_tag{std::abs(Abar.getValNoCache(event_tag))};

	// phase things
	real_t deltaC_main{PC.eval(event_main)};
	real_t deltaC_tag{PC.eval(event_tag)};
	real_t deltaD_main{Deltadelta(event_main, A, Abar)};
	real_t deltaD_tag{Deltadelta(event_tag, A, Abar)};

	return std::pow(modA_main * modAbar_tag, 2) + std::pow(modAbar_main * modA_tag, 2) - 2 * modA_main * modAbar_main * modA_tag * modAbar_tag * cos(deltaD_main + deltaC_main - (deltaD_tag + deltaC_tag));
}
