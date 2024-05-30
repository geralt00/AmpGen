#include <chrono>
#include <ctime>
#include <iostream>
#include <map>
#include <ratio>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include "AmpGen/Integrator.h"
#include "AmpGen/Chi2Estimator.h"
#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/BackgroundPdf.h"
#include "AmpGen/FlatBackgroundPdf.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/Generator.h"
#include "AmpGen/ErrorPropagator.h"
#ifdef _OPENMP
  #include <omp.h>
  #include <thread>
#endif

#if ENABLE_AVX
  #include "AmpGen/EventListSIMD.h"
  using EventList_type = AmpGen::EventListSIMD;
#else
  #include "AmpGen/EventList.h"
  using EventList_type = AmpGen::EventList; 
#endif

#include <TH1.h>
#include <TFile.h>
#include <TRandom3.h>

using namespace AmpGen;

void randomizeStartingPoint( MinuitParameterSet& mps, TRandom3& rand, double multiplier)
{
  for (auto& param : mps) {
    if ( ! param->isFree() || param->name() == "Px" || param->name() == "Py" || param->name() == "Pz" ) continue;
    double min = param->minInit();
    double max = param->maxInit();
    //double multiplier = 100;
    double new_value = rand.Uniform(param->mean()-multiplier*param->stepInit(),param->mean()+multiplier*param->stepInit());
    //double new_value = rand.Uniform(param->mean()-18*param->stepInit(),param->mean()+18*param->stepInit());
    if( min != 0 && max != 0 )
      new_value = rand.Uniform(min,max);
    param->setInit( new_value );
    param->setCurrentFitVal( new_value );
    INFO( param->name() << "  = " << param->mean() << " " << param->stepInit() );
  }
}

double calcDervivative(Minimiser mini, MinuitParameter* param, EventList current_dataset, SumPDF<EventList_type, CoherentSum&> pdf, double h, size_t j)
{
  double m = param->mean();
  param->setCurrentFitVal( m + h);
  double newLL = mini.FCN();
  double pdf_plus_variation = pdf.operator()(current_dataset[j]);
  param->setCurrentFitVal( m - h);
  newLL = mini.FCN();
  double pdf_minus_variation = pdf.operator()(current_dataset[j]);
  // Reset param to best fit value
  param->setCurrentFitVal(m);
  return (pdf_plus_variation - pdf_minus_variation) / (2*h);
}

void calcAsymptoticCorrectedCovariance(std::vector<EventList_type> data, std::vector<SumPDF<EventList_type, CoherentSum&>> pdfs, MinuitParameterSet& MPS, FitResult* fr, Minimiser mini)
{
  // Step 1 - get number of floated params and covariance matrix (reduced)
  TMatrixD covReduced = fr->getReducedCovariance();
  int nFloated = covReduced.GetNcols();
  //TMatrixD covFull = fr->cov();

  // Step 2 - initialise new cov matrix
  //TMatrixTSym<Double_t> num(nFloated);
  TMatrixD num(nFloated, nFloated);
  for (int k = 0; k < nFloated; k++) {
    for (int l = 0; l < nFloated; l++) {
      num(k, l) = 0.0;
      }
   }

  // Step 3 - get list of floated params and make sure they match with cov matrix index
  int nTot = MPS.size();
  std::vector<int> floatingIndex(nFloated,0);
  int idx = 0;
  for (int p = 0; p < nTot; p++) {
    auto param = MPS.at(p);
    if (param->isFree() == 1) {
      floatingIndex[idx] = p;
      idx += 1;
    }
  }

  // Perform calculation of D matrix (num) here
  for (size_t i = 0; i < data.size(); i++) {
    auto& current_dataset = data[i];
    for (size_t j = 0; j < current_dataset.size(); j++) {
    //for (size_t j = 0; j < 3; j++) {

      double weightSquared = current_dataset.weight(j)*current_dataset.weight(j);
      double pdf_val = pdfs[i].operator()(current_dataset[j]);
      std::vector<double> diffs(nFloated, 0.0);

      for (int k = 0; k < nFloated; k++) {
        // Here we need to populate diffs with correct derivatives
        auto parameter = MPS.at(floatingIndex[k]);

        // Use a Richardson extrapolation with two terms for derivative estimate
        // What is the best step size (h1) to use here??
        double h1 = std::sqrt( covReduced(k, k) );
        //double h1 = parameter->stepInit();
        //double h1 = 0.001;
        double h2 = h1/2;
        double nh1 = calcDervivative(mini, parameter, current_dataset, pdfs[i], h1, j);
        double nh2 = calcDervivative(mini, parameter, current_dataset, pdfs[i], h2, j);
        diffs[k] = diffs[k] + nh2 + (nh2-nh1)/3;
      }
    for (int k = 0; k < nFloated; k++) {
      for (int l = 0; l < nFloated; l++) {
        num(k, l) += weightSquared * diffs[k] * diffs[l] / (pdf_val*pdf_val);
        }
      }
    }
  }
  //num.Similarity(covFull);
  TMatrixD covNew(nFloated, nFloated);
  // Perform matrix multiplication CDC here
  for (int i = 0; i < nFloated; i++) {
    for (int j = 0; j < nFloated; j++) {
      for (int k = 0; k < nFloated; k++) {
        for (int l = 0; l < nFloated; l++) {
          covNew(i,j) += covReduced(i,k) * num(k,l) * covReduced(l,j);
        }
      }
    }
  }

  // Covariance matrix at given index
  std::cout << "Covariance matrix:" << std::endl;
  for (size_t i = 0; i < (size_t)covReduced.GetNrows(); ++i){
    for (size_t j = 0; j < (size_t)covReduced.GetNrows(); ++j) std::cout << covReduced[i][j] << " ";
   std::cout << "\n";
  }
  std::cout << "New covariance matrix" << std::endl;
  for (size_t i = 0; i < (size_t)covNew.GetNrows(); ++i){
    for (size_t j = 0; j < (size_t)covNew.GetNrows(); ++j) std::cout << covNew[i][j] << " ";
   std::cout << "\n";
  }
  //for (size_t i = 0; i < (size_t)num.GetNrows(); ++i){
  //  for (size_t j = 0; j < (size_t)num.GetNrows(); ++j) std::cout << num[i][j] << " ";
  //std::cout << "\n";
  //}

  // Errors are just sqrt of the covariance matrix diagonals!
}

int main( int argc, char* argv[] )
{
  /* The user specified options must be loaded at the beginning of the programme, 
     and these can be specified either at the command line or in an options file. */   
  OptionsParser::setArgs( argc, argv );

  /* Parameters that have been parsed can be accessed anywhere in the program 
     using the NamedParameter<T> class. The name of the parameter is the first option,
     then the default value, and then the help string that will be printed if --h is specified 
     as an option. */
  const auto datasets        = NamedParameter<std::string>("Datasets","",
      "List of data/simulated samples to fit, in the format \
      \033[3m data[0] sim[0] data[1] sim[1] ... \033[0m. \nIf a simulated sample is specified FLAT, uniformly generated phase-space events are used for integrals ").getVector();
  //std::string dataFile = NamedParameter<std::string>("DataSample", ""          , "Name of file containing data sample to fit." );
  //std::string intFile  = NamedParameter<std::string>("IntegrationSample",""    , "Name of file containing events to use for MC integration.");
  std::string logFile  = NamedParameter<std::string>("LogFile"   , "Fitter.log",     "Name of the output log file");
  std::string plotFile = NamedParameter<std::string>("Plots"     , "plots.root",     "Name of the output plot file");
  const size_t NBins   = NamedParameter<size_t>     ("nBins"     , 100         , "Number of bins used for plotting.");
  const std::string weight_branch         = NamedParameter<std::string>("WeightBranch","","Name of branch containing event weights.");
  const std::string mc_weight_branch      = NamedParameter<std::string>("MCWeightBranch","","Name of branch containing event weights.");

  std::string outOptFile = NamedParameter<std::string>("OutputOptionFile", ""  , "Name of output option file updated with the best-fit parameters");
  std::string inOptFile = NamedParameter<std::string>("InputOptionFile", argv[1] , "Name of input option file to use as template for OutputOptionFile");

  // Add named parameter to allow fixing of a given parameter
  const auto scan_string = NamedParameter<std::string>("ScanParameter", "", "Name and values for parameter to overwrite in options file");  
  
  auto bNames = NamedParameter<std::string>("Branches", std::vector<std::string>()
              ,"List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
  auto MCbNames = NamedParameter<std::string>("MCBranches", std::vector<std::string>(),
                                              "List of branch names, assumed to be \033[3m daughter1_px ... daughter1_E, daughter2_px ... \033[0m" ).getVector();
  auto pNames = NamedParameter<std::string>("EventType" , ""    
              , "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 

  auto pNames2 = NamedParameter<std::string>("EventType2" , ""    
              , "EventType to fit, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 

  double multiplier = NamedParameter<double>("Multiplier", 5, "Multiplier to apply to step size to randomize starting point");
 
  // [[maybe_unused]]
  size_t      nThreads = NamedParameter<size_t>     ("nCores"    , 8           , "Number of threads to use" );
  size_t      seed     = NamedParameter<size_t>     ("Seed"      , 0           , "Random seed used" );
   
  //if( datasets == "" ) FATAL("Must specify input with option " << italic_on << "DataSample" << italic_off );
  if( pNames.size() == 0 ) FATAL("Must specify event type with option " << italic_on << " EventType" << italic_off);

  TRandom3 rndm;
  rndm.SetSeed( seed );
  gRandom = &rndm;

  INFO("Outputs: LogFile: " << logFile << "; Plots: " << plotFile << "; Options: " << outOptFile);

#if ENABLE_AVX
#endif
  
#ifdef _OPENMP
  omp_set_num_threads( nThreads );
  INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
#endif

  MinuitParameterSet MPS;
  MPS.loadFromStream();

  if ( scan_string.size() > 1 ){
  double newval = ::atof(scan_string.getVal(1).c_str());
  MPS[scan_string.getVal(0)]->setVal(newval);
  MPS[scan_string.getVal(0)]->fix();
  }

  if ( NamedParameter<bool>("RandomizeStartingPoint",false) ){
    randomizeStartingPoint(MPS,rndm,multiplier);
    std::cout << "Starting Point Randomized" << std::endl;
  }
  EventType evtType(pNames);
//  EventType evtType2(pNames2);

  std::vector<EventList_type> events; 
  std::vector<EventList_type> eventsMC;
//  for(size_t i=0;i < datasets.size() ; i+=2 ){
//    events.emplace_back( datasets[i], evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weight_branch)  );
//    if( datasets[i+1] == "FLAT" ) eventsMC.emplace_back( Generator<>(evtType, &rndm).generate(2.5e6) );
//    else eventsMC.emplace_back( datasets[i+1], evtType, Branches(MCbNames), GetGenPdf(true), WeightBranch(mc_weight_branch) );
//  } 
  events.emplace_back( datasets[0], evtType, Branches(bNames), GetGenPdf(false), WeightBranch(weight_branch)  );
  if( datasets[1] == "FLAT" ) eventsMC.emplace_back( Generator<>(evtType, &rndm).generate(2.5e6) );
  else eventsMC.emplace_back( datasets[1], evtType, Branches(MCbNames), GetGenPdf(true), WeightBranch(mc_weight_branch) );
/*  events.emplace_back( datasets[2], evtType2, Branches(bNames), GetGenPdf(false), WeightBranch(weight_branch)  );
  eventsMC.emplace_back( datasets[3], evtType2, Branches(MCbNames), GetGenPdf(true), WeightBranch(mc_weight_branch) );
*/
  TFile* output = TFile::Open( plotFile.c_str(), "RECREATE" ); output->cd();
  bool sig_only;
  

  CoherentSum sig(evtType, MPS);
  BackgroundPdf bkg(evtType,MPS);
/*
  CoherentSum sig_fav(evtType2, MPS);
  BackgroundPdf bkg_fav(evtType2,MPS);
  FlatBackgroundPdf flat(evtType2,MPS);
*/
  sig.setMC(eventsMC[0]);
  bkg.setMC(eventsMC[0]);

/*
  sig_fav.setMC(eventsMC[1]);
  bkg_fav.setMC(eventsMC[1]);
  flat.setMC(eventsMC[1]);
*/
  bkg.prepare();
//  bkg_fav.prepare();
//  flat.prepare();


  auto LL1 = [&events, &eventsMC, &MPS, &sig, &bkg](){
	sig.prepare();

  real_t  ll_running{0};
  real_t  normalisationsig{0};
  size_t nEvents{events[0].size()};
  auto evalA = sig.amplitudeEvaluator(&events[0]);
  auto evalB = bkg.evaluator(&events[0]);
//  auto evalA_MC = sig.amplitudeEvaluator(&eventsMC[0]);
//  for (size_t i =0; i < eventsMC[0].size(); i++){
//    normalisationsig += pow(abs(evalA_MC(eventsMC[0][i])),2);
//  }
//  normalisationsig = normalisationsig/eventsMC[0].size();
  for (size_t i=0; i < nEvents; i++){
                complex_t thisA = evalA(events[0][i]);
                real_t probabilitysig = pow(abs(thisA),2);
                real_t probabilitybkg = evalB(events[0][i]);
                ll_running += log(0.9396*(probabilitysig/sig.norm())+0.0604*(probabilitybkg));
        }
        return -2*ll_running;
  };
/*
  auto LL2 = [&events, &eventsMC, &MPS,&sig,&sig_fav,&bkg_fav, &flat](){
	sig_fav.prepare();
  real_t  ll_running{0};
  real_t  normalisationsig{0};
  size_t nEvents{events[1].size()};
  auto evalA = sig_fav.amplitudeEvaluator(&events[1]);
  auto evalB = bkg_fav.evaluator(&events[1]);
  auto evalB_flat = flat.evaluator(&events[1]);

        for (size_t i=0; i < nEvents; i++){
                complex_t thisA = evalA(events[1][i]);
                real_t probabilitysig = pow(abs(thisA),2);
                real_t probabilitybkg = evalB(events[1][i]);
                real_t probabilityflat = evalB_flat(events[1][i]);
                ll_running += log(0.7377*(probabilitysig/sig_fav.norm())+0.2073*(probabilitybkg)+0.055*(probabilityflat));
        }
        return -2*ll_running;
  };
*/

auto likelihood = [&LL1] (){

  return LL1();
};

  Minimiser mini( likelihood, &MPS );
  mini.doFit();
  FitResult* fr = new FitResult(mini);
//  FitResult* fr_fav = new FitResult(mini);
  auto fitFractions = sig.fitFractions( fr->getErrorPropagator() );
//  auto fitFractions_fav = sig_fav.fitFractions( fr_fav->getErrorPropagator() );
  fr->addFractions( fitFractions );
//  fr_fav->addFractions( fitFractions_fav );
  fr->writeToFile(logFile);
//  fr_fav->writeToFile(logFile);

  /* Calculate the `fit fractions` using the signal model and the error propagator (i.e. 
     fit results + covariance matrix) of the fit result, and write them to a file. 
   */

  //fr->writeToFile( logFile.c_str() );
  if ( outOptFile != "" ) fr->writeOptions( outOptFile, inOptFile );
  
  output->Close();




}

real_t corrected_pdfsig(complex_t thisA,complex_t thisAbar,double rDG,double sigmaG,double RG,MinuitParameterSet& MPS){
        real_t Atotal;
        double radians = sigmaG * (M_PI / 180.0);
        double norm_tag = MPS["Br_CF"]->mean();
        complex_t rotation(cos(radians), sin(radians));
        complex_t forexp = rotation*complex_t(0.0, -1.0);
        real_t three = real(exp(forexp)*conj(thisA)*thisAbar);
        Atotal = (pow(abs(thisA),2) + rDG*rDG*pow(abs(thisAbar),2)-2*rDG*RG*three)*norm_tag;
        return Atotal;//|Atotal|^2 will be the pdf
}


