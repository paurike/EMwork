
#include <TClonesArray.h>
//#include "liboaanalysis/ND__TTrackerECALReconModule.h"
//#include "liboaanalysis/ND__TGlobalReconModule.h"
#include <TTree.h>
#include <TH1D.h>
#include <TFile.h>
#include "TSystem.h"
#include <iostream>
#include "TMath.h"
#include "TF1.h"
#include <math.h>
#include <TMinuit.h>
#include <vector>
#include <string>       
#include <iostream>     
#include <sstream>  
#include <TVirtualFitter.h>
#include <TMatrixD.h>



//function declarations

//Crystal Ball is an asymmetric function from high energy physics that's used for corrections; want to try and use this to fit the fractional difference 
Double_t CrystalBall(Double_t *x,Double_t *par); 


//Define possible Functions for Convolutions:
// f1:
double Gaus(double *x, double *par);

//f2:
double Landau(double *x, double *par);


//f3:
double Expo(double *x, double *par);

//Gaussian * Landau
double Func_12( double *x, double *par );

//Gaussian * Exponential
double Func_13( double *x, double *par );


//Convolutions using the two functions defined in Func_12 or Func_13

double Conv_12( double *x, double *par);
double Conv_13(double *x, double *par);


//New fitting function that takes EM-Energy/TPC-Momentum distribution and gives it a scale parameter
double Two_e_fit(double *x, double *par, TH1D *hist);


//calculates the median of a histogram
double Median(const TH1D * h1) ;

//A Superposition of a Gaussian and a Landau
double GausLandaufit(double *x, double *par);


//Functions to extract fit parameters from energy range histograms

void GaussFitEnergyRange(TH1D * hist, Double_t &Sigma, Double_t &SigmaErr, Double_t &GaussMean, Double_t &GaussMeanErr);

void GausLandauFitEnergyRange(TH1D *hist, Double_t fitMin, Double_t fitMax, Double_t &ComponentMean, Double_t &ComponentMeanErr, Double_t &ComponentSigma, Double_t &ComponentSigmaErr, Double_t &FWHM);

void GetArithmeticVars(TH1D *hist, Double_t &Mean, Double_t &MeanErr, Double_t &RMS, Double_t &RMSErr);

void GetFWHM(TF1 func, Double_t &FWHM); 

void PrintFitMatrices(TVirtualFitter * fitrp);



//__________________________________________________________________________________________________________________________

//MAIN FUNCTION

void TailReader(char * filename) {


  //make TFile and read in the file to analyse
  TFile f(filename);


  //Make TTree and set path/to/tree
  TTree * NTuple = (TTree*) f.Get("emEnergyNtuple");

 //Make Variables of interest
  Double_t TPCMomentum, MCMomentum, Energy, ReconFraction, fluxWeight=0,  MCFraction=0; 

  //Variable for containment cut, tpc cut and MC/Data distinction
  Double_t containment = 0;
  Double_t Ntpc = 0;
  Double_t IsData = 0;

  //Set Branchm Addresses
  NTuple->SetBranchAddress("emEnergyCorrected", &Energy);
  NTuple->SetBranchAddress("tpcECalEntranceMomentumMag", &TPCMomentum);
  //NTuple->SetBranchAddress("truthEcalEntranceMomentumMag", &MCMomentum);
  NTuple->SetBranchAddress("containment", &containment);
  //NTuple->SetBranchAddress("tpcNTracks", &Ntpc);
  NTuple->SetBranchAddress("isData", &IsData);
  NTuple->SetBranchAddress("fluxWeighting", &fluxWeight);
  //fluxWeight = 1;

 
  Int_t Nbins = 150;
  Int_t BinsFraction1 = Nbins/3;
  Int_t BinsFraction2 = Nbins/6;
    
  //Make Histogramms to be filled with values of interest
  TH1D * EnergyDist = new TH1D ("EMEnergy", "EM-Energy", Nbins, 0, 15e3);
  TH1D * ReconMomentumDist = new TH1D ("TPCMomentum", "TPC-Momentum", Nbins, 0, 15e3);
  TH1D * MCMomentumDist = new TH1D ("MCMomentum", "MC-Momentum", Nbins, 0, 15e3);
  TH1D * ReconFractionDist = new TH1D ("Energy-TPCMomentum", "(EMEnergy-TPCMomentum)/TPCMomentum", Nbins, -2, 10);
  TH1D * MCFractionDist = new TH1D ("Energy-MCMomentum", "(EMEnergy-MCMomentum)/MCMomentum", Nbins, -2, 10);
  TH1D * TwoTrackReconFractionDist = new TH1D ("Energy-TPCMomentum", "(EMEnergy-TPCMomentum)/TPCMomentum", Nbins, -2, 10);
  TH1D * ThreeTrackReconFractionDist = new TH1D ("Energy-TPCMomentum", "(EMEnergy-TPCMomentum)/TPCMomentum", Nbins, -2, 10);
  TH1D * FourTrackReconFractionDist = new TH1D ("Energy-TPCMomentum", "(EMEnergy-TPCMomentum)/TPCMomentum", Nbins, -2, 10);
  TH1D * Two_e_Dist = new TH1D("EMEnergy/TPCMomentum", "(EMEnergy/TPCMomentum)", Nbins, -2, 10);
 

  Double_t ranges[] = {0, 100, 150, 200, 1600};
  Int_t energyBins = 4;
  TH1D *EnergyRangeHists[energyBins];
  Double_t fitmins[] = {-2, -2, -2, -2};
  Double_t fitmaxs[] = {4, 2, 4, 4};

  for(int i=0; i< energyBins; i++){
    
    std::stringstream HistNameStream;
    HistNameStream << "Range" << i+1 << "FractionDist";
    std::string HistName = HistNameStream.str();
    const char * Name = HistName.c_str();
    EnergyRangeHists[i] = new TH1D(Name, Name, BinsFraction1, -2, 10);

  }


  //Some extra Variables for the loop 
  Int_t Events = (Int_t)NTuple->GetEntries();
  Int_t Reject = 0;

  std::cout<< "Tree contains " << Events << " Events." << std::endl;
  
  
  //loop to fill histogramms

  for(int i=0; i<Events; i++){
    
    NTuple->GetEntry(i);

    if(IsData){fluxWeight = 1;}
 
    
    EnergyDist->Fill(Energy);
    ReconMomentumDist->Fill(TPCMomentum, fluxWeight);
    //MCMomentumDist->Fill(MCMomentum);

    if(containment == 1){

      ReconFraction = (Energy-TPCMomentum)/(TPCMomentum);
      ReconFractionDist->Fill(ReconFraction, fluxWeight);
      std::cout << "Flux Weight was: " << fluxWeight << std::endl;
      Two_e_Dist->Fill(Energy/TPCMomentum);


      for(int j=0; j<energyBins; j++){

	if(ranges[j]<TPCMomentum && TPCMomentum <= ranges[j+1]){
	  
	  EnergyRangeHists[j]->Fill(ReconFraction, fluxWeight);
	  break;
	  
	}

      }

      //if(MCMomentum != 0){
      //MCFraction = (Energy-MCMomentum)/(MCMomentum);
      //MCFractionDist->Fill(MCFraction);


	  
      //} else{Reject+=1;}
     
//       if(Ntpc==2){
// 	TwoTrackReconFractionDist->Fill(ReconFraction);	  
//       }

//       if(Ntpc==3){
// 	ThreeTrackReconFractionDist->Fill(ReconFraction);	  
//       }
      
//       if(Ntpc==4){
// 	FourTrackReconFractionDist->Fill(ReconFraction);	  
//       }
   
    }
    else{Reject+=1;}
    
  }

   std::cout << "Number of rejected Events because of 0 MC-Momentum was: " << Reject << std::endl;



  //Make histogramm with only negative part of fractional difference mirrored at 0
  TH1D * ReconMirrorDist = new TH1D("ReconMirrorDist", "ReconMirrorDist", (Nbins), -2, 10);

  //Make histogramm where later the mirrored distribution gets subtracted from the total distribution leaving only the tail
  TH1D * ReconTailDist = (TH1D*)ReconFractionDist->Clone("ReconTailDist");
  
  Int_t bins = ReconMirrorDist->GetNbinsX(); //should be same as Nbins now, just as a check
  
  //loop to fill the mirrored distribution
  for(int h=0; h<(bins/6); h++){

    ReconMirrorDist->SetBinContent(h,ReconFractionDist->GetBinContent(h));
    ReconMirrorDist->SetBinContent((bins/3-h),ReconFractionDist->GetBinContent(h));
    ReconMirrorDist->SetBinContent(bins/6, ReconFractionDist->GetBinContent(Nbins/6));

  }

  //Subtract mirrored distribution from total fractional difference distribution  
  ReconTailDist->Add(ReconMirrorDist,-1);
  ReconTailDist->SetTitle("Isolating the Tail");




 
  std::cout  << "FITTING section has been reached" << std::endl;
  //___________________________________________________________________________________________________________________________________
  //FITTING


  //Create Fit-functions




  

  std::cout << "MARKER 1" << std::endl;

  //  TF1 *GausRes = new TF1("GausFitResult", "[0]*TMath::Gaus(x, [1],[2])", -2, 10);
//   GausRes->SetParameters(9.29572e2, -2.07390e-1, 2.71188e-1);
//   TF1 *LandauRes= new TF1("LandauFitResutl", "[0]*TMath::Landau(x, [1],[2])", -2, 10);
//   LandauRes->SetParameters(6.75116e3, 1.22781e-1, 1.59924e-1);

  //Crystal Ball function as possible fit function
//   //Crystal ball function for signal, parameters are 0:alpha,1:n,2:mean,3:sigma,4:normalization;
//   TF1 *crystalfit = new TF1("crystalfit",CrystalBall,-2,10,5);
//   crystalfit->SetParameters(-1, 1, ReconFractionDist->GetMean(), ReconFractionDist->GetRMS(), 1660);
//   crystalfit->SetParNames("#alpha","n","Mean","#sigma","N");
//   std::cout<< "Fitting ReconFractionDist with a Crystal Ball function" << std::endl << std::endl;
//   ReconFractionDist->Fit(crystalfit, "", "", -2, 10);
  
  std::cout << "MARKER 2" << std::endl;


  //fit total fractional difference distribution with the landau+gaus
  std::cout<< "Fitting ReconFractionDist with a Superposition of a Gauss and a Landau function" << std::endl << std::endl;
  // ReconFractionDist->Fit("GausLandaufit", "", "", -2, 2);

  std::cout << "MARKER 3" << std::endl;

  //two different ways of fitting the mirrored distribution
  std::cout<< "Fitting ReconMirrorDist with a Gauss function" << std::endl << std::endl;
  ReconMirrorDist->Fit("gaus");
  std::cout<< "Fitting ReconMirrorDist with a Gauss function just from -2 to 0" << std::endl << std::endl;
  ReconMirrorDist->Fit("gaus", "", "", -2, 0);

  std::cout << "MARKER 4" << std::endl;
  
  //fit the ReconTail Dist
  std::cout << "Fitting ReconTailDist with an exponential function" << std::endl << std::endl;
  ReconTailDist->Fit("expo", "", "", 0, 10);
  
  std::cout << "MARKER 5" << std::endl;

  //define a gaussian and give it the parameters of the best fit result 
//   TF1 *LeftHandGaus = new TF1("LeftHandGaus","[2]*TMath::Gaus(x,[0],[1])",-2, 10);
//   LeftHandGaus->SetParameters(6.83601e-2, 3.52550e-1, 1.72509e3);


  std::cout << "MARKER 6" << std::endl;

  //Try to isolate the tail by subtracting the value of LeftHandGaus from ReconFractionDist at the bin center
//   TH1D * TailDist = (TH1D*)ReconFractionDist->Clone("TailDist");
  
//   for(int k=0; k<Nbins; k++){

//     TailDist->SetBinContent(k, (TailDist->GetBinContent(k))-(LeftHandGaus->Eval(TailDist->GetBinCenter(k))));

//   }
  
  std::cout << "MARKER 7" << std::endl;


  //fit the tail (this fit regularly fails)
//   std::cout<< "Fitting TailDist with an exponential function" << std::endl << std::endl;
//   TailDist->Fit("expo", "", "", 1, 10);
//   TailDist->SetTitle("Isolating the Tail");

   
  //GausLandaufit->SetParameters(MCFractionDist->GetMean(), MCFractionDist->GetRMS(), MCFractionDist->GetMean(), MCFractionDist->GetRMS(), 1, 4);
  //MCFractionDist->Fit("GausLandaufit", "", "", -2, 2);


  std::cout << "CONVOLUTION FITS SECTION has been reached" << std::endl;
  //---------------------------------------------------------------------------------------------------------------------------------------------
  //CONVOLUTION FITS


//   TF1 *Convolution12 = new TF1("convolution", Conv_12, -2, 10, 6);
//   Convolution12->SetParNames("GausMean", "GausSigma", "GaussNorm",  "LandauMean", "LandauSigma", "LandauNorm");
//   Convolution12->SetParameters(MCFractionDist->GetMean(), MCFractionDist->GetRMS(),1800, MCFractionDist->GetMean(), MCFractionDist->GetRMS(), 180);
//   std::cout<< "Fitting ReconFractionDist with a Convolution of a Gaus and a Landau" << std::endl << std::endl;
//   ReconFractionDist->Fit(Convolution12, "", "", -2, 10);

  //create functions for Convolution fit results

  std::cout << "Creating any kind of TF1 makes this code break" << std::endl;


//   TF1 *ConvGausRes = new TF1("ConvGausFitResult", "[0]*TMath::Gaus(x, [1],[2])", -2, 10);
//   ConvGausRes->SetParameters(5.24464e2, -1.09486, 2.64800e-1);
//   TF1 *ConvLandauRes= new TF1("ConvLandauFitResutl", "[0]*TMath::Landau(x, [1],[2])", -2, 10);
//   ConvLandauRes->SetParameters(43.2089, 1.00483, 1.12558e-1);


//   TF1 *Convolution13 = new TF1("convolution", Conv_13, -2, 10, 5);
//   Convolution13->SetParNames("GausMean", "GausSigma", "GaussNorm",  "ExpoSlope", "ExpoNorm");
//   Convolution13->SetParameters(MCFractionDist->GetMean(), MCFractionDist->GetRMS(),1800, -10, 6);
//   std::cout<< "Fitting ReconFractionDist with a Convolution of a Gaus and an Exponential" << std::endl << std::endl;
//   ReconFractionDist->Fit(Convolution13, "", "", -2, 10);

  std::cout << "ENERGY RANGES section has been reached" << std::endl;
//-------------------------------------------------------------------------------------------------------------------------------------------------
//ENERGY RANGES 
 

  // Define Width and Mean Histograms to be filled with Fit Results

  TH1D *SigmaVsMom = new TH1D("SigmaVsMom", "SigmaVsMom", energyBins, ranges);
  TH1D *GaussMeanVsMom = new TH1D("GaussMeanVsMom", "GausMeanVsMom", energyBins, ranges);

  TH1D *MeanVsMom = new TH1D("MeanVsMom", "MeanVsMom", energyBins, ranges);
  TH1D *RMSVsMom = new TH1D("RMSVsMom", "RMSVsMom", energyBins, ranges);

  TH1D *ComponentSigmaVsMom = new TH1D("ComponentSigmaVsMom", "ComponentSigmaVsMom", energyBins, ranges);
  TH1D *ComponentMeanVsMom = new TH1D("ComponentMeanVsMom", "ComponentMeanVsMom", energyBins, ranges);

  TH1D *MedianVsMom = new TH1D("MeanVsMom", "MeanVsMom", energyBins, ranges);

  TH1D *FWHMVsMom = new TH1D("FWHMVsMom", "FWHMVsMom", energyBins, ranges);



// Fit every energy range with a gaussian to look at the width with respect to TPC-momentum
// Also fit every range with a GausLandau-Superposition
// Fill Histogramms for truncated arithmetic Mean and RMS of the distribution 
// Fill histogram with Median 


  std::cout << "The loop has been reached" << std::endl;

  for(int i=0; i<energyBins; i++){

    std::cout << "Simple Gauss Fit for " << ranges[i] << " MeV - " << ranges[i+1] << " Mev" << std::endl;

    Double_t Sigma=0; 
    Double_t SigmaErr=0;
    Double_t GaussMean=0; 
    Double_t GaussMeanErr=0;
    std::cout << "Variables have been initialised" << std::endl;
    GaussFitEnergyRange(EnergyRangeHists[i], Sigma, SigmaErr, GaussMean, GaussMeanErr);
    std::cout << "Parameters got set: " << Sigma << ", "<< GaussMean << ", " << SigmaErr << ", " << GaussMeanErr << std::endl;
    SigmaVsMom->SetBinContent(i+1, Sigma);
    SigmaVsMom->SetBinError(i+1, SigmaErr);
    GaussMeanVsMom->SetBinContent(i+1, GaussMean);
    GaussMeanVsMom->SetBinError(i+1, GaussMeanErr);

    std::cout << "GaussResult obtained" << std::endl;

  

    std::cout << "GaussLandau Fit for " << ranges[i] << " MeV - " << ranges[i+1] << " Mev" << std::endl;
    
    Double_t ComponentSigma=0;
    Double_t ComponentSigmaErr=0;
    Double_t ComponentMean=0;
    Double_t ComponentMeanErr=0;
    Double_t FWHM = 0;
    GausLandauFitEnergyRange(EnergyRangeHists[i], fitmins[i], fitmaxs[i], ComponentMean, ComponentMeanErr, ComponentSigma, ComponentSigmaErr, FWHM );
    ComponentMeanVsMom->SetBinContent(i+1, ComponentMean);
    ComponentMeanVsMom->SetBinError(i+1, ComponentMeanErr);
    ComponentSigmaVsMom->SetBinContent(i+1, ComponentSigma);
    ComponentSigmaVsMom->SetBinError(i+1, ComponentSigmaErr);

    FWHMVsMom->SetBinContent(i+1, FWHM);
    FWHMVsMom->SetBinError(i+1, ComponentSigmaErr);

  
    
    Double_t Mean=0;
    Double_t MeanErr=0;
    Double_t RMS=0;
    Double_t RMSErr=0;
    GetArithmeticVars(EnergyRangeHists[i], Mean, MeanErr, RMS, RMSErr);
    MeanVsMom->SetBinContent(i+1, Mean);
    MeanVsMom->SetBinError(i+1, MeanErr);
    RMSVsMom->SetBinContent(i+1, RMS);
    RMSVsMom->SetBinError(i+1, RMSErr); 


    std::cout << "The Median of this histogram is:" << Median(EnergyRangeHists[i]) << std::endl;
    MedianVsMom->SetBinContent(i+1, Median(EnergyRangeHists[i]));
    MedianVsMom->SetBinError(i+1, 1.253*Sigma/sqrt(EnergyRangeHists[i]->GetEntries()));


   

    std::cout << "Loop " << i << "completed" << std::endl;
  }


  //FIt the total ReconFractionDist
  Double_t Sigma;
  Double_t SigmaErr;
  Double_t GaussMean;
  Double_t GaussMeanErr;
  GaussFitEnergyRange(ReconFractionDist, Sigma, SigmaErr, GaussMean, GaussMeanErr);
  std::cout << "Gaussian Fit of total ReconFractionDist:" << std::endl;
  std::cout << "Sigma: " << Sigma <<" +/- " << SigmaErr << std::endl;
  std::cout << "Mean: " << GaussMean << " +/- " << GaussMeanErr << std::endl;

  TF1 * StraightGaussMean = new TF1("StraightGaussMean", "[0]", 0, 1600);
  TF1 * StraightSigma = new TF1("StraightSigma", "[0]", 0, 1600);

  StraightGaussMean->SetParameter(0, GaussMean);
  StraightGaussMean->SetParError(0, 0);
  StraightSigma->SetParameter(0, Sigma);
  StraightSigma->SetParError(0, 0);
  

  TVirtualFitter * fitrp = TVirtualFitter::GetFitter(); 

  PrintFitMatrices(fitrp);

  

  Double_t ComponentSigma;
  Double_t ComponentSigmaErr;
  Double_t ComponentMean; 
  Double_t ComponentMeanErr;
  Double_t FWHM;
  GausLandauFitEnergyRange(ReconFractionDist, fitmins[0], fitmaxs[0], ComponentMean, ComponentMeanErr, ComponentSigma, ComponentSigmaErr, FWHM );




  //Recover ReconFractionDist Fit:

  //Create Fit-functions
  //Superposition of Gaussian and Landau Distribution

  TF1 *GausLandaufit = new TF1("GausLandaufit","[4]*TMath::Gaus(x,[0],[1])+[5]*TMath::Landau(x,[2],[3])",-2, 10);
  GausLandaufit->SetParameters(ReconFractionDist->GetMean(), ReconFractionDist->GetRMS(), ReconFractionDist->GetMean(), ReconFractionDist->GetRMS(), 1, 4);
  GausLandaufit->SetParNames("GausMean", "GausRMS", "LandauMean", "LandauRMS", "GausNorm", "LandauNorm");

  TF1 *GausRes = new TF1("GausFitResult", "[0]*TMath::Gaus(x, [1],[2])", -2, 10);
  GausRes->SetParameters(9.29572e2, -2.07390e-1, 2.71188e-1);
  TF1 *LandauRes= new TF1("LandauFitResutl", "[0]*TMath::Landau(x, [1],[2])", -2, 10);
  LandauRes->SetParameters(6.75116e3, 1.22781e-1, 1.59924e-1);

std::cout<< "Fitting ReconFractionDist with a Superposition of a Gauss and a Landau function" << std::endl << std::endl;
  ReconFractionDist->Fit("GausLandaufit", "", "", -2, 2);

  



//----------------------------------------------------------------------------------------------------------------------------------------------
  std::cout << "Number of entries in:" << std::endl;
  std::cout << "Total Fractional Difference Distribution " << ReconFractionDist->GetEntries() << std::endl;
  std::cout << "0-100 MeV: " << EnergyRangeHists[0]->GetEntries() << std::endl;
  std::cout << "100-150 MeV: " << EnergyRangeHists[1]->GetEntries() << std::endl;
  std::cout << "150-200 MeV: " << EnergyRangeHists[2]->GetEntries() << std::endl;
  std::cout << "200-350 MeV: " << EnergyRangeHists[3]->GetEntries() << std::endl;
//   std::cout << "350-500 MeV: " << EnergyRangeHists[4]->GetEntries() << std::endl;
//   std::cout << "500-1600 MeV: " << EnergyRangeHists[5]->GetEntries() << std::endl;
  std::cout << "Rejected due to containment cut: " << Reject << std::endl;
//---------------------------------------------------------------------------------------------------------------------------------------------
  //DEFINE OUTPUT AND WRITE TO FILE

  TFile *outf;

  if(IsData == 1){
    outf = new TFile("output/DATA_Energy-Momentum-Fraction.root", "RECREATE");}
  else{outf = new TFile("output/MC_Energy-Momentum-Fraction.root", "RECREATE");}

  std::cout << "Output file was named " << outf->GetName() << std::endl;

  ReconFractionDist->Write("ReconFractionDist");
  MCFractionDist->Write("MCFractuionDist");
  ReconMomentumDist->Write("ReconMomentumDist");
  MCMomentumDist->Write("MCMomentumDist");
  EnergyDist->Write("EnergyDist");
  GausLandaufit->Write("GausLandaufit");
  GausRes->Write("GausFitResult");
  LandauRes->Write("LandauFitResult");
  ReconMirrorDist->Write("ReconMirrorDist");
  ReconTailDist->Write("ReconTailDist");
  //crystalfit->Write("crystalfit");
  //LeftHandGaus->Write("LeftHandGaus");
  //TailDist->Write("TailDist");
  //Convolution12->Write("Convolution12");
  //ConvGausRes->Write("ConvGausRes");
  //ConvLandauRes->Write("ConvLandauRes");
  //Convolution13->Write("Convolution13");
//   TwoTrackReconFractionDist->Write("TwoTrackReconFractionDist");
//   ThreeTrackReconFractionDist->Write("ThreeTrackReconFractionDist");
//   FourTrackReconFractionDist->Write("FourTrackReconFractionDist");
  Two_e_Dist->Write("EMEnergy/TPCMomentum");


  for (int i=0; i<energyBins; i++){
    EnergyRangeHists[i]->Write(EnergyRangeHists[i]->GetTitle()); 
  }

  SigmaVsMom->Write("SigmaVsMom");
  GaussMeanVsMom->Write("GaussMeanVsMom");
  MeanVsMom->Write("MeanVsMom");
  RMSVsMom->Write("RMSVsMom");
  ComponentMeanVsMom->Write("ComponentMeanVsMom");
  ComponentSigmaVsMom->Write("ComponentSigmaVsMom");
  MedianVsMom->Write("MedianVsMom");
  FWHMVsMom->Write("FWHMVsMom");

  StraightGaussMean->Write("StraightGaussMean");
  StraightSigma->Write("StraightSigma");

}





//-----------------------------------------------------------------------------------------------------------------------------------------------------
//FUNCTIONS
//useful functions and fit functions


//Crystal Ball is an asymmetric function from high energy physics that's used for corrections; want to try and use this to fit the fractional difference 
Double_t CrystalBall(Double_t *x,Double_t *par) {

  Double_t t = (x[0]-par[2])/par[3];
  if (par[0] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[0]);

  if (t >= -absAlpha) {
    return par[4]*exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= par[1]/absAlpha - absAlpha;
      
    return par[4]*(a/TMath::Power(b - t, par[1]));
  }
}

//Convolution of a Gaussian and a Landau
//Define possible Functions for Convolutions:

// f1:
double Gaus(double *x, double *par)
{
  return par[2]*TMath::Gaus(x[0], par[0], par[1]);
}

//f2:

double Landau(double *x, double *par)
{
  return par[2]*TMath::Landau(x[0], par[0], par[1]);
}

//f3:

double Expo(double *x, double *par)
{
  return par[1]*TMath::Exp(par[0]*x[0]); 
}

//Gaussian * Landau
double Func_12( double *x, double *par )
{
  TF1 gaus("gaus", Gaus, -2, 10, 3);
  TF1 landau("landau", Landau, -2, 10, 3);
  gaus.SetParameters(par[1], par[2], par[3]);
  landau.SetParameters(par[4], par[5], par[6]);
  return gaus.Eval(x[0])*landau.Eval(par[0] - x[0]);
}

//Gaussian * Exponential
double Func_13( double *x, double *par )
{
  TF1 gaus("gaus", Gaus, -2, 10, 3);
  gaus.SetParameters(par[1], par[2], par[3]);
  TF1 expo("exp", Expo, -2, 10, 2);
  expo.SetParameters(par[4], par[5]);
  return gaus.Eval(x[0])*expo.Eval(par[0] - x[0]);
}


//Convolutions using the two functions defined in Func_12 or Func_13
double Conv_12( double *x, double *par)
{
  TF1 f_12("f_12", Func_12, -2, 10, 7);
  f_12.SetParameters( x[0], par[0], par[1], par[2], par[3], par[4], par[5]);
  return f_12.Integral(-2, 10);
}

double Conv_13(double *x, double *par){
  TF1 f_13("f_13", Func_13, -2, 10, 6);
  f_13.SetParameters( x[0], par[0], par[1], par[2], par[3], par[4]);
  return f_13.Integral(-2, 10);
}



  //New fitting function that takes EM-Energy/TPC-Momentum distribution and gives it a scale parameter

double Two_e_fit(double *x, double *par, TH1D *hist){

  hist[0].Scale(par[1]);
  double width = hist[0].GetBinWidth(10);
  int bin = int((x[0]+0.5*width)/width);
  return hist[0].GetBinContent(bin);

  }


double Median(const TH1D * h1) {

   int n = h1->GetXaxis()->GetNbins(); 
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray();
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]);
}



double GausLandaufit(double *x, double *par){

  TF1 GausLandaufit("GausLandaufit","[4]*TMath::Gaus(x,[0],[1])+[5]*TMath::Landau(x,[2],[3])",-2, 10);
 
  GausLandaufit.SetParameters(par[0], par[1], par[2], par[3], par[4], par[5]);

  return GausLandaufit.Eval(x[0]);

 

}



//function that performs Gaussian  fits on the different range histograms and loads fit result into variables Sigma and GaussMean

void GaussFitEnergyRange( TH1D * hist, Double_t &Sigma, Double_t &SigmaErr, Double_t &GaussMean, Double_t &GaussMeanErr) {
  
  //gSystem->Load("libMinuit");

  TF1 Gaussian("Gaussian", "gaus", -0.9, 0.9);
  Gaussian.SetParameters(((hist->GetMaximum())/(hist->GetRMS())), hist->GetMean(), hist->GetRMS());
  hist->Fit("Gaussian", "", "", -0.9, 0.9);
  std::cout << "The simple Gauss fit has run properly" << std::endl;
  std::cout << "The Sigma Fitresult was: " << Gaussian.GetParameter(2) << std::endl;

  Sigma = Gaussian.GetParameter(2);
  SigmaErr = Gaussian.GetParError(2);
  GaussMean = Gaussian.GetParameter(1);
  GaussMeanErr = Gaussian.GetParError(1);

  std::cout <<"Parameters of Gaussian Fit were found as Sigma: " << Sigma << " and Mean: "<< GaussMean << std::endl;
   
 }


//function that performs the GausLandaufit on all the different range histograms and loads fit result into variables ComponentMean and ComponentSigma
void GausLandauFitEnergyRange(TH1D *hist, Double_t fitMin, Double_t fitMax, Double_t &ComponentMean, Double_t &ComponentMeanErr, Double_t &ComponentSigma, Double_t &ComponentSigmaErr, Double_t &FWHM){

  TF1 fit("fit", GausLandaufit, -2, 10, 6);
  fit.SetParNames("GausMean", "GausSigma", "LandauMean", "LandauRMS", "GausNorm", "LandauNorm");
  fit.SetParameters(hist->GetMean(), hist->GetRMS(), hist->GetMean(), hist->GetRMS(), 1, 4);
  hist->Fit("fit", "", "", fitMin, fitMax);
  ComponentMean = fit.GetParameter(0);
  ComponentMeanErr = fit.GetParError(0);
  ComponentSigma = fit.GetParameter(1);
  ComponentSigmaErr = fit.GetParError(1);

  GetFWHM(fit, FWHM);

  std::cout <<"Parameters of GausLandau Fit were found as Sigma: " << ComponentSigma << " and Mean: "<< ComponentMean << std::endl;

}

//function that gets Mean and RMS from a histogram and loads it into variables Mean and RMS

void GetArithmeticVars( TH1D *hist, Double_t &Mean, Double_t &MeanErr, Double_t &RMS, Double_t &RMSErr){


  hist->SetAxisRange(-0.9, 0.9);
  Mean = hist->GetMean();
  MeanErr = hist->GetMeanError();
  RMS = hist->GetRMS();
  RMSErr = hist->GetRMSError();
  hist->SetAxisRange(-2, 4);
  std::cout <<"Arithmetic Parameters were found as RMS: " << RMS  << " and Mean: "<< Mean << std::endl;
}

void GetFWHM (TF1 func, Double_t &FWHM){

  func.SetNpx(10000);
  TH1* hist = (TH1D*)(func.GetHistogram()->Clone());
  
  int bin1 = hist->FindFirstBinAbove((hist->GetMaximum())/2);
  int bin2 = hist->FindLastBinAbove((hist->GetMaximum())/2);
  FWHM = hist->GetBinCenter(bin2 )- hist->GetBinCenter(bin1);

  delete hist;


} 



// // code f r a g ment t o a c c e s s t he c o v a r i a n c e ma t r ix and t o c a l c u l a t e
// //       t he c o r r e l a t i o n ma t r ix a f t e r p e r f o r m i n g a f i t i n Root with t he
// // . F i t ( . . . ) method
// // g e t t he c o v a r i a n c e and c o r r e l a t i o n m a t r i c e s


void PrintFitMatrices(TVirtualFitter * fitrp) {

  int nPar = fitrp->GetNumberTotalParameters();
  TMatrixD covmat(nPar, nPar, fitrp->GetCovarianceMatrix());
 
  std::cout << "The Covariance Matrix is: " << std::endl;
  covmat.Print();

  TMatrixD cormat(covmat);
  for(int i=0; i<nPar; i++){
    for(int j=0; j<nPar; j++){
      cormat(i,j) /= sqrt(covmat(i,i)) * sqrt(covmat(j,j));
    }
  }

  std::cout << "The Correlation Matrix is: " << std::endl;
  cormat.Print();

}


// Double_t StraighLine(double *x, double *val) {
 
//   return val[0];

// }



