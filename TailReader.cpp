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



//function declarations

//Crystal Ball is an asymmetric function from high energy physics that's used for corrections; want to try and use this to fit the fractional difference 
Double_t CrystalBall(Double_t *x,Double_t *par) {


//Define possible Functions for Convolutions:
// f1:
double Gaus(double *x, double *par)

//f2:
double Landau(double *x, double *par)


//f3:
double Expo(double *x, double *par)

//Gaussian * Landau
double Func_12( double *x, double *par )

//Gaussian * Exponential
double Func_13( double *x, double *par )


//Convolutions using the two functions defined in Func_12 or Func_13

double Conv_12( double *x, double *par)
double Conv_13(double *x, double *par)


//New fitting function that takes EM-Energy/TPC-Momentum distribution and gives it a scale parameter
double Two_e_fit(double *x, double *par, TH1D *hist)


//calculates the median of a histogram
double Median(const TH1D * h1) 



//__________________________________________________________________________________________________________________________

//MAIN FUNCTION

void TailReader(char * filename) {


  //make TFile and read in the file to analyse
  TFile f(filename);


  //Make TTree and set path/to/tree
  TTree * NTuple = (TTree*) f.Get("emEnergyNtuple");

 //Make Variables of interest
  Double_t TPCMomentum, MCMomentum, Energy, ReconFraction, fluxWeight,  MCFraction=0; 

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
  NTuple->SetBranchAddress("fluxWeight", &fluxWeight);

 
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
 

  



  TH1D * E100FractionDist = new TH1D("E100FractionDist", "E100FractionDist", BinsFraction1, -2, 10);
  TH1D * Mom100FractionDist = new TH1D("Mom100FractionDist", "Mom100FractionDist", BinsFraction1, -2, 10);
  TH1D * E150FractionDist = new TH1D("E150FractionDist", "E150FractionDist",BinsFraction1 , -2, 10);
  TH1D * Mom150FractionDist = new TH1D("Mom150FractionDist", "Mom150FractionDist",BinsFraction1 , -2, 10);
  TH1D * E200FractionDist = new TH1D("E200FractionDist", "E200FractionDist",BinsFraction1 , -2, 10);
  TH1D * Mom200FractionDist = new TH1D("Mom200FractionDist", "Mom200FractionDist", BinsFraction1, -2, 10);
  TH1D * E300FractionDist = new TH1D("E300FractionDist", "E300FractionDist",BinsFraction1 , -2, 10);
  TH1D * Mom300FractionDist = new TH1D("Mom300FractionDist", "Mom300FractionDist", BinsFraction1, -2, 10);
  TH1D * E400FractionDist = new TH1D("E400FractionDist", "E400FractionDist",BinsFraction1 , -2, 10);
  TH1D * Mom400FractionDist = new TH1D("Mom400FractionDist", "Mom400FractionDist", BinsFraction1, -2, 10);
  TH1D * E800FractionDist = new TH1D("E800FractionDist", "E800FractionDist",BinsFraction1 , -2, 10);
  TH1D * Mom800FractionDist = new TH1D("Mom800FractionDist", "Mom800FractionDist", BinsFraction1, -2, 10);

  //Some extra Variables for the loop 
  Int_t Events = (Int_t)NTuple->GetEntries();
  Int_t Reject = 0;

  std::cout<< "Tree contains " << Events << " Events." << std::endl;
  
  
  //loop to fill histogramms

  for(int i=0; i<Events; i++){
    
    NTuple->GetEntry(i);

    if(IsData){fluxWeight = 1;}
 
    
    EnergyDist->Fill(Energy, fluxWeight);
    ReconMomentumDist->Fill(TPCMomentum, fluxWeight);
    //MCMomentumDist->Fill(MCMomentum);

    if(containment == 1){

      ReconFraction = (Energy-TPCMomentum)/(TPCMomentum);
      ReconFractionDist->Fill(ReconFraction, fluxWeight);
      Two_e_Dist->Fill(Energy/TPCMomentum);


      //Make Histograms for different Energy and Momentum Slices

      if(0<TPCMomentum && TPCMomentum<=100){	
	Mom100FractionDist->Fill(ReconFraction, fluxWeight);
      }

      if(0<Energy && Energy<=100){	
	E100FractionDist->Fill(ReconFraction, fluxWeight);
      }
      
      if(100<TPCMomentum && TPCMomentum<=150){	
	Mom150FractionDist->Fill(ReconFraction, fluxWeight);
      }

      if(100<Energy && Energy<=150){	
	E150FractionDist->Fill(ReconFraction, fluxWeight);
      }

      if(150<TPCMomentum && TPCMomentum<=200){	
	Mom200FractionDist->Fill(ReconFraction, fluxWeight);
      }

      if(150<Energy && Energy<=200){	
	E200FractionDist->Fill(ReconFraction, fluxWeight);
      }


       if(200<TPCMomentum && TPCMomentum<=350){	
	 Mom300FractionDist->Fill(ReconFraction, fluxWeight);
      }

      if(200<Energy && Energy<=350){	
	E300FractionDist->Fill(ReconFraction, fluxWeight);
      }

       if(350<TPCMomentum && TPCMomentum<=500){	
	 Mom400FractionDist->Fill(ReconFraction, fluxWeight);
      }

      if(350<Energy && Energy<=500){	
	E400FractionDist->Fill(ReconFraction, fluxWeight);
      }

       if(500<TPCMomentum && TPCMomentum<=1600){	
	 Mom800FractionDist->Fill(ReconFraction, fluxWeight);
      }

      if(500<Energy && Energy<=1600){	
	E800FractionDist->Fill(ReconFraction, fluxWeight);
      }
      
//        if(700<TPCMomentum && TPCMomentum<1600){	
// 	Mom1600FractionDist->Fill(ReconFraction);
//       }

//       if(700<Energy && Energy<1600){	
// 	E1600FractionDist->Fill(ReconFraction);
//       }






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




 
  
  //___________________________________________________________________________________________________________________________________
  //FITTING


  //Create Fit-functions
  //Superposition of Gaussian and Landau Distribution
  TF1 *GausLandaufit = new TF1("GausLandaufit","[4]*TMath::Gaus(x,[0],[1])+[5]*TMath::Landau(x,[2],[3])",-2, 10);
  GausLandaufit->SetParameters(ReconFractionDist->GetMean(), ReconFractionDist->GetRMS(), ReconFractionDist->GetMean(), ReconFractionDist->GetRMS(), 1, 4);
  GausLandaufit->SetParNames("GausMean", "GausRMS", "LandauMean", "LandauRMS", "GausNorm", "LandauNorm");

  TF1 *GausRes = new TF1("GausFitResult", "[0]*TMath::Gaus(x, [1],[2])", -2, 10);
  GausRes->SetParameters(9.29572e2, -2.07390e-1, 2.71188e-1);
  TF1 *LandauRes= new TF1("LandauFitResutl", "[0]*TMath::Landau(x, [1],[2])", -2, 10);
  LandauRes->SetParameters(6.75116e3, 1.22781e-1, 1.59924e-1);

  //Crystal Ball function as possible fit function
//   //Crystal ball function for signal, parameters are 0:alpha,1:n,2:mean,3:sigma,4:normalization;
//   TF1 *crystalfit = new TF1("crystalfit",CrystalBall,-2,10,5);
//   crystalfit->SetParameters(-1, 1, ReconFractionDist->GetMean(), ReconFractionDist->GetRMS(), 1660);
//   crystalfit->SetParNames("#alpha","n","Mean","#sigma","N");
//   std::cout<< "Fitting ReconFractionDist with a Crystal Ball function" << std::endl << std::endl;
//   ReconFractionDist->Fit(crystalfit, "", "", -2, 10);
  

  //fit total fractional difference distribution with the landau+gaus
  std::cout<< "Fitting ReconFractionDist with a Superposition of a Gauss and a Landau function" << std::endl << std::endl;
  ReconFractionDist->Fit("GausLandaufit", "", "", -2, 2);

  //two different ways of fitting the mirrored distribution
  std::cout<< "Fitting ReconMirrorDist with a Gauss function" << std::endl << std::endl;
  ReconMirrorDist->Fit("gaus");
  std::cout<< "Fitting ReconMirrorDist with a Gauss function just from -2 to 0" << std::endl << std::endl;
  ReconMirrorDist->Fit("gaus", "", "", -2, 0);
  
  //fit the ReconTail Dist
  std::cout << "Fitting ReconTailDist with an exponential function" << std::endl << std::endl;
  ReconTailDist->Fit("expo", "", "", 0, 10);

  //define a gaussian and give it the parameters of the best fit result 
  TF1 *LeftHandGaus = new TF1("LeftHandGaus","[2]*TMath::Gaus(x,[0],[1])",-2, 10);
  LeftHandGaus->SetParameters(6.83601e-2, 3.52550e-1, 1.72509e3);


  //Try to isolate the tail by subtracting the value of LeftHandGaus from ReconFractionDist at the bin center
  TH1D * TailDist = (TH1D*)ReconFractionDist->Clone("TailDist");
  
  for(int k=0; k<Nbins; k++){

    TailDist->SetBinContent(k, (TailDist->GetBinContent(k))-(LeftHandGaus->Eval(TailDist->GetBinCenter(k))));

  }


  //fit the tail (this fit regularly fails)
  std::cout<< "Fitting TailDist with an exponential function" << std::endl << std::endl;
  TailDist->Fit("expo", "", "", 1, 10);
  TailDist->SetTitle("Isolating the Tail");

   
  //GausLandaufit->SetParameters(MCFractionDist->GetMean(), MCFractionDist->GetRMS(), MCFractionDist->GetMean(), MCFractionDist->GetRMS(), 1, 4);
  //MCFractionDist->Fit("GausLandaufit", "", "", -2, 2);



  //---------------------------------------------------------------------------------------------------------------------------------------------
  //CONVOLUTION FITS


//   TF1 *Convolution12 = new TF1("convolution", Conv_12, -2, 10, 6);
//   Convolution12->SetParNames("GausMean", "GausSigma", "GaussNorm",  "LandauMean", "LandauSigma", "LandauNorm");
//   Convolution12->SetParameters(MCFractionDist->GetMean(), MCFractionDist->GetRMS(),1800, MCFractionDist->GetMean(), MCFractionDist->GetRMS(), 180);
//   std::cout<< "Fitting ReconFractionDist with a Convolution of a Gaus and a Landau" << std::endl << std::endl;
//   ReconFractionDist->Fit(Convolution12, "", "", -2, 10);

  //create functions for Convolution fit results
  TF1 *ConvGausRes = new TF1("ConvGausFitResult", "[0]*TMath::Gaus(x, [1],[2])", -2, 10);
  ConvGausRes->SetParameters(5.24464e2, -1.09486, 2.64800e-1);
  TF1 *ConvLandauRes= new TF1("ConvLandauFitResutl", "[0]*TMath::Landau(x, [1],[2])", -2, 10);
  ConvLandauRes->SetParameters(43.2089, 1.00483, 1.12558e-1);


//   TF1 *Convolution13 = new TF1("convolution", Conv_13, -2, 10, 5);
//   Convolution13->SetParNames("GausMean", "GausSigma", "GaussNorm",  "ExpoSlope", "ExpoNorm");
//   Convolution13->SetParameters(MCFractionDist->GetMean(), MCFractionDist->GetRMS(),1800, -10, 6);
//   std::cout<< "Fitting ReconFractionDist with a Convolution of a Gaus and an Exponential" << std::endl << std::endl;
//   ReconFractionDist->Fit(Convolution13, "", "", -2, 10);



  Double_t ranges[] = {0, 100, 150, 200, 350, 500, 1600};
  Int_t energyBins = 6;
 

  // Define Width and Mean Histograms to be filled with Fit Results

  TH1D *SigmaVsMom = new TH1D("SigmaVsMom", "SigmaVsMom", energyBins, ranges);
  TH1D *GaussMeanVsMom = new TH1D("GaussMeanVsMom", "GausMeanVsMom", energyBins, ranges);


  TH1D *MeanVsMom = new TH1D("MeanVsMom", "MeanVsMom", energyBins, ranges);
  TH1D *RMSVsMom = new TH1D("RMSVsMom", "RMSVsMom", energyBins, ranges);

  TH1D *ComponentSigmaVsMom = new TH1D("ComponentSigmaVsMom", "ComponentSigmaVsMom", energyBins, ranges);
  TH1D *ComponentMeanVsMom = new TH1D("ComponentMeanVsMom", "ComponentMeanVsMom", energyBins, ranges);

  TH1D *MedianVsMom = new TH1D("MeanVsMom", "MeanVsMom", energyBins, ranges);



// Fit every energy range with a gaussian to look at the width with respect to TPC-momentum
// Also fit every range with a GausLandau-Superposition
// Fill Histogramms for truncated arithmetic Mean and RMS of the distribution 
// Fill histogram with Median 



  TF1 *Gaussian = new TF1("Gaussian", "gaus", -0.9, 0.9);

  gSystem->Load("libMinuit");



  std::cout << "100MeV gauss fit:" << std::endl;

  Gaussian->SetParameters(((Mom100FractionDist->GetMaximum())/(Mom100FractionDist->GetRMS())), Mom100FractionDist->GetMean(), Mom100FractionDist->GetRMS());
  Mom100FractionDist->Fit("Gaussian","","", -0.9, 0.9);
  Double_t Mom100Sigma, Mom100SigmaErr;
  gMinuit->GetParameter(2, Mom100Sigma, Mom100SigmaErr);
  SigmaVsMom->SetBinContent(1, Mom100Sigma);
  SigmaVsMom->SetBinError(1, Mom100SigmaErr);
  
 //  Mom100FractionDist->Fit("GausLandaufit", "", "", -2, 4);
//   Double_t ComponentMom100Mean, ComponentMom100MeanErr;
//   Double_t ComponentMom100Sigma, ComponentMom100SigmaErr;
//   gMinuit->GetParameter(0, ComponentMom100Mean, ComponentMom100MeanErr);
//   gMinuit->GetParameter(1, ComponentMom100Sigma, ComponentMom100SigmaErr);
//   ComponentMeanVsMom->SetBinContent(1, ComponentMom100Mean);
//   ComponentMeanVsMom->SetBinError(1, ComponentMom100MeanErr);
//   ComponentSigmaVsMom->SetBinContent(1, ComponentMom100Sigma);
//   ComponentSigmaVsMom->SetBinError(1, ComponentMom100SigmaErr);
    
  Mom100FractionDist->SetAxisRange(-0.9, 0.9);
  MeanVsMom->SetBinContent(1, Mom100FractionDist->GetMean());
  MeanVsMom->SetBinError(1, Mom100FractionDist->GetMeanError());
  RMSVsMom->SetBinContent(1, Mom100FractionDist->GetRMS());
  RMSVsMom->SetBinError(1, Mom100FractionDist->GetRMSError());
  

  Mom100FractionDist->SetAxisRange(-2, 4);
  


  

  std::cout << "150MeV gauss fit:" << std::endl;

  Gaussian->SetParameters(((Mom150FractionDist->GetMaximum())/(Mom150FractionDist->GetRMS())), Mom150FractionDist->GetMean(), Mom150FractionDist->GetRMS());
  Mom150FractionDist->Fit("Gaussian", "", "", -0.9, 0.9);
  Double_t Mom150Sigma, Mom150SigmaErr;
  gMinuit->GetParameter(2, Mom150Sigma, Mom150SigmaErr);
  SigmaVsMom->SetBinContent(2, Mom150Sigma);
  SigmaVsMom->SetBinError(2, Mom150SigmaErr);

//   Mom150FractionDist->Fit("GausLandaufit", "", "", -2, 4);
//   Double_t ComponentMom150Mean, ComponentMom150MeanErr;
//   Double_t ComponentMom150Sigma, ComponentMom150SigmaErr;
//   gMinuit->GetParameter(0, ComponentMom150Mean, ComponentMom150MeanErr);
//   gMinuit->GetParameter(1, ComponentMom150Sigma, ComponentMom150SigmaErr);
//   ComponentMeanVsMom->SetBinContent(2, ComponentMom150Mean);
//   ComponentMeanVsMom->SetBinError(2, ComponentMom150MeanErr);
//   ComponentSigmaVsMom->SetBinContent(2, ComponentMom150Sigma);
//   ComponentSigmaVsMom->SetBinError(2, ComponentMom150SigmaErr);

  Mom150FractionDist->SetAxisRange(-0.9, 0.9);
  MeanVsMom->SetBinContent(2, Mom150FractionDist->GetMean());
  MeanVsMom->SetBinError(2, Mom150FractionDist->GetMeanError());
  RMSVsMom->SetBinContent(2, Mom150FractionDist->GetRMS());
  RMSVsMom->SetBinError(2, Mom150FractionDist->GetRMSError());
  Mom150FractionDist->SetAxisRange(-2, 4);





  std::cout << "200MeV gauss fit:" << std::endl;

  Gaussian->SetParameters(((Mom200FractionDist->GetMaximum())/(Mom200FractionDist->GetRMS())), Mom200FractionDist->GetMean(), Mom200FractionDist->GetRMS());
  Mom200FractionDist->Fit("Gaussian", "", "", -0.9, 0.9);
  Double_t Mom200Sigma, Mom200SigmaErr;
  gMinuit->GetParameter(2, Mom200Sigma, Mom200SigmaErr);
  SigmaVsMom->SetBinContent(3, Mom200Sigma);
  SigmaVsMom->SetBinError(3, Mom200SigmaErr);

//   Mom200FractionDist->Fit("GausLandaufit", "", "", -2, 4);
//   Double_t ComponentMom200Mean, ComponentMom200MeanErr;
//   Double_t ComponentMom200Sigma, ComponentMom200SigmaErr;
//   gMinuit->GetParameter(0, ComponentMom200Mean, ComponentMom200MeanErr);
//   gMinuit->GetParameter(1, ComponentMom200Sigma, ComponentMom200SigmaErr);
//   ComponentMeanVsMom->SetBinContent(3, ComponentMom200Mean);
//   ComponentMeanVsMom->SetBinError(3, ComponentMom200MeanErr);
//   ComponentSigmaVsMom->SetBinContent(3, ComponentMom200Sigma);
//   ComponentSigmaVsMom->SetBinError(3, ComponentMom200SigmaErr);

 
  Mom200FractionDist->SetAxisRange(-0.9, 0.9);
  MeanVsMom->SetBinContent(3, Mom200FractionDist->GetMean());
  MeanVsMom->SetBinError(3, Mom200FractionDist->GetMeanError());
  RMSVsMom->SetBinContent(3, Mom200FractionDist->GetRMS());
  RMSVsMom->SetBinError(3, Mom200FractionDist->GetRMSError());
  Mom200FractionDist->SetAxisRange(-2, 4);




  std::cout << "300MeV gauss fit:" << std::endl;

  Gaussian->SetParameters(((Mom300FractionDist->GetMaximum())/(Mom300FractionDist->GetRMS())), Mom300FractionDist->GetMean(), Mom300FractionDist->GetRMS());
  Mom300FractionDist->Fit("Gaussian", "", "", -0.9, 0.9);
  Double_t Mom300Sigma, Mom300SigmaErr;
  gMinuit->GetParameter(2, Mom300Sigma, Mom300SigmaErr);
  SigmaVsMom->SetBinContent(4, Mom300Sigma);
  SigmaVsMom->SetBinError(4, Mom300SigmaErr);

 //  Mom300FractionDist->Fit("GausLandaufit", "", "", -2, 4);
//   Double_t ComponentMom300Mean, ComponentMom300MeanErr;
//   Double_t ComponentMom300Sigma, ComponentMom300SigmaErr;
//   gMinuit->GetParameter(0, ComponentMom300Mean, ComponentMom300MeanErr);
//   gMinuit->GetParameter(1, ComponentMom300Sigma, ComponentMom300SigmaErr);
//   ComponentMeanVsMom->SetBinContent(4, ComponentMom300Mean);
//   ComponentMeanVsMom->SetBinError(4, ComponentMom300MeanErr);
//   ComponentSigmaVsMom->SetBinContent(4, ComponentMom300Sigma);
//   ComponentSigmaVsMom->SetBinError(4, ComponentMom300SigmaErr);

  Mom300FractionDist->SetAxisRange(-0.9, 0.9);
  MeanVsMom->SetBinContent(4, Mom300FractionDist->GetMean());
  MeanVsMom->SetBinError(4, Mom300FractionDist->GetMeanError());
  RMSVsMom->SetBinContent(4, Mom300FractionDist->GetRMS());
  RMSVsMom->SetBinError(4, Mom300FractionDist->GetRMSError());
  Mom300FractionDist->SetAxisRange(-2, 4);




  std::cout << "400MeV gauss fit:" << std::endl;

  Gaussian->SetParameters(((Mom400FractionDist->GetMaximum())/(Mom400FractionDist->GetRMS())), Mom400FractionDist->GetMean(), Mom400FractionDist->GetRMS());
  Mom400FractionDist->Fit("Gaussian", "", "", -0.9, 0.9);
  Double_t Mom400Sigma, Mom400SigmaErr;
  gMinuit->GetParameter(2, Mom400Sigma, Mom400SigmaErr);
  SigmaVsMom->SetBinContent(5, Mom400Sigma);
  SigmaVsMom->SetBinError(5, Mom400SigmaErr);

//   Mom400FractionDist->Fit("GausLandaufit", "", "", -1,1 );
//   Double_t ComponentMom400Mean, ComponentMom400MeanErr;
//   Double_t ComponentMom400Sigma, ComponentMom400SigmaErr;
//   gMinuit->GetParameter(0, ComponentMom400Mean, ComponentMom400MeanErr);
//   gMinuit->GetParameter(1, ComponentMom400Sigma, ComponentMom400SigmaErr);
//   ComponentMeanVsMom->SetBinContent(5, ComponentMom400Mean);
//   ComponentMeanVsMom->SetBinError(5, ComponentMom400MeanErr);
//   ComponentSigmaVsMom->SetBinContent(5, ComponentMom400Sigma);
//   ComponentSigmaVsMom->SetBinError(5, ComponentMom400SigmaErr);

  Mom400FractionDist->SetAxisRange(-0.9, 0.9);
  MeanVsMom->SetBinContent(5, Mom400FractionDist->GetMean());
  MeanVsMom->SetBinError(5, Mom400FractionDist->GetMeanError());
  RMSVsMom->SetBinContent(5, Mom400FractionDist->GetRMS());
  RMSVsMom->SetBinError(5, Mom400FractionDist->GetRMSError());
  Mom400FractionDist->SetAxisRange(-2, 4);




  std::cout << "800MeV gauss fit:" << std::endl;

  Gaussian->SetParameters(((Mom800FractionDist->GetMaximum())/(Mom800FractionDist->GetRMS())), Mom800FractionDist->GetMean(), Mom800FractionDist->GetRMS());
  Mom800FractionDist->Fit("Gaussian", "", "", -0.9, 0.9);
  Double_t Mom800Sigma, Mom800SigmaErr;
  gMinuit->GetParameter(2, Mom800Sigma, Mom800SigmaErr);
  SigmaVsMom->SetBinContent(6, Mom800Sigma);
  SigmaVsMom->SetBinError(6, Mom800SigmaErr);

//   Mom800FractionDist->Fit("GausLandaufit", "", "", -1, 1);
//   Double_t ComponentMom800Mean, ComponentMom800MeanErr;
//   Double_t ComponentMom800Sigma, ComponentMom800SigmaErr;
//   gMinuit->GetParameter(0, ComponentMom800Mean, ComponentMom800MeanErr);
//   gMinuit->GetParameter(1, ComponentMom800Sigma, ComponentMom800SigmaErr);
//   ComponentMeanVsMom->SetBinContent(6, ComponentMom800Mean);
//   ComponentMeanVsMom->SetBinError(6, ComponentMom800MeanErr);
//   ComponentSigmaVsMom->SetBinContent(6, ComponentMom800Sigma);
//   ComponentSigmaVsMom->SetBinError(6, ComponentMom800SigmaErr);

  Mom800FractionDist->SetAxisRange(-0.9, 0.9);
  MeanVsMom->SetBinContent(6, Mom800FractionDist->GetMean());
  MeanVsMom->SetBinError(6, Mom800FractionDist->GetMeanError());
  RMSVsMom->SetBinContent(6, Mom800FractionDist->GetRMS());
  RMSVsMom->SetBinError(6, Mom800FractionDist->GetRMSError());
  Mom800FractionDist->SetAxisRange(-2, 4);



//   std::cout << "1600MeV gauss fit:" << std::endl;
//   Mom1600FractionDist->Fit("gaus", "", "", -0.9, 0.9);
//   Double_t Mom1600Mean, Mom1600MeanErr;
//   gMinuit->GetParameter(2, Mom1600Mean, Mom1600MeanErr);
//   SigmaVsMom->SetBinContent(7, Mom1600Mean);
//   SigmaVsMom->SetBinError(7, Mom1600MeanErr);

//   Mom1600FractionDist->Fit("GausLandaufit", "", "", -1, 3);

//   Mom1600FractionDist->SetAxisRange(-0.9, 0.9);
//   MeanVsMom->SetBinContent(7, Mom1600FractionDist->GetMean());
//   MeanVsMom->SetBinError(7, Mom1600FractionDist->GetMeanError());
//   RMSVsMom->SetBinContent(7, Mom1600FractionDist->GetRMS());
//   RMSVsMom->SetBinError(7, Mom1600FractionDist->GetRMSError());
//   Mom1600FractionDist->SetAxisRange(-2, 4);


//----------------------------------------------------------------------------------------------------------------------------------------------
  std::cout << "Number of entries in:" << std::endl;
  std::cout << "Total Fractional Difference Distribution " << ReconFractionDist->GetEntries() << std::endl;
  std::cout << "0-100 MeV: " << Mom100FractionDist->GetEntries() << std::endl;
  std::cout << "100-150 MeV: " << Mom150FractionDist->GetEntries() << std::endl;
  std::cout << "150-200 MeV: " << Mom200FractionDist->GetEntries() << std::endl;
  std::cout << "200-350 MeV: " << Mom300FractionDist->GetEntries() << std::endl;
  std::cout << "350-500 MeV: " << Mom400FractionDist->GetEntries() << std::endl;
  std::cout << "500-1600 MeV: " << Mom800FractionDist->GetEntries() << std::endl;
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
  LeftHandGaus->Write("LeftHandGaus");
  TailDist->Write("TailDist");
  //Convolution12->Write("Convolution12");
  ConvGausRes->Write("ConvGausRes");
  ConvLandauRes->Write("ConvLandauRes");
  //Convolution13->Write("Convolution13");
//   TwoTrackReconFractionDist->Write("TwoTrackReconFractionDist");
//   ThreeTrackReconFractionDist->Write("ThreeTrackReconFractionDist");
//   FourTrackReconFractionDist->Write("FourTrackReconFractionDist");
  Two_e_Dist->Write("EMEnergy/TPCMomentum");
  E100FractionDist->Write("E100FractionDist");
  Mom100FractionDist->Write("Mom100FractionDist");
  E150FractionDist->Write("E150FractionDist");
  Mom150FractionDist->Write("Mom150FractionDist");
  E200FractionDist->Write("E200FractionDist");
  Mom200FractionDist->Write("Mom200FractionDist");
  E300FractionDist->Write("E300FractionDist");
  Mom300FractionDist->Write("Mom300FractionDist");
  E400FractionDist->Write("E400FractionDist");
  Mom400FractionDist->Write("Mom400FractionDist");
  E800FractionDist->Write("E800FractionDist");
  Mom800FractionDist->Write("Mom800FractionDist");
//E1600FractionDist->Write("E1600FractionDist");
//Mom1600FractionDist->Write("Mom1600FractionDist");
  SigmaVsMom->Write("SigmaVsMom");
  MeanVsMom->Write("MeanVsMom");
  RMSVsMom->Write("RMSVsMom");
  ComponentMeanVsMom->Write("ComponentMeanVsMom");
  ComponentSigmaVsMom->Write("ComponentSigmaVsMom");

}


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

