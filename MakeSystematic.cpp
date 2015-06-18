
#include <TClonesArray.h>
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



void MakeSystematic(Double_t MeanError, Double_t MeanError1, Double_t MeanError2, Double_t SigmaError, Double_t SigmaError1, Double_t SigmaError2){


  std::cout << "RUNNING MAKESYSTEMATIC NOW" << std::endl;


  TFile f("output/ElectronPairMC_Energy-Momentum-Fraction.root");

  TF1* CleanGaussian;
  CleanGaussian = (TF1*)f.Get("CleanGaussian")->Clone();

  TF1* CleanGaussian1;
  CleanGaussian1 = (TF1*)f.Get("CleanGaussian1")->Clone();

  TF1* CleanGaussian2;
  CleanGaussian2 = (TF1*)f.Get("CleanGaussian2")->Clone();



  TFile g("output/ElectronPairDataTailFitted.root");

  TF1* FitGaussian;
  FitGaussian = (TF1*)g.Get("Gaussian")->Clone();
  
  TF1* FitGaussian1;
  FitGaussian1 = (TF1*)g.Get("Gaussian1")->Clone();

  TF1* FitGaussian2;
  FitGaussian2 = (TF1*)g.Get("Gaussian2")->Clone();

  Double_t HighLowRegimes[] = {0, 200, 1200};
  Int_t Regimes = 2;


  TH1D* MeanSystematic = new TH1D("MeanSystematic", "MeanSystematic", Regimes, HighLowRegimes);
  TH1D* SigmaSystematic = new TH1D("SigmaSystematic", "SigmaSystematic", Regimes, HighLowRegimes);

  
    
  

  //extract Parameters for Total Fit

//   Double_t CleanNorm = CleanGaussian->GetParameter(0);
//   Double_t CleanNormErr = CleanGaussian->GetParError(0); 
  
  Double_t CleanMean = CleanGaussian->GetParameter(1);
  Double_t CleanMeanErr = CleanGaussian->GetParError(1);
  
  Double_t CleanSigma = CleanGaussian->GetParameter(2);
  Double_t CleanSigmaErr = CleanGaussian->GetParError(2);

 //  Double_t FitNorm = FitGaussian->GetParameter(0);
//   Double_t FitNormErr = FitGaussian->GetParError(0);
  
  
  Double_t FitMean = FitGaussian->GetParameter(1);
  Double_t FitMeanErr = MeanError;
  
  Double_t FitSigma = FitGaussian->GetParameter(2);
  Double_t FitSigmaErr = SigmaError;



  Double_t MeanSys = FitMean-CleanMean;
  Double_t MeanSysErr = sqrt(TMath::Power(CleanMeanErr,2)+TMath::Power(FitMeanErr, 2));

  Double_t SigmaSys = FitSigma-CleanSigma;
  Double_t SigmaSysErr =sqrt(TMath::Power(CleanSigmaErr,2)+TMath::Power(FitSigmaErr, 2));

  std::cout << "Mean Systematic is: "<< MeanSys << "+/-" << MeanSysErr << std::endl;
  std::cout << "Sigma Systematic is: " << SigmaSys << "+/-" << SigmaSysErr << std::endl;




 //  //extract Parameters for first Regime

  Double_t CleanMean1 = CleanGaussian1->GetParameter(1);
  Double_t CleanMeanErr1 = CleanGaussian1->GetParError(1);
  
  Double_t CleanSigma1 = CleanGaussian1->GetParameter(2);
  Double_t CleanSigmaErr1 = CleanGaussian1->GetParError(2);


  Double_t FitMean1 = FitGaussian1->GetParameter(1);
  Double_t FitMeanErr1 = MeanError1;
  
  Double_t FitSigma1 = FitGaussian1->GetParameter(2);
  Double_t FitSigmaErr1 = SigmaError1;


  Double_t MeanSys1 = FitMean1-CleanMean1;
  Double_t MeanSysErr1 = sqrt(TMath::Power(CleanMeanErr1,2)+TMath::Power(FitMeanErr1, 2));

  Double_t SigmaSys1 = FitSigma1-CleanSigma1;
  Double_t SigmaSysErr1 =sqrt(TMath::Power(CleanSigmaErr1,2)+TMath::Power(FitSigmaErr1, 2));


//   //extract Paramters for second Regime

  Double_t CleanMean2 = CleanGaussian2->GetParameter(1);
  Double_t CleanMeanErr2 = CleanGaussian2->GetParError(1);
  
  Double_t CleanSigma2 = CleanGaussian2->GetParameter(2);
  Double_t CleanSigmaErr2 = CleanGaussian2->GetParError(2);


  Double_t FitMean2 = FitGaussian2->GetParameter(1);
  Double_t FitMeanErr2 = MeanError2;
  
  Double_t FitSigma2 = FitGaussian2->GetParameter(2);
  Double_t FitSigmaErr2 = SigmaError2;


  Double_t MeanSys2 = FitMean2-CleanMean2;
  Double_t MeanSysErr2 = sqrt(TMath::Power(CleanMeanErr2,2)+TMath::Power(FitMeanErr2, 2));

  Double_t SigmaSys2 = FitSigma2-CleanSigma2;
  Double_t SigmaSysErr2 =sqrt(TMath::Power(CleanSigmaErr2,2)+TMath::Power(FitSigmaErr2, 2));




  MeanSystematic->SetBinContent(1, MeanSys1);
  MeanSystematic->SetBinError(1, MeanSysErr1);

  MeanSystematic->SetBinContent(2, MeanSys2);
  MeanSystematic->SetBinError(2, MeanSysErr2);


  SigmaSystematic->SetBinContent(1, SigmaSys1);
  SigmaSystematic->SetBinError(1, SigmaSysErr1);

  SigmaSystematic->SetBinContent(2, SigmaSys2);
  SigmaSystematic->SetBinError(2, SigmaSysErr2);


  TF1 * StraightMean = new TF1("StraightMean", "[0]", 0, 1200);
  TF1 * StraightSigma = new TF1("StraightSigma", "[0]", 0, 1200);


  //Error Bands
  TF1 * StraightMeanErrUp = new TF1("StraightMean", "[0]+[1]", 0, 1200);
  TF1 * StraightSigmaErrUp = new TF1("StraightSigma", "[0]+[1]", 0, 1200);

  TF1 * StraightMeanErrDown = new TF1("StraightMean", "[0]-[1]", 0, 1200);
  TF1 * StraightSigmaErrDown = new TF1("StraightSigma", "[0]-[1]", 0, 1200);

  StraightMean->SetParameter(0,MeanSys);
  StraightSigma->SetParameter(0,SigmaSys);

  StraightMeanErrUp->SetParameters(MeanSys, MeanSysErr);
  StraightMeanErrDown->SetParameters(MeanSys, MeanSysErr);

  StraightSigmaErrUp->SetParameters(SigmaSys, SigmaSysErr);
  StraightSigmaErrDown->SetParameters(SigmaSys, SigmaSysErr);


  TFile outf("output/CombinedFitSystematic.root", "RECREATE");

  SigmaSystematic->Write("SigmaSystematic");
  MeanSystematic->Write("MeanSystematic");

  StraightMean->Write("StraightMean");
  StraightSigma->Write("StraightSigma");

  StraightMeanErrUp->Write("StraightMeanErrUp");
  StraightMeanErrDown->Write("StraightMeanErrDown");

  StraightSigmaErrUp->Write("StraightSigmaErrUp");
  StraightSigmaErrDown->Write("StraightSigmaErrDown");

  outf.Close();

}
