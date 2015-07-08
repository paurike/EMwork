#include <TTree.h>
#include <TH1D.h>
#include <TFile.h>
#include "TSystem.h"
#include <iostream>
#include "TMath.h"
#include "TF1.h"
#include <math.h>
#include <TCanvas.h>

//Function Declaration
double ChiSquared(TH1D * observed, TF1* flatHypothesis);


void Systematics() {

   Double_t ranges[] = {0, 100, 300, 1600};
   Int_t energyBins = 3;

   //Double_t ranges[] = {0, 200, 250, 1600};
   //Int_t energyBins = 3;

  TFile f("output/ElectronPairMC_Energy-Momentum-Fraction.root");

  TH1D *MC_SigmaVsMom;
  TH1D *MC_GaussMeanVsMom;
  TH1D *MC_MeanVsMom;
  TH1D *MC_RMSVsMom;
  TH1D *MC_ComponentSigmaVsMom;
  TH1D *MC_ComponentMeanVsMom;
  TH1D *MC_FWHMVsMom;
  TH1D *MC_MedianVsMom;
  TF1 *MC_StraightGaussMean;
  TF1 *MC_StraightSigma;

  MC_SigmaVsMom = (TH1D*)f.Get("SigmaVsMom")->Clone();
  MC_GaussMeanVsMom =(TH1D*)f.Get("GaussMeanVsMom")->Clone();
  MC_MeanVsMom = (TH1D*)f.Get("MeanVsMom")->Clone();
  MC_RMSVsMom = (TH1D*)f.Get("RMSVsMom")->Clone();
  MC_ComponentSigmaVsMom = (TH1D*)f.Get("ComponentSigmaVsMom")->Clone();
  MC_ComponentMeanVsMom = (TH1D*)f.Get("ComponentMeanVsMom")->Clone();
  MC_FWHMVsMom = (TH1D*)f.Get("FWHMVsMom")->Clone();
  MC_MedianVsMom = (TH1D*)f.Get("MedianVsMom")->Clone();
  MC_StraightGaussMean =(TF1*)f.Get("StraightGaussMean"); 
  MC_StraightGaussMean->SetName("MC_StraightGaussMean");
  MC_StraightSigma =(TF1*)f.Get("StraightSigma");
  MC_StraightSigma->SetName("MC_StraightSigma");

  TFile g("output/ElectronPairDATA_Energy-Momentum-Fraction.root");
   
  TH1D *DATA_SigmaVsMom;
  TH1D *DATA_GaussMeanVsMom;
  TH1D *DATA_MeanVsMom;
  TH1D *DATA_RMSVsMom;
  TH1D *DATA_ComponentSigmaVsMom;
  TH1D *DATA_ComponentMeanVsMom;
  TH1D *DATA_FWHMVsMom;
  TH1D *DATA_MedianVsMom;
  TF1 *DATA_StraightGaussMean;
  TF1 *DATA_StraightSigma;

  DATA_SigmaVsMom = (TH1D*)g.Get("SigmaVsMom")->Clone();
  DATA_GaussMeanVsMom =(TH1D*)g.Get("GaussMeanVsMom")->Clone();
  DATA_MeanVsMom = (TH1D*)g.Get("MeanVsMom")->Clone();
  DATA_RMSVsMom = (TH1D*)g.Get("RMSVsMom")->Clone();
  DATA_ComponentSigmaVsMom = (TH1D*)g.Get("ComponentSigmaVsMom")->Clone();
  DATA_ComponentMeanVsMom = (TH1D*)g.Get("ComponentMeanVsMom")->Clone();
  DATA_FWHMVsMom = (TH1D*)g.Get("FWHMVsMom")->Clone();
  DATA_MedianVsMom = (TH1D*)g.Get("MedianVsMom")->Clone();
  DATA_StraightGaussMean = (TF1*)g.Get("StraightGaussMean");
  DATA_StraightGaussMean->SetName("DATA_StraightGaussMean");
  DATA_StraightSigma = (TF1*)g.Get("StraightSigma");
  DATA_StraightSigma->SetName("DATA_StraightSigma");
  

  Double_t MC=0, DATA=0, MC_error=0, DATA_error=0, Error=0;



  TH1D *SigmaVsMom = new TH1D("SigmaVsMom", "MC-DATA of Sigma of Fractional Difference", energyBins, ranges);
  TH1D *GaussMeanVsMom = new TH1D("GaussMeanVsMom", "MC-DATA of Arithmetic Mean of Fractional Difference", energyBins, ranges);
  TH1D *MeanVsMom = new TH1D("MeanVsMom", "MC-DATA of Arithmetic Mean of Fractional Difference", energyBins, ranges);
  TH1D *RMSVsMom = new TH1D("RMSVsMom", "MC-DATA of RMS of Fractional Difference", energyBins, ranges);
  TH1D *ComponentSigmaVsMom = new TH1D("ComponentSigmaVsMom", "MC-DATA of gaussian Fit Component Sigma of Fractional Difference", energyBins, ranges);
  TH1D *ComponentMeanVsMom = new TH1D("ComponentMeanVsMom", "MC-DATA of gaussian Fit Component Arithmetic Mean of Fractional Difference", energyBins, ranges);
  TH1D *FWHMVsMom = new TH1D("FWHMVsMom", "MC-DATA of FWHM of GaussLandauFit", energyBins, ranges);
  TH1D *MedianVsMom = new TH1D ("MedianVsMom", "MedianVsMom", energyBins, ranges);

  TF1 *StraightGaussMean = new TF1("StraightGaussMean", "[0]+DATA_StraightGaussMean-MC_StraightGaussMean", 0, 1600);
  StraightGaussMean->SetParameter(0, 0);
  TF1 *StraightSigma = new TF1("StraightSigma", "[0]+DATA_StraightSigma-MC_StraightSigma", 0, 1600);
  StraightSigma->SetParameter(0, 0);

  Double_t MC_GaussMeanError = MC_StraightGaussMean->GetParError(0);
  Double_t DATA_GaussMeanError = DATA_StraightGaussMean->GetParError(0);

  Double_t MC_SigmaError = MC_StraightSigma->GetParError(0);
  Double_t DATA_SigmaError = DATA_StraightSigma->GetParError(0);

  //Double_t SigmaError = sqrt(pow(MC_SigmaError, 2) + pow(DATA_SigmaError, 2));
  //Double_t GaussMeanError = sqrt(pow(MC_GaussMeanError, 2) + pow(DATA_GaussMeanError, 2));

  Double_t SigmaError = sqrt(pow(0.003, 2) + pow(0.012, 2));
  Double_t GaussMeanError = sqrt(pow(0.003, 2) + pow(0.012, 2));

  std::cout << "GaussMeanError was propogated as: " << GaussMeanError << std::endl;

  StraightGaussMean->SetParError(0, 0);
  StraightSigma->SetParError(0, 0);

  std::cout << "GaussMeanError was propogated as: " << StraightGaussMean->GetParError(0) << std::endl;

  TF1 *StraightGaussMeanError1 = new TF1("StraightGaussMeanError1", "StraightGaussMean+[0]", 0, 1600);
  StraightGaussMeanError1->SetParameter(0, GaussMeanError);
  TF1 *StraightGaussMeanError2 = new TF1("StraightGaussMeanError2", "StraightGaussMean-[0]", 0, 1600);
  StraightGaussMeanError2->SetParameter(0, GaussMeanError);

  TF1 *StraightSigmaError1 = new TF1("StraightSigmaError1", "StraightSigma+[0]", 0, 1600);
  StraightSigmaError1->SetParameter(0, SigmaError);
  TF1 *StraightSigmaError2 = new TF1("StraightSigmaError2", "StraightSigma-[0]", 0, 1600);
  StraightSigmaError2->SetParameter(0, SigmaError);
  

  for(int i=1; i<energyBins+1; i++){

    MC = MC_SigmaVsMom->GetBinContent(i);
    DATA = DATA_SigmaVsMom->GetBinContent(i);

    SigmaVsMom->SetBinContent(i, DATA-MC);

    MC_error = MC_SigmaVsMom->GetBinError(i);
    DATA_error = DATA_SigmaVsMom->GetBinError(i);

    Error = sqrt(pow(MC_error, 2) + pow(DATA_error, 2));

    SigmaVsMom->SetBinError(i, Error);


  }

  for(int i=1; i<energyBins+1; i++){

    MC = MC_GaussMeanVsMom->GetBinContent(i);
    DATA = DATA_GaussMeanVsMom->GetBinContent(i);

    GaussMeanVsMom->SetBinContent(i, DATA-MC);

    MC_error = MC_GaussMeanVsMom->GetBinError(i);
    DATA_error = DATA_GaussMeanVsMom->GetBinError(i);

    Error = sqrt(pow(MC_error, 2) + pow(DATA_error, 2));

    GaussMeanVsMom->SetBinError(i, Error);


  }


    for(int i=1; i<energyBins+1; i++){

    MC = MC_MeanVsMom->GetBinContent(i);
    DATA = DATA_MeanVsMom->GetBinContent(i);

    MeanVsMom->SetBinContent(i, DATA-MC);

    MC_error = MC_MeanVsMom->GetBinError(i);
    DATA_error = DATA_MeanVsMom->GetBinError(i);

    Error = sqrt(pow(MC_error, 2) + pow(DATA_error, 2));

    MeanVsMom->SetBinError(i, Error);


  }


  for(int i=1; i<energyBins+1; i++){

    MC = MC_RMSVsMom->GetBinContent(i);
    DATA = DATA_RMSVsMom->GetBinContent(i);

    RMSVsMom->SetBinContent(i, DATA-MC);

    MC_error = MC_RMSVsMom->GetBinError(i);
    DATA_error = DATA_RMSVsMom->GetBinError(i);

    Error = sqrt(pow(MC_error, 2) + pow(DATA_error, 2));

    RMSVsMom->SetBinError(i, Error);


  }


  for(int i=1; i<energyBins+1; i++){

    MC = MC_ComponentSigmaVsMom->GetBinContent(i);
    DATA = DATA_ComponentSigmaVsMom->GetBinContent(i);

    ComponentSigmaVsMom->SetBinContent(i, DATA-MC);

    MC_error = MC_ComponentSigmaVsMom->GetBinError(i);
    DATA_error = DATA_ComponentSigmaVsMom->GetBinError(i);

    Error = sqrt(pow(MC_error, 2) + pow(DATA_error, 2));

    ComponentSigmaVsMom->SetBinError(i, Error);


  }


  for(int i=1; i<energyBins+1; i++){

    MC = MC_ComponentMeanVsMom->GetBinContent(i);
    DATA = DATA_ComponentMeanVsMom->GetBinContent(i);

    ComponentMeanVsMom->SetBinContent(i, DATA-MC);

    MC_error = MC_ComponentMeanVsMom->GetBinError(i);
    DATA_error = DATA_ComponentMeanVsMom->GetBinError(i);

    Error = sqrt(pow(MC_error, 2) + pow(DATA_error, 2));

    ComponentMeanVsMom->SetBinError(i, Error);


  }

  for(int i=1; i<energyBins+1; i++){

    MC = MC_FWHMVsMom->GetBinContent(i);
    DATA = DATA_FWHMVsMom->GetBinContent(i);

    FWHMVsMom->SetBinContent(i, DATA-MC);

    MC_error = MC_FWHMVsMom->GetBinError(i);
    DATA_error = DATA_FWHMVsMom->GetBinError(i);

    Error = sqrt(pow(MC_error, 2) + pow(DATA_error, 2));

    FWHMVsMom->SetBinError(i, Error);


  }

 for(int i=1; i<energyBins+1; i++){

   MC = MC_MedianVsMom->GetBinContent(i);
   DATA = DATA_MedianVsMom->GetBinContent(i);

    MedianVsMom->SetBinContent(i, DATA-MC);

    MC_error = MC_MedianVsMom->GetBinError(i);
    DATA_error = DATA_MedianVsMom->GetBinError(i);

    Error = sqrt(pow(MC_error, 2) + pow(DATA_error, 2));

    MedianVsMom->SetBinError(i, Error);


  }


 //Calculating P-Value for Flat Hypothesis

 double MeanChiFlat = ChiSquared(GaussMeanVsMom,StraightGaussMean);
 double SigmaChiFlat = ChiSquared(SigmaVsMom, StraightSigma);

 std::cout << "ChiSquared for flat hypothesis for MEAN: " << MeanChiFlat << std::endl;
 std::cout << "ChiSquared for flat hypothesis for SIGMA: " << SigmaChiFlat << std::endl;




 TCanvas * canvas = new TCanvas();
 canvas->cd();
 SigmaVsMom->Draw();

 canvas->Update();
 ComponentSigmaVsMom->Draw("same");
 canvas->Update();
 

 TFile outf ("output/ElectronPair_GaussianSystematics.root", "RECREATE");

 SigmaVsMom->Write("SigmaVsMom");
 MC_SigmaVsMom->Write("MC_SigmaVsMom");
 DATA_SigmaVsMom->Write("DATA_SigmaVsMom");
 
  GaussMeanVsMom->Write("GaussMeanVsMom");
  MC_GaussMeanVsMom->Write("MC_GaussMeanVsMom");
  DATA_GaussMeanVsMom->Write("DATA_GaussMeanVsMom");

  MeanVsMom->Write("MeanVsMom");
  MC_MeanVsMom->Write("MC_MeanVsMom");
  DATA_MeanVsMom->Write("DATA_MeanVsMom");
  
  RMSVsMom->Write("RMSVsMom");
  MC_RMSVsMom->Write("MC_RMSVsMom");
  DATA_RMSVsMom->Write("DATA_RMSVsMom");

  ComponentMeanVsMom->Write("ComponentMeanVsMom");
  ComponentSigmaVsMom->Write("ComponentSigmaVsMom");

  FWHMVsMom->Write("FWHMVsMom");
  MC_FWHMVsMom->Write("MC_FWHMVsMom");
  DATA_FWHMVsMom->Write("DATA_FWHM");

  MedianVsMom->Write("MedianVsMom");
  MC_MedianVsMom->Write("MC_MedianVsMom");
  DATA_MedianVsMom->Write("DATA_MeidanVsMom");

  StraightGaussMean->Write("StraightGaussMean");
  StraightGaussMeanError1->Write("StraightGaussMeanError1");
  StraightGaussMeanError2->Write("StraightGaussMeanError2");
  StraightSigma->Write("StraightSigma");
  StraightSigmaError1->Write("StraightSigmaError1");
  StraightSigmaError2->Write("StraightSigmaError2");

  //Set up automated Drawing 
  
//   ComponentSigmaVsMom->SetLineColor(kGreen+1);
//   SigmaVsMom->SetLineColor(kRed);
//   RMSVsMom->SetLineColor(kBlue);
//   FWHMVsMom->SetLineColor(kMagenta+2);

//   ComponentSigmaVsMom->SetTitle("Component Sigma");
//   SigmaVsMom->SetTitle("Sigma");
//   RMSVsMom->SetTitle("RMS");
//   FWHMVsMom->SetTitle("FWHM");
 
  // ComponentSigmaVsMom->Draw("");
//   SigmaVsMom->Draw("same E1");

//   RMSVsMom->Draw("same E1");
//   FWHMVsMom->Draw("same E1");

}



double ChiSquared(TH1D * observed, TF1* flatHypothesis) {

  int bins = observed->GetNbinsX();
  
   

  double chiSq =0;
  double e = flatHypothesis->Eval(observed->GetBinCenter(1));
  std::cout << "Flat hypothesis is: "<< e << std:: endl;

  for (int i=1; i<bins+1; i++){

    std::cout << "Bin: "<< i << std:: endl;
    double o = observed->GetBinContent(i);
    std::cout << "Observed Value is: " << o << std::endl;
    double err = observed->GetBinError(i);
    std::cout << "Bin Error is: "<< err << std:: endl;

    chiSq = chiSq + pow((o-e)/err, 2);
    //chiSq = chiSq - (bins-1);
    
    std::cout << "Chi Squared is: " << chiSq << std::endl;


  }


  return chiSq;

}
