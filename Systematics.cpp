#include <TTree.h>
#include <TH1D.h>
#include <TFile.h>
#include "TSystem.h"
#include <iostream>
#include "TMath.h"
#include "TF1.h"
#include <math.h>


void Systematics() {

  Double_t ranges[] = {0, 100, 150, 200, 350, 500, 1600};
  Int_t energyBins = 6;

  TFile f("output/MC_Energy-Momentum-Fraction.root");

  TH1D *MC_SigmaVsMom;
  TH1D *MC_MeanVsMom;
  TH1D *MC_RMSVsMom;
  TH1D *MC_ComponentSigmaVsMom;
  TH1D *MC_ComponentMeanVsMom;
  MC_SigmaVsMom = (TH1D*)f.Get("SigmaVsMom")->Clone();
  MC_MeanVsMom = (TH1D*)f.Get("MeanVsMom")->Clone();
  MC_RMSVsMom = (TH1D*)f.Get("RMSVsMom")->Clone();
  MC_ComponentSigmaVsMom = (TH1D*)f.Get("ComponentSigmaVsMom")->Clone();
  MC_ComponentMeanVsMom = (TH1D*)f.Get("ComponentMeanVsMom")->Clone();

  TFile g("output/DATA_Energy-Momentum-Fraction.root");
   
  TH1D *DATA_SigmaVsMom;
  TH1D *DATA_MeanVsMom;
  TH1D *DATA_RMSVsMom;
  TH1D *DATA_ComponentSigmaVsMom;
  TH1D *DATA_ComponentMeanVsMom;
  DATA_SigmaVsMom = (TH1D*)g.Get("SigmaVsMom")->Clone();
  DATA_MeanVsMom = (TH1D*)g.Get("MeanVsMom")->Clone();
  DATA_RMSVsMom = (TH1D*)g.Get("RMSVsMom")->Clone();
  DATA_ComponentSigmaVsMom = (TH1D*)g.Get("ComponentSigmaVsMom")->Clone();
  DATA_ComponentMeanVsMom = (TH1D*)g.Get("ComponentMeanVsMom")->Clone();

  Double_t MC=0, DATA=0, MC_error=0, DATA_error=0, Error=0;



  TH1D *SigmaVsMom = new TH1D("SigmaVsMom", "MC-DATA of Sigma of Fractional Difference", energyBins, ranges);
  TH1D *MeanVsMom = new TH1D("MeanVsMom", "MC-DATA of Arithmetic Mean of Fractional Difference", energyBins, ranges);
  TH1D *RMSVsMom = new TH1D("RMSVsMom", "MC-DATA of RMS of Fractional Difference", energyBins, ranges);
  TH1D *ComponentSigmaVsMom = new TH1D("ComponentSigmaVsMom", "MC-DATA of gaussian Fit Component Sigma of Fractional Difference", energyBins, ranges);
  TH1D *ComponentMeanVsMom = new TH1D("ComponentMeanVsMom", "MC-DATA of gaussian Fit Component Arithmetic Mean of Fractional Difference", energyBins, ranges);

  for(int i=0; i<8; i++){

    MC = MC_SigmaVsMom->GetBinContent(i);
    DATA = DATA_SigmaVsMom->GetBinContent(i);

    SigmaVsMom->SetBinContent(i, DATA-MC);

    MC_error = MC_SigmaVsMom->GetBinError(i);
    DATA_error = DATA_SigmaVsMom->GetBinError(i);

    Error = sqrt(pow(MC_error, 2) + pow(DATA_error, 2));

    SigmaVsMom->SetBinError(i, Error);


  }

    for(int i=0; i<8; i++){

    MC = MC_MeanVsMom->GetBinContent(i);
    DATA = DATA_MeanVsMom->GetBinContent(i);

    MeanVsMom->SetBinContent(i, DATA-MC);

    MC_error = MC_MeanVsMom->GetBinError(i);
    DATA_error = DATA_MeanVsMom->GetBinError(i);

    Error = sqrt(pow(MC_error, 2) + pow(DATA_error, 2));

    MeanVsMom->SetBinError(i, Error);


  }


  for(int i=0; i<8; i++){

    MC = MC_RMSVsMom->GetBinContent(i);
    DATA = DATA_RMSVsMom->GetBinContent(i);

    RMSVsMom->SetBinContent(i, DATA-MC);

    MC_error = MC_RMSVsMom->GetBinError(i);
    DATA_error = DATA_RMSVsMom->GetBinError(i);

    Error = sqrt(pow(MC_error, 2) + pow(DATA_error, 2));

    RMSVsMom->SetBinError(i, Error);


  }


  for(int i=0; i<8; i++){

    MC = MC_ComponentSigmaVsMom->GetBinContent(i);
    DATA = DATA_ComponentSigmaVsMom->GetBinContent(i);

    ComponentSigmaVsMom->SetBinContent(i, DATA-MC);

    MC_error = MC_ComponentSigmaVsMom->GetBinError(i);
    DATA_error = DATA_ComponentSigmaVsMom->GetBinError(i);

    Error = sqrt(pow(MC_error, 2) + pow(DATA_error, 2));

    ComponentSigmaVsMom->SetBinError(i, Error);


  }


  for(int i=0; i<8; i++){

    MC = MC_ComponentMeanVsMom->GetBinContent(i);
    DATA = DATA_ComponentMeanVsMom->GetBinContent(i);

    ComponentMeanVsMom->SetBinContent(i, DATA-MC);

    MC_error = MC_ComponentMeanVsMom->GetBinError(i);
    DATA_error = DATA_ComponentMeanVsMom->GetBinError(i);

    Error = sqrt(pow(MC_error, 2) + pow(DATA_error, 2));

    ComponentMeanVsMom->SetBinError(i, Error);


  }


  TFile outf ("Sigma_Systematics.root", "RECREATE");

  SigmaVsMom->Write("SigmaVsMom");
  MC_SigmaVsMom->Write("MC_SigmaVsMom");
  DATA_SigmaVsMom->Write("DATA_SigmaVsMom");

  MeanVsMom->Write("MeanVsMom");
  MC_MeanVsMom->Write("MC_MeanVsMom");
  DATA_MeanVsMom->Write("DATA_MeanVsMom");
  
  RMSVsMom->Write("RMSVsMom");
  MC_RMSVsMom->Write("MC_RMSVsMom");
  DATA_RMSVsMom->Write("DATA_RMSVsMom");

  ComponentMeanVsMom->Write("ComponentMeanVsMom");
  ComponentSigmaVsMom->Write("ComponentSigmaVsMom");

}
