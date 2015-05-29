#include <TClonesArray.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
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
#include "MakeSystematic.cpp"




// Function Declarations

Double_t GetStdDev(TH1D * hist);

void DrawContour1D(TH1D* hist, TH1D * contour, TH1D * Min);

void DrawContour2D(TH2D* hist, TH1D* Xhist, TH1D* Yhist, TH2D* contour, TH2D * Min);

//Main Function 

void ChiSquareAnalyse(){


  TFile f("output/ElectronPairDataTailFitted.root");
  

  TH1D * NormChi((TH1D*)f.Get("NormChi"));
  TH1D * MeanChi((TH1D*)f.Get("MeanChi"));
  TH1D * SigmaChi((TH1D*)f.Get("SigmaChi"));
  TH1D * TailChi((TH1D*)f.Get("TailChi"));

  TH2D * NormTailChi((TH2D*)f.Get("NormTailChi"));
  TH2D * MeanSigmaChi((TH2D*)f.Get("MeanSigmaChi"));


  //Energy Regime1 

  TH1D * MeanChi1((TH1D*)f.Get("MeanChi1"));
  TH1D * SigmaChi1((TH1D*)f.Get("SigmaChi1"));


  //Energy Regime2

  TH1D * MeanChi2((TH1D*)f.Get("MeanChi2"));
  TH1D * SigmaChi2((TH1D*)f.Get("SigmaChi2"));


  //Calculate Errors to be passed to MakeSystematic()

  Double_t NormErr = GetStdDev(NormChi);
  Double_t MeanErr = GetStdDev(MeanChi);
  Double_t SigmaErr = GetStdDev(SigmaChi);
  Double_t TailErr = GetStdDev(TailChi);

  Double_t MeanErr1 = GetStdDev(MeanChi1);
  Double_t SigmaErr1 = GetStdDev(SigmaChi1);
  
  Double_t MeanErr2 = GetStdDev(MeanChi2);
  Double_t SigmaErr2 = GetStdDev(SigmaChi2);




  Int_t Bins = NormChi->GetNbinsX();



  //Define Histograms for Contours and Minima
  
  
  TH1D * NormContour = new TH1D("NormContour", "NormContour", Bins, NormChi->GetBinLowEdge(1), NormChi->GetBinLowEdge(Bins)+NormChi->GetBinWidth(Bins));
  TH1D * NormMin = new TH1D("NormMin", "NormMin", Bins, NormChi->GetBinLowEdge(1), NormChi->GetBinLowEdge(Bins)+NormChi->GetBinWidth(Bins));
  
  DrawContour1D(NormChi, NormContour, NormMin);

  TH1D * MeanContour = new TH1D("MeanContour", "MeanContour", Bins, MeanChi->GetBinLowEdge(1), MeanChi->GetBinLowEdge(Bins)+MeanChi->GetBinWidth(Bins));
  TH1D * MeanMin = new TH1D("MeanMin", "MeanMin", Bins, MeanChi->GetBinLowEdge(1), MeanChi->GetBinLowEdge(Bins)+MeanChi->GetBinWidth(Bins));

  DrawContour1D(MeanChi, MeanContour, MeanMin);

  TH1D * SigmaContour = new TH1D("SigmaContour", "SigmaContour", Bins, SigmaChi->GetBinLowEdge(1), SigmaChi->GetBinLowEdge(Bins)+SigmaChi->GetBinWidth(Bins));
  TH1D * SigmaMin = new TH1D("SigmaMin", "SigmaMin", Bins, SigmaChi->GetBinLowEdge(1), SigmaChi->GetBinLowEdge(Bins)+SigmaChi->GetBinWidth(Bins));

  DrawContour1D(SigmaChi, SigmaContour, SigmaMin);

  TH1D * TailContour = new TH1D("TailContour", "TailContour", Bins, TailChi->GetBinLowEdge(1), TailChi->GetBinLowEdge(Bins)+TailChi->GetBinWidth(Bins));
  TH1D * TailMin = new TH1D("TailMin", "TailMin", Bins, TailChi->GetBinLowEdge(1), TailChi->GetBinLowEdge(Bins)+TailChi->GetBinWidth(Bins));
  
  DrawContour1D(TailChi, TailContour, TailMin);




  TH2D * NormTailContour = new TH2D("NormTailContour","NormTailContour",Bins, NormChi->GetBinLowEdge(1), NormChi->GetBinLowEdge(Bins)+NormChi->GetBinWidth(Bins), Bins, TailChi->GetBinLowEdge(1), TailChi->GetBinLowEdge(Bins)+TailChi->GetBinWidth(Bins));

  TH2D * NormTailMin = new TH2D("NormTailMin","NormTailMin",Bins, NormChi->GetBinLowEdge(1), NormChi->GetBinLowEdge(Bins)+NormChi->GetBinWidth(Bins), Bins, TailChi->GetBinLowEdge(1), TailChi->GetBinLowEdge(Bins)+TailChi->GetBinWidth(Bins));

  DrawContour2D(NormTailChi, NormChi, TailChi, NormTailContour, NormTailMin);




  TH2D * MeanSigmaContour = new TH2D("MeanSigmaContour","MeanSigmaContour",Bins, MeanChi->GetBinLowEdge(1), MeanChi->GetBinLowEdge(Bins)+MeanChi->GetBinWidth(Bins), Bins, SigmaChi->GetBinLowEdge(1), SigmaChi->GetBinLowEdge(Bins)+SigmaChi->GetBinWidth(Bins));

  TH2D * MeanSigmaMin = new TH2D("MeanSigmaMin","MeanSigmaMin",Bins, MeanChi->GetBinLowEdge(1), MeanChi->GetBinLowEdge(Bins)+MeanChi->GetBinWidth(Bins), Bins, SigmaChi->GetBinLowEdge(1), SigmaChi->GetBinLowEdge(Bins)+SigmaChi->GetBinWidth(Bins));


  DrawContour2D(MeanSigmaChi, MeanChi, SigmaChi, MeanSigmaContour, MeanSigmaMin);



  std::cout << "PARAMETERS OF COMBINED FIT:" << std::endl;
  std::cout <<"____________________________________________________" << std::endl;
  std::cout << "GAUSS NORM: "<< std::endl ;
  GetStdDev(NormChi);
  std::cout <<"____________________________________________________" << std::endl;
  std::cout << "GAUSS MEAN: "<<std::endl;
  GetStdDev(MeanChi);
  std::cout <<"____________________________________________________" << std::endl;
  std::cout << "GAUSS SIGMA: "<<std::endl;
  GetStdDev(SigmaChi);
  std::cout <<"____________________________________________________" << std::endl;
  std::cout << "TAIL NORM: "<<std::endl;
  GetStdDev(TailChi);


  MakeSystematic(MeanErr, MeanErr1, MeanErr2, SigmaErr, SigmaErr1, SigmaErr2);
  
  



  TFile outf("output/ChiSquaredAnalysis.root", "RECREATE");

  NormChi->Write("NormChi");
  NormContour->Write("NormContour");
  NormMin->Write("NormMin");
  MeanChi->Write("MeanChi");
  MeanContour->Write("MeanContour");
  MeanMin->Write("MeanMin");
  SigmaChi->Write("SigmaChi");
  SigmaContour->Write("SigmaContour");
  SigmaMin->Write("SigmaMin");
  TailChi->Write("TailChi");
  TailContour->Write("TailContour");
  TailMin->Write("TailMin");
  NormTailChi->Write("NormTailChi");
  NormTailContour->Write("NormTailContour");
  NormTailMin->Write("NormTailMin");
  MeanSigmaChi->Write("MeanSigmaChi");
  MeanSigmaContour->Write("MeanSigmaContour");
  MeanSigmaMin->Write("MeanSigmaMin");
  



  



}




// Functions


Double_t GetStdDev(TH1D * hist){

  Double_t Min = hist->GetMinimum();
  Int_t MinBin = hist->GetMinimumBin();
  
  hist->Scale(-1);
  Int_t LowBin = hist->FindFirstBinAbove(-(Min+1));
  Int_t HighBin = hist->FindLastBinAbove(-(Min+1));

  Double_t StdDev1 = hist->GetBinCenter(HighBin) - hist->GetBinCenter(MinBin);
  Double_t StdDev2 = hist->GetBinCenter(MinBin) - hist->GetBinCenter(LowBin);

  Double_t StdDev;

  if(StdDev1 > StdDev2){
     StdDev = StdDev1/2;}
  else{StdDev = StdDev2/2;}

  //std::cout << "Standard Deviation of " << hist->GetName() << " is " << StdDev << std::endl;  

  std::cout <<  hist->GetBinCenter(MinBin) <<" +/- " << StdDev << std::endl;

  hist->Scale(-1);
  
  return StdDev;


}


void DrawContour1D(TH1D* hist, TH1D * contour, TH1D * Min){


  
  Int_t MinBin = hist->GetMinimumBin();
  Double_t Minimum = hist->GetMinimum();
  
  Min->SetBinContent(MinBin, Minimum);

  hist->Scale(-1);
  Int_t LowBin = hist->FindFirstBinAbove(-(Minimum+1));
  Int_t HighBin = hist->FindLastBinAbove(-(Minimum+1));
  
  for(int i=LowBin; i < HighBin+1; i++){
      contour->SetBinContent(i, Minimum+1);
    }

  hist->Scale(-1);
     
	
  
}

void DrawContour2D(TH2D* hist, TH1D* Xhist, TH1D* Yhist, TH2D* contour, TH2D * Min){

  Double_t Minimum = hist->GetMinimum();
  Int_t xMinBin= Xhist->GetMinimumBin();
  Int_t yMinBin = Yhist->GetMinimumBin();

  Min->SetBinContent(xMinBin, yMinBin, 1);
  
  Int_t Bins = hist->GetNbinsX();

  for(int i=0; i<Bins; i++){

    for(int j=0; j<Bins; j++){

      Double_t ChiSq = hist->GetBinContent(j, i);
      
      if(Minimum+0.8 < ChiSq && ChiSq < Minimum+1.2){

	contour->SetBinContent(j, i, Minimum+1);
       
      }

    }
    
  }

}
