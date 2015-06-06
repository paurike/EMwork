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


//function declarations


Double_t ChiSquared(const TH1D* Data, const TH1D * Model);

void  SetModel(Double_t norm, Double_t mean, Double_t sigma, Double_t tailNorm, const TH1D* tail, TH1D *model);

void Fitter(TH1D* Data,TH1D *MC, TH1D* CleanEvents, TH1D* TailTemplate, Double_t * Minimum, TH1D* Result,Int_t steps, TH1D * NormChi, TH1D* MeanChi, TH1D* SigmaChi, TH1D* TailChi, TH2D* NormTailChi, TH2D * MeanSigmaChi);





//Main Function
void DataTailFit(Int_t steps){
  
  TFile f("output/ElectronPairMC_Energy-Momentum-Fraction.root");
  
  TH1D* TailTemplate;
  TailTemplate = (TH1D*)f.Get("ReverseCutTailTemplate")->Clone();
  
  TH1D* CleanEvents;
  CleanEvents = (TH1D*)f.Get("SimpleCleanReconFractionDist")->Clone();

  Double_t HighLowRegimes[] = {0, 200, 1200};
  Int_t Regimes = 2;
  TH1D* RegimeHistsClean[Regimes];
  TH1D* RegimeHistsTail[Regimes];
  TH1D* RegimeHistsMC[Regimes];

  //Careful here when changing number of regimes
    RegimeHistsClean[0] = (TH1D*)f.Get("Regime1CleanSample")->Clone();
    RegimeHistsClean[1] = (TH1D*)f.Get("Regime2CleanSample")->Clone();
    RegimeHistsTail[0] = (TH1D*)f.Get("Regime1Tail")->Clone();
    RegimeHistsTail[1] = (TH1D*)f.Get("Regime2Tail")->Clone();
    RegimeHistsMC[0] = (TH1D*)f.Get("Regime1Total")->Clone();
    RegimeHistsMC[1] = (TH1D*)f.Get("Regime2Total")->Clone();

  
  TH1D* MC;
  MC = (TH1D*)f.Get("ReconFractionDist")->Clone();
  
  TFile g("output/ElectronPairDATA_Energy-Momentum-Fraction.root");
  
  TH1D* Data;
  Data = (TH1D*)g.Get("ReconFractionDist")->Clone();
  
  TH1D* RegimeHistsData[Regimes];

  RegimeHistsData[0] = (TH1D*)g.Get("Regime1Total")->Clone();
  RegimeHistsData[1] = (TH1D*)g.Get("Regime2Total")->Clone();
  
  //Int_t steps = 10;
 
 
  //fitting total distribution:

  Int_t Bins = Data->GetNbinsX();
  Double_t Minimum[] ={10000, 0, 0, 0, 0};
  TH1D *Model = new TH1D("model", "model", Bins, -2, 10);

  TH1D* NormChi = new TH1D("NormChi", "NormChi", steps, 0, 1);
  TH1D* MeanChi = new TH1D("MeanChi", "MeanChi", steps, 0, 1);
  TH1D* SigmaChi= new TH1D("SigmaChi", "SigmaChi", steps, 0,1);
  TH1D* TailChi = new TH1D("TailChi", "TailChi", steps, 0, 1);

  TH2D *NormTailChi = new TH2D("NormTailChi", "NormTailChi", steps, 0,1, steps, 0,1);
  TH2D *MeanSigmaChi = new TH2D("MeanSigmaChi", "MeanSigmaChi", steps, 0, 1, steps, 0, 1);
  
  std::cout << "everything up to the fitter has been defined" << std::endl;


  Fitter(Data, MC, CleanEvents,TailTemplate, Minimum, Model, steps, NormChi, MeanChi, SigmaChi, TailChi, NormTailChi, MeanSigmaChi);
   
  TF1 *Gaussian = new TF1("Gaussian", "[0]*TMath::Gaus(x,[1],[2])", -2, 10);
  Gaussian->SetParameters(Minimum[1], Minimum[2], Minimum[3]);

  SetModel(Minimum[1], Minimum[2], Minimum[3], Minimum[4], TailTemplate, Model);
 
  TailTemplate->Scale(Minimum[4]);


  //fitting Energy Regimes:


//   //Regime 1

  TH1D* NormChi1 = new TH1D("NormChi1", "NormChi1", steps, 0, 1);
  TH1D* MeanChi1 = new TH1D("MeanChi1", "MeanChi1", steps, 0, 1);
  TH1D* SigmaChi1= new TH1D("SigmaChi1", "SigmaChi1", steps, 0,1);
  TH1D* TailChi1 = new TH1D("TailChi1", "TailChi1", steps, 0, 1);

  TH2D *NormTailChi1 = new TH2D("NormTailChi1", "NormTailChi1", steps, 0,1, steps, 0,1);
  TH2D *MeanSigmaChi1 = new TH2D("MeanSigmaChi1", "MeanSigmaChi1", steps, 0, 1, steps, 0, 1);


  Double_t Regime1Minimum[] ={10000, 0, 0, 0, 0};
  TH1D *Regime1Model = new TH1D("Regime1model", "Regime1model", Bins, -2, 10);
 
  Fitter(RegimeHistsData[0],RegimeHistsMC[0], RegimeHistsClean[0], RegimeHistsTail[0], Regime1Minimum, Regime1Model,steps, NormChi1, MeanChi1, SigmaChi1, TailChi1, NormTailChi1, MeanSigmaChi1);

  TF1 *Gaussian1 = new TF1("Gaussian1", "[0]*TMath::Gaus(x,[1],[2])", -2, 10);
  Gaussian1->SetParameters(Regime1Minimum[1], Regime1Minimum[2], Regime1Minimum[3]);

  SetModel(Regime1Minimum[1], Regime1Minimum[2], Regime1Minimum[3], Regime1Minimum[4], RegimeHistsTail[0], Regime1Model);
 
  RegimeHistsTail[0]->Scale(Minimum[4]);



 //Regime 2


  TH1D* NormChi2 = new TH1D("NormChi2", "NormChi2", steps, 0, 1);
  TH1D* MeanChi2 = new TH1D("MeanChi2", "MeanChi2", steps, 0, 1);
  TH1D* SigmaChi2= new TH1D("SigmaChi2", "SigmaChi2", steps, 0,1);
  TH1D* TailChi2 = new TH1D("TailChi2", "TailChi2", steps, 0, 1);

  TH2D *NormTailChi2 = new TH2D("NormTailChi2", "NormTailChi2", steps, 0,1, steps, 0,1);
  TH2D *MeanSigmaChi2 = new TH2D("MeanSigmaChi2", "MeanSigmaChi2", steps, 0, 1, steps, 0, 1);
  
  Double_t Regime2Minimum[] ={10000, 0, 0, 0, 0};
  TH1D *Regime2Model = new TH1D("Regime2model", "Regime2model", Bins, -2, 10);
  
  Fitter(RegimeHistsData[1],RegimeHistsMC[1], RegimeHistsClean[1], RegimeHistsTail[1], Regime2Minimum, Regime2Model, steps, NormChi2, MeanChi2, SigmaChi2, TailChi2, NormTailChi2, MeanSigmaChi2);
  
  TF1 *Gaussian2 = new TF1("Gaussian2", "[0]*TMath::Gaus(x,[1],[2])", -2, 10);
  Gaussian2->SetParameters(Regime2Minimum[1], Regime2Minimum[2], Regime2Minimum[3]);

  SetModel(Regime2Minimum[1], Regime2Minimum[2], Regime2Minimum[3], Regime2Minimum[4], RegimeHistsTail[1], Regime2Model);
 
  RegimeHistsTail[1]->Scale(Regime2Minimum[4]);



  TFile* outf = new TFile("output/ElectronPairDataTailFitted.root", "RECREATE");

  Data->Write("Data");
  Gaussian->Write("Gaussian");
  TailTemplate->Write("TailTemplate");
  Model->Write("FittedModel");

  NormChi->Write("NormChi");
  MeanChi->Write("MeanChi");
  SigmaChi->Write("SigmaChi");
  TailChi->Write("TailChi");

  NormTailChi->Write("NormTailChi");
  MeanSigmaChi->Write("MeanSigmaChi");

  

  RegimeHistsData[0]->Write("Regime1Data");
  Gaussian1->Write("Gaussian1");
  RegimeHistsTail[0]->Write("Regime1Tail");
  Regime1Model->Write("Regime1Model");

  NormChi1->Write("NormChi1");
  MeanChi1->Write("MeanChi1");
  SigmaChi1->Write("SigmaChi1");
  TailChi1->Write("TailChi1");

  NormTailChi1->Write("NormTailChi1");
  MeanSigmaChi1->Write("MeanSigmaChi1");

  RegimeHistsData[1]->Write("Regime2Data");
  Gaussian2->Write("Gaussian2");
  RegimeHistsTail[1]->Write("Regime2Tail");
  Regime2Model->Write("Regime2Model");

  NormChi2->Write("NormChi2");
  MeanChi2->Write("MeanChi2");
  SigmaChi2->Write("SigmaChi2");
  TailChi2->Write("TailChi2");

  NormTailChi2->Write("NormTailChi2");
  MeanSigmaChi2->Write("MeanSigmaChi2");
  
  
  outf->Close();




  }



//Functions
		     




Double_t ChiSquared(const TH1D* Data, const TH1D* Model){


  Int_t bins = Data->GetNbinsX();

  if(bins != Model->GetNbinsX())
    {std::cout << "CAREFUL! DATA AND MODEL DON'T HAVE THE SAME NUMBER OF BINS" << std::endl; return 1;}
 
  Double_t val=0;
  
  for(int i=0; i<bins; i++){
    
    Double_t DATA = Data->GetBinContent(i+1);
    Double_t MC = Model->GetBinContent(i+1);   
    Double_t Diff = DATA-MC;

    val = val + (Diff*Diff)/MC;
//     std::cout << "Data is: " <<DATA << std::endl;
//     std::cout << "Model is: " << MC << std::endl;

  }

  val = val-bins;

  //std::cout << "CHI SQUARED VALUE IS: " << val << std::endl;

  return val;

}






void SetModel(Double_t norm, Double_t mean, Double_t sigma, Double_t tailNorm, const TH1D * tail, TH1D *model){


 TF1 Gaussian("Gaussian", "[0]*TMath::Gaus(x,[1],[2])", -2, 10);
 Gaussian.SetParameters(norm, mean, sigma);


 Int_t bins = tail->GetNbinsX();

 for(int i=0; i<bins; i++){

   Double_t x = tail->GetBinCenter(i+1);

   Double_t GausVal = Gaussian.Eval(x);
   Double_t TailVal = tailNorm * tail->GetBinContent(i+1);
   
   model->SetBinContent(i+1, GausVal+TailVal);
   
 }
 
}

void Fitter(TH1D* Data,TH1D *MC, TH1D* CleanEvents, TH1D* TailTemplate, Double_t * Minimum,  TH1D* Result, Int_t steps, TH1D * NormChi, TH1D* MeanChi, TH1D* SigmaChi, TH1D* TailChi, TH2D* NormTailChi, TH2D * MeanSigmaChi){

  std::cout << "marker1" << std::endl;

  Double_t ScaleFactor = Data->GetEntries()/MC->GetEntries();
  
  Double_t normStep=((CleanEvents->GetMaximum())*ScaleFactor*0.8)/steps;
  Double_t meanStep=abs((CleanEvents->GetMean())*40/steps);
  Double_t sigmaStep=(CleanEvents->GetRMS())*0.6/steps;
  Double_t tailStep=ScaleFactor*1.5/steps;
  
  Double_t normInit=(CleanEvents->GetMaximum())*0.7*ScaleFactor;
  Double_t meanInit=(CleanEvents->GetMean())*20;
  Double_t sigmaInit=(CleanEvents->GetRMS())*0.5;
  Double_t tailInit=0.5*ScaleFactor;
  
  Int_t Bins = Data->GetNbinsX();

  std::cout << "marker2" << std::endl;


  //Set ChiSquared Histogram range

  NormChi->SetBins(steps, normInit, normInit + steps*normStep);
  MeanChi->SetBins(steps, meanInit, meanInit + steps*meanStep );
  SigmaChi->SetBins(steps, sigmaInit, sigmaInit + steps*sigmaStep);
  TailChi->SetBins(steps, tailInit, tailInit + steps*tailStep);

  std::cout << "marker3" << std::endl;
  
  NormTailChi->SetBins(steps, normInit,normInit + steps*normStep, steps, tailInit,tailInit + steps*tailStep); 
  MeanSigmaChi->SetBins(steps,meanInit,meanInit + steps*meanStep, steps, sigmaInit,sigmaInit + steps*sigmaStep); 
 
  std::cout << "marker4" << std::endl;

  
  //Fill the ChiSquared Array:
  
  
  std::cout << "This is where things go wrong" << std::endl;

  Double_t ChiSquaredMap[steps][steps][steps][steps];

  Int_t MinimumBins[5] = {0,0,0,0,0};
  
  
  std::cout << "Array has been defined" << std::endl;
  
 


  TH1D *Model = new TH1D("model", "model", Bins, -2, 10);
  
  for (int norm=0; norm < steps; norm++){
    
    Double_t NormVal = normInit+normStep*norm;
    //std::cout << "-----------------------------------------------for norm: " << NormVal << std::endl;
    
    for(int mean=0; mean<steps; mean++) {

      std::cout << "Progress: " << 100.00/steps*norm + 100.00/(steps*steps) * mean <<"%"  <<std::endl;
      
      Double_t MeanVal = meanInit+meanStep*mean;
      //std::cout<< "----------------------------------------------for mean: " << MeanVal << std::endl;
      
      for(int sigma=0; sigma<steps; sigma++){
	
	Double_t SigmaVal = sigmaInit+sigmaStep*sigma;
	//std::cout << "----------------------------------------------for sigma: " << SigmaVal << std::endl;
	
	for(int tail=0; tail<steps; tail++){	 
	  
	  Double_t TailVal = tailInit+tailStep*tail;
	  //std::cout <<"for Tail norm: " << TailVal << std::endl; 
	  
	  SetModel(NormVal, MeanVal, SigmaVal, TailVal, TailTemplate, Model); 
	  
	  Double_t ChiSq = ChiSquared(Data, Model);

	  ChiSquaredMap[norm][mean][sigma][tail] = ChiSq;
	  
	  //std::cout << ChiSquaredMap[norm][mean][sigma][tail] << std::endl;
	  
	  if(ChiSq < Minimum[0]){
	    Minimum[0] = ChiSq;
	    Minimum[1] = NormVal;
	    Minimum[2] = MeanVal;
	    Minimum[3] = SigmaVal;
	    Minimum[4] = TailVal;

	    MinimumBins[1] = norm;
	    MinimumBins[2] = mean;
	    MinimumBins[3] = sigma;
	    MinimumBins[4] = tail;
	    
	    
	  }
	  
	}
	
      }
      
    }
    
  }


  //Fill ChiSquared Histograms


  for (int i=0; i<steps; i++){

    NormChi->SetBinContent(i+1,ChiSquaredMap[i][MinimumBins[2]][MinimumBins[3]][MinimumBins[4]]);
    MeanChi->SetBinContent(i+1,ChiSquaredMap[MinimumBins[1]][i][MinimumBins[3]][MinimumBins[4]]);
    SigmaChi->SetBinContent(i+1,ChiSquaredMap[MinimumBins[1]][MinimumBins[2]][i][MinimumBins[4]]);
    TailChi->SetBinContent(i+1,ChiSquaredMap[MinimumBins[1]][MinimumBins[2]][MinimumBins[3]][i]);


    for(int j=0; j<steps; j++){

      NormTailChi->SetBinContent(i+1, j+1, ChiSquaredMap[i][MinimumBins[2]][MinimumBins[3]][j]);
      MeanSigmaChi->SetBinContent(i+1, j+1, ChiSquaredMap[MinimumBins[1]][i][j][MinimumBins[4]]);

    }

  }




  std::cout << "Progress: " <<  "100%"  <<std::endl;


  std::cout << "Array has been filled" << std::endl;
    
  std::cout << "Minimum was found as: "<< std::endl;

  std::cout << "___________________________________" << std::endl;

  std::cout << "CHI SQUARED: " << Minimum[0] << " Initial Value was: "<< 10000  << std::endl;

  std::cout << "___________________________________" << std::endl;

  std::cout << "GAUSS NORM: " << Minimum[1] << " Initial Value was: "<< normInit  <<std::endl;

  std::cout << "___________________________________" << std::endl;

  std::cout << "GAUSS MEAN: " << Minimum[2] <<" Initial Value was: "<< meanInit  << std::endl;

  std::cout << "___________________________________" << std::endl;

  std::cout << "GAUSS SIGMA: " << Minimum[3] << " Initial Value was: "<< sigmaInit <<std::endl;

  std::cout << "___________________________________" << std::endl;

  std::cout << "TAIL NORM: " << Minimum[4] <<" Initial Value was: "<< tailInit  << std::endl;

  if(Minimum[0] == 10000){
    
    std::cout << "FIT HAS FAILED" << std::endl;
    
  }

  SetModel(Minimum[1], Minimum[2], Minimum[3], Minimum[4], TailTemplate, Model);
 

  for( int i=0; i<Bins; i++){

    Result->SetBinContent(i+1, Model->GetBinContent(i+1));

  }




}






 




     
