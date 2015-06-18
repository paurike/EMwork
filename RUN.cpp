#include "TailReader.cpp"
#include "DataTailFit.cpp"
#include "ChiSquareAnalyse.cpp"


void RUN(){

  TailReader("output/ntuple.magnetMC_electronPair.root", 0);
  
  TailReader("output/ntuple.DATA_electronPair.root", 0);
  
  DataTailFit(31);
  
  ChiSquareAnalyse();
  
}




