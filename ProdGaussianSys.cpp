#include "TailReader.cpp"
#include "Systematics.cpp"


void ProdGaussianSys() {


  TailReader("output/ntuple.magnetMC_electronPair.root", 0);
  TailReader("output/ntuple.DATA_electronPair.root", 0);

  Systematics();

}
