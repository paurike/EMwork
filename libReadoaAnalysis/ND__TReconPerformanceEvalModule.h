//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Feb 17 15:49:59 2015 by ROOT version 5.34/09)
//      from the StreamerInfo in file /storage/epp2/phskaj/pidSkims/prod6B/electronPair/fullSkim.electronPair.production006_B_mcp_neut_2010-02-water_magnet_run1.1of11.root
//////////////////////////////////////////////////////////


#ifndef ND__TReconPerformanceEvalModule_h
#define ND__TReconPerformanceEvalModule_h
namespace ND {
class TReconPerformanceEvalModule;
} // end of namespace.

#include "ND__TAnalysisReconModuleBase.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include <vector>
#include <utility>
#include "TObject.h"
#include <string>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "ND__TReconPerformanceEvalModule.h"

namespace ND {
class TReconPerformanceEvalModule : public ND::TAnalysisReconModuleBase {

public:
// Nested classes forward declaration.
class TGlobalReconObject;
class TGlobalTruthInfo;

public:
// Nested classes declaration.
class TGlobalTruthInfo : public TObject {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        SetOK;       //
   double      Momentum;    /// < Momentum of the true trajectory
   TVector3    Direction;    /// < Direction of the true trajectory
   double      Charge;       /// < Charge of the trajectory
   TLorentzVector Position;     /// < Initial position of the trajectory
   double         Efficiency;    /// < Efficiency of this truth matching
   double         Purity;        /// < Purity of this truth matching

   TGlobalTruthInfo();
   TGlobalTruthInfo(const TGlobalTruthInfo & );
   virtual ~TGlobalTruthInfo();

   ClassDef(TGlobalTruthInfo,3); // Generated by MakeProject.
};
class TGlobalReconObject : public TObject {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        SetOK;       //
   Int_t       NConstituents;    //
   map<string,int>* NModuleConstituents;    //
   string           SubdetectorString;      //
   string           StatusString;           //
   TClonesArray*    GlobalNodes;            //
   Int_t            NGlobalNodes;           //
   Int_t            NGlobalNodesSaved;      //
   double           Momentum;               //
   double           MomentumByRange;        //
   double           MomentumByRangeMuon;    //
   double           MomentumByRangeMuonFlip;    //
   double           MomentumByRangeElectron;    //
   double           MomentumByRangeElectronFlip;    //
   double           MomentumByRangeProton;          //
   double           MomentumByRangeProtonFlip;      //
   TLorentzVector   Position;                       //
   TVector3         Direction;                      //
   double           Charge;                         //
   double           Quality;                        //
   int              NDOF;                           //
   string           ParticleID;                     //
   double           PIDWeight;                      //
   TClonesArray*    MatchingChi2Info;               //
   int              NMatchingChi2Info;              //
   ND::TReconPerformanceEvalModule::TGlobalTruthInfo Truth;                          //
   TClonesArray*                                     Constituents;                   //
   TClonesArray*                                     DownToTrackerConstituents;      //
   Int_t                                             NDownToTrackerConstituents;     //
   UInt_t                                            UniqueID;                       //

   TGlobalReconObject();
   TGlobalReconObject(const TGlobalReconObject & );
   virtual ~TGlobalReconObject();

   ClassDef(TGlobalReconObject,3); // Generated by MakeProject.
};

public:
// Data Members.
   Int_t       fNGlobalReconObjects;    //
   TClonesArray* fGlobalReconObjects;     //
   vector<pair<string,int> > fHitInfo;                //
   bool                      fSaveAllGlobalNodes;     //

   TReconPerformanceEvalModule();
   TReconPerformanceEvalModule(const TReconPerformanceEvalModule & );
   virtual ~TReconPerformanceEvalModule();

   ClassDef(TReconPerformanceEvalModule,3); // Generated by MakeProject.
};
} // namespace
#endif
