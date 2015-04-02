//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Feb 17 15:49:59 2015 by ROOT version 5.34/09)
//      from the StreamerInfo in file /storage/epp2/phskaj/pidSkims/prod6B/electronPair/fullSkim.electronPair.production006_B_mcp_neut_2010-02-water_magnet_run1.1of11.root
//////////////////////////////////////////////////////////


#ifndef ND__TTrackerReconModule_h
#define ND__TTrackerReconModule_h
namespace ND {
class TTrackerReconModule;
} // end of namespace.

#include "ND__TAnalysisReconModuleBase.h"
#include "TClonesArray.h"
#include "TObject.h"
#include "Riostream.h"
#include <string>
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "ND__TTrackerReconModule.h"
#include "ND__TTrueVertex.h"

namespace ND {
class TTrackerReconModule : public ND::TAnalysisReconModuleBase {

public:
// Nested classes forward declaration.
class TTrackerVertex;
class TTrackerResult;
class TTPCTrack;
class TTrackOther;
class TUnusedHit;
class TTrueParticle;
class TTrackerConstituent;
class TFGDTrack;
class TTrackerNode;

public:
// Nested classes declaration.
class TTrackerNode : public TObject {

public:
// Nested classes declaration.

public:
// Data Members.
   double      Charge;      ///< The Charge (+-1)
   double      EDeposit;    ///< The Energy Deposit (number of pe's)
   TLorentzVector Position;    ///< Position 4-vector (at node) x,y,z,t in mm, ns
   TLorentzVector Variance;    ///< Position variance 4-vector (at node) var(x),var(y),var(z),var(t) in mm^2, ns^2
   TVector3       Direction;    ///< Direction vector (at node) 
   TVector3       DirectionVariance;    ///< Direction variance vector (at node) 
   double         Momentum;             ///< Track Momentum (at node) in MeV/c 
   double         MomentumError;        ///< Track Momentum uncertainty (at node) in MeV/c

   TTrackerNode();
   TTrackerNode(const TTrackerNode & );
   virtual ~TTrackerNode();

   ClassDef(TTrackerNode,3); // Generated by MakeProject.
};
class TFGDTrack : public TObject {

public:
// Nested classes declaration.

public:
// Data Members.
   UInt_t      UniqueID;    ///< Unique ID number to allow matching to Global Recon object.
   int         Detector;    ///< FGD number 1 or 2
   int         Ndof;        ///< Number of degrees of freedom in FGD track fit
   double      Chi2;        ///< Chi2 of the FGD track fit
   TLorentzVector Position;    //
   TVector3       Direction;    //
   double         EDeposit;     //

   TFGDTrack();
   TFGDTrack(const TFGDTrack & );
   virtual ~TFGDTrack();

   ClassDef(TFGDTrack,3); // Generated by MakeProject.
};
class TTrackerConstituent : public TObject {

public:
// Nested classes declaration.

public:
// Data Members.
   string      AlgorithmName;    ///< algorithm that created this object.
   int         Detectors;        ///< Detectors used
   unsigned long Status;           ///< The status for the fit.
   double        Quality;          ///< The quality of the fit.(probability)
   int           NDOF;             ///< The number of degrees of freedom.
   double        Chi2;             ///< The chi2 of the fit.
   int           NNodes;           ///< The number of nodes
   int           NHits;            ///< The number of hits.
   Int_t         NConstituents;    ///< The number of constituents this constituent is made of
   Int_t         ConstitIdx[2];    ///< Index into Constituents in TTrackerResult::Constituents of this constituent's constituents
   bool          isForward;        ///< Sense of track
   double        Charge;           ///< The Charge of this constituent (+-1)
   double        EDeposit;         ///< The deposited charge for the constituent object (number of pe's).
   TLorentzVector FrontPosition;    ///< The 4-vector position at the front of the track (x,y,z,t) in mm, ns 
   TLorentzVector BackPosition;     ///< The 4-vector position of the back of the track (x,y,z,t) in mm, ns 
   TLorentzVector FrontVariance;    ///< The 4-vector position variance at the front of the track (var(x),var(y),var(z),var(t)) in mm^2, ns^2 
   TLorentzVector BackVariance;     ///< The 4-vector position variance at the back of the track (var(x),var(y),var(z),var(t)) in mm^2, ns^2 
   TVector3       FrontDirection;    ///< The direction vector at the front of the track 
   TVector3       BackDirection;     ///< The direction vector at the back of the track
   double         FrontMomentum;     ///< the momentum at the front of the track in MeV/c 
   double         BackMomentum;      ///< the momentum at the back of the track iin MeV/c
   TLorentzVector Position;          ///< position 4-vector (x,y,z,t) in mm, ns
   TLorentzVector Variance;          ///< position variance 4-vector  (var(x),var(y),var(z),var(t)) in mm^2, ns^2
   TVector3       Direction;         ///< direction vector 
   TVector3       DirectionVariance;    ///< direction variance vector
   double         Momentum;             ///< momentum MeV/c
   double         MomentumError;        ///< uncertainty in momentum MeV/c

   TTrackerConstituent();
   TTrackerConstituent(const TTrackerConstituent & );
   virtual ~TTrackerConstituent();

   ClassDef(TTrackerConstituent,3); // Generated by MakeProject.
};
class TTrueParticle : public TObject {

public:
// Nested classes declaration.

public:
// Data Members.
   int         ID;          ///< Trajectoy  Id
   double      Pur;         ///< The purity for matching
   double      Eff;         ///< The efficiency for matching
   ND::TTrueVertex Vertex;      ///< True vertex associated to this TrueParticle

   TTrueParticle();
   TTrueParticle(const TTrueParticle & );
   virtual ~TTrueParticle();

   ClassDef(TTrueParticle,3); // Generated by MakeProject.
};
class TUnusedHit : public TObject {

public:
// Nested classes declaration.

public:
// Data Members.
   double      TotalCharge;    ///< Deposited charge (the hit EDeposit)
   TVector3    Position;       ///< The position of the hit component 0=x 1=y 2=z in mm
   TVector3    Variance;       ///< The position variance in mm
   double      Time;           ///< Time of the hit in ns
   double      TimeUnc;        ///< Time Uncertainty of hit in ns

   TUnusedHit();
   TUnusedHit(const TUnusedHit & );
   virtual ~TUnusedHit();

   ClassDef(TUnusedHit,3); // Generated by MakeProject.
};
class TTrackOther : public TObject {

public:
// Nested classes declaration.

public:
// Data Members.
   string      AlgorithmName;    ///< The name of the algorithm that created this object.
   int         Detector;         ///< Detector used (1,2,3 for TPC, or 4,5 for FGD?)
   int         NHits;            ///< The number of hits.
   TClonesArray* Hits;             ///< The hits
   double        EDeposit;         ///< The deposited charge for the object.
   TLorentzVector FrontPosition;    ///< The position of the track at its upstream-most end (x,y,z,t) in mm, ns
   TLorentzVector BackPosition;     ///< The position of the track at its downstream-most end    /// (x,y,z,t) in mm, ns
   ND::TTrackerReconModule::TUnusedHit hackHits;         ///<This is just here to fool TFile::MakeProject, not a real object. 

   TTrackOther();
   TTrackOther(const TTrackOther & );
   virtual ~TTrackOther();

   ClassDef(TTrackOther,3); // Generated by MakeProject.
};
class TTPCTrack : public TObject {

public:
// Nested classes declaration.

public:
// Data Members.
   UInt_t      UniqueID;    ///< Unique ID number to allow matching to Global Recon object.
   int         Detector;    ///< TPC number 1, 2 or 3
   int         Ndof;        ///< Number of degrees of freedom in TPC fit
   double      Chi2;        ///< TPC chi2 calculated after likelihood fit
   int         NHits;       ///< number of clusters used in TPC fit
   double      Momentum;    ///< Momentum of the TPC pid in MeV/c
   double      MomentumError;    ///< Uncertainty in the Momentum in MeV/c from the TPC pid
   double      Charge;           ///< Charge from the TPC pid (+1, or -1)
   TLorentzVector Position;         ///< Position at which kinematics are reported in mm, ns
   TLorentzVector PositionVariance;    ///< Variance in Position  in mm^2, ns^2
   TVector3       Direction;           ///< TPC pid direction vector in mm
   TVector3       DirectionVariance;    ///< TPC pid variance in vector direction in mm^2   
   double         NTrun;                ///< 70% of the number of clusters 
   double         Ccorr;                ///< Corrected truncated mean charge deposit in PID
   double         PullEle;              ///< Pull for TPC pid electron hypothesis
   double         PullMuon;             ///< Pull for TPC pid muon hypothesis
   double         PullPion;             ///< Pull for TPC pid pion hypothesis
   double         PullKaon;             ///< Pull for TPC pid kaon hypothesis
   double         PullProton;           ///< Pull for TPC pid proton hypothesis
   double         dEdxexpEle;           ///< Estimated dE/dx for electron hypothesis
   double         dEdxexpMuon;          ///< Estimated dE/dx for muon hypothesis      
   double         dEdxexpPion;          ///< Estimated dE/dx for pion hypothesis
   double         dEdxexpKaon;          ///< Estimated dE/dx for kaon hypothesis
   double         dEdxexpProton;        ///< Estimated dE/dx for proton hypothesis
   double         SigmaEle;             ///< Sigma estimated width of TPC pid electron hypothesis
   double         SigmaMuon;            ///< Sigma estimated width of TPC pid muon hypothesis
   double         SigmaPion;            ///< Sigma estimated width of TPC pid pion hypothesis
   double         SigmaKaon;            ///< Sigma estimated width of TPC pid kaon hypothesis
   double         SigmaProton;          ///< Sigma estimated width of TPC pid proton hypothesis
   double         Sigma0;               ///< TPC track diffusion sigma0 parameter
   double         Sigma1;               ///< TPC track diffusion sigma1 parameter
   double         Sigma2;               ///< TPC track diffusion sigma2 parameter
   double         MeanDrift;            ///< TPC track mean drift value used in diffusion model
   int            NConstituents;        //
   TVector3       TrDirection;          ///< track direction vector 
   TVector3       TrDirectionVar;       ///< variance in track direction vector
   double         TrCurvature;          ///< track curvature, units are 1/mm
   double         TrCurvatureVar;       ///< variance in track direction vector, units are (1/mm)^2
   bool           HasExtrapolation;     ///< extrapolation method of vertex is calculated or not
   double         ExtrapolatedVertexXX;    ///< for xbar vertex, this is x coordinate in mm
   double         ExtrapolatedVertexZX;    ///< for xbar vertex, this is z coordinate in mm
   double         ExtrapolatedVertexYY;    ///< for ybar vertex, this is y coordinate in mm
   double         ExtrapolatedVertexZY;    ///< for ybar vertex, this is z coordinate in mm
   bool           EnteringFGD;             ///< not sure

   TTPCTrack();
   TTPCTrack(const TTPCTrack & );
   virtual ~TTPCTrack();

   ClassDef(TTPCTrack,3); // Generated by MakeProject.
};
class TTrackerResult : public TObject {

public:
// Nested classes declaration.

public:
// Data Members.
   UInt_t      UniqueID;    ///< Unique ID number to allow matching to Global Recon object.
   string      AlgorithmName;    ///< The name of the algorithm that created this object.
   int         Detectors;        ///< Detectors used
   unsigned long Status;           ///< The status for the fit.
   double        Quality;          ///< The quality of the fit.(probability(chi2,ndof))
   int           NDOF;             ///< The number of degrees of freedom.
   double        Chi2;             ///< The chi2 of the fit.
   int           NHits;            ///< The number of hits.
   Int_t         NConstituents;    ///< The number of constituents (tracks and pids) used to build this track
   Int_t         ConstitIdx[2];    ///< Index into Constituents of the constituents used to build this track
   Int_t         NTotalConstituents;    ///< Number of all constituents, and constituents-constituents...
   TClonesArray* Constituents;          ///< All constituents, and constituents-constituents...
   bool          isForward;             ///< Sense of object.
   double        Charge;                ///< The Charge (+-1)
   double        EDeposit;              ///< The deposited charge for the object (number of pe's)
   double        Length;                ///< The total length of the object in mm
   int           matchingFailure_flag;    ///< Flag a object where the TPC-FGD matching failed
   vector<double> Likelihoods;             //
   vector<int>    Pids;                    ///< the PID that goes with Likelihoods
   TLorentzVector Position;                ///< track position 4-vector (x,y,z,t) in mm, ns
   TLorentzVector Variance;                ///< track position variance 4-vector var(x),var(y),var(z),var(t) in mm^2, ns^2 
   TVector3       Direction;               ///< track direction vector 
   TVector3       DirectionVariance;       ///< track direction variance 
   double         Momentum;                ///< track momentum MeV/c
   double         MomentumError;           ///< track momentum MeV/c
   ND::TTrackerReconModule::TTrueParticle TrueParticle;            ///< information about the true particle associated with this track
   Int_t                                  NTPCs;                   ///< Number of TPC tracks used to build this track
   TClonesArray*                          TPC;                     ///< Information about the TPC pids/tracks used to build this track
   Int_t                                  NFGDs;                   ///< Number of FGD Specific objects
   TClonesArray*                          FGD;                     ///< FGD objects associated with track
   int                                    NNodes;                  ///< The number of nodes (fgd hits + tpc tracks)
   TClonesArray*                          Nodes;                   ///< Kinematics of the track at each node in the track fit
   ND::TTrackerReconModule::TTrackerConstituent hackConstituentsObject;    ///<This is just here to fool TFile::MakeProject, not a real object. 
   ND::TTrackerReconModule::TTPCTrack           hackTPCTrack;              ///<This is just here to fool TFile::MakeProject, not a real object. 
   ND::TTrackerReconModule::TFGDTrack           hackFGDTrack;              ///<This is just here to fool TFile::MakeProject, not a real object. 
   ND::TTrackerReconModule::TTrackerNode        hackNodes;                 ///<This is just here to fool TFile::MakeProject, not a real object. 

   TTrackerResult();
   TTrackerResult(const TTrackerResult & );
   virtual ~TTrackerResult();

   ClassDef(TTrackerResult,3); // Generated by MakeProject.
};
class TTrackerVertex : public TObject {

public:
// Nested classes declaration.

public:
// Data Members.
   string      AlgorithmName;    ///< The name of the algorithm that created this object.
   int         Status;           ///< The status for the fit.
   double      Quality;          ///< The quality of the fit. Ie. the Prob(chi2,ndof)
   int         NDOF;             ///< The number of degrees of freedom.

   TTrackerVertex();
   TTrackerVertex(const TTrackerVertex & );
   virtual ~TTrackerVertex();

   ClassDef(TTrackerVertex,3); // Generated by MakeProject.
};

public:
// Data Members.
   Int_t       fNVertices;    ///< The number of added vertices
   TClonesArray* fVertices;     ///< The vector of trackerRecon vertices (none ever?).
   Int_t         fNTracks;      ///< The number of trackerRecon results
   TClonesArray* fTracks;       ///< The vector of overall trackerRecon results
   Int_t         fNFGDOther;    ///< The number of FGD tracks with no fit (none ever?).
   TClonesArray* fFGDOther;     ///< The vector of FGD tracks with no fit
   Int_t         fNTPCOther;    ///< The number of TPC tracks with no fit
   TClonesArray* fTPCOther;     ///< The vector of TPC tracks with no fit
   Int_t         fNTPCIso;      ///< The number of isolated TPC tracks with fits (none ever?)
   TClonesArray* fTPCIso;       ///< The vector of isolated TPC tracks with fits
   Int_t         fNTPCUnused;    ///< The number of unused TPC hits
   TClonesArray* fTPCUnused;     ///< The vector of unused TPC hits
   Int_t         fNFGDUnused;    ///< The number of unused FGD hits
   TClonesArray* fFGDUnused;     ///< The vector of unused TPC hits
   Int_t         fNTPCExtrapolation;    ///< The number of TPC tracks extrapolated into the FGD following Claudio's (2010a) method
   TClonesArray* fTPCExtrapolation;     ///< The vector of TPC tracks extrapolated into the FGD following Claudio's (2010a) method
   double        fMaxDrift;             ///< TPC maximum drift used for tpc other hits
   double        fCathodeOffset;        ///< TPC cathode offset used for tpc other hits
   double        fDriftVelocity;        ///< TPC drift velocity used for tpc other hits

   TTrackerReconModule();
   TTrackerReconModule(const TTrackerReconModule & );
   virtual ~TTrackerReconModule();

   ClassDef(TTrackerReconModule,5); // Generated by MakeProject.
};
} // namespace
#endif
