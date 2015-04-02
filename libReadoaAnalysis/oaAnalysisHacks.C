#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"

#include "libReadoaAnalysis/ND__TGlobalReconModule.h"
#include "libReadoaAnalysis/ND__TTruthTrajectoriesModule.h"
#include "libReadoaAnalysis/ND__NRooTrackerVtx.h"
#include "libReadoaAnalysis/ND__GRooTrackerVtx.h"

//#include "libReadoaAnalysis/ND__TGlobalReconModule__TGlobalPID.h"
//#include "libReadoaAnalysis/ND__TTruthTrajectoriesModule__TTruthTrajectory.h"

#ifdef __CINT__
#pragma link C++ class std::vector<int>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class std::vector<std::vector<double> >+;
#endif

//namespace ND {
//class TGlobalReconModule;
//class TGlobalReconModule::TGlobalPID;
//}

namespace warwick
{

class OaAnalysisHacks
{
    public:
        static OaAnalysisHacks& get()
        {
            static OaAnalysisHacks instance;
            return instance;
        }

    	void printDet(ND::TGlobalReconModule::TGlobalPID* obj)
    	{
    		//std::cout << obj << std::endl;
    		obj->Print();
    		std::cout << obj->DetectorUsed[0] << std::endl;
    		//std::cout << obj->EntrancePosition[0] << std::endl;
    		obj->EntrancePosition[0].Print();
    	}

    	std::vector<int> getDetectorsUsed(ND::TGlobalReconModule::TGlobalPID* obj)
    			{
    		std::vector<int> result;
    		//std::cout << obj << std::endl;
    		int size = OaAnalysisHacks::sizeOfDetectorsUsed();//obj->Detectors;
    		for(int i = 0; i < size; ++i)
    		{
    			result.push_back( obj->DetectorUsed[i] );
    		}
    		return result;
    			}


    	std::vector<TLorentzVector> getEntrancePosition(ND::TGlobalReconModule::TGlobalPID* obj)
    			{
    		std::vector<TLorentzVector> result;
    		//std::cout << obj << std::endl;
    		int size = sizeOfDetExtrapolation();//obj->Detectors;
    		for(int i = 0; i < size; ++i)
    		{
    			//TLorentzVector hlv;
    			//obj->EntrancePosition[i].Print();
    			TLorentzVector hlv = obj->EntrancePosition[i];
    			result.push_back( hlv );
    		}
    		return result;
    			}

    	std::vector<TLorentzVector> getTrajectoryEntrancePosition(ND::TTruthTrajectoriesModule::TTruthTrajectory* obj)
    			{
    		std::vector<TLorentzVector> result;
    		//std::cout << obj << std::endl;
    		int size = nSubDetectors();//obj->Detectors;
    		for(int i = 0; i < size; ++i)
    		{
    			//TLorentzVector hlv;
    			//obj->EntrancePosition[i].Print();
    			TLorentzVector hlv = obj->EntrancePosition[i];
    			result.push_back( hlv );
    		}
    		return result;
    			}

    	std::vector<TVector3> getTrajectoryEntranceMomentum(ND::TTruthTrajectoriesModule::TTruthTrajectory* obj)
    			{
    		std::vector<TVector3> result;
    		//std::cout << obj << std::endl;
    		int size = nSubDetectors();//obj->Detectors;
    		for(int i = 0; i < size; ++i)
    		{
    			//TLorentzVector hlv;
    			//obj->EntrancePosition[i].Print();
    			TVector3 hlv = obj->EntranceMomentum[i];
    			result.push_back( hlv );
    		}
    		return result;
    			}

    	 
    	std::vector<int> getTrajectoryEnteredDetector(ND::TTruthTrajectoriesModule::TTruthTrajectory* obj)
    			{
    		std::vector<int> result;
    		//std::cout << obj << std::endl;
    		int size = nSubDetectors();//obj->Detectors;
    		for(int i = 0; i < size; ++i)
    		{
    			//TLorentzVector hlv;
    			//obj->EntrancePosition[i].Print();
    			int hlv = obj->EnteredSubdetector[i];
    			result.push_back( hlv );
    		}
    		return result;
    			}
    	 

    	std::vector<TVector3> getEntranceDirection(ND::TGlobalReconModule::TGlobalPID* obj)
    			{
    		std::vector<TVector3> result;
    		//std::cout << obj << std::endl;
    		int size = sizeOfDetExtrapolation();//obj->Detectors;
    		for(int i = 0; i < size; ++i)
    		{
    			//TLorentzVector hlv;
    			//obj->EntranceDirection[i].Print();
    			TVector3 hlv = obj->EntranceDirection[i];
    			//hlv.Print();
    			result.push_back( hlv );
    		}
    		return result;
    			}

    	std::vector<Double_t> getEntranceMomentum(ND::TGlobalReconModule::TGlobalPID* obj)
    			{
    		std::vector<Double_t> result;
    		//std::cout << obj << std::endl;
    		int size = sizeOfDetExtrapolation();//obj->Detectors;
    		for(int i = 0; i < size; ++i)
    		{
    			//TLorentzVector hlv;
    			//obj->EntranceDirection[i].Print();
    			Double_t mom = obj->EntranceMomentum[i];
    			result.push_back( mom );
    		}
    		return result;
    			}

    	std::vector<std::vector<double> > getRooTrackerHepX4(ND::NRooTrackerVtx* obj)
		{
    		int nRows = obj->StdHepN;
    		double (*array)[4] = obj->StdHepX4;
    		OaAnalysisHacks& hack = OaAnalysisHacks::get();
    		return hack.convertRooTracker2DArray(array,nRows,4);
		}

    	std::vector<std::vector<double> > getRooTrackerHepP4(ND::NRooTrackerVtx* obj)
		{
    		int nRows = obj->StdHepN;
    		double (*array)[4] = obj->StdHepP4;
    		OaAnalysisHacks& hack = OaAnalysisHacks::get();
    		return hack.convertRooTracker2DArray(array,nRows,4);
		}

    	std::vector<std::vector<double> > getRooTrackerHepX4(ND::GRooTrackerVtx* obj)
		{
    		int nRows = obj->StdHepN;
    		double (*array)[4] = obj->StdHepX4;
    		OaAnalysisHacks& hack = OaAnalysisHacks::get();
    		return hack.convertRooTracker2DArray(array,nRows,4);
		}

    	std::vector<std::vector<double> > getRooTrackerHepP4(ND::GRooTrackerVtx* obj)
		{
    		int nRows = obj->StdHepN;
    		double (*array)[4] = obj->StdHepP4;
    		OaAnalysisHacks& hack = OaAnalysisHacks::get();
    		return hack.convertRooTracker2DArray(array,nRows,4);
		}

    	std::vector<std::vector<double> > convertRooTracker2DArray(double (*array)[4], int nRows, int nCols)
    	{
    	    		std::vector<std::vector<double> > result;
    	    		for(int i = 0; i < nRows; ++i)
    	    		{
    	    			std::vector<double> v;
    	    			for(int j = 0; j < nCols; ++j)
    	    			{
    	    				double d = array[i][j];
    	    				v.push_back(d);
    	    			}
    	    			result.push_back(v);
    	    		}
    	    		return result;
    	}

//    	int getTrajID(ND::TTruthTrajectoriesModule::TTruthTrajectory* traj) { return traj->TrajID; }

    	int sizeOfDetectorsUsed() { return 19; }

    	//The indices used in detector extrapolation in the Global PID
    	// are in ND::TGlobalReconModule::InitializeExtrapolationToDetectors
    	// At the time of this comment they are:
    	// 0,1,2 TPC volumes
    	// 3,4   FGD volumes
    	// 5     P0D
    	// 6     DSECal
    	// 7,8,9,10,11,12     SMRD
    	// 13,14,15,16,17,18 P0D ECal
    	// 19,20,21,22,23,24 Tracker ECal
    	// Mapping the index to the Tracker ECal:
    	// 6 = DSECal
    	// 19 = Tracker Top Left
    	// 20 = Tracker Top Right
    	// 21 = Tracker Bottom Left
    	// 22 = Tracker Bottom Right
    	// 23 = Tracker Side Left
    	// 24 = Tracker Side Right
    	int sizeOfDetExtrapolation() { return 25; }
    	int nSubDetectors() { return 13; }

    private:
        OaAnalysisHacks() { }
        OaAnalysisHacks(OaAnalysisHacks const&); //hide
        void operator=(OaAnalysisHacks const&); //hide
};

}

void oaAnalysisHacks()
{
	// nothing to do.
}

