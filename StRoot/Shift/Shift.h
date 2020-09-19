#ifndef Shift_hh
#define Shift_hh
//
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StEpdEpInfo.h"
#include "StMaker.h"

#include "StRoot/StPicoEvent/StPicoArrays.h"
#include "TClonesArray.h"
#include "StPicoEvent/StPicoConstants.h"
#include "StPicoEvent/StPicoCut.h"

#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>
#include "TNtuple.h"
//
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPileupUtil;
class TFile;
class TTree;
class TH1;
class TH2;
class TProfile;
class TNtuple;

class StPicoCut;
class TClonesArray;

#ifndef ST_NO_NAMESPACES
using std::string;
#endif
//
//  The class declaration. It innherits from StMaker.
class Shift : public StMaker {
    
public:
    Shift( const Char_t *name, StPicoDstMaker *picoMaker, const Char_t *jobid );   // constructor
    virtual ~Shift(){};                                 // destructor
    
    virtual void Clear(Option_t *option=""); // called after every event to cleanup
    virtual Int_t  Init();                   // called once at the beginning of your job
    virtual Int_t  Make();                   // invoked for every event
    virtual Int_t  Finish();                 // called once at the end
    
    
    // My functions
    Int_t  GetRunIndex( const Int_t run );
    Int_t  Centrality(Int_t gRefMult);
    bool   isGoodTrack(const StPicoTrack *ptrk);
    bool   isGoodGTrack(const StPicoTrack *ptrk, float pt, float eta, float dca);
    bool   isGoodEvent(const StPicoEvent *event);
    bool   isGoodTrigger(const StPicoEvent *event);
   
    void   setBgMode(bool nRotate)                  { mRotate = nRotate; }
    
    void   setTofPidMandantory(bool nTofPidMandantory)   { mTofPidMandantory = nTofPidMandantory; }

    void   setNSigmaDaughtersCut(float nSigmaCut)   { mNSigmaDaughters = nSigmaCut; }
    void   setMass2ProtonCut(float nMinPMass2, float nMaxPMass2)  { mMinMass2Proton = nMinPMass2; mMaxMass2Proton = nMaxPMass2;}
    void   setMass2PionCut(float nMinPiMass2, float nMaxPiMass2)  { mMinMass2Pion = nMinPiMass2;  mMaxMass2Pion = nMaxPiMass2;}
    
    void   setDcaProtonCut(float nProtonDca)        { mDcaProton = nProtonDca; }
    void   setDcaPionCut(float nPionDca)            { mDcaPion = nPionDca; }
    void   setDcaDaughtersCut(float nDcaDaughters)  { mDcaDaughters = nDcaDaughters; }
    void   setDLengthV0Cut(float nDecayLengthV0)    { mDecayLengthV0 = nDecayLengthV0; }
    void   setDcaV0Cut(float nDcaV0)                { mDcaV0 = nDcaV0; }
    
    void   setPtV0Cut(float nMinPtV0, float nMaxPtV0) { mMinPtV0 = nMinPtV0; mMaxPtV0 = nMaxPtV0; }
    
private:
    
    Int_t centnumber;
    const float C_C_LIGHT = 299792458; //(m/s);
    const double beamrapidity = 1.045;
    const double pionmass = 0.13957;
    const double protonmass = 0.93827;
    const double Lambdamass = 1.11568;
    
    //defined in StPicoEvent/StPicoConstants.h
    //const int nCentBin = 11;
    //const int nEpdBin = 6;
    //const int nPtBin = 5;
    //const int nEtaBin = 8; // -1.0, -0.6, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0
    //const int nPhiBin = 6;
    
    bool  mRotate; //default: 0
    
    bool  mTofPidMandantory; //default: false

    float mNSigmaDaughters; //default: 2.0
    float mMinMass2Proton; //default: 0.5
    float mMaxMass2Proton; //default: 1.5
    float mMinMass2Pion; //default: -0.06
    float mMaxMass2Pion; //default: 0.1
    
    float mDcaProton; //default: 0.4
    float mDcaPion; //default: 1.6
    float mDcaDaughters; //default: 0.9
    float mDecayLengthV0; //default: 4.0
    float mDcaV0; //default: 0.6
    
    float mMinPtV0; //default: 0.2
    float mMaxPtV0; //default: 5.0
    
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst   *mPicoDst;
    StPicoEvent *picoEvent;
    
    StPicoCut* mPicoCut;
    Int_t mMap2Track[nTrk];
    
    StEpdGeom *mEpdGeom;
    StEpdEpInfo *mEpdEpInfo;
    StPileupUtil* mPileupTool;
    
    // get recenter par
    TProfile *getepd_recen[2][7];
    TProfile *gettpc_A_recen[2];
    TProfile *gettpc_B_recen[2];
    // get shift par
    TProfile2D *pp_EPDshiftpar_sin[7];
    TProfile2D *pp_EPDshiftpar_cos[7];
    TProfile2D *pp_TPCshiftpar_Asin;
    TProfile2D *pp_TPCshiftpar_Acos;
    TProfile2D *pp_TPCshiftpar_Bsin;
    TProfile2D *pp_TPCshiftpar_Bcos;
    
    //Event QA plots
    TH1F *hVz;
    TH2F *hVr;
    TH1F *hRefMult_wPileup;
    TH1F *hRefMult_woPileup;
    TH1F *hCent;
    TH1F *hCentWeight;
    
    //QA plots for daughter
    //inclusive charged particle
    TH2F *h_betap;
    TH2F *h_dedxp;
    TH2F *h_mass2p;
    //proton
    TH2F *h_dedxp_proton;
    TH2F *h_mass2p_proton;
    TH2F *h_pt_y_protonplus;
    TH2F *h_pt_y_protonminus;
    TH2F *h_pt_y_protonplus_tof;
    TH2F *h_pt_y_protonminus_tof;
    //pion
    TH2F *h_dedxp_pion;
    TH2F *h_mass2p_pion;
    TH2F *h_pt_y_pionplus;
    TH2F *h_pt_y_pionminus;
    TH2F *h_pt_y_pionplus_tof;
    TH2F *h_pt_y_pionminus_tof;
    
    TFile *File;
    TString mout_shift;
    TString runNum;
    
    //11 centrality bins: 0-10%, 10-20%, 20-30%, 30-40%, 40-50%, 50-60%, 60-70%, 70-80%; 10-60%, 20-60%, 20-50%
    //Lambda QA plots
    TH1F *hLambdaMass[6];
    TH1F *hLambdaProtonDca;
    TH1F *hLambdaPionDca;
    TH1F *hLambdaDcaDaughters;
    TH1F *hLambdaDecayLength;
    TH1F *hLambdaDcaV0;
    TH1F *hLambdaPt;
    TH1F *hLambdaEta;
    TH1F *hLambdaPhi;
    TH1F *hLambdaRapidity;
    TH2F *hLambdaRapidityvsPt[4]; //daughter proton and pion has TOF info or not: (0,0); (0,1); (1,0); (1,1)
    TH2F *hProtonPhivsRFPhi;
    TH1F *hDaughterProtonTheta;
    
    //Anti-Lambda QA plots
    TH1F *hAntiLambdaMass[6];
    TH1F *hAntiLambdaProtonDca;
    TH1F *hAntiLambdaPionDca;
    TH1F *hAntiLambdaDcaDaughters;
    TH1F *hAntiLambdaDecayLength;
    TH1F *hAntiLambdaDcaV0;
    TH1F *hAntiLambdaPt;
    TH1F *hAntiLambdaEta;
    TH1F *hAntiLambdaPhi;
    TH1F *hAntiLambdaRapidity;
    TH2F *hAntiLambdaRapidityvsPt[4]; //daughter proton and pion has TOF info or not: (0,0); (0,1); (1,0); (1,1)
    TH2F *hAntiProtonPhivsRFPhi;
    TH1F *hAntiDaughterProtonTheta;
    
    
    TNtuple* EventsNtuple;
    TNtuple* LambdaTracksNtuple;
    
    ClassDef(Shift,0);
};
#endif
