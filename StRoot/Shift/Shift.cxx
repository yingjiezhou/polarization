#include <TFile.h>
#include <TTree.h>
#include <StMessMgr.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TRandom.h>
#include <StThreeVectorF.hh>
#include <StHelix.hh>
#include <TLorentzVector.h>

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StEpdEpFinder.h"
#include "StPileupUtil/StPileupUtil.h"
#include "StPicoEvent/StPicoV0.h"

#include "../run/run.h"
#include "../run/badrun.h"
#include "Shift.h"

const int epdRingIndex = 4; //for direct flow check, EPD_AB EP

ClassImp(Shift)

//__________________________________________________________________________________
Shift::Shift( const Char_t *name, StPicoDstMaker *picoMaker, const Char_t *jobid ) : StMaker(name), mRotate(0), mTofPidMandantory(0), mNSigmaDaughters(2.0), mMinMass2Proton(0.5), mMaxMass2Proton(1.5), mMinMass2Pion(-0.06), mMaxMass2Pion(0.1), mDcaProton(0.4), mDcaPion(1.6), mDcaDaughters(0.9), mDecayLengthV0(4.0), mDcaV0(0.6), mMinPtV0(0.2), mMaxPtV0(5.0) {
    mPicoDstMaker = picoMaker;
    mPicoDst = 0;
    mPicoCut = new StPicoCut();
    
    mout_shift=Form("%s_polarization.root", jobid);
    runNum=Form("%s", jobid);
}

//__________________________________________________________________________________
Int_t Shift::Init() {
    cout << "Init" << endl;
    mEpdGeom = new StEpdGeom();
    
    mPileupTool = new StPileupUtil();
    mPileupTool->init();
    
    cout << "Signal or Background mode: " << mRotate << endl;
    
    cout << "<<<<<<<<<<<<<  PID and LAMBDA TOPOLOGICAL CUT SETTING  >>>>>>>>" << endl;
    cout << "nSigmaPion/nSigmaProton cut: " << mNSigmaDaughters << endl;
    cout << "Proton Mass2 cut: " << mMinMass2Proton << ", " << mMaxMass2Proton << endl;
    cout << "Pion Mass2 cut: " << mMinMass2Pion << ", " << mMaxMass2Pion << endl;
    cout << "Daughter proton DCA cut: " << mDcaProton << endl;
    cout << "Daughter pion DCA cut: " << mDcaPion << endl;
    cout << "Daughter proton-pion DCA cut: " << mDcaDaughters << endl;
    cout << "V0 decay length cut: " << mDecayLengthV0 << endl;
    cout << "V0 DCA cut: " << mDcaV0 << endl;
    cout << "V0 pT cut: " << mMinPtV0 << " < pT < " << mMaxPtV0 << endl;
    
    //===============================
    //  Define Histograms
    //===============================
    
    File = new TFile(mout_shift.Data(),"RECREATE");
    
    // Event and track QA plots
    hVz = new TH1F("hVz","Vz distribution",500, 0, 500);
    hVr = new TH2F("hVr","Vy Vs Vx", 100, -10, 10, 100, -10, 10);
    hRefMult_wPileup = new TH1F("hRefMult_wPileup","",500,0,500);
    hRefMult_woPileup = new TH1F("hRefMult_woPileup","",500,0,500);
    hCent = new TH1F("hCentrality","",9,-0.5,8.5);
    hCentWeight = new TH1F("hCentrality_Weighted","",9,-0.5,8.5);
    hCent->Sumw2();
    hCentWeight->Sumw2();
    
    //QA plots for daughter tracks
    h_betap = new TH2F("h_betap","beta_vs_p;momentum [GeV/c];1/#beta",500,0,10,600,0,6);
    h_dedxp = new TH2F("h_dedxp",";momentum/charge;dE/dx", 500,-5.0,5.0, 500,0.0,50.0);
    h_mass2p= new TH2F("h_mass2p",";momentum/charge;mass^{2}",500,-5.0,5.0, 500,-2.0,13);
    h_betap->Sumw2();
    h_dedxp->Sumw2();
    h_mass2p->Sumw2();
    
    h_dedxp_pion = new TH2F("h_dedxp_pion",";momentum/charge;dE/dx",500,-5.0,5.0, 500,0.0,50.0);
    h_mass2p_pion= new TH2F("h_mass2p_pion",";momentum/charge;mass^{2}",500,-5.0,5.0, 500,-2.0,13.0);
    h_pt_y_pionplus = new TH2F("h_pt_y_pionplus","",300,-3.0,3.0,500,0.0,5.0);
    h_pt_y_pionminus = new TH2F("h_pt_y_pionminus","",300,-3.0,3.0,500,0.0,5.0);
    h_pt_y_pionplus_tof = new TH2F("h_pt_y_pionplus_tof","",300,-3.0,3.0,500,0.0,5.0);
    h_pt_y_pionminus_tof = new TH2F("h_pt_y_pionminus_tof","",300,-3.0,3.0,500,0.0,5.0);
    h_dedxp_pion->Sumw2();
    h_mass2p_pion->Sumw2();
    h_pt_y_pionplus->Sumw2();
    h_pt_y_pionminus->Sumw2();
    h_pt_y_pionplus_tof->Sumw2();
    h_pt_y_pionminus_tof->Sumw2();
    
    h_dedxp_proton = new TH2F("h_dedxp_proton",";momentum/charge;dE/dx",500,-5.0,5.0, 500,0.0,50.0);
    h_mass2p_proton= new TH2F("h_mass2p_proton",";momentum/charge;mass^{2}",500,-5.0,5.0, 500,-2.0,13.0);
    h_pt_y_protonplus = new TH2F("h_pt_y_protonplus","",300,-3.0,3.0,500,0.0,5.0);
    h_pt_y_protonminus = new TH2F("h_pt_y_protonminus","",300,-3.0,3.0,500,0.0,5.0);
    h_pt_y_protonplus_tof = new TH2F("h_pt_y_protonplus_tof","",300,-3.0,3.0,500,0.0,5.0);
    h_pt_y_protonminus_tof = new TH2F("h_pt_y_protonminus_tof","",300,-3.0,3.0,500,0.0,5.0);
    h_dedxp_proton->Sumw2();
    h_mass2p_proton->Sumw2();
    h_pt_y_protonplus->Sumw2();
    h_pt_y_protonminus->Sumw2();
    h_pt_y_protonplus_tof->Sumw2();
    h_pt_y_protonminus_tof->Sumw2();
    
    //QA plots for Lambda / Anti-Lambda
    for(int icut=0; icut<6; icut++) {
        hLambdaMass[icut] = new TH1F(Form("hLambdaMass_cut%d", icut),"Invariant mass of total Lambda",130,1.05,1.18);
        hAntiLambdaMass[icut] = new TH1F(Form("hAntiLambdaMass_cut%d", icut),"Invariant mass of total Anti-Lambda",130,1.05,1.18);
        hLambdaMass[icut]->Sumw2();
        hAntiLambdaMass[icut]->Sumw2();
    }
    
    hLambdaProtonDca        = new TH1F("hLambdaProtonDca","hLambdaProtonDca",500,0,10); hLambdaProtonDca->Sumw2();
    hLambdaPionDca          = new TH1F("hLambdaPionDca","hLambdaPionDca",500,0,10); hLambdaPionDca->Sumw2();
    hLambdaDcaDaughters     = new TH1F("hLambdaDcaDaughters","hLambdaDcaDaughters",500,0,10); hLambdaDcaDaughters->Sumw2();
    hLambdaDecayLength      = new TH1F("hLambdaDecayLength","hLambdaDecayLength",500,0,100); hLambdaDecayLength->Sumw2();
    hLambdaDcaV0            = new TH1F("hLambdaDcaV0","hLambdaDcaV0",500,0,10); hLambdaDcaV0->Sumw2();
    hLambdaPt               = new TH1F("LambdaPt","LambdaPt",500,0,10); hLambdaPt->Sumw2();
    hLambdaEta              = new TH1F("LambdaEta","LambdaEta",300,-1.5,1.5); hLambdaEta->Sumw2();
    hLambdaPhi              = new TH1F("LambdaPhi","LambdaPhi",314,-TMath::Pi(),TMath::Pi()); hLambdaPhi->Sumw2();
    hLambdaRapidity         = new TH1F("LambdaRapidity","LambdaRapidity",300,-1.5,1.5); hLambdaRapidity->Sumw2();
    TString tofInfoState[4] = {"p0pi0", "p0pi1", "p1pi0", "p1pi1"};
    for(int i=0; i<4; i++) {
        hLambdaRapidityvsPt[i] = new TH2F(Form("LambdaRapidityvsPt_%s", tofInfoState[i].Data()),"LambdaRapidityvsPt",300,-1.5,1.5,500,0,10);
        hLambdaRapidityvsPt[i]->Sumw2();
    }
    hProtonPhivsRFPhi       = new TH2F("hProtonPhivsRFPhi", "hProtonPhivsRFPhi;#phi_{p};#phi^{RF}_{p}", 314,-TMath::Pi(),TMath::Pi(), 314,-TMath::Pi(),TMath::Pi()); hProtonPhivsRFPhi->Sumw2();
    hDaughterProtonTheta    = new TH1F("DaughterProtonTheta","DaughterProtonTheta",314,-TMath::Pi(),TMath::Pi()); hDaughterProtonTheta->Sumw2();
    
    hAntiLambdaProtonDca    = new TH1F("hAntiLambdaProtonDca","hAntiLambdaProtonDca",500,0,10); hAntiLambdaProtonDca->Sumw2();
    hAntiLambdaPionDca      = new TH1F("hAntiLambdaPionDca","hAntiLambdaPionDca",500,0,10); hAntiLambdaPionDca->Sumw2();
    hAntiLambdaDcaDaughters = new TH1F("hAntiLambdaDcaDaughters","hAntiLambdaDcaDaughters",500,0,10); hAntiLambdaDcaDaughters->Sumw2();
    hAntiLambdaDecayLength  = new TH1F("hAntiLambdaDecayLength","hAntiLambdaDecayLength",500,0,100); hAntiLambdaDecayLength->Sumw2();
    hAntiLambdaDcaV0        = new TH1F("hAntiLambdaDcaV0","hAntiLambdaDcaV0",500,0,10); hAntiLambdaDcaV0->Sumw2();
    hAntiLambdaPt           = new TH1F("AntiLambdaPt","AntiLambdaPt",500,0,10); hAntiLambdaPt->Sumw2();
    hAntiLambdaEta          = new TH1F("AntiLambdaEta","AntiLambdaEta",300,-1.5,1.5); hAntiLambdaEta->Sumw2();
    hAntiLambdaPhi          = new TH1F("AntiLambdaPhi","AntiLambdaPhi",314,-TMath::Pi(),TMath::Pi()); hAntiLambdaPhi->Sumw2();
    hAntiLambdaRapidity     = new TH1F("AntiLambdaRapidity","AntiLambdaRapidity",300,-1.5,1.5); hAntiLambdaRapidity->Sumw2();
    for(int i=0; i<4; i++) {
        hAntiLambdaRapidityvsPt[i] = new TH2F(Form("AntiLambdaRapidityvsPt_%s", tofInfoState[i].Data()),"AntiLambdaRapidityvsPt",300,-1.5,1.5,500,0,10);
        hAntiLambdaRapidityvsPt[i]->Sumw2();
    }
    hAntiProtonPhivsRFPhi   = new TH2F("hAntiProtonPhivsRFPhi", "hAntiProtonPhivsRFPhi;#phi_{p};#phi^{RF}_{p}", 314,-TMath::Pi(),TMath::Pi(), 314,-TMath::Pi(),TMath::Pi()); hAntiProtonPhivsRFPhi->Sumw2();
    hAntiDaughterProtonTheta= new TH1F("AntiDaughterProtonTheta","AntiDaughterProtonTheta",314,-TMath::Pi(),TMath::Pi()); hAntiDaughterProtonTheta->Sumw2();
    // end of histgrom declear
    
    //NTuple containers
    const char* EventsVarlist = "runId:eventId:triggerId:nTrack:nEpdHit:refMult:nTofMult:nTofMatch:nCh:nMipSum:bField:centnumber:gweight:vtx:vty:vtz:vzVpd";
    EventsNtuple = new TNtuple("Events", "Events", EventsVarlist);
    
    const char* LambdaTracksVarlist = "runId:eventId:triggerId:nTrack:nEpdHit:refMult:bField:centnumber:gweight:vtx:vty:vtz:vzVpd:Psi1Epd0:Psi1Epd1:Psi1Epd2:Psi1Epd3:Psi1Epd4:Psi1Epd5:" //event info
    "isLambda:v0Mass:v0Pt:v0Px:v0Py:v0Pz:v0Eta:v0Phi:v0Y:v0DecLength:v0Dca:v0Dca2D:v0DcaDaughters:v0PathLen:v0x:v0y:v0z:pPhiV0RF:pThetaV0RF:" //RC Lambda (isLambda==1)/anti-Lambda (isLambda==2) info
    "pCharge:pPt:pPx:pPy:pPz:pEta:pPhi:pY:pDca:pDca2D:pNHits:pNHitsFit:pNHitsPoss:pNdedxHits:pDedx:pNSigma:pMass2:" //  RC daughter proton track
    "piCharge:piPt:piPx:piPy:piPz:piEta:piPhi:piY:piDca:piDca2D:piNHits:piNHitsFit:piNHitsPoss:piNdedxHits:piDedx:piNSigma:piMass2"; //  RC daughter pion track
    LambdaTracksNtuple = new TNtuple("LambdaTracks", "LambdaTracksVarlist", LambdaTracksVarlist);
    
    // connecting recentering and shift parameter
    TString runString = runNum;
    int found = 8;
    runString.Replace(found, runString.Length(), "");
    
    TFile *re_file = TFile::Open(Form("/star/u/ypwang/disk01/analysisFlow/polarization/FXT3_85_EPD/recenter/production_perRun/%s/%s.root",runString.Data(),runString.Data()));
    gettpc_A_recen[0]     = (TProfile*)re_file -> Get("TPCqx_A_recen");
    gettpc_A_recen[1]     = (TProfile*)re_file -> Get("TPCqy_A_recen");
    gettpc_B_recen[0]     = (TProfile*)re_file -> Get("TPCqx_B_recen");
    gettpc_B_recen[1]     = (TProfile*)re_file -> Get("TPCqy_B_recen");
    for(int i=0; i<7; i++){
        getepd_recen[0][i]     = (TProfile*)re_file -> Get(Form("EPDQx_ring%d_recen",i));
        getepd_recen[1][i]     = (TProfile*)re_file -> Get(Form("EPDQy_ring%d_recen",i));
    }
    
    TFile *shift_file = TFile::Open(Form("/star/u/ypwang/disk01/analysisFlow/polarization/FXT3_85_EPD/shiftpar/production_perRun/%s/%s.root",runString.Data(), runString.Data()));
    for(int i=0; i<7; i++){
        pp_EPDshiftpar_sin[i]  = (TProfile2D*)shift_file->Get(Form("EPDshiftpar_sin_ring%d",i));
        pp_EPDshiftpar_cos[i]  = (TProfile2D*)shift_file->Get(Form("EPDshiftpar_cos_ring%d",i));
    }
    
    pp_TPCshiftpar_Asin  = (TProfile2D*)shift_file->Get("TPCshiftpar_Asin");
    pp_TPCshiftpar_Bsin  = (TProfile2D*)shift_file->Get("TPCshiftpar_Bsin");
    pp_TPCshiftpar_Acos  = (TProfile2D*)shift_file->Get("TPCshiftpar_Acos");
    pp_TPCshiftpar_Bcos  = (TProfile2D*)shift_file->Get("TPCshiftpar_Bcos");
    
    cout << "End of Histograms" << endl;
    return kStOK;
}
//__________________________________________________________________________________
void Shift::Clear(Option_t *opt)
{
    StMaker::Clear();
}

//__________________________________________________________________________________
Int_t Shift::Finish() {
    cout << "Shift::Finish()\n";
    //===============================
    //  Write Histograms
    //===============================
    
    File->cd();
    
    //Event and track QA
    hVz->Write();
    hVr->Write();
    hRefMult_wPileup->Write();
    hRefMult_woPileup->Write();
    hCent->Write();
    hCentWeight->Write();
    
    //QA plots
    h_betap->Write();
    h_dedxp->Write();
    h_mass2p->Write();
    
    h_dedxp_pion->Write();
    h_mass2p_pion->Write();
    h_pt_y_pionplus->Write();
    h_pt_y_pionminus->Write();
    h_pt_y_pionplus_tof->Write();
    h_pt_y_pionminus_tof->Write();
    
    h_dedxp_proton->Write();
    h_mass2p_proton->Write();
    h_pt_y_protonplus->Write();
    h_pt_y_protonminus->Write();
    h_pt_y_protonplus_tof->Write();
    h_pt_y_protonminus_tof->Write();
    
    //Lambda/anti-Lambda QA plots
    for(int icut=0; icut<6; icut++) {
        hLambdaMass[icut]->Write();
        hAntiLambdaMass[icut]->Write();
    }
    
    hLambdaProtonDca->Write();
    hLambdaPionDca->Write();
    hLambdaDcaDaughters->Write();
    hLambdaDecayLength->Write();
    hLambdaDcaV0->Write();
    hLambdaPt->Write();
    hLambdaEta->Write();
    hLambdaPhi->Write();
    hLambdaRapidity->Write();
    for(int i=0; i<4; i++)
        hLambdaRapidityvsPt[i]->Write();
    hProtonPhivsRFPhi->Write();
    hDaughterProtonTheta->Write();
    
    hAntiLambdaProtonDca->Write();
    hAntiLambdaPionDca->Write();
    hAntiLambdaDcaDaughters->Write();
    hAntiLambdaDecayLength->Write();
    hAntiLambdaDcaV0->Write();
    hAntiLambdaPt->Write();
    hAntiLambdaEta->Write();
    hAntiLambdaPhi->Write();
    hAntiLambdaRapidity->Write();
    for(int i=0; i<4; i++)
        hAntiLambdaRapidityvsPt[i]->Write();
    hAntiProtonPhivsRFPhi->Write();
    hAntiDaughterProtonTheta->Write();
    
    //NTuple
    EventsNtuple->Write();
    LambdaTracksNtuple->Write();
    
    //File->Write();
    File->Close();
    
    return kStOK;
}
//__________________________________________________________________________________

Int_t Shift::GetRunIndex( const Int_t run ) {
    Int_t runindex = -999;
    for(Int_t i=0; i<nrun; i++){
        if(run==numbers[i]) runindex = i;
    }
    return runindex;
}
//---------------------------------------------------------------------------------

Int_t Shift::Centrality(int gRefMult )
{
    int centrality;
    //int centFull[9]={15,22,32,43,57,73,92,117,133};
    int centFull[9]={3,7,14,24,38,57,83,118,141};   //new final definitions by Daniel Cebra
    if      (gRefMult>=centFull[8] && gRefMult<=195) centrality=8;  //0-5%
    else if (gRefMult>=centFull[7] && gRefMult<centFull[8]) centrality=7;  //5-10%
    else if (gRefMult>=centFull[6] && gRefMult<centFull[7]) centrality=6;  //10-20%
    else if (gRefMult>=centFull[5] && gRefMult<centFull[6]) centrality=5;  //20-30%
    else if (gRefMult>=centFull[4] && gRefMult<centFull[5]) centrality=4;  //30-40%
    else if (gRefMult>=centFull[3] && gRefMult<centFull[4]) centrality=3;  //40-50%
    else if (gRefMult>=centFull[2] && gRefMult<centFull[3]) centrality=2;  //50-60%
    else if (gRefMult>=centFull[1] && gRefMult<centFull[2]) centrality=1;  //60-70%
    else if (gRefMult>=centFull[0] && gRefMult<centFull[1]) centrality=0;  //70-80%
    else centrality = 9;
    
    return centrality;
}
//---------------------------------------------------------------------------------

Int_t Shift::Make() {
    //Begining of Event loop
    
    //------------------------------------------------------------------
    if(!mPicoDstMaker) {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }
    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }
    picoEvent = (StPicoEvent*)mPicoDst->event();
    if( !picoEvent ){
        LOG_WARN << " No PicoEvent! Skip! " << endm;
        return kStWarn;
    }
    //------------------------------------------------------------------
    
    if(!isGoodTrigger(picoEvent)) return 0;  //trigger cut
    if(!isGoodEvent(picoEvent)) return 0;  //vertex cuts
    
    int countrefmult=0;
    int runnumber = picoEvent->runId();
    int runindex = GetRunIndex(runnumber);
    
    //exclude bad runs
    /*for(int ii=0; ii<24; ii++)
     {
     if(runnumber == badrun[ii])
     return 0;
     }
     */
    
    TVector3 pVertex = picoEvent->primaryVertex();
    Float_t BField = picoEvent->bField();
    const Int_t nTrack = mPicoDst->numberOfTracks();
    
    Int_t nepdHits = mPicoDst->numberOfEpdHits();
    if(nepdHits<75) return 0; 

    // primary track loop for determine refmult ----------------------------------------------
    for (Int_t itr=0;itr<nTrack;itr++) {
        const StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(itr);
        
        if(!ptrk)  continue;
        if(!ptrk->isPrimary())  continue;  // now selecting primary tracks
        if(!isGoodTrack(ptrk))  continue;
        countrefmult++;
    }
    hRefMult_wPileup->Fill(countrefmult);
    // End of track loop ---------------------------------------------
    
    //Sooraj's centrality definition and pile-up events removal tool
    mPileupTool->initEvent(mPicoDst);
    int   refMultPrim = mPileupTool->get_refMultPrim();
    int   centnumber  = mPileupTool->get_centrality9();
    double gweight    = mPileupTool->get_centralityWeight();
    int   nch         = mPileupTool->get_nch();
    float nmipsum     = mPileupTool->get_nmipsum();
    bool  isPileUp    = mPileupTool->isPileupEPD();
    
    if(isPileUp) return 0;
    if(centnumber > 8 || centnumber < 0) return kStOK;
    
    float tofMult = picoEvent->btofTrayMultiplicity();
    float tofmatch= picoEvent->nBTOFMatch();
    float refMult = refMultPrim;
    
    hVz->Fill(pVertex.Z());
    hVr->Fill(pVertex.X(), pVertex.Y());
    
    //---------------------Centrality----------------------------------
    //centnumber = Centrality(refMult);
    //double gweight = 1.0;
    //cout << "cent number = " << centnumber  << endl;
    //if(centnumber > 8 || centnumber < 0) return kStOK;
    
    hRefMult_woPileup->Fill(refMult);
    hCent->Fill(centnumber);
    hCentWeight->Fill(centnumber, gweight);
    
    //fill ntuple with qualified event
    float array2[50];
    int idx2 = 0;
    
    array2[idx2++] = picoEvent->runId(); //run id
    array2[idx2++] = picoEvent->eventId(); //event id
    array2[idx2++] = 620052; //trigger id
    // Event
    array2[idx2++] = nTrack; //track number
    array2[idx2++] = nepdHits; //epd hits
    array2[idx2++] = refMult; //qualified primary track from TPC
    array2[idx2++] = tofMult; //TOF track number
    array2[idx2++] = tofmatch; //TOF-TPC matched track number
    array2[idx2++] = nch; //number of charged particle
    array2[idx2++] = nmipsum; //EPD mip sum
    
    array2[idx2++] = BField; //b-field
    array2[idx2++] = centnumber; //centrality bin index
    array2[idx2++] = gweight; //weighting factor
    array2[idx2++] = pVertex.X();//vx
    array2[idx2++] = pVertex.Y();//vy
    array2[idx2++] = pVertex.Z();//vz
    array2[idx2++] = picoEvent->vzVpd();//vpd vz
    
    EventsNtuple->Fill(array2);
    
    //event plane estimation
    double Qx_rawep_TPC_A=0.0, Qy_rawep_TPC_A=0.0;
    double Qx_rawep_TPC_B=0.0, Qy_rawep_TPC_B=0.0;
    double Qx_recen_TPC_A=0.0, Qy_recen_TPC_A=0.0;
    double Qx_recen_TPC_B=0.0, Qy_recen_TPC_B=0.0;
    double weight_TPC_A = 0.0, weight_TPC_B = 0.0;
    
    double psi1rawep_A_TPC=0.0, psi1rawep_B_TPC=0.0;
    double psi1recen_A_TPC=0.0, psi1recen_B_TPC=0.0;
    double psi1shift_A_TPC=0.0, psi1shift_B_TPC=0.0;
    
    // primary track loop for determine refmult ----------------------------------------------
    for (Int_t itr=0;itr<nTrack;itr++) {
        const StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(itr);
        
        if(!ptrk)  continue;
        
        if(!ptrk->isPrimary())  continue;  // now selecting primary tracks
        if(!isGoodTrack(ptrk))  continue;
        
        Float_t pt  = ptrk->pMom().Perp(); // zero for global tracks
        Float_t phi = ptrk->pMom().Phi();
        Float_t eta = ptrk->pMom().PseudoRapidity();
        
        if(! (pt > 0.2 && pt < 2.0)) continue;  // using pt[0.2, 2.0]GeV/c particles to reconstruct 2nd EP
        //cout << "pt = " << pt << endl;
        
        //float trkWeight = eta + beamrapidity;
        float trkWeight = pt;
        if(eta > -2.0 && eta < -1.1)
        {
            Qx_rawep_TPC_A += trkWeight*(cos(1.0*phi));
            Qy_rawep_TPC_A += trkWeight*(sin(1.0*phi));
            Qx_recen_TPC_A += trkWeight*(cos(1.0*phi) - gettpc_A_recen[0]->GetBinContent(centnumber+1));
            Qy_recen_TPC_A += trkWeight*(sin(1.0*phi) - gettpc_A_recen[1]->GetBinContent(centnumber+1));
            weight_TPC_A++;
        }
        
        else if(eta > -1.0 && eta < 0)
        {
            Qx_rawep_TPC_B += trkWeight*(cos(1.0*phi));
            Qy_rawep_TPC_B += trkWeight*(sin(1.0*phi));
            Qx_recen_TPC_B += trkWeight*(cos(1.0*phi) - gettpc_B_recen[0]->GetBinContent(centnumber+1));
            Qy_recen_TPC_B += trkWeight*(sin(1.0*phi) - gettpc_B_recen[1]->GetBinContent(centnumber+1));
            weight_TPC_B++;
        }
        else continue;
    }
    if(Qy_rawep_TPC_A ==0. || Qx_rawep_TPC_A==0. || Qy_rawep_TPC_B==0. || Qx_rawep_TPC_B==0.) return 0;
    if(Qy_recen_TPC_A ==0. || Qx_recen_TPC_A==0. || Qy_recen_TPC_B==0. || Qx_recen_TPC_B==0.) return 0;
    if(weight_TPC_A == 0 || weight_TPC_B == 0) return 0;
    Qx_rawep_TPC_A /= weight_TPC_A;
    Qy_rawep_TPC_A /= weight_TPC_A;
    Qx_recen_TPC_A /= weight_TPC_A;
    Qy_recen_TPC_A /= weight_TPC_A;
    Qx_rawep_TPC_B /= weight_TPC_B;
    Qy_rawep_TPC_B /= weight_TPC_B;
    Qx_recen_TPC_B /= weight_TPC_B;
    Qy_recen_TPC_B /= weight_TPC_B;
    
    psi1rawep_A_TPC = (1.0/1.0) * atan2(Qy_rawep_TPC_A, Qx_rawep_TPC_A);
    psi1rawep_B_TPC = (1.0/1.0) * atan2(Qy_rawep_TPC_B, Qx_rawep_TPC_B);
    
    psi1recen_A_TPC = (1.0/1.0) * atan2(Qy_recen_TPC_A, Qx_recen_TPC_A);
    psi1recen_B_TPC = (1.0/1.0) * atan2(Qy_recen_TPC_B, Qx_recen_TPC_B);
    
    psi1shift_A_TPC = psi1recen_A_TPC;
    psi1shift_B_TPC = psi1recen_B_TPC;
    
    double shiftsin_TPC_A=0, shiftcos_TPC_A=0;
    double shiftsin_TPC_B=0, shiftcos_TPC_B=0;
    
    for(int i=0; i<20; i++)
    {
        psi1shift_A_TPC += (2.0/(i+1)) * (-pp_TPCshiftpar_Asin->GetBinContent(centnumber+1, i+1)*cos(1.0*(i+1)*psi1recen_A_TPC) + pp_TPCshiftpar_Acos->GetBinContent(centnumber+1, i+1)*sin(1.0*(i+1)*psi1recen_A_TPC));
        psi1shift_B_TPC += (2.0/(i+1)) * (-pp_TPCshiftpar_Bsin->GetBinContent(centnumber+1, i+1)*cos(1.0*(i+1)*psi1recen_B_TPC) + pp_TPCshiftpar_Bcos->GetBinContent(centnumber+1, i+1)*sin(1.0*(i+1)*psi1recen_B_TPC));
    }
    
    if(psi1shift_A_TPC < -pi)   psi1shift_A_TPC += 2*pi;
    if(psi1shift_B_TPC < -pi)   psi1shift_B_TPC += 2*pi;
    if(psi1shift_A_TPC >  pi)   psi1shift_A_TPC -= 2*pi;
    if(psi1shift_B_TPC >  pi)   psi1shift_B_TPC -= 2*pi;
    
    // TPC end
    
    //-----------------Get EPD information------------------------------
    StPicoEpdHit *epdHit;
    TVector3 StraightLine_center;
    TVector3 StraightLine_random;
    double phi_epd_center;
    double phi_epd_random;
    double eta_epd_center;
    double eta_epd_random;
    
    double mip;
    double TileWeight           = 0;
    
    double Qx_rawep_EPD[7]={0.0}, Qy_rawep_EPD[7]={0.0};
    double Qx_recen_EPD[7]={0.0}, Qy_recen_EPD[7]={0.0};
    
    double psi1rawep_EPD[7]={0.0};
    double psi1recen_EPD[7]={0.0};
    double psi1shift_EPD[7]={0.0};
    double weight_EPD[7]={0};
    
    for(Int_t iHit=0; iHit<nepdHits; iHit++){
        epdHit = mPicoDst->epdHit(iHit);
        mip = epdHit->nMIP();
        int iring = epdHit->row() -1;//(1~16)-1 -> 0-15
        if( !epdHit) continue;
        if(epdHit->id() > 0 ) continue;                  // unique tile identifier, absolute value is 100*position+tile, sign is +1/-1 for West/East
        
        int ringgroup = mEpdEpInfo->RingGroup(iring);   //0: 0-7, 1: 8-15    0-> inner most
        if(ringgroup == -1) continue;
        
        StraightLine_center = mEpdGeom->TileCenter(epdHit->id())        - picoEvent->primaryVertex();
        StraightLine_random = mEpdGeom->RandomPointOnTile(epdHit->id()) - picoEvent->primaryVertex();
        //StraightLine_center = mEpdGeom->TileCenter(epdHit->id())       ;
        //StraightLine_random = mEpdGeom->RandomPointOnTile(epdHit->id());
        
        phi_epd_center = StraightLine_center.Phi();
        eta_epd_center = StraightLine_center.Eta();
        phi_epd_random = StraightLine_random.Phi();
        eta_epd_random = StraightLine_random.Eta();
        
        if(mip < 0.3) continue;
        TileWeight = (mip > 2) ? 2 : mip;
        
        for(int i=0; i<4; i++)  // 0: EPD-A, 1: EPD-B, 2: EPD-C, 3: EPD-D
        {
            if(i==ringgroup){
                Qx_rawep_EPD[i] += TileWeight * cos(1.0*phi_epd_center);
                Qy_rawep_EPD[i] += TileWeight * sin(1.0*phi_epd_center);
                Qx_recen_EPD[i] += TileWeight * cos(1.0*phi_epd_center) - getepd_recen[0][i]->GetBinContent(centnumber+1);
                Qy_recen_EPD[i] += TileWeight * sin(1.0*phi_epd_center) - getepd_recen[1][i]->GetBinContent(centnumber+1);
                weight_EPD[i] += TileWeight;
            }
        }
        if(ringgroup == 0 || ringgroup == 1) // 4: EPD-AB
        {
            Qx_rawep_EPD[4] += TileWeight * cos(1.0*phi_epd_center);
            Qy_rawep_EPD[4] += TileWeight * sin(1.0*phi_epd_center);
            Qx_recen_EPD[4] += TileWeight * cos(1.0*phi_epd_center) - getepd_recen[0][4]->GetBinContent(centnumber+1);
            Qy_recen_EPD[4] += TileWeight * sin(1.0*phi_epd_center) - getepd_recen[1][4]->GetBinContent(centnumber+1);
            weight_EPD[4] += TileWeight;
        }
        if(ringgroup == 2 || ringgroup == 3)  // 5: EPD-CD
        {
            Qx_rawep_EPD[5] += TileWeight * cos(1.0*phi_epd_center);
            Qy_rawep_EPD[5] += TileWeight * sin(1.0*phi_epd_center);
            Qx_recen_EPD[5] += TileWeight * cos(1.0*phi_epd_center) - getepd_recen[0][5]->GetBinContent(centnumber+1);
            Qy_recen_EPD[5] += TileWeight * sin(1.0*phi_epd_center) - getepd_recen[1][5]->GetBinContent(centnumber+1);
            weight_EPD[5] += TileWeight;
        }
        if(ringgroup == 0 || ringgroup == 1 || ringgroup == 2 || ringgroup == 3)   // 6: EPD-ABCD
        {
            Qx_rawep_EPD[6] += TileWeight * cos(1.0*phi_epd_center);
            Qy_rawep_EPD[6] += TileWeight * sin(1.0*phi_epd_center);
            Qx_recen_EPD[6] += TileWeight * cos(1.0*phi_epd_center) - getepd_recen[0][6]->GetBinContent(centnumber+1);
            Qy_recen_EPD[6] += TileWeight * sin(1.0*phi_epd_center) - getepd_recen[1][6]->GetBinContent(centnumber+1);
            weight_EPD[6] += TileWeight;
        }
    }
    
    for(int i=0; i<7; i++){
        if(Qx_rawep_EPD[i] ==0. || Qy_rawep_EPD[i] == 0. || weight_EPD[i] == 0.) return 0;
        if(Qx_recen_EPD[i] ==0. || Qy_recen_EPD[i] == 0.) return 0;
        if(weight_EPD[i] == 0) return 0;
        Qx_rawep_EPD[i] /= weight_EPD[i];
        Qy_rawep_EPD[i] /= weight_EPD[i];
        Qx_recen_EPD[i] /= weight_EPD[i];
        Qy_recen_EPD[i] /= weight_EPD[i];
    }
    
    for(int i=0; i<7; i++){
        psi1rawep_EPD[i] = atan2(Qy_rawep_EPD[i], Qx_rawep_EPD[i])/1.0;
        psi1recen_EPD[i] = atan2(Qy_recen_EPD[i], Qx_recen_EPD[i])/1.0;
        psi1shift_EPD[i] = psi1recen_EPD[i];
    }
    for(int j=0; j<7; j++)
    {
        for(int i=0; i<20; i++)
        {
            psi1shift_EPD[j] += (2.0/(i+1)) * ( -pp_EPDshiftpar_sin[j]->GetBinContent(centnumber+1, i+1)*cos(1.0*(i+1)*psi1recen_EPD[j]) + pp_EPDshiftpar_cos[j]->GetBinContent(centnumber+1, i+1)*sin(1.0*(i+1)*psi1recen_EPD[j]) );
        }
    }
    for(int i=0; i<7; i++)
    {
        if(psi1rawep_EPD[i] < -pi) {psi1rawep_EPD[i] += 2*pi;}
        if(psi1rawep_EPD[i] >  pi) {psi1rawep_EPD[i] -= 2*pi;}
        if(psi1recen_EPD[i] < -pi) {psi1recen_EPD[i] += 2*pi;}
        if(psi1recen_EPD[i] >  pi) {psi1recen_EPD[i] -= 2*pi;}
        if(psi1shift_EPD[i] < -pi) {psi1shift_EPD[i] += 2*pi;}
        if(psi1shift_EPD[i] >  pi) {psi1shift_EPD[i] -= 2*pi;}
    }
    
    //calculated event plane resolution with EPD detector
    //[EPD ring index][centrality bin index]
    //Daniel centrality defiition
    //with naive pile-up removal
    /*Double_t ep_res[6][9] = {   {0.460121, 0.51597, 0.575868, 0.608141, 0.584389, 0.487711, 0.310511, 0.149618, 0.0980786}, //EPD A
     {0.343367, 0.406859, 0.490267, 0.57394, 0.628814, 0.62606, 0.523649, 0.343579, 0.260974},   //EPD B
     {0.225503, 0.260063, 0.323094, 0.414951, 0.505105, 0.555628, 0.515937, 0.372944, 0.27201},  //EPD C
     {0.137963, 0.161099, 0.213311, 0.298966, 0.399414, 0.48247, 0.495851, 0.399485, 0.309918},  //EPD D
     {0.460525, 0.55144, 0.653593, 0.728052, 0.751016, 0.704718, 0.540918, 0.309036, 0.196777},  //EPD AB
     {0.166826, 0.225428, 0.322891, 0.448159, 0.568861, 0.646506, 0.632342, 0.487365, 0.361758}};//EPD CD
     
    Double_t ep_res[6][9] = {   {0.384811, 0.477148, 0.566436, 0.606515, 0.583544, 0.487515, 0.312345, 0.154445, 0.123437}, //EPD A
        {0.323767, 0.392151, 0.485162, 0.572733, 0.628116, 0.625733, 0.525295, 0.350742, 0.297035}, //EPD B
        {0.215707, 0.25303, 0.319736, 0.413983, 0.504376, 0.554719, 0.515763, 0.375254, 0.287382}, //EPD C
        {0.138368, 0.158479, 0.211554, 0.298349, 0.398875, 0.481564, 0.495276, 0.401109, 0.321475}, //EPD D
        {0.39447, 0.513783, 0.643334, 0.726082, 0.749913, 0.703853, 0.54178, 0.313546, 0.225094}, //EPD AB
        {0.16084, 0.218058, 0.319043, 0.446904, 0.567961, 0.645357, 0.631994, 0.490368, 0.385358} }; //EPD CD
    */
    
     //Sooraj's centrality definition and pile-up method
     Double_t ep_res[6][9] = {  {0.358743, 0.435728, 0.529476, 0.594668, 0.605978, 0.554662, 0.421565, 0.263544, 0.144664}, //EPD A
                                {0.307724, 0.359375, 0.440201, 0.530430, 0.603391, 0.636488, 0.598597, 0.485208, 0.334137}, //EPD B
                                {0.207500, 0.235380, 0.283107, 0.363338, 0.456320, 0.532852, 0.553998, 0.492846, 0.357239}, //EPD C
                                {0.137093, 0.146728, 0.181296, 0.249829, 0.341917, 0.437779, 0.500766, 0.486312, 0.384718}, //EPD D
                                {0.363552, 0.459283, 0.587017, 0.692110, 0.745520, 0.742206, 0.653636, 0.487943, 0.293904}, //EPD AB
                                {0.148810, 0.188700, 0.266042, 0.379118, 0.503148, 0.608152, 0.655860, 0.613989, 0.468253}};//EPD CD
     
    //calculate direct flow for inclusive charged hadrons, pions and protons
    for (Int_t itr=0;itr<nTrack;itr++) {
        const StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(itr);
        if(!ptrk)  continue;
        if(!ptrk->isPrimary())  continue;  // now selecting primary tracks
        if(!isGoodTrack(ptrk))  continue;
        
        mMap2Track[ptrk->id()] = itr;
        
        Double_t mom = ptrk->pMom().Mag();
        
        Int_t charge = ptrk->charge();
        if(charge == 0) continue;
        
        h_dedxp->Fill(mom/charge, ptrk->dEdx(), gweight);
        
        float Mass2 = -999.0;
        //get mass2 information from TOF
        int tofIndex_trk = ptrk->bTofPidTraitsIndex();
        int   btofMatchFlag_trk =  0;
        float btofYLocal_trk    =  -999;
        float tof_trk = 0, L_trk=0, beta_trk=0.0;//mass2= 0.0;
        StPicoPhysicalHelix helix_trk = ptrk->helix(BField);
        if(tofIndex_trk>=0) {
            StPicoBTofPidTraits *tofPid = mPicoDst->btofPidTraits(tofIndex_trk);
            btofMatchFlag_trk = tofPid->btofMatchFlag();
            btofYLocal_trk    = tofPid->btofYLocal();
            if(tofPid) {
                beta_trk = tofPid->btofBeta();
                tof_trk = tofPid->btof();
                if(beta_trk<1e-4) {
                    TVector3 btofHitPos_ = tofPid->btofHitPos();
                    const StThreeVectorF *btofHitPos = new StThreeVectorF(btofHitPos_.X(),btofHitPos_.Y(),btofHitPos_.Z());
                    const StThreeVectorF *vertexPos_ = new StThreeVectorF(pVertex.X(), pVertex.Y(), pVertex.Z());
                    L_trk = tofPathLength(vertexPos_, btofHitPos, helix_trk.curvature());
                    if(tof_trk>0) beta_trk = L_trk / (tof_trk * (C_C_LIGHT/1.e9));
                    else beta_trk = -1;
                }
            }
        }
        bool isGoodTof = btofMatchFlag_trk >0 && beta_trk > 0 && fabs(btofYLocal_trk) < 1.8;
        if(isGoodTof) Mass2 = mom * mom * (1./pow(beta_trk, 2) - 1);
        else Mass2 = -999.0;
        
        if(Mass2 == -999.0) continue;
        
        if(TMath::Abs(beta_trk)>1e-5)  h_betap->Fill(mom, 1./beta_trk, gweight);
        else h_betap->Fill(mom, 0., gweight);
        
        h_mass2p->Fill(mom/charge, Mass2, gweight);
    }//end of track loop
    
    
    // Reconstruct lambda and anti-lambda
    int opt = 0;
    for(int pos=0;pos<nTrack;pos++) {
        StPicoTrack* t_pos= (StPicoTrack*) mPicoDst->track(pos);
        if(!t_pos) continue;
        
        if(t_pos->charge()!=+1) continue;
        
        //global track selection
        float gpt_pos   = t_pos->gMom().Perp();
        float geta_pos  = t_pos->gMom().PseudoRapidity();
        float gdca_pos  = t_pos->gDCA(pVertex).Mag();
        if(!isGoodGTrack(t_pos, gpt_pos, geta_pos, gdca_pos))  continue;
        //cout << "DEBUG: t_pos " << gpt_pos << ", " << geta_pos << ", " << gdca_pos << endl;
        
        bool isProtonCand = fabs(t_pos->nSigmaProton()) < 3.0;
        bool isPionCand   = fabs(t_pos->nSigmaPion()) < 3.0;
        if( !isPionCand && !isProtonCand )  continue; //nsigma PID cuts
        
        for(int neg=0;neg<nTrack;neg++) {
            StPicoTrack* t_neg= (StPicoTrack*) mPicoDst->track(neg);
            if(!t_neg) continue;
            
            if(t_neg->charge()!=-1) continue;
            
            //global track selection
            float gpt_neg   = t_neg->gMom().Perp();
            float geta_neg  = t_neg->gMom().PseudoRapidity();
            float gdca_neg  = t_neg->gDCA(pVertex).Mag();
            if(!isGoodGTrack(t_neg, gpt_neg, geta_neg, gdca_neg))  continue;
            
            //if(!mPicoCut->passV0Daughter(t_neg, picoEvent)) continue;   // nsigma, pt, dca cut
            bool isProtonCand2 = fabs(t_neg->nSigmaProton()) < 3.0;
            bool isPionCand2   = fabs(t_neg->nSigmaPion()) < 3.0;
            if( !isPionCand2 && !isProtonCand2 )  continue; //nsigma cuts
            
            opt = 0;
            if(isProtonCand && isPionCand2)         opt = 1; //lambda --> proton + pi^{-}
            else if(isPionCand && isProtonCand2)    opt = 2; //anti-lambda --> anti-proton + pi^{+}
            else continue;
            
            //lambda reconstruction
            StPicoPhysicalHelix helix_pos(t_pos->gMom(), t_pos->origin(), BField*kilogauss, t_pos->charge());
            StPicoPhysicalHelix helix_neg(t_neg->gMom(), t_neg->origin(), BField*kilogauss, t_neg->charge());
            
            if(mRotate) { //to estimate background
                TVector3 tp1_neg = helix_neg.momentum(BField*kilogauss); // momentum at origin
                TVector3 tx1_neg = helix_neg.origin();    //origin
                tp1_neg.SetX(-tp1_neg.X());
                tp1_neg.SetY(-tp1_neg.Y());
                tx1_neg.SetX(-(tx1_neg.X() - pVertex.X()) + pVertex.X());
                tx1_neg.SetY(-(tx1_neg.Y() - pVertex.Y()) + pVertex.Y());
                StPicoPhysicalHelix helixTmp_neg(tp1_neg, tx1_neg, BField*kilogauss, t_neg->charge());
                helix_neg = helixTmp_neg;
            }
            
            pair<double,double> s = helix_pos.pathLengths(helix_neg);
            TVector3 dcaP_pos = helix_pos.at(s.first);
            TVector3 dcaP_neg = helix_neg.at(s.second);
            float v0_dcaDaughters = (dcaP_pos - dcaP_neg).Mag(); //dca between daughters of lambda
            
            TVector3 v0 = (dcaP_pos+dcaP_neg)*0.5;
            TVector3 mom_pos = helix_pos.momentumAt(s.first, BField*kilogauss);
            TVector3 mom_neg = helix_neg.momentumAt(s.second, BField*kilogauss);
            TVector3 v0_mom = mom_pos + mom_neg;
            
            double rdotp = (v0 - pVertex).Dot(v0_mom); //r dot p for v0. cut on it. should be larger than 0.
            if(rdotp <= -10000.)    continue;
            
            double energy_pos = -999.;
            double energy_neg = -999.;
            if(opt==1) {
                energy_pos = sqrt(mom_pos.Mag2() + protonmass*protonmass);
                energy_neg = sqrt(mom_neg.Mag2() + pionmass*pionmass);
            }
            else {
                energy_pos = sqrt(mom_pos.Mag2() + pionmass*pionmass);
                energy_neg = sqrt(mom_neg.Mag2() + protonmass*protonmass);
            }
            
            double v0_mass = sqrt(protonmass*protonmass + pionmass*pionmass + 2.*energy_pos*energy_neg - 2.*mom_pos.Dot(mom_neg)); //lambda mass
            if(v0_mass > 1.18 || v0_mass < 1.05) continue;
            
            if(opt==1)  hLambdaMass[0]->Fill(v0_mass, gweight);
            if(opt==2)  hAntiLambdaMass[0]->Fill(v0_mass, gweight);
            
            float angle = (v0 - pVertex).Angle(v0_mom);
            float v0_dca = (v0 - pVertex).Mag()*TMath::Sin(angle); //lambda dca to primary vertex
            //double dcav0toPV = rdotp*rdotp/v0_mom.Mag2();
            //dcav0toPV = sqrt( (v0 - pVertex).Mag2() - dcav0toPV);
            
            double v0_decaylength = (v0 - pVertex).Mag(); //lambda decay-length
            
            double pxL = v0_mom.X();
            double pyL = v0_mom.Y();
            double pzL = v0_mom.Z();
        
            double v0_energy = sqrt(pxL*pxL + pyL*pyL + pzL*pzL  + v0_mass*v0_mass);
            
            TLorentzVector V0four(TVector3(pxL, pyL, pzL), v0_energy);
            float v0_phi            = V0four.Phi();
            float v0_rapidity       = V0four.Rapidity();
            float v0_pt             = V0four.Pt();
            float v0_eta            = V0four.Eta();
            
            if(v0_pt < mMinPtV0 || v0_pt > mMaxPtV0) continue;
        
            //cout << "RC Lambda: v0Mass=" << v0_mass << ",   v0_pt=" << v0_pt << ",  v0_rapidity=" << v0_rapidity << ",  v0_dca=" << v0_dca << ",    dca_pos=" << t_pos->gDCA(pVertex).Mag() << ",   dca_neg=" << t_neg->gDCA(pVertex).Mag() << endl;

            //helix of v0: straight line
            StPicoPhysicalHelix helix_v0(v0_mom, v0, 0, 0);
            double v0_pathlen = helix_v0.pathLength(v0) - helix_v0.pathLength(pVertex);
            
            //daughter tracks information
            StPicoTrack* Dau1;
            StPicoTrack* Dau2;
            if(opt==1){
                Dau1 = (StPicoTrack*) t_pos;   //+  proton
                Dau2 = (StPicoTrack*) t_neg;   //- pion
            } else if (opt==2){
                Dau1 = (StPicoTrack*) t_neg;   //-  anti-proton
                Dau2 = (StPicoTrack*) t_pos;   //+  pion
            } else continue;
            
            float px1 =  Dau1->gMom().X();  //proton
            float py1 =  Dau1->gMom().Y();
            float pz1 =  Dau1->gMom().Z();
            float px2 = Dau2->gMom().X();  //pion
            float py2 = Dau2->gMom().Y();
            float pz2 = Dau2->gMom().Z();
            
            //tof information
            float Beta1=-999., Beta2=-999;
            float M2_dau1=-999., M2_dau2=-999.;
            if(Dau1->bTofPidTraitsIndex()>=0){
                Beta1 = mPicoDst->btofPidTraits(Dau1->bTofPidTraitsIndex())->btofBeta();
            }
            if(Dau2->bTofPidTraitsIndex()>=0){
                Beta2 = mPicoDst->btofPidTraits(Dau2->bTofPidTraitsIndex())->btofBeta();
            }
            if(Beta1!=-999.&&Beta1!=0)  M2_dau1 = (Dau1->pMom().Mag2())*(1./(Beta1*Beta1)-1.);
            if(Beta2!=-999.&&Beta2!=0)  M2_dau2 = (Dau2->pMom().Mag2())*(1./(Beta2*Beta2)-1.);
            
            //daughter tracks dca to vertex
            float dcaproton =   Dau1->gDCA(pVertex).Mag();
            float dcapion   =   Dau2->gDCA(pVertex).Mag();
            
            //daughter track quality
            float etaPion           = Dau2->gMom().PseudoRapidity();
            float etaProton         = Dau1->gMom().PseudoRapidity();
            float nHitsPion         = (float) Dau2->nHitsFit();
            float nHitsProton       = (float) Dau1->nHitsFit();
            float nHitsPionMax      = (float) Dau2->nHitsMax();
            float nHitsProtonMax    = (float) Dau1->nHitsMax();
            float nSigmaProton      = Dau1->nSigmaProton();
            float nSigmaPion        = Dau2->nSigmaPion();
            
            //daughter proton QA
            Float_t mE_proton = sqrt(Dau1->gMom().Mag2() + protonmass * protonmass);
            Float_t my_proton = 0.5*log((mE_proton + pz1)/(mE_proton - pz1));   // proton's rapidity
            if(isnan(my_proton)) continue;
            Float_t rapidity_proton = my_proton + beamrapidity;
            
            h_dedxp_proton->Fill(Dau1->gMom().Mag()/Dau1->charge(), Dau1->dEdx(), gweight);
            h_mass2p_proton->Fill(Dau1->gMom().Mag()/Dau1->charge(), M2_dau1, gweight);
            if(opt==1) {
                if(M2_dau1==-999.)
                    h_pt_y_protonplus->Fill(rapidity_proton, Dau1->gMom().Perp(), gweight); //no TOF
                else {
                    h_pt_y_protonplus_tof->Fill(rapidity_proton, Dau1->gMom().Perp(), gweight);
                }
            }
            if(opt==2) {
                if(M2_dau1==-999.)
                    h_pt_y_protonminus->Fill(rapidity_proton, Dau1->gMom().Perp(), gweight); //no TOF
                else {
                    h_pt_y_protonminus_tof->Fill(rapidity_proton, Dau1->gMom().Perp(), gweight); // TOF available
                }
            }
            
            //daughter pion QA
            Float_t mE_pion = sqrt(Dau2->gMom().Mag2() + pionmass * pionmass);
            Float_t my_pion = 0.5*log((mE_pion + pz2)/(mE_pion - pz2));   // proton's rapidity
            if(isnan(my_pion)) continue;
            Float_t rapidity_pion = my_pion + beamrapidity;
            
            h_dedxp_pion->Fill(Dau2->gMom().Mag()/Dau2->charge(), Dau2->dEdx(), gweight);
            h_mass2p_pion->Fill(Dau2->gMom().Mag()/Dau2->charge(), M2_dau2, gweight);
            if(opt==1) {
                if(M2_dau2==-999.)
                    h_pt_y_pionminus->Fill(rapidity_pion, Dau2->gMom().Perp(), gweight); //no TOF
                else {
                    h_pt_y_pionminus_tof->Fill(rapidity_pion, Dau2->gMom().Perp(), gweight);
                }
            }
            if(opt==2) {
                if(M2_dau2==-999.)
                    h_pt_y_pionplus->Fill(rapidity_pion, Dau2->gMom().Perp(), gweight); //no TOF
                else {
                    h_pt_y_pionplus_tof->Fill(rapidity_pion, Dau2->gMom().Perp(), gweight); // TOF available
                }
            }
            
            //tighter cuts for daughter particles
            if(fabs(nSigmaPion) > mNSigmaDaughters)  continue;
            if(fabs(nSigmaProton) > mNSigmaDaughters)  continue;
            
            if(!mTofPidMandantory) { //TOF PID not mandantory requested
                if(M2_dau2!=-999. && (M2_dau2<mMinMass2Pion || M2_dau2>mMaxMass2Pion))   continue;
                if(M2_dau1!=-999. && (M2_dau1<mMinMass2Proton  || M2_dau1>mMaxMass2Proton))    continue;
            }
            else { //TOF PID requested
                if(M2_dau2==-999.)    continue;
                if(M2_dau1==-999.)    continue;
            }
            
            if(opt==1)  { hLambdaMass[1]->Fill(v0_mass, gweight);        hLambdaProtonDca->Fill(dcaproton, gweight);      hLambdaPionDca->Fill(dcapion, gweight);}
            if(opt==2)  { hAntiLambdaMass[1]->Fill(v0_mass, gweight);    hAntiLambdaProtonDca->Fill(dcaproton, gweight);  hAntiLambdaPionDca->Fill(dcapion, gweight);}
            
            //topological cuts for lambda or anti-lambda
            if(dcaproton < mDcaProton || dcapion < mDcaPion) continue; //tuned cut
            //if(dcaproton<0.75 || dcapion<1.5) continue; //tight cut
            //if(dcaproton<0.75 || dcapion<1.1) continue; //loose cut
            //if(dcaproton<0.6 || dcapion<2.2) continue; //Joey's cut
            if(opt==1)  {   hLambdaMass[2]->Fill(v0_mass, gweight);      hLambdaDcaDaughters->Fill(fabs(v0_dcaDaughters), gweight);  }
            if(opt==2)  {   hAntiLambdaMass[2]->Fill(v0_mass, gweight);  hAntiLambdaDcaDaughters->Fill(fabs(v0_dcaDaughters), gweight);  }
            
            if(fabs(v0_dcaDaughters) > mDcaDaughters) continue; //tuned cut
            //if(fabs(v0_dcaDaughters)>0.8) continue; //tight or loose cut
            //if(fabs(v0_dcaDaughters)>1.25) continue; //Joey's cut
            if(opt==1)  {   hLambdaMass[3]->Fill(v0_mass, gweight);      hLambdaDecayLength->Fill(v0_decaylength, gweight);   }
            if(opt==2)  {   hAntiLambdaMass[3]->Fill(v0_mass, gweight);  hAntiLambdaDecayLength->Fill(v0_decaylength, gweight);   }
            
            if(v0_decaylength < mDecayLengthV0) continue; //tuned cut
            //if(v0_decaylength<4.0) continue; //tight or loose cut
            //if(v0_decaylength<5.79) continue; //Joey's cut
            if(opt==1)  {   hLambdaMass[4]->Fill(v0_mass, gweight);      hLambdaDcaV0->Fill(fabs(v0_dca), gweight);     }
            if(opt==2)  {   hAntiLambdaMass[4]->Fill(v0_mass, gweight);  hAntiLambdaDcaV0->Fill(fabs(v0_dca), gweight); }
            
            if(fabs(v0_dca) > mDcaV0) continue; //tuned cut
            //if(fabs(v0_dca)>0.8) continue; //tight or loose cut
            //if(fabs(v0_dca)>0.75) continue; //Joey's cut
            if(opt==1)  hLambdaMass[5]->Fill(v0_mass, gweight);
            if(opt==2)  hAntiLambdaMass[5]->Fill(v0_mass, gweight);
            
            //cout << "RC Lambda: v0_mass=" << v0_mass << ",   v0_pt=" << v0_pt << ",  v0_rapidity=" << v0_rapidity << ",  v0_dca=" << v0_dca << ",    dca_pi=" << dcapion << ",   dca_p=" << dcaproton << endl;

            //lambda rest frame
            TLorentzVector protonfour(TVector3( px1, py1, pz1 ), (double)TMath::Sqrt(px1*px1 + py1*py1 + pz1*pz1 + protonmass*protonmass));
            TVector3 v0beta = V0four.BoostVector();
            protonfour.Boost( -v0beta );  //boost proton along Lambda direction
            
            float phip_V0RF = protonfour.Phi();
            float thetap_V0RF = protonfour.Theta();
            
            //test boost
            //cout << "DEBUG: v0_pT = " << v0_pt << endl;
            //cout << "DEBUG: psi1shift_EPD[epdRingIndex] - phip_V0RF= " << psi1shift_EPD[epdRingIndex] - phip_V0RF << " with Boost( -v0beta ) " << endl;
            
            //define rapidity bins
            double rapdity_lambda = v0_rapidity + beamrapidity;
            if(opt==1) {
                hLambdaPt            -> Fill(v0_pt, gweight);
                hLambdaEta           -> Fill(v0_eta, gweight);
                hLambdaPhi           -> Fill(v0_phi, gweight);
                hLambdaRapidity      -> Fill(rapdity_lambda, gweight);
                if(M2_dau1==-999. && M2_dau2==-999.)
                    hLambdaRapidityvsPt[0]  -> Fill(rapdity_lambda, v0_pt, gweight);
                else if(M2_dau1==-999. && M2_dau2!=-999.)
                    hLambdaRapidityvsPt[1]  -> Fill(rapdity_lambda, v0_pt, gweight);
                else if(M2_dau1==!-999. && M2_dau2==-999.)
                    hLambdaRapidityvsPt[2]  -> Fill(rapdity_lambda, v0_pt, gweight);
                else
                    hLambdaRapidityvsPt[3]  -> Fill(rapdity_lambda, v0_pt, gweight);
                
                hProtonPhivsRFPhi    -> Fill(protonfour.Phi(), phip_V0RF, gweight);
                hDaughterProtonTheta -> Fill(thetap_V0RF, gweight);
            }
            
            if(opt==2) {
                hAntiLambdaPt           -> Fill(v0_pt, gweight);
                hAntiLambdaEta          -> Fill(v0_eta, gweight);
                hAntiLambdaPhi          -> Fill(v0_phi, gweight);
                hAntiLambdaRapidity     -> Fill(rapdity_lambda, gweight);
                if(M2_dau1==-999. && M2_dau2==-999.)
                    hAntiLambdaRapidityvsPt[0]  -> Fill(rapdity_lambda, v0_pt, gweight);
                else if(M2_dau1==-999. && M2_dau2!=-999.)
                    hAntiLambdaRapidityvsPt[1]  -> Fill(rapdity_lambda, v0_pt, gweight);
                else if(M2_dau1==!-999. && M2_dau2==-999.)
                    hAntiLambdaRapidityvsPt[2]  -> Fill(rapdity_lambda, v0_pt, gweight);
                else
                    hAntiLambdaRapidityvsPt[3]  -> Fill(rapdity_lambda, v0_pt, gweight);
                
                hAntiProtonPhivsRFPhi   -> Fill(protonfour.Phi(), phip_V0RF, gweight);
                hAntiDaughterProtonTheta -> Fill(thetap_V0RF, gweight);
            }
            
            // fill Lambda ntuple
            float array[200];
            int idx = 0;
            
            // event info
            array[idx++] = picoEvent->runId(); //run id
            array[idx++] = picoEvent->eventId(); //event id
            array[idx++] = 620052; //trigger id
            array[idx++] = nTrack; //track number
            array[idx++] = nepdHits; //epd hits
            array[idx++] = refMult; //qualified primary track from TPC
            array[idx++] = BField; //b-field
            array[idx++] = centnumber; //centrality bin index
            array[idx++] = gweight; //weighting factor
            array[idx++] = pVertex.X(); //vx
            array[idx++] = pVertex.Y(); //vy
            array[idx++] = pVertex.Z(); //vz
            array[idx++] = picoEvent->vzVpd(); //vpd vz
            array[idx++] = psi1shift_EPD[0]; //Psi_1 of EPD A
            array[idx++] = psi1shift_EPD[1]; //Psi_1 of EPD B
            array[idx++] = psi1shift_EPD[2]; //Psi_1 of EPD C
            array[idx++] = psi1shift_EPD[3]; //Psi_1 of EPD D
            array[idx++] = psi1shift_EPD[4]; //Psi_1 of EPD AB
            array[idx++] = psi1shift_EPD[5]; //Psi_1 of EPD CD
            
            // RC lambda
            array[idx++] = opt; //isLambda: lambda == 1; anti-lambda == 2.
            array[idx++] = opt ? v0_mass : -9999.; //mass
            array[idx++] = opt ? v0_pt : -9999.; //pT
            array[idx++] = opt ? pxL : -9999.; //px
            array[idx++] = opt ? pyL : -9999.; //py
            array[idx++] = opt ? pzL : -9999.; //pz
            array[idx++] = opt ? v0_eta : -9999.; //eta
            array[idx++] = opt ? v0_phi : -9999.; //phi
            array[idx++] = opt ? v0_rapidity : -9999.; //rapidity
            array[idx++] = opt ? v0_decaylength : -9999.; //Decay length
            array[idx++] = opt ? v0_dca : -9999.; //dca
            array[idx++] = opt ? helix_v0.geometricSignedDistance(pVertex.X(), pVertex.Y()) : -9999.; //dca 2D
            array[idx++] = opt ? v0_dcaDaughters : -9999.; //dca between daughters
            array[idx++] = opt ? v0_pathlen : -9999.; //v0 path length
            array[idx++] = opt ? v0.X() : -9999.; //v0 position x
            array[idx++] = opt ? v0.Y() : -9999.; //v0 position y
            array[idx++] = opt ? v0.Z() : -9999.; //v0 position z
            array[idx++] = opt ? phip_V0RF : -9999.; //daughter proton's phi with respect to the v0 in v0 rest frame
            array[idx++] = opt ? thetap_V0RF : -9999.; //daughter proton's theta
            
            // Dau1: proton
            array[idx++] = opt ? Dau1->charge() : -9999.; //proton charge
            array[idx++] = opt ? Dau1->gMom().Perp() : -9999.; //proton pT
            array[idx++] = opt ? px1 : -9999.; //proton px
            array[idx++] = opt ? py1 : -9999.; //proton py
            array[idx++] = opt ? pz1 : -9999.; //proton pz
            array[idx++] = opt ? etaProton : -9999.; //proton eta
            array[idx++] = opt ? Dau1->gMom().Phi() : -9999.; //proton phi
            array[idx++] = opt ? my_proton : -9999.; //proton rapidity
            array[idx++] = opt ? dcaproton : -9999.; //proton dca
            array[idx++] = opt ? Dau1->gDCAxy(pVertex.X(), pVertex.Y()) : -9999.; //proton dca 2D
            array[idx++] = opt ? Dau1->nHits() : -9999.; //proton nHits
            array[idx++] = opt ? nHitsProton : -9999.; //proton nHitsFit
            array[idx++] = opt ? nHitsProtonMax : -9999.; //proton nHitsPoss
            array[idx++] = opt ? Dau1->nHitsDedx() : -9999.; //proton nHitsDedx
            array[idx++] = opt ? Dau1->dEdx() : -9999.; //proton dEdx
            array[idx++] = opt ? nSigmaProton : -9999.; //proton nSigmaProton
            array[idx++] = opt ? M2_dau1 : -9999.; //request mass2
            
            // Dau2: pion
            array[idx++] = opt ? Dau2->charge() : -9999.; //pion charge
            array[idx++] = opt ? Dau2->gMom().Perp() : -9999.; //pion pT
            array[idx++] = opt ? px2 : -9999.; //pion px
            array[idx++] = opt ? py2 : -9999.; //pion py
            array[idx++] = opt ? pz2 : -9999.; //pion pz
            array[idx++] = opt ? etaPion : -9999.; //pion eta
            array[idx++] = opt ? Dau2->gMom().Phi() : -9999.; //pion phi
            array[idx++] = opt ? my_pion : -9999.; //pion rapidity
            array[idx++] = opt ? dcapion : -9999.; //pion dca
            array[idx++] = opt ? Dau2->gDCAxy(pVertex.X(), pVertex.Y()) : -9999.; //pion dca 2D
            array[idx++] = opt ? Dau2->nHits() : -9999.; //pion nHits
            array[idx++] = opt ? nHitsPion : -9999.; //pion nHitsFit
            array[idx++] = opt ? nHitsPionMax : -9999.; //pion nHitsPoss
            array[idx++] = opt ? Dau2->nHitsDedx() : -9999.; //pion nHitsDedx
            array[idx++] = opt ? Dau2->dEdx() : -9999.; //pion dEdx (GeV/cm)
            array[idx++] = opt ? nSigmaPion : -9999.; //pion nSigmaPion
            array[idx++] = opt ? M2_dau2 : -9999.; //pion's mass2
            
            LambdaTracksNtuple->Fill(array);
            
        } //end negative track loop
    } //end positive track loop
    
    return kStOK;
}

//__________________________________________________________________________________
bool Shift::isGoodEvent(const StPicoEvent *event)
{
    Float_t vx=event->primaryVertex().X();
    Float_t vy=event->primaryVertex().Y();
    Float_t vz=event->primaryVertex().Z();
    
    if(vz < 198 || vz > 202) return false;
    if(( (vx*vx) + (vy+2.0)*(vy+2.0) ) > 4) return false;
    
    return true;
}

bool Shift::isGoodTrack(const StPicoTrack *ptrk) { //primary track selection
    const Float_t pt  = ptrk->pMom().Perp(); // zero for global tracks
    const Float_t eta = ptrk->pMom().PseudoRapidity();
    const Int_t nHits = ptrk->nHits();
    const Float_t dca = ptrk->gDCA( picoEvent->primaryVertex() ).Mag();
    const Int_t nHitsFit = ptrk->nHitsFit();
    const Int_t nHitsPoss = ptrk->nHitsMax();
    const Float_t nHitsDedx = ptrk->nHitsDedx();
    const Float_t quality = (Float_t)nHitsFit/(Float_t)nHitsPoss;
    
    if( pt < 0.15 || pt > 5.0)  return false;
    if( eta < -2.0 || 0 < eta) return false;
    if( fabs(dca)>3.0 ) return false;
    if( nHitsFit < 15 )  return false;
    if( quality < 0.52 )  return false;
    if( nHitsDedx < 5 ) return false;
    
    return true;
}

bool Shift::isGoodGTrack(const StPicoTrack *ptrk, float pt, float eta, float dca) { //global track selection
    if( pt < 0.15 || pt > 5.0)  return false;
    if( eta < -2.0 || 0 < eta) return false;
    if( fabs(dca)>3.0 ) return false;
    if( ptrk->nHitsFit() < 15 )  return false;
    if( ptrk->nHitsFit() < 0.52*ptrk->nHitsMax() )  return false;
    if( ptrk->nHitsDedx() < 5 ) return false;
    
    return true;
}

bool Shift::isGoodTrigger(const StPicoEvent* event)
{
    if(event->isTrigger(620052)) return true;
    else return false;
}

//__________________________________________________________________________________
