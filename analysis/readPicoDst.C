#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoEvent;
class Shift;
StChain *chain;
void readPicoDst(const Char_t *inputFile = "test.list",  Char_t *outputFile = "19151036")
{
  Int_t nEvents = 9999999;//10;//483;

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  //gSystem->Load("StRefMultCorr");
    gSystem->Load("StEpdUtil");
    gSystem->Load("StPicoEvent");
    gSystem->Load("StPicoDstMaker");
    gSystem->Load("StPileupUtil");
    gSystem->Load("Shift");

    chain = new StChain();

    StPicoDstMaker *picoMaker = new StPicoDstMaker(2, inputFile, "picoDst");
    Shift *anaMaker = new Shift("ana",picoMaker, outputFile);
 
    anaMaker->setBgMode(1); //signal reconstruction or background estimation

    anaMaker->setTofPidMandantory(0); //request TOF PID or not

    anaMaker->setNSigmaDaughtersCut(3.0); //default: 2.0
    anaMaker->setMass2ProtonCut(0.5, 1.5); //default: 0.5, 1.5
    anaMaker->setMass2PionCut(-0.1, 0.15); //default: -0.06, 0.1
    
    //tuned cuts (v1.0) to maximize significance in Lambda invariant mass spectrum
    anaMaker->setDcaProtonCut(0.4); //default: 0.4
    anaMaker->setDcaPionCut(1.0); //default: 1.6
    anaMaker->setDcaDaughtersCut(1.5); //default: 0.9
    anaMaker->setDLengthV0Cut(3.0); //default: 4.0
    anaMaker->setDcaV0Cut(1.0); //default: 0.6
    
    /*
    //new tuned cuts (v1.1) to balance purity and significance in Lambda invariant mass spectrum    
    anaMaker->setDcaProtonCut(0.6); 
    anaMaker->setDcaPionCut(1.8);
    anaMaker->setDcaDaughtersCut(0.7);
    anaMaker->setDLengthV0Cut(4.0);
    anaMaker->setDcaV0Cut(0.6);
    */
    //pT cut for V0 selection
    anaMaker->setPtV0Cut(0.2, 5.0);

    chain->Init();
    cout << "chain->Init();" << endl;
    int total = picoMaker->chain()->GetEntries();
    cout << " Total entries = " << total << endl;

    if (nEvents > total) nEvents = total;
    for (Int_t i = 0; i < nEvents; i++)
      {
        if (i % 1000 == 0)
	  cout << "Working on eventNumber " << i << endl;

        chain->Clear();
        int iret = chain->Make(i);

        if (iret)
	  {
            cout << "Bad return code!" << iret << endl;
            break;
	  }

        total++;

    }

    cout << "****************************************** " << endl;
    cout << "Work done... now its time to close up shop!" << endl;
    cout << "****************************************** " << endl;
    chain->Finish();
    cout << "****************************************** " << endl;
    cout << "total number of events  " << nEvents << endl;
    cout << "****************************************** " << endl;

    delete chain;

}
