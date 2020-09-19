#ifndef StPicoCut_h
#define StPicoCut_h

#include "TObject.h"

class StMuEvent;
class StMuTrack;
class StPicoV0;
class StPicoTrack;
class StPicoEvent;
class StPicoCut : public TObject {
 public:
  StPicoCut();
  ~StPicoCut();
  
  bool passEvent( StMuEvent * );
  bool passTrack( StMuTrack * );
  bool passV0Daughter( StPicoTrack * , StPicoEvent*);
  //bool passV0( StPicoV0 *, StMuEvent * );
  bool passV0( StPicoV0 *, StPicoEvent * ); //added by sss
  bool passKs( StPicoV0 * , StPicoEvent*);
  bool passLambda( StPicoV0 * , StPicoEvent*);
  bool passLbar( StPicoV0 * , StPicoEvent*);
  int  flowFlag( StMuTrack * );
  
  ClassDef(StPicoCut,1)
};

#endif
