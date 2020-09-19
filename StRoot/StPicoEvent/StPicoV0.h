#ifndef StPicoV0_hh
#define StPicoV0_hh

class StPicoTrack;
class StPicoEvent;
class StMuEvent;

#include "TObject.h"
#include "StThreeVectorF.hh"
#include "TVector3.h"

class StPicoV0 : public TObject {
 public:
  StPicoV0();
  StPicoV0(StPicoV0*);
  ~StPicoV0();
  //StPicoV0(StPicoTrack*, StPicoTrack*, StMuEvent*, Int_t*);
  StPicoV0(StPicoTrack*, StPicoTrack*, StPicoEvent*, Int_t*);
  void Clear(const Option_t* opt="");

  StPicoTrack *track(const Int_t i) const;
  Int_t       index2Track(const Int_t i) const;
  TVector3 momentum(const Int_t i) const;

  TVector3 v0Pos() const { return mV0Pos; }
  Float_t dcaDaughters() const { return (Float_t)mDcaDaughters/1000.; }
  Float_t   dca2Vertex() const { return (Float_t)mDca2Vtx/1000.; }
  Float_t            m() const { return mM; }

  Float_t     decayLength() const;
  TVector3 momentum() const;
  
  void setIndex2Track(const Int_t id_pos, const Int_t id_neg);
  void setParticleHypothesis(const Int_t ip_pos, const Int_t ip_neg);
  void rotateTrack(const Int_t i);
    
 protected:
  Short_t   mIndex2Track[2];
  TVector3 mMomentum[2];
  TVector3 mV0Pos;
  UShort_t       mDcaDaughters;
  UShort_t       mDca2Vtx;
  Float_t        mM;

  friend class StPicoDst;

  ClassDef(StPicoV0, 1)
};

#endif
