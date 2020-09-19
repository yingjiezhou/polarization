/**
 * \class StPicoBToPidTraits
 * \brief Hold information about BTOF-matched tracks
 *
 * The class stores information about tracks that matched
 * the Barrel Time-of-Flight detector
 */

#ifndef StPicoBTofPidTraits_h
#define StPicoBTofPidTraits_h

// C++ headers
#include <limits>

// ROOT headers
#include "TObject.h"
#include "TVector3.h"

//_________________
class StPicoBTofPidTraits : public TObject {

 public:
  /// Default constructor
  StPicoBTofPidTraits();
  /// Copy constructor
  StPicoBTofPidTraits(const StPicoBTofPidTraits &traits);
  /// Destructor
  virtual ~StPicoBTofPidTraits();
  /// Print TOF PID traits information
  virtual void Print(const Char_t* option = "") const;

  //
  // Getters
  //

  /// Return index of the assiciated track
  Int_t   trackIndex() const;
  /// Return BTOF cell ID (encoding = (tray-1)*192+(module-1)*6+(cell-1): -1 - no match )
  Int_t   btofCellId() const;
  /// Return matching flag (0 - no match, 1 - one-to-one, 2 - one-to-multiple)
  Int_t   btofMatchFlag() const;

  /// Return time of flight
  Float_t btof() const;
  /// Return beta (compression = beta * 20000)
  Float_t btofBeta() const;
  /// Return yLocal (compression = yLocal * 1000)
  Float_t btofYLocal() const;
  /// Return zLocal (compression = zLocal * 1000)
  Float_t btofZLocal() const;

  /// Return hit position
  TVector3 btofHitPos() const;

  /// Return x comonent of hit position
  Float_t btofHitPosX() const;
  /// Return y comonent of hit position
  Float_t btofHitPosY() const;
  /// Return z comonent of hit position
  Float_t btofHitPosZ() const;

  //
  // Setters
  //

  /// Set assiciated track index
  void setTrackIndex(Int_t inx2PicoTrack);
  /// Set TOF cell ID which encodes tray, module and cell IDs
  void setBTofCellId(Int_t tray, Int_t module, Int_t cell);
  /// Set TOF-matching flag
  void setBTofMatchFlag(UChar_t flag);
  /// Set time of flight
  void setTOF(Float_t tof);
  /// Set beta
  void setBeta(Float_t beta);
  /// Set hit position (x,y,z)
  void setHitPositionXYZ(Float_t x, Float_t y, Float_t z);
  /// Set hit position x (cm)
  void setHitPositionX(Float_t x);
  /// Set hit position y (cm)
  void setHitPositionY(Float_t y);
  /// Set hit position z (cm)
  void setHitPositionZ(Float_t z);
  /// Set yLocal
  void setYLocal(Float_t yLocal);
  /// Set zLocal
  void setZLocal(Float_t zLocal);

 private:

  /// Index to the associated picoTrack in the event
  Short_t  mTrackIndex;
  /// CellId encodes: (tray-1)*192+(module-1)*6+(cell-1): -1 - no match
  Short_t   mBTofCellId;
  /// 0 - no match, 1 - one-to-one, 2 - one-to-multiple
  UChar_t   mBTofMatchFlag;
  /// Time-Of-Flight
  Float_t mBTof;
  /// Beta * 20000
  UShort_t  mBTofBeta;
  /// ylocal * 1000
  Short_t   mBTofYLocal;
  /// zlocal * 1000
  Short_t   mBTofZLocal;
  /// Hit position projected on X plane (compression = position * 100)
  Short_t   mBTofHitPosX;
  /// Hit position projected on Y plane (compression = position * 100)
  Short_t   mBTofHitPosY;
  /// Hit position projected on Z plane (compression = position * 100)
  Short_t   mBTofHitPosZ;

  ClassDef(StPicoBTofPidTraits, 3);
};

//
// Getters
//
inline Int_t   StPicoBTofPidTraits::trackIndex() const { return mTrackIndex; }
inline Int_t   StPicoBTofPidTraits::btofCellId() const { return (Int_t)mBTofCellId; }
inline Int_t   StPicoBTofPidTraits::btofMatchFlag() const { return (Int_t)mBTofMatchFlag; }
inline Float_t StPicoBTofPidTraits::btof() const { return mBTof; }
inline Float_t StPicoBTofPidTraits::btofBeta() const { return (Float_t)mBTofBeta / 20000.f; }
inline Float_t StPicoBTofPidTraits::btofYLocal() const { return (Float_t)mBTofYLocal / 1000.; }
inline Float_t StPicoBTofPidTraits::btofZLocal() const { return (Float_t)mBTofZLocal / 1000.; }
inline TVector3 StPicoBTofPidTraits::btofHitPos() const
{ return TVector3( btofHitPosX() , btofHitPosY() , btofHitPosZ() ); }
inline Float_t StPicoBTofPidTraits::btofHitPosX() const { return (Float_t)mBTofHitPosX / 100.; }
inline Float_t StPicoBTofPidTraits::btofHitPosY() const { return (Float_t)mBTofHitPosY / 100.; }
inline Float_t StPicoBTofPidTraits::btofHitPosZ() const { return (Float_t)mBTofHitPosZ / 100.; }

//
// Setters
//
inline void StPicoBTofPidTraits::setTrackIndex(Int_t idx2PicoTrack) 
{ mTrackIndex = (idx2PicoTrack > std::numeric_limits<short>::max()) ? -1 : (Short_t)idx2PicoTrack; }
inline void StPicoBTofPidTraits::setTOF(Float_t tof) { mBTof = tof; }
inline void StPicoBTofPidTraits::setBTofCellId(Int_t tray, Int_t module, Int_t cell) 
{ mBTofCellId  = (Short_t)((tray - 1) * 192 + (module - 1) * 6 + (cell - 1)); }
inline void StPicoBTofPidTraits::setBTofMatchFlag(UChar_t flag) { mBTofMatchFlag = flag; }

#endif
