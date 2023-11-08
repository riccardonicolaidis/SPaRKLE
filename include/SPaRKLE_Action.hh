#ifndef ACTION_HH
#define ACTION_HH

#include "G4VUserActionInitialization.hh"

#include "SPaRKLE_Generator.hh"
#include "SPaRKLE_Run.hh"
#include "SPaRKLE_Stepping.hh"
#include "SPaRKLE_Event.hh"
#include "SPaRKLE_Detector.hh"
#include "SPaRKLE_Tracking.hh"


class SPaRKLE_DetectorConstruction;


class SPaRKLE_ActionInitialization : public  G4VUserActionInitialization
{
public:
  SPaRKLE_ActionInitialization(SPaRKLE_DetectorConstruction *detector);
  virtual ~SPaRKLE_ActionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;


private:
  SPaRKLE_DetectorConstruction *fDetector;
};

#endif
