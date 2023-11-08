#include "SPaRKLE_Action.hh"

SPaRKLE_ActionInitialization::SPaRKLE_ActionInitialization(SPaRKLE_DetectorConstruction *detector) 
  : G4VUserActionInitialization(),
  fDetector(detector)
{}

SPaRKLE_ActionInitialization::~SPaRKLE_ActionInitialization()
{}

void SPaRKLE_ActionInitialization::BuildForMaster() const
{
  SPaRKLE_RunAction *runAction = new SPaRKLE_RunAction(fDetector, 0);
  SetUserAction(runAction);
}


void SPaRKLE_ActionInitialization::Build() const
{
  SPaRKLE_PrimaryGenerator *primary = new SPaRKLE_PrimaryGenerator();
  SetUserAction(primary);

  SPaRKLE_RunAction *runAction = new SPaRKLE_RunAction(fDetector, primary);
  SetUserAction(runAction);

  SPaRKLE_EventAction *eventAction = new SPaRKLE_EventAction();
  SetUserAction(eventAction);

  SPaRKLE_TrackingAction* trackingAction = new SPaRKLE_TrackingAction(fDetector);
  SetUserAction(trackingAction);

  SPaRKLE_SteppingAction *steppingAction = new SPaRKLE_SteppingAction(fDetector,eventAction);
  SetUserAction(steppingAction);

}
