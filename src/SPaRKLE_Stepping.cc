#include "SPaRKLE_Stepping.hh"


SPaRKLE_SteppingAction::SPaRKLE_SteppingAction(SPaRKLE_DetectorConstruction* det, SPaRKLE_EventAction* event)
: G4UserSteppingAction(), fDetector(det), fEventAction(event)
{}

SPaRKLE_SteppingAction::~SPaRKLE_SteppingAction()
{}


void SPaRKLE_SteppingAction::UserSteppingAction(const G4Step *step)
{}