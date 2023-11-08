#ifndef STEPPING_HH
#define STEPPING_HH

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "SPaRKLE_Construction.hh"
#include "SPaRKLE_Event.hh"
#include "SPaRKLE_Analysis.hh"

class SPaRKLE_DetectorConstruction;
class SPaRKLE_EventAction;

class SPaRKLE_SteppingAction : public G4UserSteppingAction
{
public:
    SPaRKLE_SteppingAction(SPaRKLE_DetectorConstruction* ,SPaRKLE_EventAction *);
    ~SPaRKLE_SteppingAction();

    virtual void UserSteppingAction(const G4Step*);

private:
    SPaRKLE_EventAction *fEventAction;
    SPaRKLE_DetectorConstruction *fDetector;
};

#endif