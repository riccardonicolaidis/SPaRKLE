#ifndef TRACKING_HH
#define TRACKING_HH

#include "G4UserTrackingAction.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4HadronicProcessType.hh"
#include "G4SystemOfUnits.hh"

#include "SPaRKLE_Construction.hh"
#include "SPaRKLE_Analysis.hh"

class SPaRKLE_DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SPaRKLE_TrackingAction : public G4UserTrackingAction {

  public:  
    SPaRKLE_TrackingAction(SPaRKLE_DetectorConstruction*);
   ~SPaRKLE_TrackingAction();
   
    virtual void  PreUserTrackingAction(const G4Track*);   
    virtual void PostUserTrackingAction(const G4Track*);
    
  private:
    SPaRKLE_DetectorConstruction* fDetector;
};




#endif