#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "SPaRKLE_Hit.hh"
#include "SPaRKLE_Analysis.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;


class SPaRKLE_SensitiveDetector : public G4VSensitiveDetector
{
public:
    SPaRKLE_SensitiveDetector(const G4String& name);
    virtual ~SPaRKLE_SensitiveDetector();

       // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

private:

    HitsCollection* fHitsCollection;

};

#endif