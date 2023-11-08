#ifndef EVENT_HH
#define EVENT_HH

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


#include "SPaRKLE_Analysis.hh"

#include "SPaRKLE_Hit.hh"
#include "SPaRKLE_Run.hh"

class SPaRKLE_EventAction : public G4UserEventAction
{
public:
        SPaRKLE_EventAction();
        ~SPaRKLE_EventAction();

        virtual void BeginOfEventAction(const G4Event*);
        virtual void EndOfEventAction(const G4Event*);

private:


};


#endif