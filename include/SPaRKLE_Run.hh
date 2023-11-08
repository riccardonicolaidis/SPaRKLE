#ifndef RUN_HH
#define RUN_HH


#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "G4GenericMessenger.hh"


#include "SPaRKLE_Construction.hh"
#include "SPaRKLE_Generator.hh"
#include "SPaRKLE_Analysis.hh"

class SPaRKLE_DetectorConstruction;
class SPaRKLE_PrimaryGenerator;

class SPaRKLE_RunAction : public G4UserRunAction
{
public:
    SPaRKLE_RunAction(SPaRKLE_DetectorConstruction*, SPaRKLE_PrimaryGenerator*);
    ~SPaRKLE_RunAction();

    virtual void BeginOfRunAction(const G4Run* );
    virtual void EndOfRunAction(const G4Run* );

private:
    SPaRKLE_DetectorConstruction *fDetector;
    SPaRKLE_PrimaryGenerator *fPrimary;
    G4GenericMessenger *fMessenger;
    G4String TotalFileName;
    G4String TotalFileNameFinal;

};

#endif