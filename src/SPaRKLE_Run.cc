#include "SPaRKLE_Run.hh"
#include "SPaRKLE_Globals.hh"



SPaRKLE_RunAction::SPaRKLE_RunAction(SPaRKLE_DetectorConstruction* det, SPaRKLE_PrimaryGenerator* prim)
  : G4UserRunAction(),
    fDetector(det), fPrimary(prim)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man -> SetVerboseLevel(1);

    fMessenger = new G4GenericMessenger(this, "/NameOfFile/","Name of the file to save data");
    fMessenger -> DeclareProperty("NameOfFile", TotalFileName, "Name of the file to save data");


    G4int jCopyNo = 0;

    // Ntuple particle generator
    man -> CreateNtuple("Edep","Edep");
    man -> CreateNtupleDColumn("Egen"); // 0       
    man -> CreateNtupleDColumn("Xgen"); // 1
    man -> CreateNtupleDColumn("Ygen"); // 2
    man -> CreateNtupleDColumn("Zgen"); // 3
    man -> CreateNtupleDColumn("pDirX"); // 4
    man -> CreateNtupleDColumn("pDirY"); // 5
    man -> CreateNtupleDColumn("pDirZ"); // 6

    man -> CreateNtupleDColumn("Calo_A1"); // 7
    man -> CreateNtupleDColumn("Calo_A2"); // 8
    man -> CreateNtupleDColumn("Calo_B1"); // 9
    man -> CreateNtupleDColumn("Calo_B2"); // 10

    man -> CreateNtupleDColumn("SD1"); // 11
    man -> CreateNtupleDColumn("SD2"); // 12

    man -> CreateNtupleDColumn("VT"); // 13
    man -> CreateNtupleDColumn("VB"); // 14
    man -> CreateNtupleDColumn("VL0"); // 15
    man -> CreateNtupleDColumn("VL1"); // 16
    man -> CreateNtupleDColumn("VL2"); // 17
    man -> CreateNtupleDColumn("VL3"); // 18
            


}


SPaRKLE_RunAction::~SPaRKLE_RunAction()
{}

void SPaRKLE_RunAction::BeginOfRunAction(const G4Run* run)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();    


    TotalFileNameFinal = TotalFileName + ".root";

    man -> OpenFile(TotalFileNameFinal);
}

void SPaRKLE_RunAction::EndOfRunAction(const G4Run* )
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man -> Write();
    man -> CloseFile(TotalFileNameFinal);
}