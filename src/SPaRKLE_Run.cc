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

    G4int Nx = NX_SENSORS;
    G4int Ny = NY_SENSORS;

    G4int jCopyNo = 0;

    // Ntuple particle generator
    man -> CreateNtuple("Edep","Edep");
    man -> CreateNtupleDColumn("RandNumber");       // 0
    man -> CreateNtupleDColumn("RandEnergy");       // 1
    man -> CreateNtupleDColumn("Xgen");             // 2
    man -> CreateNtupleDColumn("Ygen");             // 3 
    man -> CreateNtupleDColumn("Zgen");             // 4 
    man -> CreateNtupleDColumn("pDirX");            // 5 
    man -> CreateNtupleDColumn("pDirY");            // 6 
    man -> CreateNtupleDColumn("pDirZ");            // 7

    
            
    G4String NameTupleColumn;
    G4String ThickName;
    G4String ThinName;



    for(int i = 0; i < (N_PL_SCINT_NO_VETO); ++i)
    {
        NameTupleColumn = "Ed_Calo" + std::to_string(i);
        man -> CreateNtupleDColumn(NameTupleColumn);          // 8
    }
    man -> CreateNtupleDColumn("Ed_DrilledVeto");   // 9
    man -> CreateNtupleDColumn("Ed_BottomVeto");    // 10


    man -> CreateNtupleDColumn("NumberID");         // 11


    for(G4int ix = 0; ix < Nx ; ++ix)
    {
        for(G4int iy = 0; iy < Ny ; ++iy)
        {           
            NameTupleColumn = "x"+std::to_string(ix)+"_y"+std::to_string(iy)+"_ID"+std::to_string(jCopyNo);
            ThinName = "Thin_"+NameTupleColumn;            
            man -> CreateNtupleDColumn(ThinName); 
            ++jCopyNo;
        }
    }
}


SPaRKLE_RunAction::~SPaRKLE_RunAction()
{}

void SPaRKLE_RunAction::BeginOfRunAction(const G4Run* run)
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();    
    G4int TkThick = TK_THICK/um;


    TotalFileNameFinal = TotalFileName + "_" + std::to_string(TkThick) +"_um.root";

    man -> OpenFile(TotalFileNameFinal);
}

void SPaRKLE_RunAction::EndOfRunAction(const G4Run* )
{
    G4AnalysisManager *man = G4AnalysisManager::Instance();
    man -> Write();
    man -> CloseFile(TotalFileNameFinal);
}