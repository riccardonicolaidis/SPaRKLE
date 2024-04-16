#include "SPaRKLE_Event.hh"
#include "SPaRKLE_Globals.hh"

SPaRKLE_EventAction::SPaRKLE_EventAction()
{
  G4RunManager::GetRunManager()->SetPrintProgress(1e6);
}

SPaRKLE_EventAction::~SPaRKLE_EventAction()
{}



void SPaRKLE_EventAction::BeginOfEventAction(const G4Event*)
{}


void SPaRKLE_EventAction::EndOfEventAction(const G4Event* event)
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  // DEFINITIONS OF THE HISTOGRAMNS WITH COUNTS
  G4VHitsCollection* hc_Thin        = event->GetHCofThisEvent()->GetHC(0);
  G4VHitsCollection* hc_Calo        = event->GetHCofThisEvent()->GetHC(1);
  G4VHitsCollection* hc_VetoDrilled = event->GetHCofThisEvent()->GetHC(2);
  G4VHitsCollection* hc_Bottom      = event->GetHCofThisEvent()->GetHC(3);

  // DEFINITIONS OF THE VECTORS
  G4int Nx = NX_SENSORS;
  G4int Ny = NY_SENSORS;


  G4double Ed_Thin [N_TOTAL_SENSORS];
  G4double Ed_Calo[N_PL_SCINT_NO_VETO];
  G4double Ed_VetoDrilled = 0.;
  G4double Ed_VetoBottom = 0.;

  G4double EdepInside = 0.;

  G4int jCopyNo = 0;

  for( G4int ix = 0; ix < Nx; ++ix)
  {
    for( G4int iy = 0; iy < Ny; ++iy)
    {
      Ed_Thin [jCopyNo] = 0.;
      ++jCopyNo;
    }
  }

  for(G4int i = 0; i < (N_PL_SCINT_NO_VETO); ++i)
  {
    Ed_Calo[i] = 0.;
  }

  // ENERGY MEASUREMENTS

  G4int ReplicaNo;
  
  for(G4int iH=0;iH < (hc_Thin->GetSize());++iH){
    HitClass* hit=static_cast<HitClass*>(hc_Thin->GetHit(iH));
    ReplicaNo = hit-> GetReplicaNb();
    Ed_Thin[ReplicaNo]+=hit->GetEdep(); //adding the energies of the steps inside each detector, identified with chamber number
    EdepInside += hit->GetEdep();
  } 


  for(G4int iH=0;iH < (hc_Calo->GetSize());++iH){
    HitClass* hit=static_cast<HitClass*>(hc_Calo->GetHit(iH));
    ReplicaNo = hit-> GetReplicaNb();
    Ed_Calo[ReplicaNo]+=hit->GetEdep(); //adding the energies of the steps inside each detector, identified with chamber number
  }

  
  for(G4int iH=0;iH < (hc_VetoDrilled->GetSize());++iH){
    HitClass* hit=static_cast<HitClass*>(hc_VetoDrilled->GetHit(iH));
    Ed_VetoDrilled+=hit->GetEdep(); //adding the energies of the steps inside each detector, identified with chamber number
  }

  for(G4int iH=0;iH < (hc_Bottom->GetSize());++iH){
    HitClass* hit=static_cast<HitClass*>(hc_Bottom->GetHit(iH));
    Ed_VetoBottom+=hit->GetEdep(); //adding the energies of the steps inside each detector, identified with chamber number
  }

  G4bool HasEdep = false; // set it to a default value of true to save all the events
  
  HasEdep = (Ed_VetoBottom > 0. || Ed_VetoDrilled > 0.) ? true : false;
  
  if(!HasEdep) {
    for(G4int counter = 0; counter < N_TOTAL_SENSORS; counter++ ){
      if(Ed_Thin[counter]== 0.) continue;
      HasEdep = true;
      break;
    }
  }

  if(!HasEdep) {
    for(G4int counter = 0; counter < N_PL_SCINT_NO_VETO; counter++ ){
      if(Ed_Calo[counter] == 0.) continue;
      HasEdep = true;
      break;
    }
  }

  // HasEdep = true;

  if(HasEdep){

    // FILLING THE NTUPLE
    G4PrimaryVertex* primaryV = event->GetPrimaryVertex();
    G4PrimaryParticle* primary = event->GetPrimaryVertex()->GetPrimary();

    man -> FillNtupleDColumn(0, 0, 0);
    man -> FillNtupleDColumn(0, 1, primary->GetKineticEnergy());  // Energy of the primary particle
    man -> FillNtupleDColumn(0, 2, (primaryV->GetPosition())[0]); // x component of the position
    man -> FillNtupleDColumn(0, 3, (primaryV->GetPosition())[1]); // y component of the position
    man -> FillNtupleDColumn(0, 4, (primaryV->GetPosition())[2]); // z component of the position
    man -> FillNtupleDColumn(0, 5, (primary->GetMomentum())[0]);  // x component of the momentum
    man -> FillNtupleDColumn(0, 6, (primary->GetMomentum())[1]);  // y component of the momentum
    man -> FillNtupleDColumn(0, 7, (primary->GetMomentum())[2]);  // z component of the momentum

    G4int CurrentTuple = 8;
    for(G4int i = 0; i < (N_PL_SCINT_NO_VETO);++i) {
      man -> FillNtupleDColumn(0, CurrentTuple, Ed_Calo[i]);
      ++CurrentTuple;
    }
    
    man -> FillNtupleDColumn(0, CurrentTuple++, Ed_VetoDrilled);
    man -> FillNtupleDColumn(0, CurrentTuple++, Ed_VetoBottom);

    //  THRESHOLDS
    //G4double Eth_Veto = 1*keV;
    //G4double Eth_VetoDrilled = 1*keV;
        
    G4bool AtLeastOne = false;

    // nTupleID 0 = particle gun
    jCopyNo = 0;
    for( G4int ix = 0; ix < Nx; ++ix){
      for( G4int iy = 0; iy < Ny; ++iy){
        man -> FillNtupleDColumn(0, (CurrentTuple+1)+jCopyNo, Ed_Thin[jCopyNo]);

        if((Ed_Thin[jCopyNo]>0)){
          AtLeastOne = true;
          man -> FillNtupleDColumn(0, CurrentTuple, jCopyNo);
        }
        else man -> FillNtupleDColumn(0, CurrentTuple, -1);
        ++jCopyNo;
      }
    }

    man -> AddNtupleRow(0);
  }

  

  G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if (printModulo==0 || event->GetEventID() % printModulo != 0) return;
  G4cout << ">>> Event " << event->GetEventID() << G4endl;


  G4cout << G4endl
    << ">>> Event " << event->GetEventID() << " >>> Simulation truth : "
    << event->GetPrimaryVertex()->GetPrimary()->GetG4code()->GetParticleName()
    << " PDG code " << event->GetPrimaryVertex()->GetPrimary()->GetPDGcode()
    << " energy " << event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy()
    << " momentum vector " << event->GetPrimaryVertex()->GetPrimary()->GetMomentum()
	  << " position " << event->GetPrimaryVertex(0)->GetPosition()
	<< G4endl;
    
}
