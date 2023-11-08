#include "SPaRKLE_Event.hh"
#include "SPaRKLE_Globals.hh"

SPaRKLE_EventAction::SPaRKLE_EventAction()
{}

SPaRKLE_EventAction::~SPaRKLE_EventAction()
{}



void SPaRKLE_EventAction::BeginOfEventAction(const G4Event*)
{}


void SPaRKLE_EventAction::EndOfEventAction(const G4Event* event)
{
  G4AnalysisManager *man = G4AnalysisManager::Instance();
  //man -> FillNtupleDColumn(0,);

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


  G4int CurrentTuple = 8;
  for(G4int i = 0; i < (N_PL_SCINT_NO_VETO);++i)
  {
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
  for( G4int ix = 0; ix < Nx; ++ix)
  {
    for( G4int iy = 0; iy < Ny; ++iy)
    {
      man -> FillNtupleDColumn(0, (CurrentTuple+1)+jCopyNo, Ed_Thin[jCopyNo]);

      if((Ed_Thin[jCopyNo]>0))
      {
        AtLeastOne = true;
        man -> FillNtupleDColumn(0, CurrentTuple, jCopyNo);
      }
      else
      {
        man -> FillNtupleDColumn(0, CurrentTuple, -1);
      }
      ++jCopyNo;
    }
  }




  // && (Ed_Veto<=Eth_Veto) && (Ed_VetoDrilled<=Eth_VetoDrilled)
  //if(AtLeastOne)
  //{
  //    man -> AddNtupleRow(0);
  //}

  man -> AddNtupleRow(0);

 

}
