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



  G4VHitsCollection* hc_Calo  = event->GetHCofThisEvent()->GetHC(0);
  G4VHitsCollection* hc_SD    = event->GetHCofThisEvent()->GetHC(1);
  G4VHitsCollection* hc_VT    = event->GetHCofThisEvent()->GetHC(2);
  G4VHitsCollection* hc_VB    = event->GetHCofThisEvent()->GetHC(3);
  G4VHitsCollection* hc_VL    = event->GetHCofThisEvent()->GetHC(4);
  

  G4double Ed_Calo_A1 = 0.;
  G4double Ed_Calo_A2 = 0.;
  G4double Ed_Calo_B1 = 0.;
  G4double Ed_Calo_B2 = 0.;

  G4double Ed_SD1 = 0.;
  G4double Ed_SD2 = 0.;

  G4double Ed_VT = 0.;
  G4double Ed_VB = 0.;

  G4double Ed_VL0 = 0.;
  G4double Ed_VL1 = 0.;
  G4double Ed_VL2 = 0.;
  G4double Ed_VL3 = 0.;


  G4int ReplicaNo;
  
  for(G4int iH=0;iH < (hc_Calo->GetSize());++iH){
    HitClass* hit=static_cast<HitClass*>(hc_Calo->GetHit(iH));
    ReplicaNo = hit-> GetReplicaNb();

    if(ReplicaNo==0)
      Ed_Calo_A1 += hit->GetEdep();
    
    if(ReplicaNo==1)
      Ed_Calo_A2 += hit->GetEdep();
    
    if(ReplicaNo==2)
      Ed_Calo_B1 += hit->GetEdep();
    
    if(ReplicaNo==3)
      Ed_Calo_B2 += hit->GetEdep();
  } 

  for(G4int iH=0;iH < (hc_SD->GetSize());++iH){
    HitClass* hit=static_cast<HitClass*>(hc_SD->GetHit(iH));
    ReplicaNo = hit-> GetReplicaNb();

    if(ReplicaNo==0)
      Ed_SD1 += hit->GetEdep();
  
    if(ReplicaNo==1)
      Ed_SD2 += hit->GetEdep();
  } 

  for(G4int iH=0;iH < (hc_VT->GetSize());++iH){
    HitClass* hit=static_cast<HitClass*>(hc_VT->GetHit(iH));
    Ed_VT += hit->GetEdep();    
  } 

  for(G4int iH=0;iH < (hc_VB->GetSize());++iH){
    HitClass* hit=static_cast<HitClass*>(hc_VB->GetHit(iH));
    Ed_VB += hit->GetEdep();    
  }

  for(G4int iH=0;iH < (hc_VL->GetSize());++iH){
    HitClass* hit=static_cast<HitClass*>(hc_VL->GetHit(iH));
    ReplicaNo = hit-> GetReplicaNb();

    if(ReplicaNo==0)
      Ed_VL0 += hit->GetEdep();
    
    if(ReplicaNo==1)
      Ed_VL1 += hit->GetEdep();
    
    if(ReplicaNo==2)
      Ed_VL2 += hit->GetEdep();
    
    if(ReplicaNo==3)
      Ed_VL3 += hit->GetEdep();
  }










  
  man -> FillNtupleDColumn(0, 7, Ed_Calo_A1);
  man -> FillNtupleDColumn(0, 8, Ed_Calo_A2);
  man -> FillNtupleDColumn(0, 9, Ed_Calo_B1);
  man -> FillNtupleDColumn(0, 10, Ed_Calo_B2);

  man -> FillNtupleDColumn(0, 11, Ed_SD1);
  man -> FillNtupleDColumn(0, 12, Ed_SD2);

  man -> FillNtupleDColumn(0, 13, Ed_VT);
  man -> FillNtupleDColumn(0, 14, Ed_VB);

  man -> FillNtupleDColumn(0, 15, Ed_VL0);
  man -> FillNtupleDColumn(0, 16, Ed_VL1);
  man -> FillNtupleDColumn(0, 17, Ed_VL2);
  man -> FillNtupleDColumn(0, 18, Ed_VL3);




  // && (Ed_Veto<=Eth_Veto) && (Ed_VetoDrilled<=Eth_VetoDrilled)
  //if(AtLeastOne)
  //{
  //    man -> AddNtupleRow(0);
  //}

  man -> AddNtupleRow(0);

 

}
