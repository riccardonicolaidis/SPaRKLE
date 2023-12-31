#include "SPaRKLE_Detector.hh"

SPaRKLE_SensitiveDetector::SPaRKLE_SensitiveDetector(const G4String& name) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(name);
}

SPaRKLE_SensitiveDetector::~SPaRKLE_SensitiveDetector()
{}

void SPaRKLE_SensitiveDetector::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new HitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}



G4bool SPaRKLE_SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist)
{
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep==0.) return false;

  HitClass* newHit = new HitClass();

  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
                                               ->GetCopyNumber());
  newHit->SetEdep(edep);
  newHit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
  newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());

  fHitsCollection->insert( newHit );

  //newHit->Print();

  return true;
}



void SPaRKLE_SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     G4int nofHits = fHitsCollection->entries();
     G4cout
       << G4endl 
       << "-------->Hits Collection: in this event they are " << nofHits 
       << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

