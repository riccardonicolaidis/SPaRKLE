#include "SPaRKLE_Tracking.hh"

SPaRKLE_TrackingAction::SPaRKLE_TrackingAction(SPaRKLE_DetectorConstruction* det)
:G4UserTrackingAction(), fDetector(det)
{}

SPaRKLE_TrackingAction::~SPaRKLE_TrackingAction()
{}


void SPaRKLE_TrackingAction::PreUserTrackingAction(const G4Track* track)
{}


void SPaRKLE_TrackingAction::PostUserTrackingAction(const G4Track* )
{}