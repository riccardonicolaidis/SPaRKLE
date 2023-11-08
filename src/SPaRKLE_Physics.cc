#include "SPaRKLE_Physics.hh"

SPaRKLE_PhysicsList::SPaRKLE_PhysicsList()
{
  RegisterPhysics (new G4EmStandardPhysics());
  RegisterPhysics (new G4OpticalPhysics());  
  RegisterPhysics (new G4DecayPhysics());
  RegisterPhysics (new G4RadioactiveDecayPhysics());
}

SPaRKLE_PhysicsList::~SPaRKLE_PhysicsList()
{}
