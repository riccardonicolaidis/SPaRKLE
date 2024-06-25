#include "SPaRKLE_Generator.hh"
#include "SPaRKLE_Globals.hh"
#include <math.h>

using namespace std;

SPaRKLE_PrimaryGenerator::SPaRKLE_PrimaryGenerator(): G4VUserPrimaryGeneratorAction(),
   fParticleSource(0)
{ 
  // PARTICLE GUN DEFINITION
  fParticleSource = new G4GeneralParticleSource();
  
  // COUNTER
  ParticleNumber = 0;

  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4String          particleName = "geantino";
  G4ParticleDefinition *particle = particleTable -> FindParticle(particleName);
  
}

SPaRKLE_PrimaryGenerator::~SPaRKLE_PrimaryGenerator()
{
  delete fParticleSource;
}

void SPaRKLE_PrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
  // PI DEFINITION
  //G4double  pi  = 3.14159265358979323846;

  // FILL ALL THE DETAILS REQUIRED FOR THE ANALYSIS
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  
  

    fParticleSource -> GeneratePrimaryVertex(anEvent);
    G4ThreeVector ParticlePosition = fParticleSource -> GetParticlePosition();
    G4ThreeVector ParticleMomentumDirection = fParticleSource -> GetParticleMomentumDirection();
    man -> FillNtupleDColumn(0, 0, 0);
    man -> FillNtupleDColumn(0, 0, fParticleSource -> GetParticleEnergy());
    man -> FillNtupleDColumn(0, 1, ParticlePosition.getX());
    man -> FillNtupleDColumn(0, 2, ParticlePosition.getY());
    man -> FillNtupleDColumn(0, 3, ParticlePosition.getZ());
    man -> FillNtupleDColumn(0, 4, ParticleMomentumDirection.getX());
    man -> FillNtupleDColumn(0, 5, ParticleMomentumDirection.getY());
    man -> FillNtupleDColumn(0, 6, ParticleMomentumDirection.getZ());
  
  // COUNTER
  if(ParticleNumber%10000 == 0)
  {
    G4cout << ParticleNumber/10000 << "e4" << G4endl;
  }

  ++ParticleNumber;
  
}
