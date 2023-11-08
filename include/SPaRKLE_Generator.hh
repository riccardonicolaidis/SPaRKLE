#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4Geantino.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"


#include "SPaRKLE_Analysis.hh"


class SPaRKLE_PrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
  SPaRKLE_PrimaryGenerator();
  ~SPaRKLE_PrimaryGenerator();

  virtual void GeneratePrimaries(G4Event*);

private:
  G4GeneralParticleSource *fParticleSource;
  G4ParticleGun *fParticleGun;
  G4int ParticleNumber;
  G4GenericMessenger *fMessenger_1;

  G4double Ek_min, Ek_max, RndNum, Ek_random; 
  G4bool IsExponential;

  G4int Ix, Iy;

  G4bool DebuggingModeIsOn;

  G4double LActive;

  G4double xgen, ygen, zgen;
  G4double cos2theta;
  G4double costheta;
  G4double sintheta;
  G4double phi;
  G4ThreeVector myParticlePosition;


  G4int NxHoles;
  G4int NyHoles;

  G4double LThin;
  G4double xDelta;
  G4double yDelta;
  G4double xCurrentHole;
  G4double yCurrentHole;
  G4double thetaMax;
  G4double RMax;
  G4double distanceR;
  G4double theta;
  G4ThreeVector holePosition;
  G4ThreeVector holeDirection;
  G4ThreeVector momentumDirection;

};

#endif
