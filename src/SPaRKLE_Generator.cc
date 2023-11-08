#include "SPaRKLE_Generator.hh"
#include "SPaRKLE_Globals.hh"
#include <math.h>

using namespace std;

SPaRKLE_PrimaryGenerator::SPaRKLE_PrimaryGenerator(): G4VUserPrimaryGeneratorAction(),
   fParticleSource(0)
{ 
  // PARTICLE GUN DEFINITION
  fParticleSource = new G4GeneralParticleSource();
  fParticleGun = new G4ParticleGun();
  
  // COUNTER
  ParticleNumber = 0;

  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4String          particleName = "geantino";
  G4ParticleDefinition *particle = particleTable -> FindParticle(particleName);


  fParticleGun -> SetParticlePosition(G4ThreeVector(0.,0.,0.));
  fParticleGun -> SetParticleDefinition(particle);
  
  // DEFINE THE USER DEFINED MESSENGERS
  fMessenger_1 = new G4GenericMessenger(this, "/particleEnergyRandom/","Range Energy Random");
  fMessenger_1 -> DeclareProperty("Ek_min", Ek_min, "Minimum kinetic energy in MeV");
  fMessenger_1 -> DeclareProperty("Ek_max", Ek_max, "Maximum kinetic energy in MeV");


  // INDEXING FOR THE DETECTPR
  // ONLY IF DEBUGGING MODE IS ON
  Ix = 0;
  Iy = 0;

  // DEBUGGING MODE 
  DebuggingModeIsOn = false;
  fMessenger_1 -> DeclareProperty("DebbuggingModeIsOn", DebuggingModeIsOn, "Activate the debugging mode");
  IsExponential = false;
  fMessenger_1 -> DeclareProperty("IsExponential", IsExponential, "Activate the exponential distribution for Energy (Works only if the Debugging mode is on)");
  
}

SPaRKLE_PrimaryGenerator::~SPaRKLE_PrimaryGenerator()
{
  delete fParticleSource;
  delete fParticleGun;
}

void SPaRKLE_PrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
  // PI DEFINITION
  //G4double  pi  = 3.14159265358979323846;

  // FILL ALL THE DETAILS REQUIRED FOR THE ANALYSIS
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  
  
  // TEST CONFIGURATION comment if necessary
  if(DebuggingModeIsOn)
  {
    LActive            = 14.*cm;
    zgen               =  -2.5*cm;
    NxHoles            = NX_SENSORS;
    NyHoles            = NY_SENSORS;
    LThin              = 5.8 * cm;
    xDelta             = LThin / (NxHoles);
    yDelta             = LThin / (NyHoles);
    xCurrentHole       = -(LThin/2.) + (xDelta/2.) + Ix * xDelta;
    yCurrentHole       = -(LThin/2.) + (yDelta/2.) + Iy * yDelta;

    thetaMax = 35.*deg;
    RMax = std::sqrt(2) * 1.5;

    G4double Delta_ii_x = NxHoles/2. -0.5;
    G4double Delta_ii_y = NyHoles/2. -0.5;


    
    phi = std::atan2((Iy-Delta_ii_y),(Ix-Delta_ii_x));
    distanceR = std::sqrt((Iy-Delta_ii_y)*(Iy-Delta_ii_y) + (Ix-Delta_ii_x)*(Ix-Delta_ii_x));
    theta = distanceR * thetaMax / RMax;

    //xCurrentHole = xCurrentHole - (abs(zgen))*std::tan(theta) * std::cos(phi);
    //yCurrentHole = yCurrentHole - (abs(zgen))*std::tan(theta) * std::sin(phi);


    //phi                = std::atan2((Iy-1.5),(Ix-1.5));
    //distanceR          = std::sqrt((Iy-1.5)*(Iy-1.5) + (Ix-1.5)*(Ix-1.5));
    //theta              = distanceR * thetaMax / RMax;
    sintheta           = std::sin(theta);
    costheta           = std::cos(theta);
    holePosition       = G4ThreeVector(xCurrentHole, yCurrentHole, 0.);
    holeDirection      = G4ThreeVector(std::cos(phi)*std::sin(theta),std::sin(phi)*std::sin(theta), std::cos(theta));
    myParticlePosition = holePosition - (zgen/std::cos(theta)) * holeDirection;
    xgen               = myParticlePosition.getX();
    ygen               = myParticlePosition.getY();
    myParticlePosition.setZ(-myParticlePosition.getZ());
    momentumDirection = G4ThreeVector(-cos(phi)*sintheta,-sin(phi)*sintheta,costheta);

    // SET MOMENTUM DIRECTION AND POSITION
    fParticleGun->SetParticleMomentumDirection(momentumDirection);
    fParticleGun->SetParticlePosition(myParticlePosition);
  
    // TAKE A RANDOM NUMBER AND GENERATE THE ENERGY
    RndNum = G4UniformRand();
    Ek_random = Ek_min * MeV + (Ek_max * MeV - Ek_min * MeV) * RndNum;

    if(IsExponential)
    {
      G4double logEkMin = std::log10(Ek_min * MeV);
      G4double logEkMax = std::log10(Ek_max * MeV);
      Ek_random = std::pow(10.,logEkMin + (logEkMax - logEkMin) * RndNum);

    }


    man -> FillNtupleDColumn(0, 0, RndNum);
    man -> FillNtupleDColumn(0, 1, Ek_random);
    man -> FillNtupleDColumn(0, 2, xgen);
    man -> FillNtupleDColumn(0, 3, ygen);
    man -> FillNtupleDColumn(0, 4, zgen);
    man -> FillNtupleDColumn(0, 5, -cos(phi)*sintheta);
    man -> FillNtupleDColumn(0, 6, -sin(phi)*sintheta);
    man -> FillNtupleDColumn(0, 7, costheta);
    // SET PARTICLE ENERGY
    fParticleGun -> SetParticleEnergy(Ek_random);
    fParticleGun  -> GeneratePrimaryVertex(anEvent);

    if(Iy < (NyHoles-1))
    {
      ++Iy;
    } else
    {
      Iy = 0;
      if(Ix < (NxHoles-1))
      {
        ++Ix;
      } 
      else 
      {
        Ix = 0;
      } 
    }
  } 
  else
  {
    fParticleSource -> GeneratePrimaryVertex(anEvent);
    G4ThreeVector ParticlePosition = fParticleSource -> GetParticlePosition();
    G4ThreeVector ParticleMomentumDirection = fParticleSource -> GetParticleMomentumDirection();
    man -> FillNtupleDColumn(0, 0, 0);
    man -> FillNtupleDColumn(0, 1, fParticleSource -> GetParticleEnergy());
    man -> FillNtupleDColumn(0, 2, ParticlePosition.getX());
    man -> FillNtupleDColumn(0, 3, ParticlePosition.getY());
    man -> FillNtupleDColumn(0, 4, ParticlePosition.getZ());
    man -> FillNtupleDColumn(0, 5, ParticleMomentumDirection.getX());
    man -> FillNtupleDColumn(0, 6, ParticleMomentumDirection.getY());
    man -> FillNtupleDColumn(0, 7, ParticleMomentumDirection.getZ());
  }
  
  // COUNTER
  if(ParticleNumber%20000 == 0)
  {
    G4cout << ParticleNumber/1000 << "e3" << G4endl;
  }
  ++ParticleNumber;
  
}
