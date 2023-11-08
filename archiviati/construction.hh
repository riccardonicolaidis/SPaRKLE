#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "myglobals.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4GenericMessenger.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4GenericMessenger.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4GDMLParser.hh"
#include "G4MaterialTable.hh"
#include <cmath>
#include "G4VisAttributes.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"


#include "Analysis.hh"
#include "detector.hh"


#if NEW_GEOMETRY == 0



class LEM_DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  LEM_DetectorConstruction();
  ~LEM_DetectorConstruction();

  virtual G4VPhysicalVolume *Construct();

private:

  virtual void ConstructSDandField();


  // World
  G4Box             *solidWorld;
  G4LogicalVolume   *logicWorld;
  G4VPhysicalVolume *physWorld;
  G4VisAttributes   *visWorld;

  // Container
  G4Box             *solidContainer;
  G4LogicalVolume   *logicContainer;
  G4VPhysicalVolume *physContainer;
  G4VisAttributes   *visContainer;

  // Thin detector
  G4Tubs            *solidSiDetThin;
  G4LogicalVolume   *logicSiDetThin;
  G4VPhysicalVolume *physSiDetThin;
  G4VisAttributes   *visSiDetThin;

  // SILICON FRAME ALUMINIUM
  G4Tubs            *solidFrameCircle;
  G4Tubs            *solidConnector;
  G4UnionSolid      *solidFrame[2][2]; // 4 possible orientations
  G4LogicalVolume   *logicFrame[2][2];
  G4VPhysicalVolume *physFrame;
  G4VisAttributes   *visFrame[2][2];
  
  // Thick detector
  G4Tubs            *solidSiDetThick;
  G4LogicalVolume   *logicSiDetThick;
  G4VPhysicalVolume *physSiDetThick;
  G4VisAttributes   *visSiDetThick;

  // Scintillator Veto
  G4Box             *solidScintVeto;
  G4LogicalVolume   *logicScintVeto;
  G4VPhysicalVolume *physScintVeto; 
  G4VisAttributes   *visScintVeto;
  G4LogicalSkinSurface *logicSurfScintVeto[(N_PL_SCINT_NO_VETO+1)];

  // DRILLED VETO
  // A: Front Side Drilled
  G4Box              *solidDrilledVeto0;
  G4Box              *solidDrilledVeto1Ext; // External
  G4Box              *solidDrilledVeto1Int; // Internal
  G4SubtractionSolid *solidDrilledVeto1;    // Subtraction
  G4UnionSolid       *solidDrilledVeto;
  G4LogicalVolume    *logicDrilledVeto;
  G4VPhysicalVolume  *physDrilledVeto;
  G4SubtractionSolid *solidFinalDrilledVeto;
  G4VisAttributes    *visDrilledVeto;
  G4LogicalSkinSurface *logicSurfDrilledVeto;

  G4OpticalSurface   *opsurfVeto;


  // DRILLED ALUMINIUM
  G4Box              *solidDrilledAl0;
  G4Box              *solidDrilledAl1Ext; // External
  G4Box              *solidDrilledAl1Int; // Internal
  G4SubtractionSolid *solidDrilledAl1;    // Subtraction
  G4UnionSolid       *solidDrilledAl;
  G4LogicalVolume    *logicDrilledAl;
  G4VPhysicalVolume  *physDrilledAl;
  G4SubtractionSolid *solidFinalDrilledAl;
  G4VisAttributes    *visDrilledAl;


  // HOLE 
  G4Tubs             *solidHole;
  G4Tubs             *solidHoleLong;

  void DefineMaterials();

  // Materials
  G4Material *worldMat, *EJ200, *SiMat, *CdZnTe, *Al, *bachelite;
  G4Element *Si;
  
  // World
  G4double xWorld, yWorld, zWorld;
  // Container
  G4double xContainer, yContainer, zContainer;
  G4double LxContainer, LyContainer, LzContainer; 

  // Thin detector 
  // Active layer zone
  G4double xThin, yThin, zThin; // Position of the active layer
  G4double TkThin, LThin;       // Thickness (100 um, 300 um) and size of acyive layer
  // Circular detector DATA
  G4double rExtAlDet, rIntAlDet, TkAlDet; // Aluminum frame details
  G4double rThinDet;                      // Radius of the SILICON detector
  G4double LConnector;                    // Length of the connector
  G4double LTranslationConnector;         // Length of the connector translation
  G4double rConnector;                    // Radius of the connector
  
  // Thick Detector
  G4double xThick, yThick, zThick; // Position of the active layer
  G4double TkThick, LThick; // Size and thickness of the active layer (500 um 2 mm)
  G4double lBoxThick;       // Size of the box for CZT square

  // Scintillator Veto
  G4double xScintVeto, yScintVeto, zScintVeto;
  G4double TkScintVeto, LScintVeto;

  // Drilled Veto
  // 0 is the front drilled face
  // 1 is the external box
  G4double LxDrilledVeto0, LyDrilledVeto0;
  G4double xDrilledVeto, yDrilledVeto, zDrilledVeto;
  G4double TkDrilledVeto0;
  G4double LxDrilledVeto1, LyDrilledVeto1;
  G4double LzDrilledVeto1;
  G4double TkCompenetration;


  // Drilled Al
  // 0 is the front drilled face
  // 1 is the external box
  G4double LxDrilledAl0, LyDrilledAl0;
  G4double xDrilledAl, yDrilledAl, zDrilledAl;
  G4double TkDrilledAl0;
  G4double LxDrilledAl1, LyDrilledAl1;
  G4double LzDrilledAl1;

  
  // USEFUL QUANTITIES 
  G4double xCurrentHole, yCurrentHole;
  G4double xCurrentHoleThick, yCurrentHoleThick;
  G4double xDelta, yDelta;

  G4int NxHoles;
  G4int NyHoles;

  G4double rHoles;
  G4double LengthHole;

  G4double theta;
  G4double phi;
  G4double distanceR;
  G4double thetaMax;
  G4double RMax;
  G4int CopyNo;

  G4double TotalMass;





};

#endif
#endif
