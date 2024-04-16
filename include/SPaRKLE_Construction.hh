#ifndef CONSTRUCTION_NEW
#define CONSTRUCTION_NEW

#include "SPaRKLE_Globals.hh"

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
// #include "G4GDMLParser.hh"
#include "G4MaterialTable.hh"
#include <cmath>
#include "G4VisAttributes.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"


#include "SPaRKLE_Analysis.hh"
#include "SPaRKLE_Detector.hh"


#if NEW_GEOMETRY == 1



class SPaRKLE_DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  SPaRKLE_DetectorConstruction();
  ~SPaRKLE_DetectorConstruction();

  virtual G4VPhysicalVolume *Construct();

  void ClearSD();

  G4UnionSolid * SolidoVaschetta(G4double hInternal, G4double TkBottom, G4double xIn, G4double yIn, G4double xOut, G4double);

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


  // Scintillator Veto
  G4Box             *solidScintVeto;
  G4LogicalVolume   *logicScintVeto;
  G4VPhysicalVolume *physScintVeto; 
  G4VisAttributes   *visScintVeto;

  // DRILLED VETO
  // A: Front Side Drilled
  G4Box              *solidDrilledVeto;
  G4LogicalVolume    *logicDrilledVeto;
  G4VPhysicalVolume  *physDrilledVeto;
  G4SubtractionSolid *solidFinalDrilledVeto;
  G4VisAttributes    *visDrilledVeto;



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
  G4Material *worldMat, *CsI, *EJ200, *SiMat, *Al, *bachelite;
  G4Element *Si;
  
  G4Material *GAGG, *BGO, *KaptonMat;


  // World
  G4double xWorld, yWorld, zWorld;
  // Container
  G4double xContainer, yContainer, zContainer;
  G4double LxContainer, LyContainer, LzContainer; 

  // Thin detector 
  // Active layer zone
  G4double xThin, yThin, zThin; // Position of the active layer
  G4double TkThin, LThin;       // Thickness (100 um, 300 um) and size of acyive layer
  

  // Define the size of the active layer. 
  // The size of the square of the Thick detectors is a constraint!!!
  G4double LActive, LSquareCentersThin;




  // Circular detector DATA
  G4double rThinDet;                      // Radius of the SILICON detector

  // Scintillator Veto
  G4double xScintVeto, yScintVeto, zScintVeto;
  G4double TkScintVeto, LScintVeto;

  // Drilled Veto
  // 0 is the front drilled face
  // 1 is the external box
  G4double LxDrilledVeto, LyDrilledVeto, LzDrilledVeto;
  G4double xDrilledVeto, yDrilledVeto, zDrilledVeto;
  G4double TkDrilledVeto;
  G4double TkCompenetration;


  // Drilled Al
  // 0 is the front drilled face
  // 1 is the external box
  G4double LxDrilledAl0, LyDrilledAl0;
  G4double xDrilledAl, yDrilledAl, zDrilledAl;
  G4double TkDrilledAl0;
  G4double LxDrilledAl1, LyDrilledAl1;
  G4double LzDrilledAl1;


  // PCB
  G4double TkPCBThin;
  G4double LxPCB, LyPCB;


  // PCB THIN DETECTOR
  G4Box              *solidPCBThin;
  G4SubtractionSolid *solidPCBThinFinal;
  G4LogicalVolume    *logicPCBThin;
  G4VPhysicalVolume  *physPCBThin;
  G4VisAttributes    *visPCBThin;


  // HOLES IN THE PCB
  G4Tubs             *solidHolePCB;

  // LATERAL VETO
  G4Box              *solidLateralVeto;
  G4LogicalVolume    *logicLateralVeto;
  G4VPhysicalVolume  *physLateralVeto;
  G4VisAttributes    *visLateralVeto;


  G4UnionSolid *VaschettaVeto;
  G4LogicalVolume *logicVaschettaVeto;
  G4VPhysicalVolume *physVaschettaVeto;
  G4VisAttributes *visVaschettaVeto;


  
  // USEFUL QUANTITIES 
  G4double xCurrentHole, yCurrentHole;
  G4double xDelta, yDelta;


  G4double LxLateralVeto, LzLateralVeto, TkLateralVeto;

  G4double ExtraSpacing;
  G4double SpacingPCBComponents;

  G4int NxHoles;
  G4int NyHoles;

  G4double rHoles;
  G4double rHolesPCB;
  G4double LengthHole;

  G4double theta;
  G4double phi;
  G4double distanceR;
  G4double thetaMax;
  G4double RMax;
  G4int CopyNo;

  G4double TotalMass;



  // Data for sensitive detector
  std::vector<G4String> names;
  std::vector<G4LogicalVolume*> volumes;
  std::vector<SPaRKLE_SensitiveDetector*> sensitives; 

  G4SDManager * SDman = NULL;




};











#endif
#endif