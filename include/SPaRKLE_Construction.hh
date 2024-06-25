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
#include "G4GDMLParser.hh"
#include "G4MaterialTable.hh"
#include <cmath>
#include "G4VisAttributes.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"


#include "SPaRKLE_Analysis.hh"
#include "SPaRKLE_Detector.hh"





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
  G4Tubs            *solidSiDet;
  G4LogicalVolume   *logicSiDet;
  G4VPhysicalVolume *physSiDet;
  G4VisAttributes   *visSiDet;

  // Al Frame 
  G4Tubs            *solidAlFrame;
  G4LogicalVolume   *logicAlFrame;
  G4VPhysicalVolume *physAlFrame;
  G4VisAttributes   *visAlFrame;


  // Scintillator Veto
  G4Box             *solidScintVeto;
  G4LogicalVolume   *logicScintVeto;
  G4VPhysicalVolume *physScintVeto; 
  G4VisAttributes   *visScintVeto;

  //Calo A1
  G4Box             *solidCalo;
  G4LogicalVolume   *logicCalo;
  G4VPhysicalVolume *physCalo; 
  G4VisAttributes   *visCalo;
  

  // DRILLED VETO
  // A: Front Side Drilled
  G4Box              *solidDrilledVeto;
  G4LogicalVolume    *logicDrilledVeto;
  G4VPhysicalVolume  *physDrilledVeto;
  G4SubtractionSolid *solidFinalDrilledVeto;
  G4VisAttributes    *visDrilledVeto;


  // VETO LATERAL
  G4Box              *solidLateralVeto;
  G4LogicalVolume    *logicLateralVeto;
  G4VPhysicalVolume  *physLateralVeto;
  G4VisAttributes    *visLateralVeto;

  // VETO BOTTOM
  G4Box              *solidBottomVeto;
  G4LogicalVolume    *logicBottomVeto;
  G4VPhysicalVolume  *physBottomVeto;
  G4VisAttributes    *visBottomVeto;


  //ALUMINIUM
  G4Box              *solidAlTop;
  G4SubtractionSolid *solidAlTopFinal;
  G4LogicalVolume    *logicAlTop;
  G4VPhysicalVolume  *physAlTop;
  G4VisAttributes    *visAl;

  G4Box              *solidAlLateral1;
  G4LogicalVolume    *logicAlLateral1;
  G4VPhysicalVolume  *physAlLateral1;

  G4Box              *solidAlLateral2;
  G4LogicalVolume    *logicAlLateral2;
  G4VPhysicalVolume  *physAlLateral2;




  // HOLE 
  G4Tubs             *solidHole;
  G4Tubs             *solidHoleLong;

  void DefineMaterials();

  // Materials
  G4Material *worldMat, *CsI, *EJ200, *SiMat, *Al, *bachelite;
  G4Element *Si;
  
  G4Material *GAGG, *BGO;


  // World
  G4double xWorld, yWorld, zWorld;
  // Container
  G4double xContainer, yContainer, zContainer;
  G4double LxContainer, LyContainer, LzContainer; 

  // Thin detector 
  // Active layer zone
  




  // Data for sensitive detector
  std::vector<G4String> names;
  std::vector<G4LogicalVolume*> volumes;
  std::vector<SPaRKLE_SensitiveDetector*> sensitives; 

  G4SDManager * SDman = NULL;




};











#endif