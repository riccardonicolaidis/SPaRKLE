#include "SPaRKLE_Construction.hh"
#include "SPaRKLE_Globals.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "G4SDManager.hh"
#include "G4SDStructure.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Region.hh" 
#include "G4UserLimits.hh"


using namespace std;



SPaRKLE_DetectorConstruction::SPaRKLE_DetectorConstruction()
{
  // Define the materials only once
  DefineMaterials();  
}

SPaRKLE_DetectorConstruction::~SPaRKLE_DetectorConstruction()
{
}

void SPaRKLE_DetectorConstruction::DefineMaterials()
{
  G4NistManager *nist = G4NistManager::Instance();

  GAGG = new G4Material("GAGG", 6.63*g/cm3, 4);
  GAGG -> AddElement(nist -> FindOrBuildElement("Gd"), 3);
  GAGG -> AddElement(nist -> FindOrBuildElement("Al"), 2);
  GAGG -> AddElement(nist -> FindOrBuildElement("Ga"), 3);
  GAGG -> AddElement(nist -> FindOrBuildElement("O"), 12);



  BGO = nist -> FindOrBuildMaterial("G4_BGO");



  Si    = nist -> FindOrBuildElement("Si");
  SiMat = new G4Material("Si", 2.328*g/cm3, 1);
  SiMat -> AddElement(Si,1);

  EJ200     = nist -> FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  worldMat  = nist -> FindOrBuildMaterial("G4_Galactic");
  Al        = nist -> FindOrBuildMaterial("G4_Al");
  bachelite = nist -> FindOrBuildMaterial("G4_BAKELITE");
  CsI       = nist -> FindOrBuildMaterial("G4_CESIUM_IODIDE");


  std::vector<G4double> RIndexWorld   = {1.                , 1.                 , 1.};
  std::vector<G4double> RIndexSi      = {3.88163           , 3.88163            , 3.88163};
  std::vector<G4double> RIndexAl      = {3.88163           , 3.88163            , 3.88163};
  std::vector<G4double> RIndexEJ200   = {1.58              , 1.58               , 1.58};
 

  std::vector<G4double> Wavelength0 = { 500 * nm, 425*nm, 400*nm};
  std::vector<G4double> Energy;

  for (int i = 0; i < Wavelength0.size(); i++)
  {
    Energy.push_back(1240*nm/Wavelength0[i]);
  }


  std::vector<G4double> ABSL        = { 380. * cm, 380. * cm , 380. * cm  };
  std::vector<G4double> ABSL_Al     = { 0.1 * um, 0.1 * um , 0.1 * um  };
  std::vector<G4double> ABSL_Si     = { 0.1 * um, 0.1 * um , 0.1 * um  };






  G4MaterialPropertiesTable *mptWorld   = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable *mptSi      = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable *mptEJ200   = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable *mptAl      = new G4MaterialPropertiesTable();

  

  mptWorld -> AddProperty("RINDEX", Energy, RIndexWorld, 3);
  worldMat -> SetMaterialPropertiesTable(mptWorld);

  mptSi    -> AddProperty("RINDEX", Energy, RIndexSi, 3);
  mptSi    -> AddProperty("ABSLENGTH", Energy, ABSL_Si, 3);
  SiMat    -> SetMaterialPropertiesTable(mptSi);

  mptEJ200 -> AddProperty("RINDEX", Energy, RIndexEJ200, 3);
  mptEJ200 -> AddProperty("ABSLENGTH", Energy, ABSL, 3);
  EJ200    -> SetMaterialPropertiesTable(mptEJ200);
  
  mptAl       -> AddProperty("RINDEX", Energy, RIndexAl, 3);
  mptAl       -> AddProperty("ABSLENGTH", Energy, ABSL_Al, 3);
  Al          -> SetMaterialPropertiesTable(mptAl);


}

G4VPhysicalVolume *SPaRKLE_DetectorConstruction::Construct()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Define a file to save the angles of the holes
  ofstream myfile;
  myfile.open ("../output/GeometryReport.txt");
  myfile << "Report of the geometry of the detector" << endl;
  myfile << "======================================" << endl;
  // Print the date and time at the execution of the file
  time_t now = time(0);
  char* dt = ctime(&now);

  myfile << "\nLocal Date and Time is: " << dt << endl;
  myfile << "Compilation Date and Time: " << __DATE__ << " " << __TIME__ << endl << endl;


  

  TotalMass = 0;

  // .oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo
  // .oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo
  //
  //        NUMERICAL DATA - PARAMETRIZATION
  //
  // .oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo
  // .oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo
  
  // World dimensions
  xWorld = 3*m;
  yWorld = 3*m;
  zWorld = 3*m;

  myfile << "Dimensions are in cm unless otherwise specified" << endl << endl;
  myfile << "World dimensions: " << 2*xWorld/cm << " x " << 2*yWorld/cm << " x " << 2*zWorld/cm << " cm" << endl;


  // Container position
  xContainer = 0. *cm;
  yContainer = 0. *cm;
  zContainer = 0. *cm; 

  myfile << "Container position: " << xContainer/cm << " x " << yContainer/cm << " x " << zContainer/cm << " cm" << endl;

  // Container dimensions 
  LxContainer = 2 *m;
  LyContainer = 2 *m;
  LzContainer = 2 *m; 

  myfile << "Container dimensions: " << LxContainer/cm << " x " << LyContainer/cm << " x " << LzContainer/cm << " cm" << endl;

  // Dimensions of the detectors
  // ACTIVE AREA OF THE DETECTOR

  LSquareCentersThin = 5.8 * cm;
  LActive = 8 * cm;

  // THICKNESS OF THE MATERIAL
  TkThin        = TK_THIN;

  myfile << "Thin detector thickness: " << TkThin/um << " um" << endl;
  
  // DATA OF THE GEOMETRY OF THE DETECTORS
  // SILICON DETECTORS

  rThinDet = 5. * mm;

  // HOLES IN THE PCB
  rHolesPCB = rThinDet + 0.3 * mm ;

  // THICKNESS OF THE PCB BOARDS
  TkPCBThin     = 2 * mm;

  // Extra thickness of the detector
  ExtraSpacing         = 100 * um;
  SpacingPCBComponents = 4 * mm;

  // Dimensions of the PCB
  LxPCB = LActive - 1.2*cm;
  LyPCB = LActive - 1.2*cm;

  // Dimensions of the Veto
  LScintVeto  = LActive - 0*cm;
  TkScintVeto = TK_GAGG;

  // Dimensions of the GAGG calorimeters
  LCalo = 4*cm;
  TkCalo = TK_GAGG;

  myfile << "Plastic scintillator thickness: " << TkScintVeto/cm << " cm" << endl;


  // ACTIVE LAYER POSITION
  xThin = 0.;
  yThin = 0.;
  zThin = 0.;

  myfile << "Thin detector plane position: " << xThin/cm << " x " << yThin/cm << " x " << zThin/cm << " cm" << endl;


  xScintVeto = 0.;
  yScintVeto = 0.;
  zScintVeto = zThin + TkScintVeto/2. + TkPCBThin/2. + ExtraSpacing + SpacingPCBComponents/2.;

 //myfile << "First plastic scintillator position : " << xScintVeto/cm << " x " << yScintVeto/cm << " x " << zScintVeto/cm << " cm" << endl;

 
  xCalo_A1 = -2.05*cm;
  yCalo_A1= 2.05*cm;

  xCalo_A2 = 2.05*cm;
  yCalo_A2= -2.05*cm;

  xCalo_B1 = 2.05*cm;
  yCalo_B1= 2.05*cm;

  xCalo_B2 = -2.05*cm;
  yCalo_B2= -2.05*cm;

  zCalo = zThin + TkScintVeto/2. + TkPCBThin/2. + ExtraSpacing + SpacingPCBComponents/2.;

  myfile << "Calo A1 position : " << xCalo_A1/cm << " x " << yCalo_A1/cm << " x " << zCalo/cm << " cm3" << endl;

  myfile << "Calo A2 position : " << xCalo_A2/cm << " x " << yCalo_A2/cm << " x " << zCalo/cm << " cm3" << endl;

  myfile << "Calo B1 position : " << xCalo_B1/cm << " x " << yCalo_B1/cm << " x " << zCalo/cm << " cm3" << endl;

  myfile << "Calo B2 position : " << xCalo_B2/cm << " x " << yCalo_B2/cm << " x " << zCalo/cm << " cm3" << endl;


  // Construction of the Drilled Veto
  // BIGGER than the Veto on the back of the instrument
  LxDrilledVeto = LActive- 0.3*cm;
  LyDrilledVeto = LActive - 0.3 * cm;
  TkDrilledVeto = 0.8*cm;

  myfile << "Drilled Veto thickness: " << TkDrilledVeto/cm << " cm" << endl;

  
  xDrilledVeto = 0.;
  yDrilledVeto = 0.;
  zDrilledVeto = zThin  - TkDrilledVeto/2. - TkPCBThin/2. - ExtraSpacing - SpacingPCBComponents;

  myfile << "Drilled Veto position: " << xDrilledVeto/cm << " x " << yDrilledVeto/cm << " x " << zDrilledVeto/cm << " cm" << endl;

  // Construction of the Drilled Aluminium
  // BIGGER than the Veto Drilled
  LxDrilledAl0 = 10 * cm;
  LyDrilledAl0 = 10 * cm;
  TkDrilledAl0 = 0.8 * cm;
  LxDrilledAl1 = (10 - 0.3*2 )* cm;
  LyDrilledAl1 = (10 - 0.3*2 )* cm;
  TkCompenetration = 0.1 * cm;

  xDrilledAl = 0.;
  yDrilledAl = 0.;
  zDrilledAl = zDrilledVeto -TkDrilledVeto/2. - TkDrilledAl0/2. - ExtraSpacing;

  // LENGTH OF THE SIDE WALLS OF THE ALUMINIUM AND THE DRILLED VETO
  LzDrilledVeto = (zScintVeto +  (N_PL_SCINT_NO_VETO + 1)*TkScintVeto + (ExtraSpacing)*N_PL_SCINT_NO_VETO) - zDrilledVeto - TkDrilledVeto/2.;
  LzDrilledAl1 =  (zScintVeto + (N_PL_SCINT_NO_VETO + 1)*TkScintVeto + (ExtraSpacing)*N_PL_SCINT_NO_VETO) - zDrilledAl - TkDrilledAl0/2. + TkCompenetration;

  // HOLES DATA
  NxHoles = NX_SENSORS;
  NyHoles = NY_SENSORS;
  rHoles = rThinDet;
  LengthHole = 2*cm;

 
  solidWorld = new G4Box("solidWorld",xWorld*0.5, yWorld*0.5, zWorld*0.5);  
  logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
  visWorld   = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visWorld   -> SetVisibility(false);
  logicWorld -> SetVisAttributes(visWorld);
  physWorld  = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicWorld, "physWorld", 0, false, 0, true);
    
  solidContainer = new G4Box("solidContainer",LxContainer*0.5, LyContainer*0.5, LzContainer*0.5);  
  logicContainer = new G4LogicalVolume(solidContainer, worldMat, "logicContainer");
  visContainer   = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visContainer   -> SetVisibility(false);
  logicContainer -> SetVisAttributes(visContainer);


  solidScintVeto = new G4Box("solidScintVeto", LScintVeto*0.5, LScintVeto*0.5, TkCalo*0.5);
  logicScintVeto = new G4LogicalVolume(solidScintVeto, GAGG, "logicScintVeto");
  visScintVeto   = new G4VisAttributes(G4Colour(1.,1.,0.));
  logicScintVeto -> SetVisAttributes(visScintVeto);


  /*
  myfile << "Plastic scintillators \n\n" << endl;

  G4double zBottomVeto =0.;

  for(G4int i = 0; i < (N_PL_SCINT_NO_VETO); ++ i)
  {
    if(i == 0)// Add skin surface
    {
      physScintVeto  = new G4PVPlacement(0,G4ThreeVector(xScintVeto, yScintVeto, zScintVeto),logicScintVeto, "physScintVeto", logicContainer, false, i, true );
      zBottomVeto = zScintVeto;
    }
    else
    {
      physScintVeto  = new G4PVPlacement(0,G4ThreeVector(xScintVeto, yScintVeto, (zScintVeto + i*(0.1* mm + TkScintVeto + ExtraSpacing))),logicScintVeto, "physScintVeto", logicContainer, false, i, true );
      zBottomVeto = zScintVeto + i*(0.1* mm + TkScintVeto + ExtraSpacing);
    }
    myfile << "Scintillator " << i << " placed at " << zScintVeto + i*(0.1* mm + TkScintVeto + ExtraSpacing) << " mm" << endl;
    TotalMass += (logicScintVeto->GetMass())/kg;
  }

  zBottomVeto += 0.1*mm + TkScintVeto + ExtraSpacing + 1* mm;
*/


  solidCalo_A1= new G4Box("solidCalo_A1", LCalo*0.5, LCalo*0.5, TkCalo*0.5);
  logicCalo_A1 = new G4LogicalVolume(solidCalo_A1, GAGG, "logicCalo_A1");
  visCalo_A1  = new G4VisAttributes(G4Colour(1.,1.,0.));
  logicCalo_A1 -> SetVisAttributes(visCalo_A1);

  solidCalo_A2= new G4Box("solidCalo_A2", LCalo*0.5, LCalo*0.5, TkCalo*0.5);
  logicCalo_A2 = new G4LogicalVolume(solidCalo_A2, GAGG, "logicCalo_A2");
  visCalo_A2  = new G4VisAttributes(G4Colour(1.,1.,0.));
  logicCalo_A2 -> SetVisAttributes(visCalo_A2);

  solidCalo_B1= new G4Box("solidCalo_B1", LCalo*0.5, LCalo*0.5, TkCalo*0.5);
  logicCalo_B1 = new G4LogicalVolume(solidCalo_B1, GAGG, "logicCalo_B1");
  visCalo_B1  = new G4VisAttributes(G4Colour(1.,1.,0.));
  logicCalo_B1 -> SetVisAttributes(visCalo_B1);

  solidCalo_B2= new G4Box("solidCalo_B2", LCalo*0.5, LCalo*0.5, TkCalo*0.5);
  logicCalo_B2 = new G4LogicalVolume(solidCalo_B2, GAGG, "logicCalo_B2");
  visCalo_B2  = new G4VisAttributes(G4Colour(1.,1.,0.));
  logicCalo_B2 -> SetVisAttributes(visCalo_B2);
  physCalo_B2 = new G4PVPlacement(0, G4ThreeVector(xCalo_B2, yCalo_B2, zCalo), logicCalo_B2, "physCalo_B2", logicContainer, false, 0, true);


  G4double zBottomVeto =0.;

  for(G4int i = 0; i < (N_PL_SCINT_NO_VETO); ++ i)
  {
    if(i == 0)// Add skin surface
    {
      physCalo_A1 = new G4PVPlacement(0, G4ThreeVector(xCalo_A1, yCalo_A1, zCalo), logicCalo_A1, "physCalo_A1", logicContainer, false, i, true);
      physCalo_A2 = new G4PVPlacement(0, G4ThreeVector(xCalo_A2, yCalo_A2, zCalo), logicCalo_A2, "physCalo_A2", logicContainer, false, i, true);
      physCalo_B1 = new G4PVPlacement(0, G4ThreeVector(xCalo_B1, yCalo_B1, zCalo), logicCalo_B1, "physCalo_B1", logicContainer, false, i, true);
      physCalo_B2 = new G4PVPlacement(0, G4ThreeVector(xCalo_B2, yCalo_B2, zCalo), logicCalo_B2, "physCalo_B2", logicContainer, false, i, true);

      zBottomVeto = zCalo;
    }
    else
    {
      physCalo_A1 = new G4PVPlacement(0, G4ThreeVector(xCalo_A1, yCalo_A1, zCalo + i*(0.1* mm + TkCalo + ExtraSpacing)), logicCalo_A1, "physCalo_A1", logicContainer, false, i, true);
      physCalo_A2 = new G4PVPlacement(0, G4ThreeVector(xCalo_A2, yCalo_A2, zCalo + i*(0.1* mm + TkCalo + ExtraSpacing)), logicCalo_A2, "physCalo_A2", logicContainer, false, i, true);
      physCalo_B1 = new G4PVPlacement(0, G4ThreeVector(xCalo_B1, yCalo_B1, zCalo + i*(0.1* mm + TkCalo + ExtraSpacing)), logicCalo_B1, "physCalo_B1", logicContainer, false, i, true);
      physCalo_B2 = new G4PVPlacement(0, G4ThreeVector(xCalo_B2, yCalo_B2, zCalo + i*(0.1* mm + TkCalo + ExtraSpacing)), logicCalo_B2, "physCalo_B2", logicContainer, false, i, true);

      zBottomVeto = zCalo + i*(0.1* mm + TkCalo + ExtraSpacing);
    }
    myfile << "Layer " << i << " of calorimeters placed at " << zCalo + i*(0.1* mm + TkCalo + ExtraSpacing) << " mm" << endl;
    
    TotalMass += (logicCalo_A1->GetMass())/kg;
    TotalMass += (logicCalo_A2->GetMass())/kg;
    TotalMass += (logicCalo_B1->GetMass())/kg;
    TotalMass += (logicCalo_B2->GetMass())/kg;
  }

  zBottomVeto += zCalo + 0.1*mm + TkCalo + ExtraSpacing + 1* mm;


  VaschettaVeto = SPaRKLE_DetectorConstruction::SolidoVaschetta((zBottomVeto-zDrilledVeto - TkScintVeto/2.), 1*cm, LScintVeto+3*mm, LScintVeto+3*mm, LScintVeto+11*mm, LScintVeto+11*mm);
  logicVaschettaVeto = new G4LogicalVolume(VaschettaVeto, EJ200, "logicVaschettaVeto");
  visVaschettaVeto = new G4VisAttributes(G4Colour(0.3,0.3,1.0));
  logicVaschettaVeto -> SetVisAttributes(visVaschettaVeto);
  physVaschettaVeto = new G4PVPlacement(0,G4ThreeVector(xScintVeto, yScintVeto, zBottomVeto ),logicVaschettaVeto, "physVaschettaVeto", logicContainer, false, 0, true );

  // ...oooOOOooo...oooOOOooo...oooOOOooo...oooOOOooo...oooOOOooo...oooOOOooo...
  // ...oooOOOooo...oooOOOooo...oooOOOooo...oooOOOooo...oooOOOooo...oooOOOooo...

  // SILICON DETECTORS DELTA E - E PADS 
  // LOGIC AND SOLID DEFINITIONS

  solidSiDetThin  = new G4Tubs("solidSiDetThin", 0.,rThinDet, TkThin*0.5, 0., 360. *deg);
  

  logicSiDetThin = new G4LogicalVolume(solidSiDetThin, SiMat, "logicSiDetThin");
  TotalMass += (logicSiDetThin->GetMass())*16/kg;
  
  
  // Drilled Veto Solid without holes
  solidDrilledVeto    = new G4Box("solidDrilledVeto",LxDrilledVeto*0.5, LyDrilledVeto*0.5, TkDrilledVeto*0.5);  
  
  // Drilled Aluminium without holes
  solidDrilledAl0    = new G4Box("solidDrilledAl0",LxDrilledAl0*0.5, LyDrilledAl0*0.5, TkDrilledAl0*0.5);  
  solidDrilledAl1Ext = new G4Box("solidDrilledAl1Ext", LxDrilledAl0*0.5, LyDrilledAl0*0.5, (LzDrilledAl1+TkCompenetration)*0.5);
  solidDrilledAl1Int = new G4Box("solidDrilledAl1Int", LxDrilledAl1*0.5, LyDrilledAl1*0.5, (LzDrilledAl1+TkCompenetration + 10*um)*0.5);
  solidDrilledAl1 = new G4SubtractionSolid("solidDrilledAl1", solidDrilledAl1Ext, solidDrilledAl1Int);
  

  // Objects for the transformations 
  G4RotationMatrix Identity;
  G4ThreeVector position1;
  G4Transform3D tr1;

  // Drilled Veto  
  Identity = G4RotationMatrix();
  position1 = G4ThreeVector(0.,0.,(TkDrilledVeto+LzDrilledVeto)/2. -TkCompenetration);
  tr1 = G4Transform3D(Identity, position1);
  
  // Drilled Aluminium  
  Identity = G4RotationMatrix();
  position1 = G4ThreeVector(0.,0.,(TkDrilledAl0+LzDrilledAl1)/2. -TkCompenetration);
  tr1 = G4Transform3D(Identity, position1);
  solidDrilledAl = new G4UnionSolid("solidDrilledAl", solidDrilledAl0, solidDrilledAl1, tr1);
  
  xDelta = LSquareCentersThin / (NxHoles);
  yDelta = LSquareCentersThin / (NyHoles);

  G4RotationMatrix HoleRotation;
  G4ThreeVector HolePositionVeto;
  G4Transform3D HoleTransformVeto;
  G4ThreeVector HolePositionAl;
  G4Transform3D HoleTransformAl;
  G4ThreeVector RotationAxis;
  G4ThreeVector HolePositionPCB;
  G4RotationMatrix HoleRotationPCB;
  G4Transform3D HoleTransformPCB;


  thetaMax = 40*deg;
  RMax = std::sqrt(2) * 1.5;

  /* -------------------------------------------------------------------------- */
  /*                           DEFINITION OF THE FRAME                          */
  /* -------------------------------------------------------------------------- */



  // DEFINITION OF THE HOLES 
  solidHole     = new G4Tubs("solidHole", 0., rHoles, LengthHole+1*cm, 0.*deg, 360.*deg);
  solidHoleLong = new G4Tubs("solidHoleLong", 0., rHoles, LengthHole+3*cm, 0.*deg, 360.*deg);
  solidHolePCB = new G4Tubs("solidHolePCB", 0., rHolesPCB, 1*cm, 0.*deg, 360.*deg);

  G4double LHoleVisualization = 20 * cm;
  G4Tubs *solidHoleVisualization = new G4Tubs("solidHoleVisualization", 0., rHoles, LHoleVisualization/2., 0.*deg, 360.*deg);
  G4LogicalVolume *logicHoleVisualization = new G4LogicalVolume(solidHoleVisualization, worldMat, "logicHoleVisualization");
  G4VisAttributes *visHoleVisualization = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  logicHoleVisualization -> SetVisAttributes(visHoleVisualization);




  // DEFINITION OF THE PCB BOARDS
  solidPCBThin = new G4Box("solidPCBThin", LxPCB*0.5, LyPCB*0.5, TkPCBThin*0.5);
  
  myfile << "Angles Direction" << endl;
  myfile << "================================================" << endl;
  myfile << "ix \t iy \t phi \t theta \t x \t y \t x' \t y'" << endl;

  G4int CopyNoMatrix[NX_SENSORS][NY_SENSORS];
  G4int jCopyNo = 0;
  for(G4int iy = 0; iy < NyHoles; ++iy)
  {
    for(G4int ix = 0; ix <NxHoles; ++ix)
    {
      CopyNoMatrix[ix][iy] = (jCopyNo++);

      //xCurrentHole = -(LSquareCentersThin/2.) + (xDelta/2.) + ix * xDelta;
      //yCurrentHole = -(LSquareCentersThin/2.) + (yDelta/2.) + iy * yDelta;

      xCurrentHole = -LCalo*0.5 + iy*LCalo;
      yCurrentHole = LCalo*0.5 - iy*LCalo;

      G4double Delta_ii_x = NxHoles/2. -0.5;
      G4double Delta_ii_y = NyHoles/2. -0.5;


      HoleRotation = G4RotationMatrix();
      phi = std::atan2((iy-Delta_ii_y),(ix-Delta_ii_x));
      distanceR = std::sqrt((iy-Delta_ii_y)*(iy-Delta_ii_y) + (ix-Delta_ii_x)*(ix-Delta_ii_x));
      theta = distanceR * thetaMax / RMax;

      myfile << ix << "\t" << iy << "\t" << phi/deg << "\t" << theta/deg << "\t" << xCurrentHole/cm << "\t" << yCurrentHole/cm  << endl;

      if((ix == 0) && (iy == 0) )
      {
        HolePositionPCB = G4ThreeVector(xCurrentHole, yCurrentHole, zThin);
        HoleRotationPCB = G4RotationMatrix();
        HoleTransformPCB = G4Transform3D(HoleRotationPCB, HolePositionPCB);
        solidPCBThinFinal = new G4SubtractionSolid("solidPCBThinFinal", solidPCBThin, solidHolePCB, HoleTransformPCB);
      }
      else
      {
        HolePositionPCB = G4ThreeVector(xCurrentHole, yCurrentHole, zThin);
        HoleRotationPCB = G4RotationMatrix();
        HoleTransformPCB = G4Transform3D(HoleRotationPCB, HolePositionPCB);
        solidPCBThinFinal = new G4SubtractionSolid("solidPCBThinFinal", solidPCBThinFinal, solidHolePCB, HoleTransformPCB);
      }



      RotationAxis = G4ThreeVector(-std::sin(phi)*std::sin(theta), std::cos(phi)*std::sin(theta), 0.);
      HoleRotation.rotate(theta, RotationAxis);
      HoleRotation.invert(); 
      
      HolePositionVeto = G4ThreeVector(xCurrentHole, yCurrentHole, std::abs(zDrilledVeto));
      HoleTransformVeto = G4Transform3D(HoleRotation, HolePositionVeto);

      HolePositionAl = G4ThreeVector(xCurrentHole, yCurrentHole, std::abs(zDrilledAl));
      HoleTransformAl = G4Transform3D(HoleRotation, HolePositionAl);



      if((ix == 0) && (iy == 0) )
      {
        solidFinalDrilledVeto = new G4SubtractionSolid("solidFinalDrilledVeto", solidDrilledVeto, solidHole, HoleTransformVeto);
        solidFinalDrilledAl   = new G4SubtractionSolid("solidFinalDrilledAl", solidDrilledAl, solidHoleLong, HoleTransformAl);
      } else 
      {
        solidFinalDrilledVeto = new G4SubtractionSolid("solidFinalDrilledVeto", solidFinalDrilledVeto, solidHole, HoleTransformVeto);
        solidFinalDrilledAl   = new G4SubtractionSolid("solidFinalDrilledAl", solidFinalDrilledAl, solidHoleLong, HoleTransformAl);
      }

      G4double translationL = LHoleVisualization*0.5 + 1 * cm;
      HolePositionVeto = G4ThreeVector(xCurrentHole + translationL*std::sin(theta)*std::cos(phi), yCurrentHole + translationL*std::sin(theta)*std::sin(phi), zThin - translationL*std::cos(theta));
      HoleTransformVeto = G4Transform3D(HoleRotation, HolePositionVeto);


      //new G4PVPlacement(HoleTransformVeto, logicHoleVisualization, "physHoleVisualization", logicContainer, false, CopyNoMatrix[ix][iy], true);

      physSiDetThin   = new G4PVPlacement(0, G4ThreeVector(xCurrentHole,yCurrentHole,zThin) , logicSiDetThin, "physSiDetThin", logicContainer, false, CopyNoMatrix[ix][iy], true);
      
      //physSiDetThick  = new G4PVPlacement(0, G4ThreeVector(xCurrentHoleThick,yCurrentHoleThick,zThick) , logicSiDetThick, "physSiDetThick", logicContainer, false, CopyNoMatrix[ix][iy], true);
    }
  }

/*
  

  physSiDetThin   = new G4PVPlacement(0, G4ThreeVector(-LCalo*0.5, LCalo*0.5,zThin) , logicSiDetThin, "physSiDetThin", logicContainer, false, 0, true);

  physSiDetThin   = new G4PVPlacement(0, G4ThreeVector(LCalo*0.5, -LCalo*0.5,zThin) , logicSiDetThin, "physSiDetThin", logicContainer, false, 1, true);
*/  
  
  // LOGIC DRILLED VETO
  logicDrilledVeto = new G4LogicalVolume(solidFinalDrilledVeto, EJ200, "logicDrilledVeto");
  TotalMass += (logicDrilledVeto->GetMass())/kg;
  visDrilledVeto = new G4VisAttributes(G4Colour(0.3,0.3,1.0));
  logicDrilledVeto->SetVisAttributes(visDrilledVeto);
  
  // LOGIC DRILLED ALUMINIUM
  logicDrilledAl = new G4LogicalVolume(solidFinalDrilledAl, Al, "logicDrilledAl");
  TotalMass += (logicDrilledAl->GetMass())/kg;
  visDrilledAl = new G4VisAttributes(G4Colour(0.6,0.6,1.));
  logicDrilledAl->SetVisAttributes(visDrilledAl);

  // LOGIC PCB
  logicPCBThin = new G4LogicalVolume(solidPCBThinFinal, bachelite , "logicPCBThin");
  TotalMass += (logicPCBThin->GetMass())/kg;
  visPCBThin = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
  logicPCBThin->SetVisAttributes(visPCBThin);


  // PHYSICAL PLACEMENT DRILLED ALUMINIUM
  physDrilledAl  = new G4PVPlacement(0, G4ThreeVector(xDrilledAl,yDrilledAl,zDrilledAl) , logicDrilledAl, "physDrilledAl", logicContainer, false, 0, true);

  // PHYSICAL PLACEMENT DRILLED VETO
  physDrilledVeto  = new G4PVPlacement(0, G4ThreeVector(xDrilledVeto,yDrilledVeto,zDrilledVeto) , logicDrilledVeto, "physDrilledVeto", logicContainer, false, 0, true);
  
  
  // PHYSICAL PLACEMENT PCB
  physPCBThin  = new G4PVPlacement(0, G4ThreeVector(xThin,yThin,zThin) , logicPCBThin, "physPCBThin", logicContainer, false, 0, true);
  

  // Quantities for the lateral veto

  G4double z1 = zThin - SpacingPCBComponents/2. - TkPCBThin/2.;
  G4double z2 = zScintVeto + N_PL_SCINT_NO_VETO *(TkScintVeto + 100*um + ExtraSpacing) + TkScintVeto/2. + 4 * mm;


  G4double xLateralVeto = 0;
  G4double yLateralVeto = LActive/2. + TkLateralVeto/2. + 2 * mm;
  G4double zLateralVeto = (z1+z2)/2.;

  LzLateralVeto = abs(z2 - z1);
  LxLateralVeto = LScintVeto + 2 * mm;
  TkLateralVeto =  TkScintVeto;



  solidLateralVeto = new G4Box("solidLateralVeto", LxLateralVeto*0.5, TkLateralVeto*0.5, LzLateralVeto*0.5);
  logicLateralVeto = new G4LogicalVolume(solidLateralVeto, EJ200, "logicLateralVeto");
  TotalMass += (logicLateralVeto->GetMass())/kg;
  visLateralVeto = new G4VisAttributes(G4Colour(0.3,0.3,1.0));
  logicLateralVeto->SetVisAttributes(visLateralVeto);


  G4ThreeVector positionLateralVeto = G4ThreeVector(xLateralVeto,yLateralVeto,zLateralVeto);
  //new G4PVPlacement(0, positionLateralVeto, logicLateralVeto, "physLateralVeto", logicContainer, false, 0, true);
  positionLateralVeto = G4ThreeVector(xLateralVeto,-yLateralVeto,zLateralVeto);
  //new G4PVPlacement(0, positionLateralVeto, logicLateralVeto, "physLateralVeto", logicContainer, false, 1, true);



  G4RotationMatrix* rotLateralVeto = new G4RotationMatrix();
  rotLateralVeto->rotateZ(90*deg);
  positionLateralVeto = G4ThreeVector(yLateralVeto,xLateralVeto,zLateralVeto);
  //new G4PVPlacement(rotLateralVeto, positionLateralVeto, logicLateralVeto, "physLateralVeto", logicContainer, false, 2, true);
  positionLateralVeto = G4ThreeVector(-yLateralVeto,xLateralVeto,zLateralVeto);
  //new G4PVPlacement(rotLateralVeto, positionLateralVeto, logicLateralVeto, "physLateralVeto", logicContainer, false, 3, true);


  
  physContainer  = new G4PVPlacement(0, G4ThreeVector(xContainer,yContainer,zContainer), logicContainer, "physContainer", logicWorld, false, 0, true);
  


  // CONTROL
  
  for(G4int M = 0; M <= 10; ++M)
  {
    G4cout << "#####################" << G4endl;
  }

  G4cout << "LThin =  " << TkThin << G4endl;
  G4cout << "Total mass = " << TotalMass << " kg" << G4endl;

  myfile << "LThin =  " << TkThin/cm << G4endl;
  myfile << "Total mass = " << TotalMass << " kg" << G4endl;


  for(G4int M = 0; M <= 10; ++M)
  {
    G4cout << "#####################" << G4endl;
  }


  
  if(GENERATE_GDML)
  {
    G4GDMLParser parser;

    // Create a string with date and time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime (buffer,80,"%d-%m-%Y_%H-%M-%S",timeinfo);
    std::string str(buffer);

    // Create a string with the name of the file + date and time
    std::string str1 = "../CadGeometry/SPaRKLE_geometry_";
    std::string str2 = ".gdml";
    std::string str3 = str1 + str + str2;

    parser.Write(str3, physContainer);


    str1 = "../CadGeometry/DrilledVeto_";
    str3 = str1 + str + str2;

    parser.Write(str3, physDrilledVeto);

    myfile << "GDML file generated" << G4endl;
    myfile << "GDML file name: " << str3 << G4endl;
    myfile.close();
  }
  

  
  return physWorld;  
}


void SPaRKLE_DetectorConstruction::ConstructSDandField()
{
  if(SDman == NULL)
  {

  
    SDman = G4SDManager::GetSDMpointer();


    names.push_back("SiThin");
    volumes.push_back(logicSiDetThin);

    names.push_back("Calo");
    volumes.push_back(logicScintVeto);

    names.push_back("DrilledVeto");
    volumes.push_back(logicDrilledVeto);

    names.push_back("BottomVeto");
    volumes.push_back(logicVaschettaVeto);

    for(G4int i = 0; i < (G4int)names.size(); ++i)
    {
      SPaRKLE_SensitiveDetector* aSD = new SPaRKLE_SensitiveDetector(names[i]);
      sensitives.push_back(aSD);

      SDman->AddNewDetector(sensitives[i]);
      SetSensitiveDetector( volumes[i], sensitives[i]);

    }
  }


}


void SPaRKLE_DetectorConstruction::ClearSD()
{

  return;
}


G4UnionSolid * SPaRKLE_DetectorConstruction::SolidoVaschetta(G4double hInternal, G4double TkBottom, G4double xIn, G4double yIn, G4double xOut, G4double yOut)
{
  G4UnionSolid *solidoVaschetta;

  G4double ExtraSpacing2 = 0.1 * mm; // Only for union to be well defined

  G4Box *extEnvelope = new G4Box("extEnvelope", xOut/2., yOut/2., (hInternal + ExtraSpacing2)/2.);
  G4Box *intEnvelope = new G4Box("intEnvelope", xIn/2., yIn/2., (hInternal + 2*ExtraSpacing2)/2.); // The 2* is for the subtraction to be well defined
  G4SubtractionSolid *LateralEnvelope = new G4SubtractionSolid("LateralEnvelope", extEnvelope, intEnvelope);


  G4Box *bottom = new G4Box("bottom", xOut/2., yOut/2., TkBottom/2.);

  // Ora voglio unire le pareti laterali con il fondo
  // Ma voglio che le pareti laterali siano spostate in modo tale da formare una vaschetta
  // La traslazione deve essere lungo z

  G4double zTranslation = (hInternal - ExtraSpacing2)/2. + TkBottom/2.;
  G4ThreeVector translation(0.,0.,-zTranslation);

  G4Transform3D tr = G4Transform3D(G4RotationMatrix(), translation);

  solidoVaschetta = new G4UnionSolid("solidoVaschetta", bottom, LateralEnvelope,  tr);


  return (G4UnionSolid*) solidoVaschetta;
}


