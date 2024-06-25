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


  




  // WORLD DEFINITION
 
  solidWorld = new G4Box("solidWorld",xWorld*0.5, yWorld*0.5, zWorld*0.5);  
  logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
  visWorld   = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visWorld   -> SetVisibility(false);
  logicWorld -> SetVisAttributes(visWorld);
  physWorld  = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicWorld, "physWorld", 0, false, 0, true);
    
  // CONTAINER DEFINITION

  solidContainer = new G4Box("solidContainer",LxContainer*0.5, LyContainer*0.5, LzContainer*0.5);  
  logicContainer = new G4LogicalVolume(solidContainer, worldMat, "logicContainer");
  visContainer   = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visContainer   -> SetVisibility(false);
  logicContainer -> SetVisAttributes(visContainer);


  // GAGG DIMENSIONS
  // Dimensions of the GAGG calorimeters
  G4double LCalo = 4*cm;
  G4double TkCalo = 1*cm;
  G4double Displacing_Calo = LCalo/2. + 0.05*cm;
  G4double xCalo_A1 = -Displacing_Calo;
  G4double yCalo_A1= Displacing_Calo;
  G4double xCalo_A2 = Displacing_Calo;
  G4double yCalo_A2= -Displacing_Calo;
  G4double xCalo_B1 = Displacing_Calo;
  G4double yCalo_B1= Displacing_Calo;
  G4double xCalo_B2 = -Displacing_Calo;
  G4double yCalo_B2= -Displacing_Calo;
  G4double zCalo = 0.;


  myfile << "Calo A1 position : " << xCalo_A1/cm << " x " << yCalo_A1/cm << " x " << zCalo/cm << " cm3" << endl;
  myfile << "Calo A2 position : " << xCalo_A2/cm << " x " << yCalo_A2/cm << " x " << zCalo/cm << " cm3" << endl;
  myfile << "Calo B1 position : " << xCalo_B1/cm << " x " << yCalo_B1/cm << " x " << zCalo/cm << " cm3" << endl;
  myfile << "Calo B2 position : " << xCalo_B2/cm << " x " << yCalo_B2/cm << " x " << zCalo/cm << " cm3" << endl;

  solidCalo= new G4Box("solidCalo", LCalo*0.5, LCalo*0.5, TkCalo*0.5);
  logicCalo = new G4LogicalVolume(solidCalo, GAGG, "logicCalo");
  visCalo  = new G4VisAttributes(G4Colour(1.,1.,0.));
  logicCalo -> SetVisAttributes(visCalo);




  new G4PVPlacement(0, G4ThreeVector(xCalo_A1, yCalo_A1, zCalo), logicCalo, "physCalo_A1", logicContainer, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(xCalo_A2, yCalo_A2, zCalo), logicCalo, "physCalo_A2", logicContainer, false, 1, true);
  new G4PVPlacement(0, G4ThreeVector(xCalo_B1, yCalo_B1, zCalo), logicCalo, "physCalo_B1", logicContainer, false, 2, true);
  new G4PVPlacement(0, G4ThreeVector(xCalo_B2, yCalo_B2, zCalo), logicCalo, "physCalo_B2", logicContainer, false, 3, true);




  // DATA OF THE GEOMETRY OF THE DETECTORS
  // SILICON DETECTORS
  G4double Area_Silicon = 150 * mm * mm;
  G4double rDet = sqrt(Area_Silicon/3.14);
  G4double SiDetTk = 100 * um;

  G4double AlFrameRadiusExt = 1*cm;
  G4double AlFrameRadiusInt = rDet+0.2*mm;
  G4double AlFrameTk = 5*mm;

  G4double Displacement_Silicon = LCalo/2.;
  G4double xSi_1 = -Displacement_Silicon;
  G4double ySi_1 = Displacement_Silicon;
  G4double xSi_2 = Displacement_Silicon;
  G4double ySi_2 = -Displacement_Silicon;

  G4double zSi = -zCalo - TkCalo/2. - AlFrameTk/2. - 1*mm;


  solidSiDet = new G4Tubs("SiDet", 0., rDet, SiDetTk, 0., 360.*deg);
  logicSiDet = new G4LogicalVolume(solidSiDet, SiMat, "logicSiDet");

  solidAlFrame = new G4Tubs("AlFrame", AlFrameRadiusInt, AlFrameRadiusExt, AlFrameTk/2., 0., 360.*deg);
  logicAlFrame = new G4LogicalVolume(solidAlFrame, Al, "logicAlFrame");
  visAlFrame  = new G4VisAttributes(G4Colour(1.,0.,0.));
  logicAlFrame -> SetVisAttributes(visAlFrame);

  new G4PVPlacement(0, G4ThreeVector(xSi_1, ySi_1,zSi), logicSiDet, "physSD1", logicContainer, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(xSi_2, ySi_2,zSi), logicSiDet, "physSD2", logicContainer, false, 1, true);

  new G4PVPlacement(0, G4ThreeVector(xSi_1, ySi_1,zSi), logicAlFrame, "physAlFrame1", logicContainer, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(xSi_2, ySi_2,zSi), logicAlFrame, "physAlFrame2", logicContainer, false, 1, true);

  

  // DRILLED SCINTILLATOR TOP
  G4double VetoTopTk = 1*cm;
  G4double zVetoTop = zSi - AlFrameTk/2. - VetoTopTk/2. - 1*mm;
  solidDrilledVeto = new G4Box("solidVetoTop",Displacing_Calo*2,Displacing_Calo*2, VetoTopTk/2.);


  // HOLES TO THE DRILLED VETO
  G4double rHole = rDet;
  G4double LHole = 10*cm;

  G4Tubs *solidHole = new G4Tubs("solidHole", 0, rHole, LHole, 0., 360.*deg);

  G4double theta = 30*deg;

  G4ThreeVector normal_SD1 = G4ThreeVector(-1, -1, 0);
  G4ThreeVector *center_SD1 = new G4ThreeVector(xSi_1, ySi_1,-zVetoTop+zSi);

  G4ThreeVector normal_SD2 = G4ThreeVector(+1, +1, 0);
  G4ThreeVector *center_SD2 = new G4ThreeVector(xSi_2, ySi_2,-zVetoTop+zSi);

  G4RotationMatrix *rot = new G4RotationMatrix();
  rot -> rotate(theta, normal_SD1);

  solidFinalDrilledVeto = new G4SubtractionSolid("solidDrilledvetoFinal",solidDrilledVeto, solidHole, rot, *center_SD1);

  rot = new G4RotationMatrix();
  rot -> rotate(theta, normal_SD2);
  solidFinalDrilledVeto = new G4SubtractionSolid("solidDrilledvetoFinal",solidFinalDrilledVeto, solidHole, rot, *center_SD2);

  logicDrilledVeto = new G4LogicalVolume(solidFinalDrilledVeto, EJ200, "logicDrilledVeto");
  visDrilledVeto = new G4VisAttributes(G4Color(0.,0.,1.,0.3));
  logicDrilledVeto -> SetVisAttributes(visDrilledVeto);

  new G4PVPlacement(0, G4ThreeVector(0,0,zVetoTop), logicDrilledVeto, "physDrilledVeto", logicContainer, false,0,true);





  // BOTTOM VETO SCINTILLATOR
  G4double VetoBottomTk = 1*cm;
  G4double zVetoBottom = zCalo + TkCalo/2. + VetoBottomTk/2. + 1*mm;
  solidBottomVeto = new G4Box("solidVetoBottom",Displacing_Calo*2,Displacing_Calo*2, VetoTopTk/2.);
  logicBottomVeto = new G4LogicalVolume(solidBottomVeto, EJ200, "logicBottomVeto");
  visBottomVeto = new G4VisAttributes(G4Color(0.,0.,1.,0.3));
  logicBottomVeto -> SetVisAttributes(visBottomVeto);

  new G4PVPlacement(0, G4ThreeVector(0,0,zVetoBottom), logicBottomVeto, "physBottomVeto", logicContainer, false,0,true);

  // LATERAL VETO SCINTILLATOR
  G4double VetoLateralTk = 1*cm;
  G4double xyVetoLateral = LCalo + 2*mm + VetoLateralTk/2.;
  G4double zVetoLateral = ((zVetoTop-VetoTopTk/2.) + (zVetoBottom+VetoBottomTk/2.))/2.;
  G4double VetoLateralHeight = - (zVetoTop-VetoTopTk/2.) + (zVetoBottom+VetoBottomTk/2.);

  solidLateralVeto = new G4Box("solidLateralVeto", 2*Displacing_Calo, VetoLateralTk/2., VetoLateralHeight/2.);
  logicLateralVeto = new G4LogicalVolume(solidLateralVeto, EJ200, "logicLateralVeto");
  visLateralVeto = new G4VisAttributes(G4Color(0.,0.,1.,0.3));
  logicLateralVeto -> SetVisAttributes(visLateralVeto);

  new G4PVPlacement(0, G4ThreeVector(0, xyVetoLateral, zVetoLateral), logicLateralVeto, "physLateralVeto", logicContainer, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(0, -xyVetoLateral, zVetoLateral), logicLateralVeto, "physLateralVeto", logicContainer, false, 1, true);

  G4RotationMatrix *rotation = new G4RotationMatrix();
  rotation -> rotateZ(90.*deg);

  new G4PVPlacement(rotation, G4ThreeVector(xyVetoLateral, 0, zVetoLateral), logicLateralVeto, "physLateralVeto", logicContainer, false, 2, true);
  new G4PVPlacement(rotation, G4ThreeVector(-xyVetoLateral,0,  zVetoLateral), logicLateralVeto, "physLateralVeto", logicContainer, false, 3, true);




  // ALUMINIUM TOP
  G4double AlTopTk = 7*mm;
  G4double zAlTop = zVetoTop - VetoTopTk/2. - AlTopTk/2. - 1*mm;
  G4double LxyAlTop = 2*LCalo + 2*VetoLateralTk + 3*mm;
  solidAlTop = new G4Box("solidAlTop",LxyAlTop/2., LxyAlTop/2., AlTopTk/2.);

  normal_SD1 = G4ThreeVector(-1, -1, 0);
  center_SD1 = new G4ThreeVector(xSi_1, ySi_1,-zAlTop+zSi);

  normal_SD2 = G4ThreeVector(+1, +1, 0);
  center_SD2 = new G4ThreeVector(xSi_2, ySi_2,-zAlTop+zSi);

  rot = new G4RotationMatrix();
  rot -> rotate(theta, normal_SD1);

  solidAlTopFinal = new G4SubtractionSolid("solidAlTopFinal",solidAlTop, solidHole, rot, *center_SD1);

  rot = new G4RotationMatrix();
  rot -> rotate(theta, normal_SD2);
  solidAlTopFinal = new G4SubtractionSolid("solidAlTopFinal",solidAlTopFinal, solidHole, rot, *center_SD2);


  logicAlTop = new G4LogicalVolume(solidAlTopFinal, Al, "logicAlTop");
  visAl = new G4VisAttributes(G4Color(0.4, 0.4, 0.4, 0.3));
  logicAlTop -> SetVisAttributes(visAl);
  new G4PVPlacement(0, G4ThreeVector(0.,0.,zAlTop), logicAlTop, "physAlTop", logicContainer, false, 0, true);


  // LATERAL ALLUMINIUM 1
  G4double TkLateral = 3*mm;
  G4double LxLateral_1 = LxyAlTop;
  G4double LxLateral_2 = LxyAlTop + 2*TkLateral;
  G4double LzAlLateral = (zVetoBottom + VetoBottomTk/2.) - (zAlTop - AlTopTk/2.);
  G4double zAlLateral = ((zVetoBottom + VetoBottomTk/2.) + (zAlTop - AlTopTk/2.))/2.;
  G4double yAlLateral_1 = LxyAlTop/2. + 1*mm + TkLateral/2.;
  G4double xAlLateral_2 = LxyAlTop/2. + 1*mm + TkLateral/2.;




  solidAlLateral1 = new G4Box("solidLateralAl1", LxLateral_1/2., TkLateral/2., LzAlLateral/2.);
  logicAlLateral1 = new G4LogicalVolume(solidAlLateral1, Al, "logicAlLateral1");
  logicAlLateral1 -> SetVisAttributes(visAl);
  new G4PVPlacement(0, G4ThreeVector(0., -yAlLateral_1, zAlLateral), logicAlLateral1, "physAlLateral_1", logicContainer, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(0., +yAlLateral_1, zAlLateral), logicAlLateral1, "physAlLateral_1", logicContainer, false, 1, true);

  solidAlLateral2 = new G4Box("solidLateralAl1", TkLateral/2., LxLateral_2/2., LzAlLateral/2.);
  logicAlLateral2 = new G4LogicalVolume(solidAlLateral2, Al, "logicAlLateral2");
  logicAlLateral2 -> SetVisAttributes(visAl);
  new G4PVPlacement(0, G4ThreeVector(-xAlLateral_2, 0.,zAlLateral), logicAlLateral2, "physAlLateral_2", logicContainer, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(+xAlLateral_2, 0., zAlLateral), logicAlLateral2, "physAlLateral_2", logicContainer, false, 1, true);



      // HoleRotation = G4RotationMatrix();
      // phi = std::atan2((iy-Delta_ii_y),(ix-Delta_ii_x));
      // distanceR = std::sqrt((iy-Delta_ii_y)*(iy-Delta_ii_y) + (ix-Delta_ii_x)*(ix-Delta_ii_x));
      // theta = distanceR * thetaMax / RMax;

        // HolePositionPCB = G4ThreeVector(xCurrentHole, yCurrentHole, zThin);
        // HoleRotationPCB = G4RotationMatrix();
        // HoleTransformPCB = G4Transform3D(HoleRotationPCB, HolePositionPCB);
        // solidPCBThinFinal = new G4SubtractionSolid("solidPCBThinFinal", solidPCBThin, solidHolePCB, HoleTransformPCB);
  

  new G4PVPlacement(0, G4ThreeVector(0,0,0), logicContainer, "Container", logicWorld, false, 0, true);

  return physWorld;  
}


void SPaRKLE_DetectorConstruction::ConstructSDandField()
{
  if(SDman == NULL)
  {

  
    SDman = G4SDManager::GetSDMpointer();

    names.push_back("Calo");
    volumes.push_back(logicCalo);

    names.push_back("Si");
    volumes.push_back(logicSiDet);

    names.push_back("VetoTop");
    volumes.push_back(logicDrilledVeto);

    names.push_back("VetoBottom");
    volumes.push_back(logicBottomVeto);

    names.push_back("VetoLateral");
    volumes.push_back(logicLateralVeto);
    

    
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


