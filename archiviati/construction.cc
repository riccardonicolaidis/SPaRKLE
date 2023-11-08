#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "construction.hh"
#include "myglobals.hh"


#if NEW_GEOMETRY == 0

using namespace std;

LEM_DetectorConstruction::LEM_DetectorConstruction()
{
  // Define the materials only once
  DefineMaterials();  
}

LEM_DetectorConstruction::~LEM_DetectorConstruction()
{}

void LEM_DetectorConstruction::DefineMaterials()
{
  G4NistManager *nist = G4NistManager::Instance();

  

  Si    = nist -> FindOrBuildElement("Si");
  SiMat = new G4Material("Si", 2.328*g/cm3, 1);
  SiMat -> AddElement(Si,1);

  EJ200 = nist -> FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  CdZnTe    = nist -> FindOrBuildMaterial("G4_CADMIUM_TELLURIDE");
  worldMat  = nist -> FindOrBuildMaterial("G4_Galactic");
  Al        = nist -> FindOrBuildMaterial("G4_Al");
  bachelite = nist -> FindOrBuildMaterial("G4_BAKELITE");

  
  std::vector<G4double> RIndexWorld   = {1.                , 1.                 , 1.};
  std::vector<G4double> RIndexSi      = {3.88163           , 3.88163            , 3.88163};
  std::vector<G4double> RIndexAl      = {3.88163           , 3.88163            , 3.88163};
  std::vector<G4double> RIndexEJ200   = {1.58              , 1.58               , 1.58};
  std::vector<G4double> RIndexCdZnTe  = {3.09              , 3.09               , 3.09};
 

  std::vector<G4double> Energy      = { 7.0 * eV, 7.07 * eV, 7.14 * eV };
  std::vector<G4double> SCINT       = { 0.1     , 1.0      , 0.1       };
  std::vector<G4double> ABSL        = { 35. * cm, 35. * cm , 35. * cm  };
  std::vector<G4double> ABSL_Al     = { 0.1 * um, 0.1 * um , 0.1 * um  };
  std::vector<G4double> ABSL_Si     = { 0.1 * um, 0.1 * um , 0.1 * um  };
  std::vector<G4double> ABSL_CdZnTe = { 10 * um , 10 * um  , 10 * um   };



  G4MaterialPropertiesTable *mptWorld   = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable *mptSi      = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable *mptEJ200   = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable *mptCdZnTe  = new G4MaterialPropertiesTable();
  G4MaterialPropertiesTable *mptAl      = new G4MaterialPropertiesTable();

  

  mptWorld -> AddProperty("RINDEX", Energy, RIndexWorld, 3);
  worldMat -> SetMaterialPropertiesTable(mptWorld);

  mptSi    -> AddProperty("RINDEX", Energy, RIndexSi, 3);
  mptSi    -> AddProperty("ABSLENGTH", Energy, ABSL_Si, 3);
  SiMat    -> SetMaterialPropertiesTable(mptSi);

  mptEJ200    -> AddProperty("RINDEX", Energy, RIndexEJ200, 3);
  mptEJ200    -> AddProperty("ABSLENGTH", Energy, ABSL, 3);
  
  if(OPTICAL_PROCESSES == 1)
  {
  mptEJ200 -> AddProperty("SCINTILLATIONCOMPONENT1", Energy, SCINT);
  mptEJ200 -> AddProperty("SCINTILLATIONCOMPONENT2", Energy, SCINT);
  mptEJ200 -> AddProperty("RINDEX", Energy, RIndexEJ200);
  mptEJ200 -> AddProperty("ABSLENGTH", Energy, ABSL);
  mptEJ200 -> AddConstProperty("SCINTILLATIONYIELD", 10000. / MeV);
  mptEJ200 -> AddConstProperty("RESOLUTIONSCALE", 1.0);
  mptEJ200 -> AddConstProperty("SCINTILLATIONTIMECONSTANT1", 20. * ns);
  mptEJ200 -> AddConstProperty("SCINTILLATIONTIMECONSTANT2", 45. * ns);
  mptEJ200 -> AddConstProperty("SCINTILLATIONYIELD1", 1.0);
  mptEJ200 -> AddConstProperty("SCINTILLATIONYIELD2", 0.0);
  EJ200    -> SetMaterialPropertiesTable(mptEJ200);
  EJ200    -> GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
  }
  else
  {
    EJ200    -> SetMaterialPropertiesTable(mptEJ200);
  }
  
  

  mptCdZnTe   -> AddProperty("RINDEX", Energy, RIndexCdZnTe, 3);
  mptCdZnTe   -> AddProperty("ABSLENGTH", Energy, ABSL_CdZnTe, 3);
  CdZnTe      -> SetMaterialPropertiesTable(mptCdZnTe);

  mptAl       -> AddProperty("RINDEX", Energy, RIndexAl, 3);
  mptAl       -> AddProperty("ABSLENGTH", Energy, ABSL_Al, 3);
  Al          -> SetMaterialPropertiesTable(mptAl);


}

G4VPhysicalVolume *LEM_DetectorConstruction::Construct()
{

  // Define a file to save the angles of the holes
  ofstream myfile;
  myfile.open ("../GeometryReport.txt");
  myfile << "Report of the geometry of the detector" << endl;
  myfile << "======================================" << endl;
  // Print the date and time at the execution of the file
  time_t now = time(0);
  char* dt = ctime(&now);

  myfile << "\nLocal Date and Time is: " << dt << endl;
  myfile << "Compilation Date and Time: " << __DATE__ << " " << __TIME__ << endl << endl;
  



  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  TotalMass = 0;

  // .oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo
  // .oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo
  //
  //        NUMERICAL DATA - PARAMETRIZATION
  //
  // .oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo
  // .oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo
  
  // World dimensions
  xWorld = .2*m;
  yWorld = .2*m;
  zWorld = .2*m;

  myfile << "World dimensions: " << 2*xWorld << " " << 2*yWorld << " " << 2*zWorld << endl;

  // Container position
  xContainer = 0. *cm;
  yContainer = 0. *cm;
  zContainer = 0. *cm; 

  myfile << "Container position: " << xContainer << " " << yContainer << " " << zContainer << endl;

  // Container dimensions 
  LxContainer = .19 *m;
  LyContainer = .19 *m;
  LzContainer = .19 *m; 

  myfile << "Container dimensions: " << 2*LxContainer << " " << 2*LyContainer << " " << 2*LzContainer << endl;

  // Dimensions of the detectors
  // ACTIVE AREA OF THE DETECTOR
  LThin         = 9*cm;
  LThick        = 9*cm;

  // THICKNESS OF THE MATERIAL
  TkThin        = TK_THIN;
  TkThick       = TK_THICK;

  myfile << "Thin detector thickness: " << TkThin << endl;
  myfile << "Thick detector thickness: " << TkThick << endl;
  


  // DATA OF THE GEOMETRY OF THE DETECTORS
  // SILICON DETECTORS

  rThinDet = 5. * mm;
  rExtAlDet = 9.7 * mm;
  TkAlDet = 7.9 *mm;
  rIntAlDet = rThinDet + 0.5*mm;

  myfile << "Silicon detector radius: " << rThinDet << endl;
  myfile << "Aluminum detector radius: " << rExtAlDet << endl;
  myfile << "Aluminum detector thickness: " << TkAlDet << endl;
  myfile << "Aluminum detector inner radius: " << rIntAlDet << endl;

  // ACTIVE LAYER POSITION
  xThin = 0.;
  yThin = 0.;
  zThin = 0.;

  myfile << "Thin silicon detector plane position: " << xThin << " " << yThin << " " << zThin << endl;

  xThick = 0.;
  yThick = 0.;
  zThick = zThin + TkThin/2. + TkAlDet/2. + TkAlDet/2. + 100*um;

  myfile << "Thick silicon detector plane position: " << xThick << " " << yThick << " " << zThick << endl;


  // Dimensions of the Veto
  LScintVeto  = 8*cm;
  TkScintVeto = 0.5*cm;

  myfile << "Veto detector size: " << 2*LScintVeto << endl;
  myfile << "Veto detector thickness: " << TkScintVeto << endl;


  xScintVeto = 0.;
  yScintVeto = 0.;
  zScintVeto = zThick + TkAlDet/2. + TkScintVeto/2. + 0.1*mm;

  myfile << "Veto detector plane position: " << xScintVeto << " " << yScintVeto << " " << zScintVeto << endl;

  // Construction of the Drilled Veto
  // BIGGER than the Veto on the back of the instrument
  LxDrilledVeto0 = LThin + 0.5*cm;
  LyDrilledVeto0 = LThin + 0.5*cm;
  TkDrilledVeto0 = 0.8*cm;
  LxDrilledVeto1 = LThin + 0.1*cm;
  LyDrilledVeto1 = LThin + 0.1*cm;
  TkCompenetration = 2* mm;

  xDrilledVeto = 0.;
  yDrilledVeto = 0.;
  zDrilledVeto = zThin - TkAlDet/2. - TkDrilledVeto0/2. - 0.1*mm;

  // Construction of the Drilled Aluminium
  // BIGGER than the Veto Drilled
  LxDrilledAl0 = LxDrilledVeto0 + 0.5 * cm;
  LyDrilledAl0 = LyDrilledVeto0 + 0.5 * cm;
  TkDrilledAl0 = 0.8 *cm;
  LxDrilledAl1 = LxDrilledVeto0 + 0.1 * cm;
  LyDrilledAl1 = LyDrilledVeto0 + 0.1 * cm;

  xDrilledAl = 0.;
  yDrilledAl = 0.;
  zDrilledAl = zDrilledVeto -TkDrilledVeto0/2. - TkDrilledAl0/2. - 0.1*mm ;

  // LENGTH OF THE SIDE WALLS OF THE ALUMINIUM AND THE DRILLED VETO
  LzDrilledVeto1 = (zScintVeto +  (N_PL_SCINT_NO_VETO + 1)*TkScintVeto) - zDrilledVeto - TkDrilledVeto0/2. + TkCompenetration;
  LzDrilledAl1 =  (zScintVeto + (N_PL_SCINT_NO_VETO + 1)*TkScintVeto) - zDrilledAl - TkDrilledAl0/2. + TkCompenetration;

  // HOLES DATA
  NxHoles = NxSensors;
  NyHoles = NySensors;
  rHoles = rThinDet;
  LengthHole = 3*cm;

 
  solidWorld = new G4Box("solidWorld",xWorld*0.5, yWorld*0.5, zWorld*0.5);  
  logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
  visWorld   = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visWorld->SetVisibility(false);
  logicWorld->SetVisAttributes(visWorld);
  physWorld  = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicWorld, "physWorld", 0, false, 0, true);
    
  solidContainer = new G4Box("solidContainer",LxContainer*0.5, LyContainer*0.5, LzContainer*0.5);  
  logicContainer = new G4LogicalVolume(solidContainer, worldMat, "logicContainer");
  visContainer   = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visContainer->SetVisibility(false);
  logicContainer->SetVisAttributes(visContainer);

  solidScintVeto = new G4Box("solidScintVeto", LScintVeto*0.5, LScintVeto*0.5, TkScintVeto*0.5);
  logicScintVeto = new G4LogicalVolume(solidScintVeto, EJ200, "logicScintVeto");
  visScintVeto   = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  logicScintVeto -> SetVisAttributes(visScintVeto);
  
  
  // Optical surface for the Drilled Veto
  G4MaterialPropertiesTable *mptSurface = new G4MaterialPropertiesTable();
  std::vector<G4double> energy        = {0.01*eV     ,      8.0 * eV    ,    1000* eV};
  std::vector<G4double> REFLECTIVITY  = { 1.              , 1.0              , 1.0 };
  mptSurface->AddProperty("REFLECTIVITY", energy , REFLECTIVITY, 3);
  opsurfVeto  = new G4OpticalSurface("Surface");
  opsurfVeto -> SetType(dielectric_metal);
  opsurfVeto -> SetFinish(polished);
  opsurfVeto -> SetModel(glisur);
  opsurfVeto -> SetMaterialPropertiesTable(mptSurface);

  
  myfile << "Plastic scintillator position" << endl;

  for(G4int i = 0; i < (N_PL_SCINT_NO_VETO + 1); ++ i)
  {
    if(i == 0)
    {
      physScintVeto  = new G4PVPlacement(0,G4ThreeVector(xScintVeto, yScintVeto, zScintVeto),logicScintVeto, "physScintVeto", logicContainer, false, i, true );
      if(OPTICAL_PROCESSES == 1)
      {
        G4String name = "SurfScintVeto" + std::to_string(0);
        logicSurfScintVeto[0] = new G4LogicalSkinSurface(name, logicScintVeto, opsurfVeto);
      }
      
    }
    else
    {
      physScintVeto  = new G4PVPlacement(0,G4ThreeVector(xScintVeto, yScintVeto, (zScintVeto + i*(2*um + TkScintVeto))),logicScintVeto, "physScintVeto", logicContainer, false, i, true );
      if(OPTICAL_PROCESSES == 1)
      {
        G4String name = "SurfScintVeto" + std::to_string(i);
        logicSurfScintVeto[i] = new G4LogicalSkinSurface(name, logicScintVeto, opsurfVeto);
      }
    }
    myfile << "Plastic scintillator position " << i << " " << zScintVeto + i*(2*um + TkScintVeto) << endl;
    TotalMass += (logicScintVeto->GetMass())/kg;
  }
  

  // ...oooOOOooo...oooOOOooo...oooOOOooo...oooOOOooo...oooOOOooo...oooOOOooo...
  // ...oooOOOooo...oooOOOooo...oooOOOooo...oooOOOooo...oooOOOooo...oooOOOooo...

  // SILICON DETECTORS DELTA E - E PADS 
  // LOGIC AND SOLID DEFINITIONS

  solidSiDetThin  = new G4Tubs("solidSiDetThin", 0.,rThinDet, TkThin*0.5, 0., 360. *deg);
  solidSiDetThick = new G4Tubs("solidSiDetThick", 0.,rThinDet, TkThick*0.5, 0., 360. * deg);
  

  logicSiDetThin = new G4LogicalVolume(solidSiDetThin, SiMat, "logicSiDetThin");
  TotalMass += (logicSiDetThin->GetMass())*16/kg;
  logicSiDetThick = new G4LogicalVolume(solidSiDetThick, CdZnTe, "logicSiDetThick");
  TotalMass += (logicSiDetThick->GetMass())*16/kg;
  visSiDetThin = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visSiDetThick = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  logicSiDetThin -> SetVisAttributes(visSiDetThin);
  logicSiDetThick -> SetVisAttributes(visSiDetThick);

  // Drilled Veto Solid without holes
  solidDrilledVeto0    = new G4Box("solidDrilledVeto0",LxDrilledVeto0*0.5, LyDrilledVeto0*0.5, TkDrilledVeto0*0.5);  
  solidDrilledVeto1Ext = new G4Box("solidDrilledVeto1Ext", LxDrilledVeto0*0.5, LyDrilledVeto0*0.5, (LzDrilledVeto1+TkCompenetration)*0.5);
  solidDrilledVeto1Int = new G4Box("solidDrilledVeto1Int", LxDrilledVeto1*0.5, LyDrilledVeto1*0.5, (LzDrilledVeto1+TkCompenetration + 10*um)*0.5);
  solidDrilledVeto1 = new G4SubtractionSolid("solidDrilledVeto1", solidDrilledVeto1Ext, solidDrilledVeto1Int);
  
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
  position1 = G4ThreeVector(0.,0.,(TkDrilledVeto0+LzDrilledVeto1)/2. -TkCompenetration);
  tr1 = G4Transform3D(Identity, position1);
  solidDrilledVeto = new G4UnionSolid("solidDrilledVeto", solidDrilledVeto0, solidDrilledVeto1, tr1);
  
  // Drilled Aluminium  
  Identity = G4RotationMatrix();
  position1 = G4ThreeVector(0.,0.,(TkDrilledAl0+LzDrilledAl1)/2. -TkCompenetration);
  tr1 = G4Transform3D(Identity, position1);
  solidDrilledAl = new G4UnionSolid("solidDrilledAl", solidDrilledAl0, solidDrilledAl1, tr1);
  

  xDelta = LThin / (NxHoles);
  yDelta = LThin / (NyHoles);

  G4RotationMatrix HoleRotation;
  G4ThreeVector HolePositionVeto;
  G4Transform3D HoleTransformVeto;
  G4ThreeVector HolePositionAl;
  G4Transform3D HoleTransformAl;
  G4ThreeVector RotationAxis;
  
  thetaMax = 35.*deg;
  RMax = std::sqrt(2) * 1.5;

  /* -------------------------------------------------------------------------- */
  /*                           DEFINITION OF THE FRAME                          */
  /* -------------------------------------------------------------------------- */

  
  rConnector = 3.8*mm;
  LConnector = 3.7 * mm;
  LTranslationConnector =  rExtAlDet + LConnector/3.;

  solidFrameCircle  = new G4Tubs("solidFrameCircle", rIntAlDet,rExtAlDet, TkAlDet*0.5, 0., 360. *deg);
  solidConnector    = new G4Tubs("solidConnector", 0., rConnector, LConnector*0.5, 0., 360. *deg); 
  
  
  int Index1[4][4] = {{1,1,1,1},{0,0,1,1},{0,0,1,1},{0,0,0,0}};
  int Index2[4][4] = {{0,0,0,1},{0,0,0,1},{0,1,1,1},{0,1,1,1}};
  
  // Define the transformation for the connector
  G4RotationMatrix RotationConnector = G4RotationMatrix();
  G4ThreeVector RotationAxisConnector = G4ThreeVector(-1.,1.,0.);
  RotationConnector.rotate(-90.*deg, RotationAxisConnector);
  RotationConnector.invert();
  G4ThreeVector PositionConnector = G4ThreeVector(LTranslationConnector/sqrt(2.),LTranslationConnector/sqrt(2.),0.);
  G4Transform3D TransformConnector = G4Transform3D(RotationConnector, PositionConnector);
  solidFrame[0][0] = new G4UnionSolid("solidFrame00", solidFrameCircle, solidConnector, TransformConnector);


  RotationConnector = G4RotationMatrix();
  RotationAxisConnector = G4ThreeVector(-1.,1.,0.);
  RotationConnector.rotate(90.*deg, RotationAxisConnector);
  RotationConnector.invert();
  PositionConnector = G4ThreeVector(-LTranslationConnector/sqrt(2.),-LTranslationConnector/sqrt(2.),0.);
  TransformConnector = G4Transform3D(RotationConnector, PositionConnector);
  solidFrame[1][1] = new G4UnionSolid("solidFrame11", solidFrameCircle, solidConnector, TransformConnector);

  RotationConnector = G4RotationMatrix();
  RotationAxisConnector = G4ThreeVector(1.,1.,0.);
  RotationConnector.rotate(90.*deg, RotationAxisConnector);
  RotationConnector.invert();
  PositionConnector = G4ThreeVector(-LTranslationConnector/sqrt(2.),LTranslationConnector/sqrt(2.),0.);
  TransformConnector = G4Transform3D(RotationConnector, PositionConnector);
  solidFrame[1][0] = new G4UnionSolid("solidFrame10", solidFrameCircle, solidConnector, TransformConnector);

  RotationConnector = G4RotationMatrix();
  RotationAxisConnector = G4ThreeVector(1.,1.,0.);
  RotationConnector.rotate(-90.*deg, RotationAxisConnector);
  RotationConnector.invert();
  PositionConnector = G4ThreeVector(LTranslationConnector/sqrt(2.),-LTranslationConnector/sqrt(2.),0.);
  TransformConnector = G4Transform3D(RotationConnector, PositionConnector);
  solidFrame[0][1] = new G4UnionSolid("solidFrame01", solidFrameCircle, solidConnector, TransformConnector);

  logicFrame[0][0] = new G4LogicalVolume(solidFrame[1][1], Al, "logicFrame00");
  TotalMass += (logicFrame[0][0]->GetMass())*16/kg;
  logicFrame[1][1] = new G4LogicalVolume(solidFrame[0][0], Al, "logicFrame11");
  logicFrame[1][0] = new G4LogicalVolume(solidFrame[1][0], Al, "logicFrame10");
  logicFrame[0][1] = new G4LogicalVolume(solidFrame[0][1], Al, "logicFrame01");

  visFrame[0][0] = new G4VisAttributes(G4Colour(0.8,0.0,0.1));
  visFrame[1][1] = new G4VisAttributes(G4Colour(0.8,0.2,0.1));
  visFrame[1][0] = new G4VisAttributes(G4Colour(0.8,0.5,0.1));
  visFrame[0][1] = new G4VisAttributes(G4Colour(0.8,0.7,0.1));

  logicFrame[0][0]->SetVisAttributes(visFrame[0][0]);
  logicFrame[1][1]->SetVisAttributes(visFrame[1][1]);
  logicFrame[1][0]->SetVisAttributes(visFrame[1][0]);
  logicFrame[0][1]->SetVisAttributes(visFrame[0][1]);




  // DEFINITION OF THE HOLE
  solidHole     = new G4Tubs("solidHole", 0., rHoles, LengthHole+1*cm, 0.*deg, 360.*deg);
  solidHoleLong = new G4Tubs("solidHoleLong", 0., rHoles, LengthHole+3*cm, 0.*deg, 360.*deg);
  

  // Print the directions of the holes in the myfile fstream
  myfile << "Angles of the holes" << endl;
  myfile << "kx\tky\ttheta\tphi" << endl;


  G4int CopyNoMatrix[(G4int)NxHoles][(G4int)NyHoles];
  G4int jCopyNo = 0;
  for(G4int iy = 0; iy < NyHoles; ++iy)
  {
    for(G4int ix = 0; ix <NxHoles; ++ix)
    {
      CopyNoMatrix[ix][iy] = (jCopyNo++);
      xCurrentHole = -(LThin/2.) + (xDelta/2.) + ix * xDelta;
      yCurrentHole = -(LThin/2.) + (yDelta/2.) + iy * yDelta;

      HoleRotation = G4RotationMatrix();
      phi = std::atan2((iy-1.5),(ix-1.5));
      distanceR = std::sqrt((iy-1.5)*(iy-1.5) + (ix-1.5)*(ix-1.5));
      theta = distanceR * thetaMax / RMax;

      G4cout << ix << "\t" << iy << "\t" << (theta*180/3.141592) << "\t" << (phi*180/3.141592) << G4endl;
      myfile << ix << "\t" << iy << "\t" << (theta*180/3.141592) << "\t" << (phi*180/3.141592) << endl;

      xCurrentHoleThick = xCurrentHole - (zThick-zThin)*std::tan(theta) * std::cos(phi);
      yCurrentHoleThick = yCurrentHole - (zThick-zThin)*std::tan(theta) * std::sin(phi);

      

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

      physSiDetThin   = new G4PVPlacement(0, G4ThreeVector(xCurrentHole,yCurrentHole,zThin) , logicSiDetThin, "physSiDetThin", logicContainer, false, CopyNoMatrix[ix][iy], true);
      physSiDetThick  = new G4PVPlacement(0, G4ThreeVector(xCurrentHoleThick,yCurrentHoleThick,zThick) , logicSiDetThick, "physSiDetThick", logicContainer, false, CopyNoMatrix[ix][iy], true);
      physFrame   = new G4PVPlacement(0, G4ThreeVector(xCurrentHole,yCurrentHole,zThin) , logicFrame[Index1[ix][iy]][Index2[ix][iy]], "physFrameThin", logicContainer, false, 2*CopyNoMatrix[ix][iy], true);
      physFrame   = new G4PVPlacement(0, G4ThreeVector(xCurrentHoleThick,yCurrentHoleThick,zThick) , logicFrame[Index1[ix][iy]][Index2[ix][iy]], "physFrame", logicContainer, false, 2*CopyNoMatrix[ix][iy]+1, true);
    
      //G4cout << "    ix = " + ix << "  " << "iy = " << iy << G4endl;
    }
  }


  
  
  // LOGIC DRILLED VETO
  logicDrilledVeto = new G4LogicalVolume(solidFinalDrilledVeto, EJ200, "logicDrilledVeto");
  TotalMass += (logicDrilledVeto->GetMass())/kg;
  visDrilledVeto = new G4VisAttributes(G4Colour(0.3,0.3,1.0));
  logicDrilledVeto->SetVisAttributes(visDrilledVeto);
  // Add skin surface
  //G4LogicalSkinSurface* skinDrilledVeto = new G4LogicalSkinSurface("skinDrilledVeto", logicDrilledVeto, opsurfVeto);

  // LOGIC DRILLED ALUMINIUM
  logicDrilledAl = new G4LogicalVolume(solidFinalDrilledAl, Al, "logicDrilledAl");
  TotalMass += (logicDrilledAl->GetMass())/kg;
  visDrilledAl = new G4VisAttributes(G4Colour(0.6,0.6,0.6));
  logicDrilledAl->SetVisAttributes(visDrilledAl);


  // PHYSICAL PLACEMENT DRILLED ALUMINIUM
  physDrilledAl  = new G4PVPlacement(0, G4ThreeVector(xDrilledAl,yDrilledAl,zDrilledAl) , logicDrilledAl, "physDrilledAl", logicContainer, false, 0, true);
  // PHYSICAL PLACEMENT DRILLED VETO
  physDrilledVeto  = new G4PVPlacement(0, G4ThreeVector(xDrilledVeto,yDrilledVeto,zDrilledVeto) , logicDrilledVeto, "physDrilledVeto", logicContainer, false, 0, true);
  if(OPTICAL_PROCESSES == 1)
  {
    logicSurfDrilledVeto = new G4LogicalSkinSurface("LogicSurfDrilledVeto", logicDrilledVeto, opsurfVeto);
  }
  

  physContainer  = new G4PVPlacement(0, G4ThreeVector(xContainer,yContainer,zContainer), logicContainer, "physContainer", logicWorld, false, 0, true);
  


  // CONTROL
  
  for(G4int M = 0; M <= 10; ++M)
  {
    G4cout << "#####################" << G4endl;
  }

  G4cout << "LThick =  " << TkThick << G4endl;
  G4cout << "LThin =  " << TkThin << G4endl;
  G4cout << "Total mass = " << TotalMass << " kg" << G4endl;

  for(G4int M = 0; M <= 10; ++M)
  {
    G4cout << "#####################" << G4endl;
  }


  myfile << "Total mass of the detector " << TotalMass << " kg" << G4endl; 

  if(GENERATE_GDML && (OPTICAL_PROCESSES == 0))
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
    std::string str1 = "../CadGeometry/LEM_geometry_";
    std::string str2 = ".gdml";
    std::string str3 = str1 + str + str2;

    parser.Write(str3, physContainer);


    str1 = "../CadGeometry/DrilledVeto_";
    str3 = str1 + str + str2;

    parser.Write(str3, physDrilledVeto);

  }
  

  
  myfile.close();

  return physWorld;  
}


void LEM_DetectorConstruction::ConstructSDandField()
{
  // Thin Silicon 100 um or 300 um 
  G4String ThinDetectorSD = "SiThin";
  LEM_SensitiveDetector* aThinDetectorSD = new LEM_SensitiveDetector(ThinDetectorSD, "SiThin");
  G4SDManager::GetSDMpointer()->AddNewDetector(aThinDetectorSD);
  SetSensitiveDetector( logicSiDetThin, aThinDetectorSD ); 
  
  // Thick detector 500 Si or 2 mm CZT
  G4String ThickDetectorSD = "SiThick";
  LEM_SensitiveDetector* aThickDetectorSD = new LEM_SensitiveDetector(ThickDetectorSD, "SiThick");
  G4SDManager::GetSDMpointer()->AddNewDetector(aThickDetectorSD);
  SetSensitiveDetector( logicSiDetThick, aThickDetectorSD ); 

  // Plastic scintillator VETO BACK
  G4String VetoSD = "Veto";
  LEM_SensitiveDetector* aVetoSD = new LEM_SensitiveDetector(VetoSD, "Veto");
  G4SDManager::GetSDMpointer()->AddNewDetector(aVetoSD);
  SetSensitiveDetector( logicScintVeto, aVetoSD ); 

  // Plastic scintillator VETO DRILLED
  G4String VetoDrilledSD = "VetoDrilled";
  LEM_SensitiveDetector* aVetoDrilledSD = new LEM_SensitiveDetector(VetoDrilledSD, "VetoDrilled");
  G4SDManager::GetSDMpointer()->AddNewDetector(aVetoDrilledSD);
  SetSensitiveDetector( logicDrilledVeto, aVetoDrilledSD ); 

}


#endif