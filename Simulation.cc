#include <iostream>
#include <string>
#include <string_view>

#include "G4Types.hh"
#include "G4GenericPhysicsList.hh"
#include "G4VModularPhysicsList.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4EmLowEPPhysics.hh"


#include "SPaRKLE_Construction.hh"
#include "SPaRKLE_Action.hh"
#include "SPaRKLE_Globals.hh"

int main(int argc, char **argv)
{
  // Construct the default run manager
  G4UIExecutive *ui = nullptr;
  // Single Thread mode
  // auto* runManager = G4RunManagerFactory::CreateRunManager();

  G4RunManager *runManager = new G4RunManager();

  G4VModularPhysicsList *physicsList = new FTFP_BERT;
  //physicsList   -> ReplacePhysics(new G4EmLowEPPhysics());
  physicsList   -> ReplacePhysics(new G4EmStandardPhysics_option4());
  
  

  // Detector geometry and materials
  SPaRKLE_DetectorConstruction *det = new SPaRKLE_DetectorConstruction;

  runManager->SetUserInitialization(det);
  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new SPaRKLE_ActionInitialization(det));

  if (argc == 1)
  {
    ui = new G4UIExecutive(argc, argv);
  }

  G4VisManager *visManager = new G4VisExecutive();
  visManager->Initialize();

  G4UImanager *UImanager = G4UImanager::GetUIpointer();

  //runManager->Initialize();

  if (ui)
  {
    UImanager->ApplyCommand("/control/execute default_generator.mac");
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
  }
  else
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand("/control/execute default_generator.mac");
    UImanager->ApplyCommand(command + fileName);
  }

  delete visManager;
  delete runManager;

  return 0;
}
