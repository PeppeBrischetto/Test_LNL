//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// This is the modified version of the example /basic/B3 specifically tailored for a Geant4 school
// This code simulated a naked NaI detector irradiated by a point-like radioactive source



#include "G4UImanager.hh"
#include "G4RunManagerFactory.hh"

#include "Randomize.hh"

#include "TRDetectorConstruction.hh"
#include "TRPhysicsList.hh"
#include "TRActionInitialization.hh"
#include "TRAnalysisManager.hh"
#include "Shielding.hh"
#include "G4Timer.hh"
#include <time.h>
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  //Instantiate the G4Timer object, to monitor the CPU time spent for 
  //the entire execution
  G4Timer* theTimer = new G4Timer();
  //Start the benchmark
  theTimer->Start();

  //
  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4Random::setTheSeed(time(0));
    
  // Construct the default run manager. Pick the proper run 
  // manager depending if the multi-threading option is 
  // active or not.
  //
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  
  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(new TRDetectorConstruction);
  //
  runManager->SetUserInitialization(new Shielding()); 
  //runManager->SetUserInitialization(new TRPhysicsList);
    
  // Set user action initialization
  //
  runManager->SetUserInitialization(new TRActionInitialization());  
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
  
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // if an argument is given after the name of the executable 
  // (i.e. argc > 1), then take the argument as a Geant4 macro 
  // and execute it
  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  // otherwise (only the executable is given), start a user 
  // interface session. An initialization macro is executed 
  // by default. The macro which is executed depends on the 
  // activation (or not) of the visualization
  else
    {  // interactive mode : define UI session
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      if (ui->IsGUI())
       UImanager->ApplyCommand("/control/execute init_vis.mac"); 
      else
       UImanager->ApplyCommand("/control/execute init.mac"); 

      // start the session here: make the Geant4 prompt Idle>
      // available to the user
      ui->SessionStart();
      delete ui;
    }

  //Stop the benchmark here
  theTimer->Stop();

  G4cout << "The simulation took: " << theTimer->GetRealElapsed() << " s to run (real time)" 
	 << G4endl;

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  //save histograms
  TRAnalysisManager* man = TRAnalysisManager::getInstance();
  man->CloseFile();

  delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
