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
// $Id: TRRunAction.cc 71323 2013-06-13 16:54:23Z gcosmo $
//
/// \file TRRunAction.cc
/// \brief Implementation of the TRRunAction class
#include "TRAnalysisManager.hh"
#include "TRRunAction.hh"
#include "TRPrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "TRAnalysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TRRunAction::TRRunAction()
 : G4UserRunAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TRRunAction::~TRRunAction()
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  //delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TRRunAction::BeginOfRunAction(const G4Run* run)
{ 
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  

  // Create analysis manager
  // Notice: it must be done the same way in master and workers
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  //TRAnalysisManager*  analysisManager = TRAnalysisManager::getInstance();
  //analysisManager->Book();

  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstNtupleId(1);
  
  //Uncomment these lines in task 4b
  //Create a one-column ntuple
  analysisManager->CreateNtuple("tree", "Energy and position");
  // 1) total energy released in the crystals (double), MeV
  // Ricordati che devi aggiungere una colonna con l'EventID
  analysisManager->CreateNtupleIColumn("evtID");
  analysisManager->CreateNtupleIColumn("logicalVolume");
  analysisManager->CreateNtupleIColumn("copyNumber");  
  analysisManager->CreateNtupleDColumn("eDep");
  analysisManager->CreateNtupleDColumn("x");
  analysisManager->CreateNtupleDColumn("y");
  analysisManager->CreateNtupleDColumn("z");
  analysisManager->CreateNtupleDColumn("cosThetaIn");
  analysisManager->CreateNtupleDColumn("KinE");
  analysisManager->CreateNtupleDColumn("xf");
  analysisManager->CreateNtupleDColumn("yf");
  analysisManager->CreateNtupleDColumn("tf");
  analysisManager->CreateNtupleDColumn("phf");
  analysisManager->CreateNtupleDColumn("tpol");
  analysisManager->CreateNtupleDColumn("phpol");
  analysisManager->CreateNtupleDColumn("tof");
  analysisManager->CreateNtupleDColumn("time");
  analysisManager->CreateNtupleDColumn("stepLen");
 
  //ok, done
  analysisManager->FinishNtuple();
  
  /*
  analysisManager->CreateH1("h1","Energy",2000,0.,2000.);
  analysisManager->CreateH1("h2","Energy",2000,0.,2000.);
  analysisManager->CreateH1("h3","Energy",2000,0.,40000.);
  */

  // Create a new output file
  analysisManager->OpenFile("TR");
  

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TRRunAction::EndOfRunAction(const G4Run* run)
{
  //retrieve the number of events produced in the run
  G4int nofEvents = run->GetNumberOfEvent();

  //do nothing, if no events were processed
  if (nofEvents == 0) return;
  
  // Run conditions
  // This retrieves the UserPrimaryGeneratorAction object: it is retrieved through the 
  // G4RunManager. 
  //
  // Following the TRActionInitialization, the UserPrimaryGeneratorAction 
  // exists for all workers (-> Build()) but not for the master (-> BuildForMaster()). The 
  // TRRunAction instead exists for the master and for the worker. So, when the function is 
  // executed by the master, no pointer is found for the primary generator and 
  // generatorAction = NULL

  const TRPrimaryGeneratorAction* generatorAction = static_cast<const TRPrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
 

  G4String partName;
  if (generatorAction) 
  {
      // The GetParticleGun() is defined inside the TRPrimaryGeneratorAction.hh file and
      // it returns a pointer to the current concrete implementation of the G4VPrimaryGeneratorAction
      // Note that GetParticleDefinition return a 'const' type in the case of the use of the G4ParticleGun but not
      // in the case of the G4ParticleGeneralSource. Hence the GetParticleGun in the TRPrimaryGeneratorAction.hh
      // must not be defined 'const' in the case of the GeneralParticleSource
      //
    G4ParticleDefinition* particle = generatorAction -> GetParticleGun() -> GetParticleDefinition();
    partName = particle->GetParticleName();
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
