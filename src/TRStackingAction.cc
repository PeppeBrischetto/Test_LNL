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
// $Id: TRStackingAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
// 
/// \file TRStackingAction.cc
/// \brief Implementation of the TRStackingAction class

#include "TRStackingAction.hh"

#include "G4Track.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh" // GB 2021-05-10

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TRStackingAction::TRStackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TRStackingAction::~TRStackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
TRStackingAction::ClassifyNewTrack(const G4Track* track)
{
  //kill secondary neutrino
  if ( (track->GetDefinition() == G4NeutrinoE::NeutrinoE() || track->GetDefinition() == G4AntiNeutrinoE::AntiNeutrinoE() ) && 
      track->GetParentID()>0) // GB 2021-05-10
    return fKill;
  //otherwise, return what Geant4 would have returned by itself
  else 
    return G4UserStackingAction::ClassifyNewTrack(track);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
