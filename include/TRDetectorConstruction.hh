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
// $Id: TRDetectorConstruction.hh 71323 2013-06-13 16:54:23Z gcosmo $
//
/// \file TRDetectorConstruction.hh
/// @brief Definition of the TRDetectorConstruction class (Mandatory)

#ifndef TRDetectorConstruction_h
#define TRDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4UserLimits;

/// Detector construction class to define materials (with their physical properties) and detector geometry.
///

class TRDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    /// constructor
    TRDetectorConstruction();
    /// destructor
    virtual ~TRDetectorConstruction();

  public:
    /// Defines the detector geometry and returns a pointer to the physical World Volume
    virtual G4VPhysicalVolume* Construct();
    /// Register some of the detector's volumes as "sensitive"
    virtual void ConstructSDandField();
               
  private:
    /// Defines all the materials the detector is made of.
    void DefineMaterials();

    G4bool  fCheckOverlaps;
	 G4UserLimits* fStepLimit; // pointer to user step limits
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
