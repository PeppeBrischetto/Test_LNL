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
// $Id: TRPrimaryGeneratorAction.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file TRPrimaryGeneratorAction.hh
/// \brief Definition of the TRPrimaryGeneratorAction class

#ifndef TRPrimaryGeneratorAction_h
#define TRPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

/// The primary generator action class with particle gun.
///
/// It defines (as default) a Cs-137 nucleus, at rest, from a point-like
/// source at (0,0,-6) cm. All parameters (particle type, energy, direction, 
/// position) can be changed via the UI commands /gun/...
///

class TRPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    TRPrimaryGeneratorAction();    
    virtual ~TRPrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);         

    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
    G4double GetCosTheta() const {return fCosTheta;};
    G4double GetKinE() const {return fKinE;};
    G4double Getxf() const {return fxf;};
    G4double Getyf() const {return fyf;};
    G4double Gettf() const {return ftf;};
    G4double Getphf() const {return fphf;};
    G4double Gettpol() const {return ftpol;};
    G4double Getphpol() const {return fphpol;};
    G4double Gettof() const {return ftof;};

  
  private:
    G4ParticleGun*  fParticleGun;
    G4double fCosTheta;
    G4double fxf, fyf, ftf, fphf, ftpol, fphpol;
    G4double fKinE;
    G4double ftof;
    std::vector<G4double> fBufferxf;
    std::vector<G4double> fBufferyf;
    std::vector<G4double> fBuffertf;
    std::vector<G4double> fBufferphf;
    std::vector<G4double> fBufferKinEnergy;
    std::vector<G4double> fBuffertof;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


