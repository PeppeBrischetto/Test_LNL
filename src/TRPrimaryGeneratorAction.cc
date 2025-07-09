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
// $Id: TRPrimaryGeneratorAction.cc 73766 2013-09-10 12:57:13Z gcosmo $
//
/// \file TRPrimaryGeneratorAction.cc
/// \brief Implementation of the TRPrimaryGeneratorAction class

#include "TRPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Proton.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <fstream>
#define twopi ((6.2831853072))

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TRPrimaryGeneratorAction::TRPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("chargedgeantino");
  fParticleGun->SetParticleDefinition(particle);
  // fixed position
      G4double x0 = 0*cm, y0= 0*cm;
      G4double z0 = 46.*cm;
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  
  G4double energy = 100. * MeV;
  fParticleGun->SetParticleEnergy(energy);
  
  // The default direction is the z-axis (i.e. towards the detector). 
  // However, if the primary particle is an unstable nucleus, Geant4 
  // will take care of the production of the final decay state, and the 
  // products will be emitted isotropically.
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  
  
  fCosTheta = 0.;
  fKinE = 0.;
  fxf = 0., fyf = 0., ftf = 0., fphf = 0., ftpol = 0., fphpol = 0., ftof = 0.;
 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TRPrimaryGeneratorAction::~TRPrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TRPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{  
  G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
 
  // If the primary particle is defined to be a charged geantino (default), 
  // a Oxygen 18 nucleus is generated instead. The primary particle can be 
  // overridden at run time by the command /gun/particle
  //
  
  if (particle == G4ChargedGeantino::ChargedGeantino()) {
    // Oxygen-18 or Neon-18
    G4int Z = 10, A = 18;
    
    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 0.*keV;
    
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);

    fParticleGun->SetParticleDefinition(ion);
    // fParticleGun->SetParticleEnergy(35.*MeV); //at rest
    fParticleGun->SetParticleCharge(ionCharge);
  }
  
     //posizione di sparo
  G4double zOfSource = -46.0*cm;
  G4double xOfSource = -0.*mm;
  G4double yOfSource = -0.*mm;
  
  fParticleGun -> SetParticlePosition(G4ThreeVector(xOfSource, yOfSource, zOfSource));
  
  /* Randomizzazione del punto di generazione della particella: 
         ogni particella viene generata in un piano adiacente alla faccia d'ingresso dei rivelatori.
         Il piano deve essere grande tanto quanto la faccia del muro */
  
  //G4double xOfSource = (-5. + (5. + 5.)*G4UniformRand())*mm;    // Rivelatore 1cm*1cm 
  //G4double yOfSource = (-5. + (5. + 5.)*G4UniformRand())*mm;    // Rivelatore 1cm*1cm
  //G4double xOfSource = (17. + (7. - 17.)*G4UniformRand())*mm;
  //G4double xOfSource = (-7.5 + (7.5 + 7.5)*G4UniformRand())*mm;   // Rivelatore 1.5cm*1.5cm
  //G4double yOfSource = (-7.5 + (7.5 + 7.5)*G4UniformRand())*mm;   // Rivelatore 1.5cm*1.5cm
  //G4double xOfSource = (-10. + (10. + 10.)*G4UniformRand())*mm;   // Rivelatore 2cm*2cm
  //G4double yOfSource = (-10. + (10. + 10.)*G4UniformRand())*mm;   // Rivelatore 2cm*2cm



    
  
  //fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(0,0,1));
  
  //Randomizzazione della direzione della particella
  
  //Proposta Lombardo
  // 3. Direzione casuale entro il cono
    G4double radius = 3.4435 * cm;    // raggio effettivo
    G4double distance = 46. * cm;     // distanza sorgente-rivelatore
    G4double theta_max = std::atan(radius / distance);  

    // Genera theta e phi casuali nel cono solido
    G4double cosThetaMin = std::cos(theta_max);
    G4double cosTheta = cosThetaMin + (1.0 - cosThetaMin) * G4UniformRand();  // uniforme in cos(θ)
    G4double theta = std::acos(cosTheta);
    G4double phi = 2. * CLHEP::pi * G4UniformRand();

    // Conversione in vettore di direzione (rispetto all’asse Z)
    G4double dirX = std::sin(theta) * std::cos(phi);
    G4double dirY = std::sin(theta) * std::sin(phi);
    G4double dirZ = std::cos(theta);

    G4ThreeVector direction(dirX, dirY, dirZ);
    fParticleGun->SetParticleMomentumDirection(direction);

    // 4. Imposta energia della particella
    fParticleGun->SetParticleEnergy(250.*MeV);  // o qualsiasi valore ti serva

    // 5. Genera il vertice primario
    //fParticleGun->GeneratePrimaryVertex(anEvent);
    
    //Usato da Brischetto in passato
  //G4double cosTheta = G4UniformRand(); // Questo per una direzione da 0° a 90°
  //G4double cosTheta = 0.939693 + (1. - 0.939693)*G4UniformRand();
  //G4double phi = twopi*G4UniformRand();
  //G4double sinTheta = sqrt(1-cosTheta*cosTheta);
  //these are the cosines for an isotropic direction
  //fCosTheta = cosTheta;
  //fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(sinTheta*cos(phi), sinTheta*sin(phi),cosTheta));
  // Fino qui 
  
  //fCosTheta = 1.;
  //fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));


  //G4double energy = 250. * MeV;
  //fParticleGun->SetParticleEnergy(energy);

  // Adesso facciamo variare l'energia di 1 GeV dello 0.5%
  //G4double energy = (995. + (1005. - 995.)*G4UniformRand())*MeV;
  //G4double energy = (900.673 + (910.220 - 900.673)*G4UniformRand())*MeV;
  //G4double energy = (500. + (1000. - 500.)*G4UniformRand())*MeV;
  //G4double energy = 800.*MeV;

  //fParticleGun -> SetParticleEnergy(energy);


  
  G4int evID = anEvent->GetEventID();


  //create vertex
  fParticleGun->GeneratePrimaryVertex(anEvent);

}
