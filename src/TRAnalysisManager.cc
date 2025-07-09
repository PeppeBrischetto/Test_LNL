#include "TRAnalysisManager.hh"
#include "G4AutoLock.hh"
#include "TError.h"
#include "G4Version.hh"
#include "G4SystemOfUnits.hh"


TRAnalysisManager* TRAnalysisManager::instance = 0;

namespace { 
  //Mutex to acquire access to singleton instance check/creation
  G4Mutex instanceMutex = G4MUTEX_INITIALIZER;
  //Mutex to acquire accss to histograms creation/access
  //It is also used to control all operations related to histos 
  //File writing and check analysis
  G4Mutex dataManipulationMutex = G4MUTEX_INITIALIZER;
}

TRAnalysisManager::TRAnalysisManager() : 
  fFile(0),fTree(0)
{;} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
TRAnalysisManager::~TRAnalysisManager()
{
  //No need to mutex, this is a real singleton.
  //loop over all histograms 
  if (fTree)
    delete fTree; 
  if (fFile) 
    delete fFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TRAnalysisManager* TRAnalysisManager::getInstance()
{
 
G4AutoLock l(&instanceMutex);
  if (instance == 0) 
    instance = new TRAnalysisManager();
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TRAnalysisManager::Book()
{
  //Booking of histograms has to be protected.
  //In addition there are issues with ROOT that is 
  //heavily non thread-safe. In particular I/O related operations
  //are not thread safe. To avoid problems let's mutex everything
  //here
  G4AutoLock l(&dataManipulationMutex);
  if (!fFile)
    {
      //Use the version name for the ROOT file
      TString filename = "ROOToutput.root"; 
      fFile = new TFile(filename,"RECREATE");
    }

  if (!fTree)
    {
      fTree = new TTree("tree","Global results");
      fTree->Branch("NCrystalsSiC",&fNCrystalsSiC,"NCrystalsSiC/I");
      fTree->Branch("EnergySiC",fEnergySiC,"EnergySiC[NCrystalsSiC]/D");
      fTree->Branch("IDSiC",fIDSiC,"IDSiC[NCrystalsSiC]/I");
      fTree->Branch("NCrystalsDL",&fNCrystalsDL,"NCrystalsDL/I");
      fTree->Branch("EnergyDL",fEnergyDL,"EnergyDL[NCrystalsDL]/D");
      fTree->Branch("IDDL",fIDDL,"IDDL[NCrystalsDL]/I");
      fTree->Branch("NCrystalsCsI",&fNCrystalsCsI,"NCrystalsCsI/I");
      fTree->Branch("EnergyCsI",fEnergyCsI,"EnergyCsI[NCrystalsCSI]/D");
      fTree->Branch("IDCsI",fIDCsI,"IDCsI[NCrystalsCsI]/I");

    }
  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TRAnalysisManager::AddEvent(std::vector<G4double> ene1,std::vector<G4int> ID1,
				 std::vector<G4double> ene2,std::vector<G4int> ID2,
				 std::vector<G4double> ene3,std::vector<G4int> ID3)
{
  G4AutoLock l(&dataManipulationMutex);
  fNCrystalsSiC = (Int_t) ene1.size();
  for (size_t i=0;i<ene1.size();i++)    
    {
      fEnergySiC[i] = ene1.at(i);
      fIDSiC[i] = ID1.at(i);
    }
  fNCrystalsDL = (Int_t) ene2.size();
  for (size_t i=0;i<ene2.size();i++)    
    {
      fEnergyDL[i] = ene2.at(i);
      fIDDL[i] = ID2.at(i);
    }

  fNCrystalsCsI = (Int_t) ene3.size();
  for (size_t i=0;i<ene3.size();i++)    
    {
      fEnergyCsI[i] = ene3.at(i);
      fIDCsI[i] = ID3.at(i);
    }


  fTree->Fill();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TRAnalysisManager::CloseFile()
{
  G4AutoLock l(&dataManipulationMutex);
  if (!fFile) //file not created at all: e.g. for a vis-only execution
    return;
  if (!fFile->IsOpen())
    {
      G4Exception("TRAnalysisManager::CloseFile()","tst67_02",FatalException,
                  "Trying to close a ROOT file which is not open");
      return;
    }
  fFile->cd(); 
  if (fTree)
   fTree->Write(fTree->GetName());
  fFile->Close();
}

