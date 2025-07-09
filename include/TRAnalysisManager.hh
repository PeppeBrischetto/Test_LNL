#ifndef TRAnalysisManager_h
#define TRAnalysisManager_h 1

#include "globals.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>

class TRAnalysisManager
{
public:
  virtual ~TRAnalysisManager();
 

  ///method to call to create an instance of this class
  static TRAnalysisManager* getInstance();

  void Book();
  void AddEvent(std::vector<G4double>,std::vector<G4int>,
		std::vector<G4double>,std::vector<G4int>,
		std::vector<G4double>,std::vector<G4int>);
  void CloseFile();


private:
  static const Int_t MAX_DET=100;

  ///private constructor in order to create a singleton
  TRAnalysisManager();
  static TRAnalysisManager* instance; 

  TFile* fFile;
  TTree* fTree;

  ///Store energy of the individual crystals
  Double_t fEnergySiC[MAX_DET];

  ///How many crystals are fired
  Int_t fNCrystalsSiC;

  ///The ID's of the detector which have energy deposit
  Int_t fIDSiC[MAX_DET];

  ///Store energy of the individual crystals
  Double_t fEnergyDL[MAX_DET];
  Int_t fNCrystalsDL;
  Int_t fIDDL[MAX_DET];
  
  ///Store energy of the individual crystals
  Double_t fEnergyCsI[MAX_DET];
  Int_t fNCrystalsCsI;
  Int_t fIDCsI[MAX_DET];

};

#endif

