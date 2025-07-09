/**
* @file:TRSteppingAction.cc
* @brief:
*
* 
*/

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Ions.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "TRSteppingAction.hh"
#include "TRPrimaryGeneratorAction.hh"
#include "TRAnalysis.hh"
#include "Randomize.hh"

TRSteppingAction::TRSteppingAction(TRPrimaryGeneratorAction* pri) : 
  G4UserSteppingAction(), fNavigator(nullptr), fGen(pri)
{;}
				       
  
void TRSteppingAction::UserSteppingAction (const G4Step *aStep)
{
  if (!aStep)
    return;

  /*
  G4double lPreEnergy = 0;
  G4double lPostEnergy = 0;
  G4double lLength = 0;
      			
  lPreEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
  lPostEnergy = aStep->GetPostStepPoint()->GetKineticEnergy();
  lLength = aStep->GetStepLength();
  */
  G4double edep = aStep->GetTotalEnergyDeposit();
  //if (!edep)  // GB 2021-05-17
     //return;  // GB 2021-05-17

  G4ThreeVector theGlobalPoint = aStep->GetTrack()->GetPosition();
  if (aStep->GetTrack()->GetVolume() == nullptr)
    return;

  if (!fNavigator)
    fNavigator = 
      G4TransportationManager::GetTransportationManager()->
      GetNavigatorForTracking();

  G4ThreeVector theLocalPoint = fNavigator->
    GetGlobalToLocalTransform().
    TransformPoint(theGlobalPoint);
  
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  G4String volumeName = aStep->GetTrack()->GetVolume()->GetName();
  G4LogicalVolume* targetLV = aStep->GetTrack()->GetVolume()->GetLogicalVolume();
  // G4LogicalVolume* targetLV = 0;

  G4int eID = 0;
  const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
  if (evt)
    eID = evt->GetEventID();

  G4int copyNo = -1;
  G4int depth = -1;
  G4int moduleCopyNo = -1;
  G4int detectorCopyNo = -1;
  //G4cout << "depth = " << depth << G4endl;

  //moduleCopyNo = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetVolume(1)->GetCopyNo();
  //detectorCopyNo = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetVolume(0)->GetCopyNo();
  depth = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistoryDepth();
  
  
  //if (targetLV->GetName()=="windowLV")
  //G4cout << "depth = " << depth << G4endl;
  //G4int detectorCopyNo = aStep->GetTrack()->GetTouchable()->GetCopyNumber(0);
  //G4int moduleCopyNo = aStep->GetTrack()->GetTouchable()->GetCopyNumber(1);
  //G4int copyNo = 1000*moduleCopyNo + detectorCopyNo;

  for (G4int iv=0; iv<=aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistoryDepth(); iv++) {
      //G4cout << "Level " << iv << " -> " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()
	//							->GetVolume(iv)->GetName() << G4endl;
  }

  if ( depth == 2 ) {
     moduleCopyNo = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetVolume(1)->GetCopyNo();
     detectorCopyNo = aStep->GetTrack()->GetVolume()->GetCopyNo();
     //G4cout << "Depth=2  Module: " << moduleCopyNo << " - detector: " << detectorCopyNo << G4endl;
  }

  if ( depth == 3 ) {
     moduleCopyNo = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetVolume(1)->GetCopyNo();
     detectorCopyNo = aStep->GetPreStepPoint()->GetTouchableHandle()->GetHistory()->GetVolume(2)->GetCopyNo();
     //G4cout << "Depth=3  Module: " << moduleCopyNo << " - detector: " << detectorCopyNo << G4endl;
  }


  //if (targetLV->GetName()=="windowLV") 
  //	G4cout << "depth = " << depth << "\t " << moduleCopyNo << "\t" << detectorCopyNo << G4endl;

  

  copyNo = 1000*moduleCopyNo + detectorCopyNo;

  G4StepPoint* p1 = aStep->GetPreStepPoint();
  G4StepPoint* p2 = aStep->GetPostStepPoint();

  G4ThreeVector coord1 = p1->GetPosition();
  G4ThreeVector coord2 = p2->GetPosition();

  const G4AffineTransform transformation = p1->GetTouchable()->
                                           GetHistory()->GetTopTransform();

  G4ThreeVector localPosition1 = transformation.TransformPoint(coord1);
  G4ThreeVector localPosition2 = transformation.TransformPoint(coord2);

  // This part is used for generating a random point 
  // between the PreStepPoint and the PostStepPoint
  G4ThreeVector PreStepPoint = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector PostStepPoint = aStep->GetPostStepPoint()->GetPosition();
  
  G4double r = G4UniformRand();

  G4ThreeVector RandomPoint = PreStepPoint + r*( PostStepPoint - PreStepPoint );    

  G4ThreeVector randomPoint = localPosition1 + r*( localPosition2 - localPosition1 );
  
  G4double KinE = aStep->GetTrack()->GetKineticEnergy();

  G4double time = -10.; // RP

  G4double cosThetaIn = fGen->GetCosTheta();
  //G4double KinE = fGen->GetKinE();
  G4double xf = fGen->Getxf();
  G4double yf = fGen->Getyf();
  G4double tf = fGen->Gettf();
  G4double phf = fGen->Getphf();
  G4double tpol = fGen->Gettpol();
  G4double phpol = fGen->Getphpol();
  G4double tof = fGen->Gettof();

  xf = xf - 19.*mm*tan(tf) + 620.*mm;  // RP 
  yf = yf + 19.*mm*tan(phf);           // RP

  time =  aStep->GetPreStepPoint()->GetGlobalTime(); // RP
  G4double stepLen = aStep->GetStepLength();                  // RP

  //G4double xCheck = -1000.*m;
  //G4double yCheck = -1000.*m;

  G4String parName = aStep->GetTrack()->GetDefinition()->GetParticleName(); 

  //if (volumeName != "World" ) // GB 2021-05-18
  {
 
     //if(edep == 0 && targetLV->GetName() == "SicLV" && parName != "gamma" && parName != "neutron")     // GB 2021-05-10
      //G4cout << "AAAAAAAAAAAA " << parName 
      //       << "\t step length " << aStep->GetStepLength() 
      //       << "\t kin energy " << aStep->GetTrack()->GetKineticEnergy() << G4endl;   // GB 2021-05-10

    /*
    G4cout << " Pre = ( " << PreStepPoint.getX()/um << ", " << PreStepPoint.getY()/um << ", " << PreStepPoint.getZ()/um<< " ) \n"
           << " Post = ( " << PostStepPoint.getX()/um << ", " << PostStepPoint.getY()/um << ", " << PostStepPoint.getZ()/um << " ) \n"
           << " Rand = ( " << RandomPoint.getX()/um << ", " << RandomPoint.getY()/um << ", " << RandomPoint.getZ()/um << " )" 
           << G4endl;
   */
    //G4cout << "EvID " << eID << " : eDep in " << targetLV->GetName() << " (module: " << moduleCopyNo << " det: " 
    //       << detectorCopyNo << "; copyNo = " << copyNo <<") in pos ("   << randomPoint.getX() / um << ", " << 
    //          randomPoint.getY() / um << ", " << randomPoint.getZ() / um << ") is " << edep / keV << " keV" << G4endl;
   
    // Adesso provo a riempire il tree. 
    analysis->FillNtupleIColumn(0, eID);

    //  Poiche' il metodo FillNtupleDColumn accetta un int 
    // e un double, ho pensato di risolvere il problema facendogli restituire 1 per il SiC,
    // 2 per il deadlayer e 3 per il CsI
    G4String namePV = targetLV->GetName();
    if ( namePV == "moduleLV") 
      analysis->FillNtupleIColumn(1, 0);
    else if ( namePV == "SicLV") 
      analysis->FillNtupleIColumn(1, 1);
    else if (namePV == "LaydeaLV") 
      analysis->FillNtupleIColumn(1, 2);
    else if (namePV == "CsILV") 
      analysis->FillNtupleIColumn(1, 3);
    else if (namePV == "gridLV") 
      analysis->FillNtupleIColumn(1, 4);
    else if (namePV == "World") 
      analysis->FillNtupleIColumn(1, 7);

//    else if (namePV == "resLV"){			// GB 2022-11-11 Il numero 10 è usato nell'EventAction
//      analysis->FillNtupleIColumn(1, 98);		// GB 2022-11-11
      //G4cout <<targetLV->GetName() <<G4endl;	// GB 2022-11-11
//    }	
    else
      analysis->FillNtupleIColumn(1, 15);  //Ordine numerico crescente, stesso ordine del TRrunaction

    analysis->FillNtupleIColumn(2, copyNo);
    analysis->FillNtupleDColumn(3, edep/keV);
    analysis->FillNtupleDColumn(4, randomPoint.getX() / um);
    analysis->FillNtupleDColumn(5, randomPoint.getY() / um);
    analysis->FillNtupleDColumn(6, randomPoint.getZ() / um);
    analysis->FillNtupleDColumn(7, KinE / MeV);
    analysis->FillNtupleDColumn(8, stepLen/um);    
/*    
    analysis->FillNtupleDColumn(7, cosThetaIn);
    analysis->FillNtupleDColumn(9, xf / m);
    analysis->FillNtupleDColumn(10, yf / m);
    analysis->FillNtupleDColumn(11, tf);
    analysis->FillNtupleDColumn(12, phf);
    analysis->FillNtupleDColumn(13, tpol);
    analysis->FillNtupleDColumn(14, phpol);
    analysis->FillNtupleDColumn(15, tof);
    analysis->FillNtupleDColumn(16, time); // RP
    analysis->FillNtupleDColumn(17, stepLen/um);
*/    
    analysis->AddNtupleRow();
  }

  
  //G4cout << aStep->GetTrack()->GetVolume()->GetName() << " Glob" << theGlobalPoint/mm << "; theLocalPoint " 
  //       << theLocalPoint/mm << "; localPosition2 " << localPosition2/mm << G4endl;
  
  //G4cout << "Cos(theta) is " << cosThetaIn << "; so theta is " << std::acos(cosThetaIn)/rad << G4endl;

  
  //This Stepping action kills long-lived nuclei (they do not decay)
  G4String particleType =aStep->GetTrack()->GetDefinition()->GetParticleType(); 
  //if (particleType == "nucleus" && aStep->GetTrack()->GetParentID()>0) // GB 2021-05-17
  if (particleType == "nucleus")
    {
      G4double energy = aStep->GetTrack()->GetKineticEnergy();
      if (energy < 0.1*keV)
	{
	  G4Ions* ion = (G4Ions*) aStep->GetTrack()->GetDefinition();
	  G4double lifetime = ion->GetPDGLifeTime();
	  G4double excitationEnergy = ion->GetExcitationEnergy();
	  // stable and excited nuclei --> track them as usual
	  if (lifetime < 0 || excitationEnergy > 0) return;
	  if (lifetime > 1.0*microsecond) //kill long-lived nuclei
	    {
	      G4String particleName = aStep->GetTrack()->GetDefinition()->GetParticleName(); 
	      aStep->GetTrack()->SetTrackStatus(fStopAndKill);
	    }   //stable nuclei are unaffected 
	}
    }
  
}



