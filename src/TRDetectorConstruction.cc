
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
// $Id: TRDetectorConstruction.cc 71323 2013-06-13 16:54:23Z gcosmo $
//
/// \file TRDetectorConstruction.cc
/// \brief Implementation of the TRDetectorConstruction class

#include "TRDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"             // RP
#include "G4PSNofSecondary.hh"            // RP
#include "G4SDChargedFilter.hh"           // RP
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh" // RP
#include "G4Region.hh"
#include "G4RegionStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TRDetectorConstruction::TRDetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true)
{
  // **Material definition**
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TRDetectorConstruction::~TRDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TRDetectorConstruction::DefineMaterials()
{

  // **Materials from the NIST database**
  G4NistManager* man = G4NistManager::Instance();
   
  G4bool isotopes = false;
 
  // Element definition
  // Retrieve materials required to build the material
  G4Element* Si = man->FindOrBuildElement("Si" , isotopes); 
  G4Element*  C = man->FindOrBuildElement("C", isotopes);
  G4Element* O = man->FindOrBuildElement("O", isotopes);
  G4Element* H = man->FindOrBuildElement("H", isotopes);

  // Declare the material with its density and number of components
  G4Material* SiC = new G4Material("SiC", //its name
				   3.22*g/cm3, //its density
				   2); //number of components
//resina EPOXY resin EPOTEK 301-1
 G4Material* res = new G4Material ("res",
 						1.16*g/cm3,
 						3);
 res->AddElement(C,19);
 res->AddElement(H,20);
 res->AddElement(O,4);
 
  //Add Element for Material "NaI" specifiyng the number of each element
  SiC->AddElement(Si,1);
  SiC->AddElement(C,1);
  


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* TRDetectorConstruction::Construct()
{  
 
  //
  G4NistManager* nist = G4NistManager::Instance();

  // **Retrieve Nist Materials** 
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* cryst_semi   = nist->FindOrBuildMaterial("SiC");      
  G4Material* cryst_scin   = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  G4Material* grid_mat = nist->FindOrBuildMaterial("G4_BRASS");
  G4Material* res_material = nist->FindOrBuildMaterial("res");
  G4Material* window_mat = nist->FindOrBuildMaterial("Mylar"); 	// GB 2022-11-11
  G4Material* poly_frame   = nist->FindOrBuildMaterial("G4_NYLON-6/6"); // GB 2025-07-15  
  
  //     
  // G4cout << "Pointer to material: " << cryst_scin << G4endl;
  //
  // Create the world volume as a box. This is big enough to contain 
  // all the detector setup
  //
  G4double world_sizeXYZ = 2.0*m;
  
  G4Box* solidWorld =    
    new G4Box("World",                                      //its name
	      world_sizeXYZ, world_sizeXYZ, world_sizeXYZ); //its size
  
  // World Logical Volume definition    
  // RP changed world material to isobutane (from default_mat)
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        default_mat,             //its material
                        "World");            //its name
                      
  // World Physical Volume Placement at (0,0,0)              
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps 



  //
  //  ***** Module envelope ****

  G4double moduleX = 3.08002*cm;//3.12002*cm;
  G4double moduleY = 15.44002*cm;
  G4double moduleZ = 1.11201*cm;
  G4Box* moduleSolid = new G4Box("moduleSolid", moduleX/2., moduleY/2., moduleZ/2.);

  // Module Logical Volume definition    
  G4LogicalVolume* moduleLogical = new G4LogicalVolume(moduleSolid,          //its solid
						    default_mat,        //its m
						    "moduleLV");          //its name
  
  G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
  G4double alpha = 0.*deg; //35.*deg;
  //G4cout << "pi" << pi << "\t sin(alpha)= " << sin(alpha) << "\t cos(alpha)= " << cos(alpha) << G4endl;
  rotationMatrix->rotateY(-1*alpha);

  G4int nModuleColumns = 1; // il numero è 36 a 35deg, 39 a 0deg
  G4int nModuleRows = 1;
  // G4double moduleDistanceY = 0.5*mm; // This value is for 2 rows of modules
  G4double moduleDistanceY = 0.*mm;
  G4double moduleDistanceX = 4.*mm;  // Prima era 5 mm; il valore per 35deg è 4*mm
  G4double manualDistanceX = 0.*cm;  // Prima era 25 cm

  G4double ky=0.;

  for (G4int iColumn=0; iColumn<nModuleRows; iColumn++) 
  {
            for (G4int iRow=0; iRow<nModuleColumns; iRow++)
            {
                G4int kx = iRow - nModuleColumns/2;
                if (nModuleRows%2==0)
                   ky = iColumn - 0.5;
                else
                   ky = iColumn - nModuleRows/2.;
                //G4double y = (ky)*(moduleY+moduleDistanceY); // This value is for 2 rows of modules
                G4double y = 0.;
                G4double x = (kx)*(moduleX+moduleDistanceX) + manualDistanceX;
                char volumename[100];
                sprintf(volumename, "module_row%d_col%d", iColumn, iRow);
                G4int copyNo = iRow + iColumn*nModuleColumns;
	         new G4PVPlacement(rotationMatrix,                //no rotation
                            G4ThreeVector(x,y,0.),  //at (0,0,0)
			       moduleLogical,             //its logical volume
                            volumename,             //its name
                            logicWorld,             //its mother  volume
                            false,                  //no boolean operation
                            copyNo,                 //copy number
                            fCheckOverlaps);        // checking overlaps 

                G4cout << "Sono qui: iColumn=" << iColumn << "\t kx=" << kx << "\t ky=" << ky << " \t x=" << x/cm << "\t y=" << y/cm 
                       << "\t copyNo=" << copyNo << G4endl;
            }
  }


/*  //This was for a single module
  new G4PVPlacement(rotationMatrix2,                //no rotation
                            G4ThreeVector(0.,0.,0.),  //at (0,0,0)
			       moduleLogical,             //its logical volume
                            "moduleName",             //its name
                            logicWorld,             //its mother  volume
                            false,                  //no boolean operation
                            0,                 //copy number
                            fCheckOverlaps);        // checking overlaps 
*/


  //     
  // ***** SiC semiconductor *****
  //
 
  // Solid
  // SiC Solid

  G4double sicXY;
  G4double sicZ;

  std::ifstream fileSiCLength("sic_geometry.txt");
  fileSiCLength >> sicXY;
  fileSiCLength >> sicZ;  
  fileSiCLength.close();


    
  std::cout << "SiC Active Volume width = " << sicXY << " mm, SiC Active Volume thickness = " << sicZ << "mm" << std::endl;


  G4Box* sicSolid = new G4Box("SiCSolid",sicXY/2.,sicXY/2.,sicZ/2.);

  // SiC Logical Volume definition    
  G4LogicalVolume* sicLogical = new G4LogicalVolume(sicSolid,          //its solid
						    cryst_semi,        //its m
						    "SicLV");          //its name
  G4int nColumns = 2;
  G4int nRows = 10;
  G4double crystalDistanceY = 10.*nm; // Valore messo per sicXY=15.4 mm
  G4double crystalDistanceX = 10.*nm; // Valore messo per sicXY=15.4 mm
  //G4double crystalDistanceY = 0.2*mm; // Valore che avevamo messo quando siCXY=15.2 mm
  //G4double crystalDistanceX = 0.2*mm; // Valore che avevamo messo quando siCXY=15.2 mm


  //G4RotationMatrix* rotationMatrix2 = new G4RotationMatrix();
  //G4double alpha2 = 35.*deg;
  //G4cout << "pi" << pi << "\t sin(alpha)= " << sin(alpha) << "\t cos(alpha)= " << cos(alpha) << G4endl;
  //rotationMatrix2->rotateY(-1*alpha2);

  // create column & row

  for (G4int iColumn=0;iColumn<nRows; iColumn++) 
  {
            for (G4int iRow=0;iRow<nColumns;iRow++)
            {
                G4int jx = iRow - nColumns/2;
                G4int jy = iColumn - nRows/2;
                // G4double y = (jy)*(sicXY+crystalDistanceY); // This value is for 5 rows of telescopes
                G4double y = (jy)*(sicXY+crystalDistanceY) + sicXY/2. + crystalDistanceX/2.;
                G4double x = (jx)*(sicXY+crystalDistanceX) + sicXY/2. + crystalDistanceX/2.;
                char volumename[100];
                sprintf(volumename, "SiC_row%d_col%d", iColumn, iRow);
                G4int copyNo = iRow + iColumn*nColumns;
	         new G4PVPlacement(0,                //no rotation
                            G4ThreeVector(x,y,0.),  //at (0,0,0)
			       sicLogical,             //its logical volume
                            volumename,             //its name
                            moduleLogical,             //its mother  volume
                            false,                  //no boolean operation
                            copyNo,                 //copy number
                            fCheckOverlaps);        // checking overlaps 

                G4cout << "Sono qui: iColumn " << iColumn << "\t jx " << jx << "\t jy " << jy << " \t x " << x/cm << "\t y " << y/cm 
                       << "\t copyNo " << copyNo << G4endl;
            }
  }
  
  // ******************************************************************************************************** //
  //
  // *****   Polyamide frame *****//
  
  G4double poly_ext_XY = 14.70 * mm;
  G4double poly_width = 0.03 * mm;
  G4double poly_int_XY = poly_ext_XY - 2*poly_width;
  G4double polyZ = 0.02 * mm;
  G4Box* poly_ext_Solid = new G4Box("polyExtSolid",poly_ext_XY/2.,poly_ext_XY/2.,polyZ/2.);
  G4Box* poly_int_Solid = new G4Box("polyIntSolid",poly_int_XY/2.,poly_int_XY/2.,polyZ/2.);
  G4VSolid* poly_solid = new G4SubtractionSolid("polySolid",poly_ext_Solid,poly_int_Solid);
  
  // polyamide frame Logical Volume definition    
  G4LogicalVolume* polyLogical = new G4LogicalVolume(poly_solid,          //its solid
						       poly_frame,          //its material
						       "polyLV");           //its name


  for (G4int iColumn=0;iColumn<nRows; iColumn++) 
  {
            for (G4int iRow=0;iRow<nColumns;iRow++)
            {
                G4int jx = iRow - nColumns/2;
                G4int jy = iColumn - nRows/2;
                // G4double y = (jy)*(sicXY+crystalDistanceY); // This value is for 5 rows of telescopes
                G4double y = (jy)*(sicXY+crystalDistanceY) + sicXY/2. + crystalDistanceX/2.;
                G4double x = (jx)*(sicXY+crystalDistanceX) + sicXY/2. + crystalDistanceX/2.;
                G4double z = -sicZ/2. - polyZ/2. - 2.*nm;
                char volumename[100];
                sprintf(volumename, "poly_row%d_col%d", iColumn, iRow);
                G4int copyNo = iRow + iColumn*nColumns;
	         new G4PVPlacement(0,                //no rotation
                            G4ThreeVector(x,y,z),  //at (0,0,0)
			     polyLogical,             //its logical volume
                            volumename,             //its name
                            moduleLogical,             //its mother  volume
                            false,                  //no boolean operation
                            copyNo,                 //copy number
                            fCheckOverlaps);        // checking overlaps 

                G4cout << "Sono qui: iColumn " << iColumn << "\t jx " << jx << "\t jy " << jy << " \t x " << x/cm << "\t y " << y/cm 
                       << "\t copyNo " << copyNo << G4endl;
            }
  }
						       
						       

  //     
  // ***** strato morto SiC *****
  //
  G4double laydeaXY = sicXY;
  std::cout << "SiC deadlayer width = " << laydeaXY << " mm" << std::endl;

  //G4double laydeaZ = 350*um;
  G4double laydeaZ = 10.*um;

 
  // layerdead solid
  G4Box* laydeaSolid = new G4Box("LayDeaSolid",laydeaXY/2.,laydeaXY/2.,laydeaZ/2.);   
  // Layer Logical  Volume
  G4LogicalVolume* laydeaLogical = new G4LogicalVolume(laydeaSolid,
    						       cryst_semi,
  						       "LaydeaLV");

  //Placement: its upper surface must be +0.02 cm from the lower surface of the Sic
  //Lower surface of dead layer: -sicZ
  G4double zLaydea = sicZ/2. + laydeaZ/2. + 10*nm;
  for (G4int iColumn=0;iColumn<nRows; iColumn++)
  {
            for (G4int iRow=0;iRow<nColumns;iRow++)                    
	    {
		G4int jx = iRow - nColumns/2;
		G4int jy = iColumn - nRows/2;
              // G4double y = (jy)*(sicXY+crystalDistanceY); // This value is for 5 rows of telescopes
              G4double y = (jy)*(sicXY+crystalDistanceY) + sicXY/2. + crystalDistanceX/2.;
		//G4double x = (jx)*(sicXY+crystalDistanceX); // Original: telescope not rotated
		G4double x = (jx)*(sicXY+crystalDistanceX) + laydeaXY/2. + crystalDistanceX/2.;   // Telescope in the module
              G4double z = zLaydea;                         // Original: telescope not rotated, good also for the telescope in module
		//G4double x = (jx)*(sicXY+crystalDistanceX) + zLaydea*sin(alpha);
              //G4double z = zLaydea*cos(alpha);
		char volumename[100];                                                            
		sprintf(volumename, "Laydea_row%d_col%d", iColumn, iRow);    
		G4int copyNo = iRow + iColumn*nColumns;
		new G4PVPlacement(0,                           // no rotation
				  G4ThreeVector(x,y,z),  // at (0,0,0)
				  laydeaLogical,               // its logical volume
				  volumename,                  // its name
				  moduleLogical,                  // its mother  volume
				  false,                       // no boolean operation
				  copyNo,                      // copy number
				  fCheckOverlaps);             // checking overlaps 
		
                G4cout << "Sono qui: iColumn " << iColumn << "\t jx " << jx << "\t jy " << jy << " \t x " << x/cm << "\t y " << y/cm 
                       << "\t z " << z/cm << "\t zLaydea " << zLaydea/cm << "\t copyNo " << copyNo << G4endl;
	    } 
  }


// Copper grid -- New implementation

  // ***************************Solid dimensions*************************
  //External grid dimensions
  G4double grid_X = 15.4*mm;
  G4double grid_Y = grid_X;
  G4double grid_Z = 0.5*mm;

  //Grid-hole dimensions
  G4double gridHole_X = 15.*mm; 
  G4double gridHole_Y = gridHole_X;

  // **************************Solid definition*****************************
  //Grid without hole
  G4VSolid* externalSolidGrid = new G4Box("externalSolidGrid",grid_X/2.,grid_Y/2.,grid_Z/2.);
  //Grid-hole solid
  G4VSolid* gridHole = new G4Box("gridhole", gridHole_X/2., gridHole_Y/2.,grid_Z/2.);
    
  // *****************************Boolean solid composition**************************
  //Subtraction of the 10 frameHole
  G4VSolid* gridUnitSolid = new G4SubtractionSolid("grid_cell",externalSolidGrid,
  gridHole);

  // *************************Solid placement**********************

  //G4double X_hole2 =  interframe/2. +frameHole_X/2.;   // 7.7 mm
  //G4double X_hole1 = -X_hole2;  //le file sono simmetriche rispetto al centro
  
  
  G4LogicalVolume* gridLogical = new G4LogicalVolume(gridUnitSolid,
    						       grid_mat,
     						       "gridLV");


  G4double Zgrid = zLaydea + laydeaZ/2. + grid_Z/2. + 10*nm;
    
  
    for (G4int iColumn=0;iColumn<nRows; iColumn++)
  {
            for (G4int iRow=0;iRow<nColumns;iRow++)                    
	    {
		G4int jx = iRow - nColumns/2;
		G4int jy = iColumn - nRows/2;
              // G4double y = (jy)*(sicXY+crystalDistanceY); // This value is for 5 rows of telescopes
              G4double y = (jy)*(sicXY+crystalDistanceY) + sicXY/2. + crystalDistanceX/2.;
		//G4double x = (jx)*(sicXY+crystalDistanceX); // Original: telescope not rotated
		G4double x = (jx)*(sicXY+crystalDistanceX) + laydeaXY/2. + crystalDistanceX/2.;   // Telescope in the module
              G4double z = Zgrid;                         // Original: telescope not rotated, good also for the telescope in module
		//G4double x = (jx)*(sicXY+crystalDistanceX) + zLaydea*sin(alpha);
              //G4double z = zLaydea*cos(alpha);
		char volumename[100];                                                            
		sprintf(volumename, "grid_row%d_col%d", iColumn, iRow);    
		G4int copyNo = iRow + iColumn*nColumns;
		new G4PVPlacement(0,                           // no rotation
				  G4ThreeVector(x,y,z),  // at (0,0,0)
				  gridLogical,               // its logical volume
				  volumename,                  // its name
				  moduleLogical,                  // its mother  volume
				  false,                       // no boolean operation
				  copyNo,                      // copy number
				  fCheckOverlaps);             // checking overlaps 
		
                G4cout << "Sono qui: iColumn " << iColumn << "\t jx " << jx << "\t jy " << jy << " \t x " << x/cm << "\t y " << y/cm 
                       << "\t z " << z/cm << "\t zLaydea " << zLaydea/cm << "\t copyNo " << copyNo << G4endl;
	    } 
  }
  


  // ***RESINA + cesio interno***  

  G4double csiXY = 1.5*cm; 
  G4double csiZ = 0.5*cm;
  G4double resZ = csiZ + 5.*nm;
  G4double resXY = csiXY + 0.4*mm;
  
  G4Box* resSolid = new G4Box("resSolid",resXY/2.,resXY/2.,resZ/2.);

  // resina Logical Volume     

  G4LogicalVolume* resLogical = new G4LogicalVolume(resSolid,        //its solid
						      res_material,    //its material
						      "resLV");        //its name


  crystalDistanceX = 10.*nm;
  crystalDistanceY = 10.*nm;
  G4double zres = Zgrid + 0.5*grid_Z + 0.5*csiZ + 10*nm;
  //G4double zres = frameZ + 0.5*frame_Z + 0.5*csiZ + 10*nm;
  //G4double zCsi = 0.56002*mm + 0.5*csiZ + 10*nm; // 0.56002 è l'equivalente di frameZ+0.5*frame_Z

  for (G4int iColumn=0;iColumn<nRows; iColumn++)
  {
        for (G4int iRow=0;iRow<nColumns;iRow++)
        {
              G4int jx = iRow - nColumns/2;
              G4int jy = iColumn - nRows/2;
              // G4double y = (jy)*(resXY+crystalDistanceY); // This value is for 5 rows of telescopes
              G4double y = (jy)*(resXY+crystalDistanceY) + resXY/2. + crystalDistanceX/2.;
		//G4double x = (jx)*(csiXY+crystalDistanceX);   // Original: telescope not rotated
		G4double x = (jx)*(resXY+crystalDistanceX) + resXY/2. + crystalDistanceX/2.;   // Telescope in the module
              G4double z = zres;                            // Original: telescope not rotated, good also for the telescope in module
		//G4double x = (jx)*(sicXY+crystalDistanceX) + zCsi*sin(alpha);
              //G4double z = zCsi*cos(alpha);
              char volumename[100]; 
              sprintf(volumename, "res_row%d_col%d", iColumn, iRow);
              G4int copyNo = iRow + iColumn*nColumns;
              new G4PVPlacement(0,                       //no rotation
                                G4ThreeVector(x,y,z), //at (0,0,0)
                                resLogical,              //its logical volume
                                volumename,              //its name
                                moduleLogical,              //its mother  volume
                                false,                   //no boolean operation
                                copyNo,                  //copy number
                                fCheckOverlaps);         // checking overlaps 

                G4cout << "Sono qui: iColumn " << iColumn << "\t jx " << jx << "\t jy " << jy << " \t x " << x/cm << "\t y " << y/cm 
                       << "\t z " << z/cm << "\t zres " << zres/cm << "\t copyNo " << copyNo << G4endl;
        }
  }
  
  //
  // ***** CsI *****
  //



  // CsI solid
  G4Box* csiSolid = new G4Box("CsISolid",csiXY/2.,csiXY/2.,csiZ/2.);

  // CsI Logical Volume     

  G4LogicalVolume* csiLogical = new G4LogicalVolume(csiSolid,        //its solid
						    cryst_scin,      //its material
						    "CsILV");        //its name

  // Placement: its upper surface must be -0.055 cm from the lower surface of the NaI
  //  //Lower surface of CsI: layer dead
  //

  //crystalDistanceX = 0.4*mm;
  G4double zCsi = Zgrid + 0.5*grid_Z + 0.5*csiZ + 10*nm;
//  G4double zCsi = frameZ + 0.5*frame_Z + 0.5*csiZ + 10*nm;
  //G4double zCsi = 0.56002*mm + 0.5*csiZ + 10*nm; // 0.56002 è l'equivalente di frameZ+0.5*frame_Z

//  for (G4int iColumn=0;iColumn<nRows; iColumn++)
//  {
//        for (G4int iRow=0;iRow<nColumns;iRow++)
 //       {
//              G4int jx = iRow - nColumns/2;
//              G4int jy = iColumn - nRows/2;
              //G4double y = (jy)*(csiXY/*+crystalDistanceY*/);
		//G4double x = (jx)*(csiXY+crystalDistanceX);   // Original: telescope not rotated
	//	G4double x = (jx)*(csiXY+crystalDistanceX) + csiXY/2. + crystalDistanceX/2.;   // Telescope in the module
             // G4double z = zCsi;                            // Original: telescope not rotated, good also for the telescope in module
		//G4double x = (jx)*(sicXY+crystalDistanceX) + zCsi*sin(alpha);
              //G4double z = zCsi*cos(alpha);
//              char volumename[100]; 
 //             sprintf(volumename, "CsI_row%d_col%d", iColumn, iRow);
              G4int copyNoCsI = 0;
              new G4PVPlacement(0,                       //no rotation
                                G4ThreeVector(0,0,0), //at (0,0,0)
                                csiLogical,              //its logical volume
                                "csi",              //its name
                                resLogical,              //its mother  volume
                                false,                   //no boolean operation
                                copyNoCsI,                  //copy number
                                fCheckOverlaps);         // checking overlaps 

                //G4cout << "Sono qui: iColumn " << iColumn << "\t jx " << jx << "\t jy " << jy << " \t x " << x/cm << "\t y " << y/cm << "\t z " << z/cm << "\t zCsi " << zCsi/cm << "\t copyNoCsI " << copyNoCsI << G4endl;
        //}
  //}

	      /* commented out by RP
  // Now we put two Check-Volumes to evaluate if the particle generated by the Simul_diana hit the telescope
  // G4double checkvolX = 100.*cm;
  // G4double checkvolY = 0.1*cm;
  // G4double checkvolZ = 0.5*mm;
  // G4Box* checkvolSolid = new G4Box("checkvolSolid",checkvolX/2.,
  //                                  checkvolY/2.,checkvolZ/2.);
  // G4LogicalVolume* checkvolLogical = new G4LogicalVolume(checkvolSolid,      //its solid
  // 						               default_mat,        //its material
  // 						               "checkvolLV");      //its name

  // G4double zcheckvol = -133.*mm;
  // G4double xcheckvol = 0.*cm;
  // G4double ycheckvol = 10.*cm;
  // new G4PVPlacement(0, G4ThreeVector(xcheckvol,ycheckvol,zcheckvol), 
  //                   checkvolLogical, "CheckVolume", logicWorld, false, 0, fCheckOverlaps);
  // G4Box* checkvolSolid2 = new G4Box("checkvolSolid2",checkvolX/2.,checkvolY/2.,checkvolZ/2.);
  // G4LogicalVolume* checkvolLogical2 = new G4LogicalVolume(checkvolSolid2,      //its solid
  // 						               default_mat,        //its material
  // 						               "checkvolLV2");      //its name
  // G4double zcheckvol2 = -51.*mm;
  // new G4PVPlacement(0, G4ThreeVector(xcheckvol,ycheckvol,zcheckvol2), checkvolLogical2, "CheckVolume2", logicWorld, false, 0, fCheckOverlaps);
  Commented out by RP */

  /* ****** Now we place an "extra" telescope to evaluate showering from upstream detectors ******  */
  /*
  G4double zCsiAhead = 1*mm + csiZ/2.;
  G4double zLaydeaAhead = zCsiAhead + csiZ/2. + laydeaZ/2.;
  G4double zSicAhead = zLaydeaAhead + laydeaZ/2. + 10*nm + sicZ/2.;
  

  new G4PVPlacement(0,                     
                    G4ThreeVector(0*mm, 0*mm, -zSicAhead),  
                    sicLogical,            
                    "SiC_ahead",             
                    logicWorld,                
                    false,                
                    20,                   
                    fCheckOverlaps);      
  
  new G4PVPlacement(0,                  
                    G4ThreeVector(0*mm, 0*mm, -zLaydeaAhead),     
                    laydeaLogical,          
                    "Laydea_ahead",           
                    logicWorld,                 
                    false,                
                    20,                    
                    fCheckOverlaps);       


  new G4PVPlacement(0,                 
                    G4ThreeVector(0*mm, 0*mm, -zCsiAhead),    
                    csiLogical,          
                    "CsI_ahead",          
                    logicWorld,              
                    false,              
                    20,                  
                    fCheckOverlaps);      

  */
  
  
 G4Region* towersRegion = new G4Region("towers_region");
 towersRegion->AddRootLogicalVolume(moduleLogical);

  // Visualization attributes
  //
  G4Colour red (1., 0., 0.) ;  // red
  G4Colour green (0.,1.,0.);   //green
  G4Colour blue (0.,0.,1.);    //blue
  G4Colour yellow (0.5, 0.5, 0.);
  G4Colour lightgray (0.8, 0.8, 0.8);
  G4Colour lightgray2 (0.6, 0.6, 0.6);
  G4Colour vergadoro (.855,.647,.125);
  G4Colour  gold  (.75, .55, 0.0) ;
  G4Colour  orange  (1., .5, 0.1) ;
  G4Colour white (1., 1., 1.);
  logicWorld->SetVisAttributes(new G4VisAttributes(false));
  moduleLogical->SetVisAttributes(new G4VisAttributes(false));
  //moduleLogical->SetVisAttributes(new G4VisAttributes(false));
  //sicLogical_test->SetVisAttributes(new G4VisAttributes(red));
  //laydeaLogical_test->SetVisAttributes(new G4VisAttributes(green));
  //csiLogical_test->SetVisAttributes(new G4VisAttributes(blue));
  //resLogical->SetVisAttributes(G4VisAttributes(false));
  resLogical->SetVisAttributes(new G4VisAttributes(lightgray2));
  sicLogical->SetVisAttributes(new G4VisAttributes(red));
  polyLogical->SetVisAttributes(new G4VisAttributes(G4Color::Cyan()));  //GB 2025-07-15
  laydeaLogical->SetVisAttributes(new G4VisAttributes(green));
  gridLogical->SetVisAttributes(new G4VisAttributes(vergadoro));
  //frameLogical->SetVisAttributes(G4VisAttributes(false));
  //resLogical->SetVisAttributes(new G4VisAttributes(lightgray));
  csiLogical->SetVisAttributes(new G4VisAttributes(blue));
  //checkvolLogical->SetVisAttributes(G4VisAttributes(false));
  //checkvolLogical2->SetVisAttributes(G4VisAttributes(false));
  /* commented out by RP
  checkvolLogical->SetVisAttributes(new G4VisAttributes(G4Color::Gray()));
  checkvolLogical2->SetVisAttributes(new G4VisAttributes(G4Color::Gray()));
  commented out by RP */
  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl; 

  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TRDetectorConstruction::ConstructSDandField()
{
  //
  // Register some of the volumes as "sensitive" and decide the 
  // type of sensitivity that they have
  //
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  
  // declare SiC semiconduttore as a MultiFunctionalDetector scorer
  //  
  
  // Create a new scorer (G4MultiFunctionalDetector) and set its 
  // "capability" to G4PSEnergyDeposit (will score total energy deposit)
  G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("SiC");
  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  cryst->RegisterPrimitive(primitiv1);
  // Attach the scorer to the logical volume
  SetSensitiveDetector("SicLV",cryst);
  //
  G4MultiFunctionalDetector* cryst2 = new G4MultiFunctionalDetector("DeadLayer");
  G4VPrimitiveScorer* primitiv2 = new G4PSEnergyDeposit("edep");
  cryst2->RegisterPrimitive(primitiv2);
  // Attach the scorer to the logical volume
  SetSensitiveDetector("LaydeaLV",cryst2);
  // 
  G4MultiFunctionalDetector* cryst3 = new G4MultiFunctionalDetector("CsI"); 
  G4VPrimitiveScorer* primitiv3 = new G4PSEnergyDeposit("edep");
  cryst3->RegisterPrimitive(primitiv3); 
  // Attach the scorer to the logical volume
  SetSensitiveDetector("CsILV",cryst3);

  G4SDChargedFilter* chargedFilter = new G4SDChargedFilter("chargedFilter");  // RP

  
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst2);
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst3);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
