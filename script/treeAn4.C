
// This macro reads the files produced by the Geant4 simulation, called TR_1.root, TR_2.root, etc., and creates a root file (simul.root). This file contains a tree (simul_tree) with the relevant information: ejectile identity, inital energy, angles and position, energy deposit in SiC and in CsI, tof.


#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <map>
#include <iterator>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF2.h"
#include "TF3.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLegend.h"
#include <TString.h>
#include "TLine.h"
#include "TApplication.h"  // GB 2021-05-10



// This function is the linear response function which parameterizes
// the dead and the partially active region of the SiC detector
Double_t lin_f(Double_t *x, Double_t* par) {

   Double_t SiCX;

   // GB 2021-05-10 - begin
   Double_t SiCZ;
   Double_t part_alive_layerX, dead_layerX;
   Double_t part_alive_layerZ, dead_layerZ;

   // We read the SiC width and thickness together with 
   // its partially active region and dead region extensions 
   // from the "sic_geometry.txt", where these quantities are expressed in mm
   TString fileSiCGeometryName = "sic_geometry.txt";
   std::ifstream fileSiCGeometry(fileSiCGeometryName);
   if (!fileSiCGeometry) {
      std::cerr << "File " << fileSiCGeometryName << " does not exist.\n" << std::endl;
      gApplication->Terminate();
      //exit(1);
   }
   //else 
      //std::cout << "File " << fileSiCGeometryName << " correctly opened!\n" << std::endl;

   fileSiCGeometry >> SiCX >> SiCZ >> part_alive_layerX >> dead_layerX >> part_alive_layerZ >> dead_layerZ;
   fileSiCGeometry.close();
   
   //std::cout << "SiC length = " << SiCX << " mm\t" 
   //          << "SiC thickness = " << SiCZ << " mm\t"
   //          << "SiC part_alive_layerX = " << part_alive_layerX << " mm\t"
   //          << "SiC dead_layerX = " << dead_layerX << " mm\t"
   //          << "SiC part_alive_layerZ = " << part_alive_layerZ << " mm\t"
   //          << "SiC dead_layerZ = " << dead_layerZ << " mm"
   //          << std::endl;



   // Convertiamo la lunghezza in um e ne prendiamo la metà (per la simmetria del sistema è più facile così)
   SiCX *= 1000./2.;
   SiCZ *= 1000./2.;
   part_alive_layerX *= 1000.;
   dead_layerX *= 1000.;
   part_alive_layerZ *= 1000.;
   dead_layerZ *= 1000.;

   //std::cout << "Half SiC length = " << SiCX << " um\t" 
   //          << "Half SiC thickness = " << SiCZ << " um" << std::endl;


   //Double_t part_alive_layer = 20.;  // Corrisponde a 20um -> Suggerito da Tudisco
   //Double_t part_alive_layerX = 50.;   // Ricordati di mettere 50!! Corrisponde a 50um -> lo metto per evidenziare effetti di bordo
   //Double_t part_alive_layerX = 100.;	//**prova strato parzialmente attivo da 100 
   //Double_t dead_layerX = 200.;        // Ricordati di mettere 200!! Corrisponde a 200um -> Confermato da Tudisco
   //Double_t dead_layerX = 100.; 		//**prova strato morto da 100
   //Double_t part_alive_layerZ = 1.;    // Corrisponde a 1um -> lo metto per evidenziare effetti di bordo
   //Double_t dead_layerZ = 1.;          // Corrisponde a 1um -> Confermato da Tudisco

   // GB 2021-05-10 - end


   Double_t b = SiCX - dead_layerX;
   Double_t a = b - part_alive_layerX;
   Double_t d = SiCZ - dead_layerZ;
   Double_t c = d - part_alive_layerZ;

   Double_t result;


   x[0] = TMath::Abs(x[0]);
   x[1] = TMath::Abs(x[1]);
   x[2] = TMath::Abs(x[2]);
   
   if ( x[0] <= a && x[1] <= a && x[2] <= c ) 
      result = 1;
   
   if ( x[0] > a && x[0] < b && x[1] <= x[0] && x[2] <= ( (d-c)*x[0]/(b-a) - (a*d-b*c)/(b-a) ) )
      result = x[0]/( a - b ) - b/( a - b );

   if ( x[0] < x[1] && x[1] > a && x[1] < b && x[2] < ( (d-c)*x[1]/(b-a) - (a*d-b*c)/(b-a) ) )
      result = x[1]/( a - b ) - b/( a - b );

   if ( x[0] < ( (b-a)*x[2]/(d-c) + (a*d-b*c)/(d-c) ) && x[1] < ( (b-a)*x[2]/(d-c) + (a*d-b*c)/(d-c) ) && x[2] > c && x[2] < d )
      result = x[2]/( c - d ) - d/( c - d );

   /*
   if ( x[0] < a && x[1] > a && x[1] < b && x[2] < ( (c-d)*x[1]/(a-b) + (a*d-b*c)/(a-b) ) )
      return x[2]/( c - d ) - d/( c - d );
   */

   if ( x[0] > b || x[1] > b || x[2] > d ) 
      result = 0;
  
   return result;

}


//n_phot = norm*(d2/atom_number + d3 + d4*atom_number)*( ECsI/mass_number + (d1*atom_number)*(TMath::Exp(-ECsI/(mass_number*d1*atom_number)) -1) ); // Mastinu, NIMA 338 (1994)
//n_phot_meas = ((kk1/(kk2)*kk3))*((cc1+(cc2/atom_number))*ECsI + (bb1 + (bb2*atom_number))*(TMath::Exp(-(aa1+aa2/atom_number)*ECsI)-1)); // BIRKS + BECCHETTI Sperimentale (Finocchiaro)

Double_t f_n_phot_Mastinu(Int_t atom_number, Int_t mass_number, Double_t ECsI){  // Mastinu, NIMA 338 (1994)
	
	 Double_t d1= 12.9;
 	 Double_t d2= 8.17;
 	 Double_t d3=0.89;
 	 Double_t d4= 0.0025;
  	 Double_t norm= 46255.5;
  	
  	return norm*(d2/atom_number + d3 + d4*atom_number)*( ECsI/mass_number + (d1*atom_number)*(TMath::Exp(-ECsI/(mass_number*d1*atom_number)) -1) );
  	
	}
	
Double_t f_n_phot_Finocchiaro(Int_t atom_number, Double_t ECsI){  // 
	
	  Double_t aa1 = 0.01; 
	  Double_t aa2 = 0.13; 
	  Double_t bb1 = 32.;
	  Double_t bb2 = 2.6;
	  Double_t cc1 = 0.58;
	  Double_t cc2 =3.87;
	  Double_t kk1 = 32000.;
	  Double_t kk2 = 18.;
	  Double_t kk3 = 7.25; 
	  	
  	return ((kk1/(kk2)*kk3))*((cc1+(cc2/atom_number))*ECsI + (bb1 + (bb2*atom_number))*(TMath::Exp(-(aa1+aa2/atom_number)*ECsI)-1));
  	
	}
	
Double_t f_n_phot_Horn(Int_t atom_number, Int_t mass_number, Double_t ECsI){  // 
	
	  Double_t a1 = 8.145; 
	  Double_t a2 = 0.326;
	  Double_t k1 = 32000./8.; //Paolo Finocchiaro
 
	  	
  	return (k1)*(a1)*(ECsI - a2*mass_number*TMath::Power(atom_number,2)*TMath::Log(TMath::Abs((ECsI+a2*mass_number*TMath::Power(atom_number,2))/(a2*mass_number*TMath::Power(atom_number,2)))));
  	
	}


void treeAn4() {
  
   char fileName[50];

   const char* treeInName = "tree";

   const char* fileOutName = "prova2.root";
   const char* treeOutName = "simul_tree";


   string PrintOrNot;

   //gStyle->SetOptFit(0111);
   //gStyle->SetOptStat(0000);

   std::cout << "Do you want to see the energy depositions in each event? [Y]es or [N]o " << std::endl;
   cin >> PrintOrNot;


   Int_t evtID=0, evtIDpost, dump=-1, n_evt=0, logVol, copyNb, max_copyNb=36010;
   Int_t copy_number=-1;
   Int_t nEntries=0;
   Int_t atom_number=0, mass_number=0;
   Double_t eDep, x, y, z;
   Double_t xSiC=-100., ySiC=-100.;	// These are the impact point on the SiC detector
   Double_t ESiC, ECsI, sigma_SiC=0., sigma_CsI=0.;
   Double_t cosTheta, cosThetaPre=0.;
   Double_t KinEnergy, KinEnergyPre=0.;
   Double_t xf, xfPre=0.;
   Double_t yf, yfPre=0.;
   Double_t tf, tfPre=0.;
   Double_t phf, phfPre=0.;
   Double_t tpol, tpolPre=0.;
   Double_t phpol, phpolPre=0.;
   Double_t tof, tofPre=0.;
   Double_t thetafoc_meas=0.;		// The variable thetafoc_meas mimics the measurement process of the angle thetafoc.
   Double_t threshold=0.1;
   Double_t n_phot_Mastinu=0.;
   Double_t n_phot_Finocchiaro=0.;
   // Double_t n_phot_Horn=0.;                  // RP
   Double_t dEtrk = -10.;                       // RP
   Int_t sic_multiplicity=0;
   Int_t csi_multiplicity=0;
   Int_t sic_mult = 0;                          // RP
   Double_t aESiC[100] = {0.};                   // RP 
   Int_t    sicID[100] = {0};                    // RP
   Int_t    xSic[100] = {0};                     // RP
   Int_t    ySic[100] = {0};                     // RP
   Double_t tSiC[100]= {-1.};                    // RP
   Int_t csi_mult = 0;                          // RP
   Double_t aECsI[100] = {0.};                   // RP 
   Int_t    csiID[100] = {0};                    // RP
   Double_t n_phot_Horn[100] = {0.};             // RP
   Double_t tCsI[100]= {-1.};                    // RP
   //k=1,2,3,4,5,6,7,8,9
   Int_t Zk[9] = {8,8,9,9,9,10,10,8,8};         // RP
   Int_t Ak[9] = {20,19,19,20,21,20,21,16,19};  // RP
   
   Double_t thetafoc_resol=0.00425531;	// From simul.out we get the exact value of thetafoc, called tf.
						// The variable thetafoc_resol represents the thetafoc resolution in sigma,
						// which in FWHM is 10 mrad.
   Double_t tSiCfake=-100.;   // This variable represents the arrival time in the SiC 

   //Double_t k1 = 0.35, k2 = 0.065;  // With these values we get an energy resolution for CsI of 2%
   //Double_t k3 = 0.001, k4 = 0.001; // With these values we get an energy resolution for Sic of 6.05%
   Double_t k1 = 0.35, k2 = 0.05;  // With these values we get an energy resolution for CsI of 2%
   Double_t k3 = 0.0025, k4 = 0.00014; // With these values we get an energy resolution for Sic of 6.05%

  //sigma_CsI = TMath::Sqrt( k1 + k2*ECsI );   // GB 2021-05-10
  //sigma_SiC = TMath::Sqrt( k3 + k4*ESiC );   // GB 2021-05-10


   TTree* t;

   std::map <std::pair<Int_t,Int_t>, Double_t> EnergyMap;
   map <pair<Int_t,Int_t>, Double_t> MaxEnergyMap;          // RP
   map <pair<Int_t,Int_t>, Double_t> xSiCMap;       // RP
   map <pair<Int_t,Int_t>, Double_t> ySiCMap;       // RP
   map <pair<Int_t,Int_t>, Double_t> tSiCMap;       // RP
   map <pair<Int_t,Int_t>, Double_t> tCsIMap;       // RP
   

   // GB 2021-05-10 - begin
   Double_t SiCX, SiCZ;


   TString fileSiCGeometryName = "sic_geometry.txt";
   std::ifstream fileSiCGeometry(fileSiCGeometryName);
   if (!fileSiCGeometry) {
      std::cerr << "File " << fileSiCGeometryName << " does not exist.\n" << std::endl;
      gApplication->Terminate();
      //exit(1);
   }
   else 
      std::cout << "File " << fileSiCGeometryName << " correctly opened!\n" << std::endl;

   fileSiCGeometry >> SiCX >> SiCZ;
   fileSiCGeometry.close();

   // Convertiamo la lunghezza in um
   SiCX *= 1000.;
   SiCZ *= 1000.;

   TF3 *resp_funct = new TF3("resp_funct", lin_f, 0., SiCX, 0., SiCX, 0., SiCZ);

   // GB 2021-05-10 - end

   // TF3 *resp_funct = new TF3("resp_funct", lin_f,                                                    // GB 2021-05-10
   //                           0., 15200., 0., 15200., 0., 50.); //con le nuove misure 1.48+100+100    // GB 2021-05-10

                              //0., 15000., 0., 15000., 0., 50.);
                              //4600., 5200., 4600., 5200., 0., 50.);
   //resp_funct->Draw("SURF");



   TFile* fileOut = new TFile(fileOutName,"recreate");
   TTree* simul_tree = new TTree("simul_tree", "Tree from simulation");
   simul_tree->Branch("copy_number", &copy_number, "copy_number/I" );



   const char* fileOutName2 = "test2.root";
   const char* treeOutName2 = "tree";
   TFile* fileOut2 = new TFile(fileOutName2,"recreate");
   TTree* simul_tree2 = new TTree(treeOutName2, "Tree from simulation for event building");
   simul_tree2->Branch("evtID", &evtID, "evtID/I" );
   simul_tree2->Branch("atom_number", &atom_number, "atomic_number/I" );
   simul_tree2->Branch("mass_number", &mass_number, "mass_number/I" );
   simul_tree2->Branch("KinE", &KinEnergy, "KinEnergy/D" );
   simul_tree2->Branch("xf", &xf, "xf/D" );
   simul_tree2->Branch("yf", &yf, "yf/D" );
   simul_tree2->Branch("tf", &tf, "tf/D" );
   simul_tree2->Branch("phf", &phf, "phf/D" );
   simul_tree2->Branch("dEtrk", &dEtrk, "dEtrk/D" );           // RP
   simul_tree2->Branch("sic_mult", &sic_mult, "sic_mult/I" );  // RP
   simul_tree2->Branch("aESiC", aESiC, "aESiC[sic_mult]/D" );  // RP
   simul_tree2->Branch("sicID",sicID,"sicID[sic_mult]/I");     // RP
   simul_tree2->Branch("xSic",xSic,"xSic[sic_mult]/I");        // RP
   simul_tree2->Branch("ySic",ySic,"ySic[sic_mult]/I");        // RP
   simul_tree2->Branch("tSiC", tSiC, "tSiC[sic_mult]/D");      // RP
   simul_tree2->Branch("csi_mult", &csi_mult, "csi_mult/I" );  // RP
   simul_tree2->Branch("aECsI", aECsI, "aECsI[csi_mult]/D" );  // RP
   simul_tree2->Branch("csiID",csiID,"csiID[csi_mult]/I");     // RP
   simul_tree2->Branch("tCsI", tCsI, "tCsI[csi_mult]/D");      // RP
   simul_tree2->Branch("n_phot_Horn",n_phot_Horn,"n_phot_Horn[csi_mult]/D");//RP

   // la coppia individua i copyNb del SiC e del CsI, 
   // non c'è bisogno del logicVolume


   TRandom3 *rnd = new TRandom3();
   

   for (Int_t k=1; k<2; k++) { // This is a cycle on the TR.root files

       sprintf(fileName, "TR_%d.root", k);

       TFile *f = new TFile(fileName);
       if (!f) {
          std::cerr << "File " << fileName << " does not exist." << std::endl;
          exit(1);
       }
       else std::cout << "File \"" << fileName << "\" correctly opened!\n" << std::endl;

       t =(TTree*)f->Get(treeInName);

       t->SetBranchAddress("evtID", &evtIDpost);
       t->SetBranchAddress("logicalVolume", &logVol);
       t->SetBranchAddress("copyNumber", &copyNb);
       t->SetBranchAddress("eDep", &eDep);
       t->SetBranchAddress("x", &x);
       t->SetBranchAddress("y", &y);
       t->SetBranchAddress("z", &z);
       t->SetBranchAddress("cosThetaIn", &cosThetaPre);
       t->SetBranchAddress("KinE", &KinEnergyPre);
       t->SetBranchAddress("xf", &xfPre);
       t->SetBranchAddress("yf", &yfPre);
       t->SetBranchAddress("tf", &tfPre);
       t->SetBranchAddress("phf", &phfPre);
       t->SetBranchAddress("tpol", &tpolPre);
       t->SetBranchAddress("phpol", &phpolPre);
       t->SetBranchAddress("tof", &tofPre);
       t->SetBranchAddress("tSiC", &tSiCfake);


       nEntries = t->GetEntries();
       std::cout << "The number of entries is " << nEntries << std::endl;

       dump=-1; n_evt=0;

       //for (Int_t i=0; i<150; i++) {
       for (Int_t i=0; i<nEntries; i++) {

           ESiC = ECsI = -100.;
           sic_multiplicity=0;
           csi_multiplicity=0;

           t->GetEntry(i);
           
           //if(logVol==99)
             //cout << "copyNb" << copyNb  << endl;
           
	    //if(evtIDpost%10000==0)
	      //cout <<"Analyzing entry number: " << evtIDpost<< endl;
           // if (tSiCfake != -10.)
           //    tSiC = tSiCfake;

	   //cout <<evtIDpost<<"  "<<logVol<<"  "<<copyNb<<endl;

           if (( dump != -1 && evtIDpost != dump ) || i == nEntries - 1 ) {
              if ( PrintOrNot == "Y" || PrintOrNot == "y" || PrintOrNot == "Yes" ) {
             
                 std::cout << "In the event " << evtID << " this is what happened: " << std::endl;

                 for ( auto it = EnergyMap.begin(); it != EnergyMap.end(); ++it ) {
                     std::cout << "****** LogVol = " << (it->first).first << "; copyNb = " << (it->first).second 
                               << "\t eDep = " << it->second << " MeV" << std::endl;
                 }
              }
/*
              for (Int_t k=0; k<max_copyNb; k++) {

                  if ( EnergyMap.count(std::make_pair(1,k)) && EnergyMap.at(std::make_pair(1,k))/1000.>threshold )
                     sic_multiplicity++;

                  if ( EnergyMap.count(std::make_pair(3,k)) && EnergyMap.at(std::make_pair(3,k))/1000.>threshold )
                     csi_multiplicity++;

              }
*/
// RP - begin
	      for (auto it=EnergyMap.cbegin(); it!=EnergyMap.cend(); ++it){
		// if((it->first).first==1)
		//   cout <<(it->first).first<< "  "<<(it->first).second<<"  "
		//        <<it->second<<endl;
		double edep = it->second;
		switch((it->first).first)
		  {
		  case 1:
		    // Facciamo il sampling gaussiano dei soli eventi in cui
		    // il SiC è stato colpito
                  sigma_SiC =  TMath::Sqrt( k3 + k4*edep);  // GB 2021-05-10
		    if(edep > 0.)   
		      edep = rnd->Gaus(edep, sigma_SiC);
		    // Se il SiC è stato colpito, ma con il sampling gaussiano
		    // otteniamo un valore negativo, diciamo che l'energia
		    // depositata è stata di 0.000001 MeV
		    aESiC[sic_mult] = (edep < 0.)? 1e-06 : edep ;
		    sicID[sic_mult] = (it->first).second;
		    xSic[sic_mult] = xSiCMap[make_pair((it->first).first,
						       (it->first).second)];
		    ySic[sic_mult] = ySiCMap[make_pair((it->first).first,
						       (it->first).second)];
		    tSiC[sic_mult] = tSiCMap[make_pair((it->first).first,
						       (it->first).second)];
		    ++sic_mult;

		    break;
		  case 3:		    
		    // Facciamo il sampling gaussiano dei soli eventi in cui
		    // il CsI è stato colpito
		    if(edep > 0.)
                    sigma_CsI =  TMath::Sqrt( k1 + k2*edep);  // GB 2021-05-10   
		      edep = rnd->Gaus(edep, sigma_CsI);
		    // Se il CsI è stato colpito, ma con il sampling gaussiano
		    // otteniamo un valore negativo, diciamo che l'energia
		    // depositata è stata di 0.000001 MeV
		    aECsI[csi_mult] = (edep < 0.)? 1e-06 : edep ;
		    csiID[csi_mult] = (it->first).second;
		    tCsI[csi_mult] = tSiCfake;
		    tCsI[sic_mult] = tCsIMap[make_pair((it->first).first,
						       (it->first).second)];
		    n_phot_Horn[csi_mult] = f_n_phot_Horn(Zk[k-1],
							  Ak[k-1],
							  aECsI[csi_mult]);
		    ++csi_mult;
		    break;
		  default:
		    break;
		  }		
	      }

	      // thetafoc_meas = rnd->Gaus(tf, thetafoc_resol);
	      // atom_number = Zk[k-1];
	      // mass_number = Ak[k-1];
	      
	      // n_phot_Mastinu = f_n_phot_Mastinu(Zk[k-1], Ak[k-1], ECsI);
	      // n_phot_Finocchiaro = f_n_phot_Finocchiaro(Zk[k-1], ECsI);
	      // n_phot_Horn = f_n_phot_Horn(Zk[k-1], Ak[k-1], ECsI);
		atom_number=Zk[k-1];
		mass_number=Ak[k-1];
	      simul_tree->Fill();
	      simul_tree2->Fill();	      
	      
             EnergyMap.clear();
	      MaxEnergyMap.clear();
	      xSiCMap.clear();
	      ySiCMap.clear();
	      tSiCMap.clear();
	      tCsIMap.clear();
	      
	      for(int a=0; a<sic_mult; ++a){
		aESiC[a] = 0.;
		sicID[a] = 0;
		tSiC[a] = 0.;
	      }
	      sic_mult = 0;
	      
	      for(int a=0; a<csi_mult; ++a){
		aECsI[a] = 0.;
		csiID[a] = 0;
		n_phot_Horn[a] = 0.;
		tCsI[a] = 0.;
	      }
	      csi_mult = 0;
	      
	      // RP - end
	      
           } //chiudo l'if(dump!=-1 && evtIDpost!=dump)       


	   // RP
	   if(logVol==1){
	     if(MaxEnergyMap.count(make_pair(logVol,copyNb))){
	       if(eDep>MaxEnergyMap[make_pair(logVol,copyNb)]){
		 MaxEnergyMap[make_pair(logVol,copyNb)] = eDep;
		 xSiCMap[make_pair(logVol,copyNb)]=x;
		 ySiCMap[make_pair(logVol,copyNb)]=y;
		 tSiCMap[make_pair(logVol,copyNb)]=tSiCfake;
	       }
	     }
	     else{
	       xSiCMap[make_pair(logVol,copyNb)]=x;
	       ySiCMap[make_pair(logVol,copyNb)]=y;
	       tSiCMap[make_pair(logVol,copyNb)]=tSiCfake;
	     }
	     eDep = eDep*(resp_funct->Eval(x,y,z));
	   }
	   else if(logVol==3){
	     if(MaxEnergyMap.count(make_pair(logVol,copyNb))){
	       if(eDep>MaxEnergyMap[make_pair(logVol,copyNb)]){
		 MaxEnergyMap[make_pair(logVol,copyNb)] = eDep;
		 tCsIMap[make_pair(logVol,copyNb)]=tSiCfake;
		 //////// tCsIMap[make_pair(logVol,copyNb)]=tCsIfake; CAMBIARE
	       }
	     }
	     else
	       tCsIMap[make_pair(logVol,copyNb)]=tSiCfake;
	   }
	    // converto L'energia in MeV 
	   eDep *= 0.001;
           if ( EnergyMap.count( std::make_pair(logVol,copyNb) ) )
	     EnergyMap[std::make_pair(logVol,copyNb)] += eDep;
	   else{
	     EnergyMap.insert( std::pair <std::pair<Int_t,Int_t>, Double_t> 
			       ( std::make_pair(logVol,copyNb), eDep ) );
	     MaxEnergyMap.insert(std::pair <std::pair<Int_t,Int_t>, Double_t>
				 (make_pair(logVol,copyNb), eDep));
	   }

	   /*	   // commented by RP
           if ( EnergyMap.count( std::make_pair(logVol,copyNb) ) ) {
	     if ( logVol == 1 ) {
	       EnergyMap[std::make_pair(logVol,copyNb)] += eDep*(resp_funct->Eval(x,y,z));
	     }
	     else
	       EnergyMap[std::make_pair(logVol,copyNb)] += eDep;
	   }
	   
	   else {
	     if ( logVol == 1 ) {
	       EnergyMap.insert( std::pair <std::pair<Int_t,Int_t>, Double_t> 
				 ( std::make_pair(logVol,copyNb), eDep*(resp_funct->Eval(x,y,z)) ) );
	       xSiC = x;
	       ySiC = y;
	     }
	     else
	       EnergyMap.insert( std::pair <std::pair<Int_t,Int_t>, Double_t> 
				 ( std::make_pair(logVol,copyNb), eDep ) );
	   }
	   */

	   
       
           if ( evtIDpost != dump ) {
              dump = evtIDpost;
              evtID = evtIDpost;
              cosTheta = cosThetaPre;
              KinEnergy = KinEnergyPre;
              xf = xfPre;
              yf = yfPre;
              tf = tfPre;
	      // da scommentare nel momento in cui avremo l'energia persa nel
	      // tracker direttamente da simulazione MonteCarlo
	      //dEtrk = dEtrkPre;   // RP
              phf = phfPre;
              tpol = tpolPre;
              phpol = phpolPre;
              tof = tofPre;
              n_evt++;
           }
       
       } // chiudo il ciclo for sulle entries

	//f->Close();
	EnergyMap.clear();
	MaxEnergyMap.clear();
   } 

   fileOut->Write();
   fileOut->Close();
  
   std::cout << "\n\nFile " << fileOutName << " has been created!" << std::endl;

   fileOut2->Write();
   fileOut2->Close();
  
   std::cout << "\n\nFile " << fileOutName2 << " has been created!" << std::endl;
   
}
