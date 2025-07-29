void macro()
{

	//TFile *cut = new TFile("cut.root");

	//TCutG *prova = (TCutG*)cut->Get("prova");

	//TCut prova_cut = "prova"; 

        //cut->Close();
        
        


  int atom_number, mass_number;
  int sic_mult = -1;
  double aESiC[50] = {0.};
  int sicID[50] = {0};
  int csi_mult = -1;
  double aECsI[50] = {0.};
  int csiID[50] = {0};
  int deadsic_mult = -1;
  double aEdeadSiC[50] = {0.};
  int deadsicID[50] = {0};  
  
  
  double sumE = 0.;
  //TH2D* histo2 = new TH2D("histo","histo",100,0.,100.,1000,0.,1000.);
  TH2D* histo2 = new TH2D("histo2","histo2",300,0.,300.,200,0.,200.);
  TH2D* histo3 = new TH2D("histo3","histo3",1000,0.,1000.,150,0.,240.);
  TH2D* histo4 = new TH2D("histo4","histo4",1000,0.,1000.,150,0.,150.);
  TH2D* histo5 = new TH2D("histo5","histo5",1000,0.,1000.,150,0.,150.);
  TH2D* histo6 = new TH2D("histo6","histo6",1000,0.,1000.,150,0.,150.);
  TH1D* histo1 = new TH1D("histo1","histo1",40,240.,260.);
  
  TFile* file = new TFile("output_postprocessing.root");
  TTree* tree = (TTree*)file->Get("tree");
  tree->SetBranchAddress("atom_number",&atom_number);
  tree->SetBranchAddress("mass_number",&mass_number);
  tree->SetBranchAddress("sic_mult",&sic_mult);
  tree->SetBranchAddress("aESiC",&aESiC);
  tree->SetBranchAddress("sicID",&sicID);
  tree->SetBranchAddress("csi_mult",&csi_mult);
  tree->SetBranchAddress("aECsI",&aECsI);
  tree->SetBranchAddress("csiID",&csiID);
  tree->SetBranchAddress("deadsic_mult",&deadsic_mult);
  tree->SetBranchAddress("aEdeadSiC",&aEdeadSiC);
  tree->SetBranchAddress("deadsicID",&deadsicID);  

  for(int a=0; a<tree->GetEntries(); ++a){
    tree->GetEntry(a);
    if (sic_mult == csi_mult == 1)
        histo1->Fill(aESiC[0]+aECsI[0]+aEdeadSiC[0]);
    
    for(int b=0; b<sic_mult; ++b)
      if(aESiC[b]!=0)
	for(int c=0; c<csi_mult; ++c) {
	  if(sicID[b]==csiID[c] && aESiC[b] > 90. )
	    histo2->Fill(aECsI[c],aESiC[b]); // versione originale Rino: histo2->Fill(aESiC[b],aECsI[c]);
	  
//	  if(sicID[b]==25009&&csiID[b]==25009)
//  	    histo3->Fill(aECsI[c],aESiC[b]);
  	    
//	  if(sicID[b]==25010&&csiID[b]==25010)
//  	    histo4->Fill(aECsI[c],aESiC[b]);
  	    
  	    
//	  if((sicID[b]==1010&&csiID[b]==1010) && prova->IsInside(aECsI[c],aESiC[b]) )
//	  if((sicID[b]==25010&&csiID[b]==25010) && atom_number==8 && mass_number==20 )
// 	    histo5->Fill(aECsI[c],aESiC[b]);  
  	    
//    	  if((sicID[b]==1010&&csiID[b]==1010) && prova->IsInside(aECsI[c],aESiC[b]) )
//	  if((sicID[b]==25010&&csiID[b]==25010) && atom_number==8 && mass_number==19 )
//  	    histo6->Fill(aECsI[c],aESiC[b]);  	    
	}
  }
  histo2->Draw();
  //histo1->Draw();
//  TCanvas* c2 = new TCanvas("c2","c2");
//  histo3->Draw();
  
  
//  TCanvas* c3 = new TCanvas("c3","c3");
//  histo4->Draw();
/*  
  TCanvas* c4 = new TCanvas("c4","c4");
  histo5->Draw(); 
  histo6->SetMarkerColor(kRed); 
  histo6->Draw("same");  
*/
}
