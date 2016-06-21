#include "TText.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveStats.h"
/*
#ifndef ROOT_Gtypes
#include "Gtypes.h"
#endif
#ifndef ROOT_TString
#include "TString.h"
#endif
*/

//spin : 0 spin0, 2 spin2.
//ch : 0 4mu, 1 4e, 2 2e2mu.
short local_ZZCandType; //1 for merged jet (J), 2 for two resolved jets (jj)
int exclude; //0 for exclude nothing; 1 exclude events with one jet type; 2 exclude events with both jet types

TGraphErrors* makegr(int spin=0, int ch=0, Color_t color=2, int marker=20, int line=1){
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetOptStat(0000);
  gStyle->SetTitleFontSize(0.05);

//get the right directory automatically
  TFile* f;
  char dest[PATH_MAX];
  sprintf(dest, "%s", gSystem->pwd());
  char * pchar = strstr(dest, "prod");
  strcpy(pchar, "prod/\0");
  TString inputDir = dest;
 
  int inputfiles_ggH[]={1000, 2000, 200, 250, 350, 400, 450, 500, 600, 700};
  
  int Nfiles_ggH=sizeof(inputfiles_ggH)/sizeof(*inputfiles_ggH);
  char inputfile[1000];
  vector<TString> files_ggH;

   for (int i=0; i<Nfiles_ggH; i++) {
	 sprintf(inputfile,"2016_2l2q/mytree/PT13TeV/ggH_zz2l2q_M%d/ZZ2l2qAnalysis.root",inputfiles_ggH[i]);
	 files_ggH.push_back(inputfile);
   }

  TChain *candTree = new TChain("ZZTree/candTree");

  for (vector<TString>::const_iterator file = files_ggH.begin(); file!=files_ggH.end(); ++file) 
  candTree->Add(inputDir+(*file));

// draw raw histogram with 1GeV binning, 0GeV to 3050GeV.
// const Int_t m=3050;
 const Int_t m=4500;
 int zzflav=28561;
 TH1F *hgen = new TH1F("hgen","hgen",m,0,m); //histogram
 TH1F *hreco = new TH1F("hreco","hreco",m,0,m);

 switch (ch){
 case 0: zzflav=28561;  break; //4mu
 case 1: zzflav=14641;  break; //4e
 case 2: zzflav=20449;  break; //2e2mu
 default:  cout<<"channel unknown\t"<<ch<<endl;
 }

	vector<short> *ZZsel=0,*ZZCandType=0;
	float GenHMass=0;
	vector<float> *ZZMass=0;
	short genFinalState;
	float PUWeight,genHEPMCweight;

	candTree->SetBranchAddress("ZZCandType",&ZZCandType);
	candTree->SetBranchAddress("ZZsel",&ZZsel);
	candTree->SetBranchAddress("ZZMass",&ZZMass);
	candTree->SetBranchAddress("GenHMass",&GenHMass);
	candTree->SetBranchAddress("genFinalState",&genFinalState);
	candTree->SetBranchAddress("PUWeight",&PUWeight);
	candTree->SetBranchAddress("genHEPMCweight",&genHEPMCweight);

  for (int i=0; i<candTree->GetEntries(); i++) {
      candTree->GetEntry(i);
 
	  switch(exclude) {
	case 0:
		if((ZZCandType->size() == 1 || ZZCandType->size() == 2) && (genFinalState==ch) && ((ZZCandType->size()==1)? ((ZZCandType->at(0) == local_ZZCandType)) : (((ZZCandType->at(0) == local_ZZCandType) || (ZZCandType->at(1) == local_ZZCandType))) )) hgen->Fill(GenHMass,(genHEPMCweight*PUWeight));

		if((ZZCandType->size() == 1 || ZZCandType->size() == 2) && (genFinalState==ch) && ((ZZCandType->size()==1)? ((ZZCandType->at(0) == local_ZZCandType)) : ((((ZZCandType->at(0) == local_ZZCandType) && (ZZsel->at(0)>=100)) || ((ZZCandType->at(1) == local_ZZCandType) && (ZZsel->at(1)>=100)))) )) hreco->Fill(GenHMass,(genHEPMCweight*PUWeight));
		break;
	case 1:
		if((ZZCandType->size() == 2) && ((ZZCandType->at(0) == local_ZZCandType) || (ZZCandType->at(1) == local_ZZCandType)) && (genFinalState==ch)) hgen->Fill(GenHMass,(genHEPMCweight*PUWeight));

		if((ZZCandType->size() == 2) && (((ZZCandType->at(0) == local_ZZCandType) && (ZZsel->at(0)>=100)) || ((ZZCandType->at(1) == local_ZZCandType) && (ZZsel->at(1)>=100))) && (genFinalState==ch)) hreco->Fill(GenHMass,(genHEPMCweight*PUWeight));
		break;
	case 2:
		if((ZZCandType->size() == 1) && (ZZCandType->at(0) == local_ZZCandType) && (genFinalState==ch)) hgen->Fill(GenHMass,(genHEPMCweight*PUWeight));

		if((ZZCandType->size() == 1) && (ZZCandType->at(0) == local_ZZCandType) && (genFinalState==ch) && (ZZsel->at(0)>=100))  hreco->Fill(GenHMass,(genHEPMCweight*PUWeight));
		break;
		}
 }



// bin contents of raw histograms
 double M_raw[m]={0};
 double gen_raw[m]={0};
 double reco_raw[m]={0};

// bin contents of merged histograms
 const Int_t n=14;
 double M[n]={1000,1200, 1400,1600,1800, 2000, 200, 250, 350, 400, 450, 500, 600, 700};

 double massE[n]={0};
 double gen[n]={0};
 double genE[n]={0};
 double reco[n]={0};
 double recoE[n]={0};
 double eff[n]={0};
 double effE[n]={0};

 double mass,width;

 for (int bin=1;bin<=m;bin++){
   M_raw[bin-1] = hreco->GetXaxis()->GetBinCenter(bin);
   gen_raw[bin-1] = hgen->GetBinContent(bin);
   reco_raw[bin-1] = hreco->GetBinContent(bin);
 }

 cout<<"spin="<<spin<<",ch="<<ch<<endl;

 for (int i=0;i<n;i++){
  mass = M[i];
  if(mass<=180) width=1;
  else if (mass>180 && mass<=300) width =5;
  else if (mass>300&&mass<600) width = 10;
  else if (mass>=600&&mass<1000) width = 25;
  else if (mass>=1000&&mass<2000) width = 100;
  else if (mass>=2000&&mass<2300) width=150;
  else width=250;

//merge
  for (int j=(mass-width);j<(mass+width);j++){
  //cout<<"mass="<<mass<<endl;
	//cout<<"gen_bin"<<i<<"="<<gen_raw[j]<<endl;
	//cout<<"reco_bin"<<i<<"="<<reco_raw[j]<<endl;
		gen[i]+=gen_raw[j];
		reco[i]+=reco_raw[j];
  }

  genE[i]=sqrt(gen[i]);
  recoE[i]=sqrt(reco[i]);

  // calculate efficiency and efficiency error
  double g,l,r;
  g=gen[i];
  r=reco[i];
  l=(g-r);
  if(g!=0){
	eff[i]=r/g;
	effE[i]=sqrt(l*l*r+r*r*l)/(g*g);
  }
  else {
	eff[i]=0;
	effE[i]=0;
  }
  cout<<eff[i]<<",";
}

  TF1 *polyFunctot= new TF1("polyFunctot","([0]+[1]*TMath::Erf( (x-[2])/[3] ))*([4]+[5]*x+[6]*x*x+[10]*x*x*x)+[7]*TMath::Gaus(x,[8],[9])", 110., 4100);
  polyFunctot->SetParameters(-4.42749e+00,4.61212e+0,-6.21611e+01,1.13168e+02,2.14321e+00,1.04083e-03,4.89570e-07, 0.03, 200, 30,0);
  polyFunctot->SetParLimits(7,0,0.2);
  polyFunctot->SetParLimits(8,160,210);
  polyFunctot->SetParLimits(9,10,70);
  polyFunctot->SetLineColor(color);
 
  if(spin==0)   polyFunctot->SetLineStyle(1);
  else if(spin==2) {
	polyFunctot->SetParLimits(7,0,0);
	polyFunctot->SetParLimits(8,0,0);
	polyFunctot->SetParLimits(9,0,0); 
	polyFunctot->SetLineStyle(2);
  }

  TGraphErrors *gr = new TGraphErrors (n,M,eff,massE,effE);

  gr->SetMarkerColor(color);
  gr->SetMarkerStyle(marker);
  gr->SetMarkerSize(1);
  cout<<endl<<endl<<"parameters:"<<endl;
  gr->Fit("polyFunctot");
  gr->SetLineColor(color);
  //gr->SetLineStyle(line);
  //TCanvas* c = new TCanvas("c", "c", 1000, 10, 1400, 800);
  //c->cd();
  //gr->Draw("AP"); 
  //c->SaveAs("/afs/cern.ch/user/c/cayou/www/HighMass/160309/test.png");
  delete hreco;
  return gr;
}

void ploteff_sig_spin2_80X_2(){
  TCanvas* c2 = new TCanvas("c2", "c2", 1000, 10, 1400, 800);
  //c2->SetLogx(); 
  c2->SetFillColor(0);
  c2->SetRightMargin(0.13);
  TMultiGraph *mg = new TMultiGraph();

  //TGraphErrors* ggH_4mu = makegr(0,0,kRed,20,1);
  //TGraphErrors* ggH_4e = makegr(0,1,kGreen,20,1);
  //TGraphErrors* ggH_2e2mu = makegr(0,2,kBlue,20,1);
  //TGraphErrors* spin2_4mu = makegr(2,0,kRed+2,24,2);
  //TGraphErrors* spin2_4e = makegr(2,1,kGreen+2,24,2);
  //TGraphErrors* spin2_2e2mu = makegr(2,2,kBlue+2,24,2);
  TGraphErrors* ggH_2l2q = makegr(2,99,kViolet,24,2);

  //mg->Add(ggH_4e);
  //mg->Add(ggH_4mu);
  //mg->Add(ggH_2e2mu);
  mg->Add(ggH_2l2q);
  //mg->Add(spin2_4e);
  //mg->Add(spin2_4mu);
  //mg->Add(spin2_2e2mu);
  mg->Draw("AP");
  mg->GetXaxis()->SetTitle("genHMass [GeV]");
  mg->GetYaxis()->SetTitle("efficiency*acceptance");
  mg->GetYaxis()->SetRangeUser(0,1);
  mg->GetXaxis()->SetRangeUser(100,4100);

  TLegend* leg3 = new TLegend(.9,0.3,0.99,.85);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  //leg3->AddEntry(ggH_4mu,"ggH->2b+ 4mu","pl");
  //leg3->AddEntry(ggH_2e2mu,"ggH->2b+ 2e2mu","pl");
  //leg3->AddEntry(ggH_4e,"ggH->2b+ 4e","pl");
  leg3->AddEntry(ggH_2l2q,"ggh->2b+ 2l2q","pl");

  //leg3->AddEntry(ggH_4mu,"ggH->0+ 4mu","pl");
  //leg3->AddEntry(ggH_2e2mu,"ggH->0+ 2e2mu","pl");
  //leg3->AddEntry(ggH_4e,"ggH->0+ 4e","pl");  
  //leg3->AddEntry(spin2_4mu,"ggH->2b+ 4mu","pl");
  //leg3->AddEntry(spin2_2e2mu,"ggH->2b+ 2e2mu","pl");
  //leg3->AddEntry(spin2_4e,"ggH->2b+ 4e","pl");
  leg3->Draw();


//auto resolve directory name
   char dest[PATH_MAX];  
   sprintf(dest, "%s", gSystem->pwd());

   char * pchar = strstr(dest, "/w/wqin");
   if(pchar != 0)
   strcpy(pchar, "/w/wqin/www/\0");

   pchar = strstr(dest, "/r/rbarr");
   if(pchar != 0)
   strcpy(pchar, "/r/rbarr/www/2l2q/\0");

   pchar = strstr(dest, "/c/cayou");
   if(pchar != 0)
   strcpy(pchar, "/c/cayou/www/HighMass/\0");
   
   sprintf(dest+strlen(dest), "efficiency_%smerged_exclude_%d", (local_ZZCandType==2)?"un":"", exclude);

   TString resolve(dest);

   c2->Update();
   c2->SaveAs(resolve + ".png");
   c2->SaveAs(resolve + ".pdf");
}

void ploteff_sig_spin2_80X() {
	for(exclude=0;exclude<3; ++exclude)
	for(local_ZZCandType=1;local_ZZCandType<3;++local_ZZCandType)
	ploteff_sig_spin2_80X_2();
}
