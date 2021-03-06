#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "TChain.h"
#include "TCanvas.h"
#include <vector>

using namespace RooFit;

void angularDistributions_spin2_mh750_2e2mu(int plotIndex=1, int isqq=0, int is2mp=0,TString decay="all", int binning=40 ){

  gROOT->ProcessLine(".L  ./RooSpinTwo_7D.cxx+");  
 // gROOT->ProcessLine(".L  ./AngularPdfFactory_HWW.cc+");
//  gROOT->ProcessLine(".L  ./TensorPdfFactory_HWW.cc+");
  gROOT->ProcessLine(".L  ./TensorPdfFactory.h");
  gROOT->ProcessLine(".x tdrstyle.cc");


 TString mod ="2mplus";
 TString prod="gg_";
 if(isqq)
	 prod="qqbar_";
	if(!is2mp)
               mod ="2bplus";
 TString linkname = prod+mod;


  RooRealVar* z1mass = new RooRealVar("Gen_mZ1","m_{1} [GeV]",1e-9,120);
  RooRealVar* z2mass = new RooRealVar("Gen_mZ2","m_{2} [GeV]",1e-9,120);
  RooRealVar* hs = new RooRealVar("Gen_costhetastar","cos#theta*",-1,1);
  RooRealVar* h1 = new RooRealVar("Gen_costheta1","cos#theta_{1}",-1,1);
  RooRealVar* h2 = new RooRealVar("Gen_costheta2","cos#theta_{2}",-1,1);
  RooRealVar* Phi = new RooRealVar("Gen_phi","#Phi",-TMath::Pi(),TMath::Pi());
  RooRealVar* Phi1 = new RooRealVar("Gen_phistar1","#Phi_{1}",-TMath::Pi(),TMath::Pi());

/*
  RooRealVar* z1mass = new RooRealVar("Z1Mass","m_{1} [GeV]",1e-9,120);
  RooRealVar* z2mass = new RooRealVar("Z2Mass","m_{2} [GeV]",1e-9,120);
  RooRealVar* hs = new RooRealVar("costhetastar","cos#theta*",-1,1);
  RooRealVar* h1 = new RooRealVar("helcosthetaZ1","cos#theta_{1}",-1,1);
  RooRealVar* h2 = new RooRealVar("helcosthetaZ2","cos#theta_{2}",-1,1);
  RooRealVar* Phi = new RooRealVar("helphi","#Phi",-TMath::Pi(),TMath::Pi());
  RooRealVar* Phi1 = new RooRealVar("phistarZ1","#Phi_{1}",-TMath::Pi(),TMath::Pi());
*/
    //RooRealVar* z1mass = new RooRealVar("wplusmass","m_{l#nu} [GeV]",1e-09,120);
    //RooRealVar* z2mass = new RooRealVar("wminusmass","m(W-)",1e-09,120);
    //RooRealVar* hs = new RooRealVar("costhetastar","cos#theta*",-1,1); 
    //  RooRealVar* h1 = new RooRealVar("costheta1","cos#theta_{1}",-1,1);
    //RooRealVar* h1 = new RooRealVar("costheta1","cos#theta_{1} or cos#theta_{2}",-1,1);
    //RooRealVar* h2 = new RooRealVar("costheta2","cos#theta_{2}",-1,1);
    //RooRealVar* Phi = new RooRealVar("phi","#Phi",-TMath::Pi(),TMath::Pi());
    //RooRealVar* Phi1 = new RooRealVar("phistar1","#Phi_{1}",-TMath::Pi(),TMath::Pi());


  vector<RooRealVar*> measureables;
  measureables.push_back(z1mass);
  measureables.push_back(z2mass);
  measureables.push_back(hs);
  measureables.push_back(h1);
  measureables.push_back(h2);
  measureables.push_back(Phi);
  measureables.push_back(Phi1);

  RooRealVar* mzz = new RooRealVar("mzz","mzz",750,100,1000);

  TensorPdfFactory* Grav = new TensorPdfFactory(z1mass,z2mass,hs,h1,h2,Phi,Phi1,mzz);
  //Grav->make2MH10();
  //Grav->makeMinGrav();
//  Grav->make2BP();
  if(isqq)
	  Grav->makeQQB();
	else
	  Grav->makeGG();
	if(is2mp)
	  Grav->makeMinGrav();
	else
  	Grav->make2bPlus();

//  TChain* treeGrav = new TChain("SelectedTree");
  TChain* treeGrav = new TChain("candTree");
//  treeGrav->Add("../../data/ggH_2bp_zz4l_M750_Ga247_all.root");
  treeGrav->Add("../test_m750_w0.root");
  if(treeGrav->GetEntries()<=0){ cout << "couldn't load minGrav data" << endl; return;}
  RooDataSet* dataGrav = new RooDataSet("dataGrav","dataGrav",treeGrav,RooArgSet(*z1mass,*z2mass,*hs,*h1,*h2,*Phi,*Phi1));

  RooPlot* plot = measureables[plotIndex]->frame(binning);
  plot->GetXaxis()->CenterTitle();
  plot->GetYaxis()->CenterTitle();
  plot->GetYaxis()->SetTitle(" ");
  plot->GetXaxis()->SetNdivisions(-505);

 dataGrav->plotOn(plot,MarkerColor(kRed),MarkerStyle(3),MarkerSize(1.5),LineWidth(0),XErrorSize(0),DataError(RooAbsData::None)); 
 Grav->PDF->plotOn(plot,LineColor(kRed),LineWidth(2));

  TGaxis::SetMaxDigits(3);

  TCanvas* can =new TCanvas("can","can",600,600);

  gStyle->SetPadLeftMargin(0.05);

  char temp[150];
  sprintf(temp,"%s>>minGrav_histo(%i,%i,%i)",measureables[plotIndex]->GetName(),binning,(int)measureables[plotIndex]->getMin(),(int)measureables[plotIndex]->getMin());
  treeGrav->Draw(temp);
  TH1F* minGrav_histo = (TH1F*) gDirectory->Get("minGrav_histo");


//plot->GetYaxis()->SetRangeUser(0,minGrav_histo->GetMaximum()*2.0/1000.);  
plot->Draw();
  
  sprintf(temp,"/afs/cern.ch/user/c/cayou/www/HighMass/test/%s_mH750_%s_%s.eps",measureables[plotIndex]->GetName(),linkname.Data(),decay.Data());
  can->SaveAs(temp);
  sprintf(temp,"/afs/cern.ch/user/c/cayou/www/HighMass/test/%s_mH750_%s_%s.png",measureables[plotIndex]->GetName(),linkname.Data(),decay.Data());
  can->SaveAs(temp);

  delete Grav;

}
