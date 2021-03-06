#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "TChain.h"
#include "TCanvas.h"
#include <vector>

using namespace RooFit;

void angularDistributions_spin2_mh750_4l(int plotIndex=1, int isqq=0, int is2mp=0,int binning=40 ){

  gROOT->ProcessLine(".L  ./RooSpinTwo_7D.cxx+");  
 // gROOT->ProcessLine(".L  ./AngularPdfFactory_HWW.cc+");
//  gROOT->ProcessLine(".L  ./TensorPdfFactory_HWW.cc+");
  gROOT->ProcessLine(".L  ./TensorPdfFactory.h");
	  gROOT->ProcessLine(".x tdrstyle.cc");
//  gROOT->ProcessLine(".L  ./tdrstyle.cc");
//  setTDRStyle();

 TString outname ="2mplus"; 
 TString finals="gg";
 if(isqq)
	 finals="qq";
	if(!is2mp)
		outname="2bplus";
 TString linkname = outname+"_"+finals;

/* 
  RooRealVar* z1mass = new RooRealVar("Z1_m","m_{1} [GeV]",1e-9,120);
  RooRealVar* z2mass = new RooRealVar("Z2_m","m_{2} [GeV]",1e-09,120);
  RooRealVar* hs = new RooRealVar("costhetastar_2","cos#theta*",-1,1);
  RooRealVar* h1 = new RooRealVar("helcosthetaZ1_2","cos#theta_{1}",-1,1);
  RooRealVar* h2 = new RooRealVar("helcosthetaZ2_2","cos#theta_{2}",-1,1);
  RooRealVar* Phi = new RooRealVar("helphi_2","#Phi",-TMath::Pi(),TMath::Pi());
  RooRealVar* Phi1 = new RooRealVar("phistarZ1_2","#Phi_{1}",-TMath::Pi(),TMath::Pi());
  RooRealVar* channel= new RooRealVar("channel","channel",0,3);
*/

  RooRealVar* z1mass = new RooRealVar("Z1Mass","m_{1} [GeV]",1e-9,120);
  RooRealVar* z2mass = new RooRealVar("Z2Mass","m_{2} [GeV]",1e-9,120);
  RooRealVar* hs = new RooRealVar("costhetastar","cos#theta*",-1,1);
  RooRealVar* h1 = new RooRealVar("helcosthetaZ1","cos#theta_{1}",-1,1);
  RooRealVar* h2 = new RooRealVar("helcosthetaZ2","cos#theta_{2}",-1,1);
  RooRealVar* Phi = new RooRealVar("helphi","#Phi",-TMath::Pi(),TMath::Pi());
  RooRealVar* Phi1 = new RooRealVar("phistarZ1","#Phi_{1}",-TMath::Pi(),TMath::Pi());
  RooRealVar* channel= new RooRealVar("flavortype","flavortype",0,3); 

    //RooRealVar* z1mass = new RooRealVar("wplusmass","m_{l#nu} [GeV]",1e-09,120);
    //RooRealVar* z2mass = new RooRealVar("wminusmass","m(W-)",1e-09,120);
    //RooRealVar* hs = new RooRealVar("costhetastar","cos#theta*",-1,1); 
//  //  RooRealVar* h1 = new RooRealVar("costheta1","cos#theta_{1}",-1,1);
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

  TChain* treeGrav = new TChain("newTree");
//  treeGrav->Add("/afs/cern.ch/work/c/cayou/JHUGen/projection_ZZ/inputFile/Graviton2BPToZZTo4L_M-750_13TeV-JHUGenV6_2e2mu.root");
  treeGrav->Add("/afs/cern.ch/work/c/cayou/SampleFromHeshy/mH750_"+linkname+"_2e2mu.root");
  if(treeGrav->GetEntries()<=0){ cout << "couldn't load minGrav data" << endl; return;}
  RooDataSet* dataGrav = new RooDataSet("dataGrav","dataGrav",treeGrav,RooArgSet(*z1mass,*z2mass,*hs,*h1,*h2,*Phi,*Phi1));
//  RooDataSet* dataGrav = new RooDataSet("dataGrav","dataGrav",treeGrav,RooArgSet(*z1mass,*z2mass,*hs,*h1,*h2,*Phi,*Phi1,*channel),"channel==1");

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
  
  char temp[150];
  sprintf(temp,"/afs/cern.ch/user/c/cayou/www/HighMass/151218/%s_mH750_%s_4l.eps",measureables[plotIndex]->GetName(),linkname.Data());
  can->SaveAs(temp);
  sprintf(temp,"/afs/cern.ch/user/c/cayou/www/HighMass/151218/%s_mH750_%s_4l.png",measureables[plotIndex]->GetName(),linkname.Data());
  can->SaveAs(temp);

  delete Grav;

}
