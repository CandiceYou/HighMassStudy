#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "TChain.h"
#include "TCanvas.h"
#include <vector>

using namespace RooFit;

void angularDistributions_spin2_2bplus(int plotIndex=0, int binning=80){

  //gROOT->ProcessLine(".L  ../PDFs/RooXZsZs_5D.cxx+");
//  gROOT->ProcessLine(".L  ../PDFs/RooSpinZero_7DComplex.cc+");
//  gROOT->ProcessLine(".L  ../PDFs/RooSpinOne_7D.cxx+");  
//  gROOT->ProcessLine(".L  ../PDFs/RooSpinTwo_7D.cxx+");  
//  gROOT->ProcessLine(".L  ../src/AngularPdfFactory.cc+");
//  gROOT->ProcessLine(".L  ../src/ScalarPdfFactory.cc+");
//  gROOT->ProcessLine(".L  ../src/VectorPdfFactory.cc+");
//  gROOT->ProcessLine(".L  ../src/TensorPdfFactory.cc+");
//  gROOT->ProcessLine(".L  ~/tdrstyle.C");
//  setTDRStyle();

  gROOT->ProcessLine(".L  ./RooSpinTwo_7D.cxx+");
  gROOT->ProcessLine(".L  ./TensorPdfFactory.h");
          gROOT->ProcessLine(".x tdrstyle.cc");
  // observables
 /* 
  RooRealVar* z1mass = new RooRealVar("z1mass","m_{1} [GeV]",40,110);
  RooRealVar* z2mass = new RooRealVar("z2mass","m_{2} [GeV]",1e-09,65);
  RooRealVar* hs = new RooRealVar("costhetastar","cos#theta*",-1,1);
  RooRealVar* h1 = new RooRealVar("costheta1","cos#theta_{1}",-1,1);
  RooRealVar* h2 = new RooRealVar("costheta2","cos#theta_{2}",-1,1);
  RooRealVar* Phi = new RooRealVar("phi","#Phi",-TMath::Pi(),TMath::Pi());
  RooRealVar* Phi1 = new RooRealVar("phistar1","#Phi_{1}",-TMath::Pi(),TMath::Pi());
 */
 
  RooRealVar* z1mass = new RooRealVar("Z1Mass","m_{1} [GeV]",1e-9,110);
  RooRealVar* z2mass = new RooRealVar("Z2Mass","m_{2} [GeV]",1e-9,110);
  RooRealVar* hs = new RooRealVar("costhetastar","cos#theta*",-1,1);
  RooRealVar* h1 = new RooRealVar("helcosthetaZ1","cos#theta_{1}",-1,1);
  RooRealVar* h2 = new RooRealVar("helcosthetaZ2","cos#theta_{2}",-1,1);
  RooRealVar* Phi = new RooRealVar("helphi","#Phi",-TMath::Pi(),TMath::Pi());
  RooRealVar* Phi1 = new RooRealVar("phistarZ1","#Phi_{1}",-TMath::Pi(),TMath::Pi());
 

  vector<RooRealVar*> measureables;
  measureables.push_back(z1mass);
  measureables.push_back(z2mass);
  measureables.push_back(hs);
  measureables.push_back(h1);
  measureables.push_back(h2);
  measureables.push_back(Phi);
  measureables.push_back(Phi1);

  RooRealVar* mzz = new RooRealVar("mzz","mzz",750,100,1000);

  TensorPdfFactory* MinGrav = new TensorPdfFactory(z1mass,z2mass,hs,h1,h2,Phi,Phi1,mzz);
  MinGrav->makeMinGrav(); 
  MinGrav->g1Val->setVal(0.0); 
  MinGrav->makeParamsConst(true);

 // TensorPdfFactory* TwohPlus = new TensorPdfFactory(z1mass,z2mass,hs,h1,h2,Phi,Phi1,mzz);
 // TwohPlus->make2hPlus();
 // TwohPlus->makeParamsConst(true);

 // TensorPdfFactory* TwohMinus = new TensorPdfFactory(z1mass,z2mass,hs,h1,h2,Phi,Phi1,mzz);
 // TwohMinus->make2hMinus();
 // TwohMinus->makeParamsConst(true);

 // TChain* treeMinGrav = new TChain("angles");
 // treeMinGrav->Add("/scratch0/hep/whitbeck/OLDHOME/4lHelicity/generatorJHU_V02-01-00/minGrav_store/minGrav_125GeV.root");
  TChain* treeMinGrav = new TChain("newTree");
//  treeMinGrav->Add("../inputFile/Graviton2BPToZZTo4L_M-750_13TeV-JHUGenV6_2e2mu.root");
  treeMinGrav->Add("/afs/cern.ch/work/c/cayou/JHUGen/projection_ZZ/inputFile/Graviton2BPToZZTo4L_M-750_13TeV-JHUGenV6_2e2mu.root");
  if(treeMinGrav->GetEntries()<=0){ cout << "couldn't load minGrav data" << endl; return;}
  RooDataSet* dataSM = new RooDataSet("dataMinGrav","dataMinGrav",treeMinGrav,RooArgSet(*z1mass,*z2mass,*hs,*h1,*h2,*Phi,*Phi1));

//  TChain* tree2hPlus = new TChain("angles");
//  tree2hPlus->Add("/scratch0/hep/whitbeck/OLDHOME/4lHelicity/generatorJHU_V02-01-00/2hPlus_store/2hPlus_125GeV.root");
//  if(tree2hPlus->GetEntries()<=0){ cout << "couldn't load 2hPlus data" << endl; return;}
//  RooDataSet* data = new RooDataSet("data2hPlus","data2hPlus",tree2hPlus,RooArgSet(*z1mass,*z2mass,*hs,*h1,*h2,*Phi,*Phi1));

//  TChain* tree2hMinus = new TChain("angles");
//  tree2hMinus->Add("/scratch0/hep/whitbeck/OLDHOME/4lHelicity/generatorJHU_V02-01-00/2hMinus_store/2hMinus_125GeV.root");
//  if(tree2hMinus->GetEntries()<=0){ cout << "couldn't load 2hMinus data" << endl; return;}
//  RooDataSet* data2hMinus = new RooDataSet("data2hMinus","data2hMinus",tree2hMinus,RooArgSet(*z1mass,*z2mass,*hs,*h1,*h2,*Phi,*Phi1));

  RooPlot* plot = measureables[plotIndex]->frame(binning);
  plot->GetXaxis()->CenterTitle();
  plot->GetYaxis()->CenterTitle();
  plot->GetYaxis()->SetTitle(" ");
  plot->GetXaxis()->SetNdivisions(-505);

  dataMinGrav->plotOn(plot,MarkerColor(kRed),MarkerStyle(4),MarkerSize(1.5),LineWidth(0),XErrorSize(0),Rescale((3./15.)*.001),DataError(RooAbsData::None)); 
  MinGrav->PDF->plotOn(plot,LineColor(kRed),LineWidth(2),Normalization( (3./15.)*0.001 ));

//  data2hPlus->plotOn(plot,MarkerColor(kBlue),MarkerStyle(27),MarkerSize(1.9),XErrorSize(0),Rescale(.001),DataError(RooAbsData::None)); 
//  TwohPlus->PDF->plotOn(plot,LineColor(kBlue),LineWidth(2),Normalization( 0.001 ));

//  data2hMinus->plotOn(plot,MarkerColor(kGreen+3),MarkerStyle(25),MarkerSize(1.5),XErrorSize(0),Rescale(.001),DataError(RooAbsData::None)); 
//  TwohMinus->PDF->plotOn(plot,LineColor(kGreen+3),LineWidth(2),Normalization( 0.001 ));

  TGaxis::SetMaxDigits(3);

  TCanvas* can =new TCanvas("can","can",600,600);

  gStyle->SetPadLeftMargin(0.05);

  char temp[150];
  sprintf(temp,"%s>>minGrav_histo(%i,%i,%i)",measureables[plotIndex]->GetName(),binning,(int)measureables[plotIndex]->getMin(),(int)measureables[plotIndex]->getMin());
  treeMinGrav->Draw(temp);
  TH1F* minGrav_histo = (TH1F*) gDirectory->Get("minGrav_histo");
//  sprintf(temp,"%s>>TwohPlus_histo(%i,%i,%i)",measureables[plotIndex]->GetName(),binning,(int)measureables[plotIndex]->getMin(),(int)measureables[plotIndex]->getMin());
//  tree2hPlus->Draw(temp);
//  TH1F* TwohPlus_histo = (TH1F*) gDirectory->Get("TwohPlus_histo");
//  sprintf(temp,"%s>>TwohMinus_histo(%i,%i,%i)",measureables[plotIndex]->GetName(),binning,(int)measureables[plotIndex]->getMin(),(int)measureables[plotIndex]->getMin());
//  tree2hMinus->Draw(temp);
//  TH1F* TwohMinus_histo = (TH1F*) gDirectory->Get("TwohMinus_histo");
  
//  plot->GetYaxis()->SetRangeUser(0,max(max(TwohMinus_histo->GetMaximum(),TwohPlus_histo->GetMaximum()),(2./15.)*minGrav_histo->GetMaximum())*1.3/1000.);
//treeSM->Add("../../rootfiles/300k/2e2mu/SMHiggsToZZTo4L_M-125_8TeV_POWHEG-JHUgenV3-pythia6_false_2e2mu.root");
plot->GetYaxis()->SetRangeUser(0,(2./15.)*minGrav_histo->GetMaximum()*1.7/10000.);  
plot->Draw();
  
  char temp[150];
  sprintf(temp,"../plots/%s_13TeV_750GeV_2bplus.eps",measureables[plotIndex]->GetName());
  can->SaveAs(temp);
  sprintf(temp,"../plots/%s_13TeV_750GeV_2bplus.png",measureables[plotIndex]->GetName());
  can->SaveAs(temp);

  delete MinGrav;
//  delete TwohPlus;
//  delete TwohMinus;

}
