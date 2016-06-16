#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "THStack.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TCut.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "HZZ2L2QRooPdfs.h"
#include "HZZ2L2QRooPdfs.cc+"
//#include "ZZAnalysis/AnalysisStep/interface/Category.h"
//#include "Math/GenVector/PtEtaPhiM4D.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

//#include "Math/GenVector/LorentzVector.h"
//#include "Math/GenVector/PtEtaPhiM4D.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooSimultaneous.h"
//#include "RooCategory.h"
#include "RooWorkspace.h"
#include "RooHist.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "Math/MinimizerOptions.h"
#include <iomanip>
#include "RooAbsCollection.h"
#include "RooWorkspace.h"

using namespace RooFit;


#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
    void IndiFit(char* channel ="4e", int massBin[]={}, int maxMassBin=0, int inputfiles[]={}, int Nfiles=0)                            
 {

 // ------ root settings ---------
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kFALSE);
  gStyle->SetPadGridY(kFALSE);
  gStyle->SetOptStat("iourme");
  gStyle->SetOptFit(11);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------

  ROOT::Math::MinimizerOptions::SetDefaultTolerance( 1.E-7);

//  float m4l,genM,pm;
//  Short_t z1flav, z2flav,zzsel;
//  float m4l,genM,pm,p_0p,p_2bp,Z1Mass,Z2Mass;
  vector<float> *m4l=0,*p_0p=0,*p_2bp=0,*Z1Mass=0,*Z2Mass=0;
  float genM;
//  double wt_2bp;
  vector<float> *wt=0;
  vector<Short_t> *z1flav=0,*z2flav=0,*zzsel=0;
  int ZZ_pass_ID,ZZ_pass_SIP,ZZ_pass_ISO;
  float weight = 0,PUWeight=0,genHEPMCweight=0;
  char tempmass[100];
  char tempmass2[100];
  double width[100];
  double xMinF=-2000;
  double xMaxF=1500;
  double xMin[100];
  double xMax[100];
//  cout << "Fit range: [" << xMin << " , " << xMax << "]." << endl;

  RooRealVar x("reso","m_{reco}-m_{true}",0.,xMinF,xMaxF,"GeV");
//  x.setBins(100);
  RooRealVar w("myW","myW",1.0,-2000.,1500.);
  RooCategory massrc("massrc","massrc");

  for (int i=0; i<maxMassBin; i++) {
    sprintf(tempmass,"mh%d",massBin[i]);
    massrc.defineType(tempmass,massBin[i]);

  if(channel=="4e") width[i] = 1.9891+0.00554202*(massBin[i])+3.83558e-07*(massBin[i])*(massBin[i]);
  else if(channel=="4mu") width[i] = -4.58023+0.0191778*(massBin[i])+3.74327e-06*(massBin[i])*(massBin[i]);
  else if(channel=="2e2mu") width[i] = -3.28297+0.0153095*(massBin[i])+2.09897e-06*(massBin[i])*(massBin[i]);
  else if(channel=="2l2q") width[i] = -4.58023+0.0191778*(massBin[i])+3.74327e-06*(massBin[i])*(massBin[i]); //FIXME by figuring out what the width should be
//  if(channel=="4e") width[i] = 5.92172+0.00639677*(massBin[i]-600)-1.06763e-06*(massBin[i]-600)*(massBin[i]-600)+4.10566e-09*(massBin[i]-600)*(massBin[i]-600)*(massBin[i]-600);
//  else if(channel=="4mu") width[i] = 8.8625+0.023313*(massBin[i]-600)+0.000014677*(massBin[i]-600)*(massBin[i]-600);
//  else if(channel=="2e2mu") width[i] = 7.5026+0.0156385*(massBin[i]-600)+7.8748e-06*(massBin[i]-600)*(massBin[i]-600)+2.35478e-09*(massBin[i]-600)*(massBin[i]-600)*(massBin[i]-600);
//  xMin[i] = width[i]*(-15);
//  xMax[i] = width[i]*(10);
  xMin[i] = width[i]*(-30);
  xMax[i] = width[i]*(25);
  }

  RooArgSet ntupleVarSet(x,w,massrc);
  RooDataSet dataset("resoM","resoM",ntupleVarSet,WeightVar("myW"));


	char dest[PATH_MAX];
	sprintf(dest, "%s", gSystem->pwd());
	char * pchar = strstr(dest, "prod");
	strcpy(pchar, "prod/\0");
	TString inputDir = dest;

  vector<TString> files;

  char inputfile[100];
  for (int i=0; i<Nfiles; i++) {
//    sprintf(inputfile,"gen_ggH_noCut/mytree/PT13TeV/ggH%d/ZZ4lAnalysis.root",inputfiles[i]);
//    sprintf(inputfile,"ggH_spin0_rw/mytree/PT13TeV/ggH%d/ZZ4lAnalysis.root",inputfiles[i]);
     sprintf(inputfile,"wqin_test/my_tree/PT13TeV/ggH_zz2l2q_M%d/ZZ2l2qAnalysis.root",inputfiles[i]);
    files.push_back(inputfile);
  }

  TChain *candTree = new TChain("ZZTree/candTree");
//  TChain *candTree = new TChain("candTree");

  for (vector<TString>::const_iterator file = files.begin(); file!=files.end(); ++file) {
   candTree->Add(inputDir+(*file));
   }

    int  nentries = candTree->GetEntries();

    //--- ggTree part
    candTree->SetBranchAddress("ZZMass",&m4l);
    candTree->SetBranchAddress("GenHMass",&genM);
    candTree->SetBranchAddress("Z1Flav",&z1flav);
    candTree->SetBranchAddress("Z2Flav",&z2flav);
//    candTree->SetBranchAddress("overallEventWeight",&weight);
    candTree->SetBranchAddress("PUWeight",&PUWeight);
    candTree->SetBranchAddress("genHEPMCweight",&genHEPMCweight);
//    candTree->SetBranchAddress("PoleMass",&pm);
    candTree->SetBranchAddress("ZZsel",&zzsel);
    candTree->SetBranchAddress("p2bplus_VAJHU",&p_2bp);
//    candTree->SetBranchAddress("reweightingweights",&wt);
    candTree->SetBranchAddress("p0plus_VAJHU",&p_0p);
    //candTree->SetBranchAddress("ZZ_pass_ID",&ZZ_pass_ID);
    //candTree->SetBranchAddress("ZZ_pass_ISO",&ZZ_pass_ISO);
    //candTree->SetBranchAddress("ZZ_pass_SIP",&ZZ_pass_SIP);
    candTree->SetBranchAddress("Z1Mass",&Z1Mass);
    candTree->SetBranchAddress("Z2Mass",&Z2Mass);

    for(int k=0; k<nentries; k++){
//    for(int k=0; k<100000; k++){
      candTree->GetEvent(k);

      if(channel=="4mu" && (z1flav->at(0))*(z2flav->at(0)) != 28561) continue;
      if(channel=="4e" && (z1flav->at(0))*(z2flav->at(0)) != 14641) continue;
      if(channel=="2e2mu" && (z1flav->at(0))*(z2flav->at(0)) != 20449) continue;
      if (PUWeight*genHEPMCweight <= 0 ) cout << "Warning! Negative weight events" << endl;

     for (int i=0; i<maxMassBin; i++) {
      ntupleVarSet.setCatIndex("massrc",massBin[i]);
      ntupleVarSet.setRealValue("reso",(m4l->at(0))-genM);
      ntupleVarSet.setRealValue("myW",PUWeight*genHEPMCweight);
//      ntupleVarSet.setRealValue("myW",weight*(wt->at(7)));
      if(((zzsel->at(0))>=100 && x.getVal()>xMin[i] && x.getVal()<xMax[i])&&((massBin[i]<1000&&genM>(massBin[i]-5)&&genM<(massBin[i]+5))||(massBin[i]>=1000&&genM>(massBin[i]*0.95)&&genM<(massBin[i]*1.05))))
//      if((zzsel>=100 && Z1Mass>40 && Z1Mass<120 && Z2Mass>4 && Z2Mass<120 && m4l>100 && x.getVal()>xMin[i] && x.getVal()<xMax[i])&&((massBin[i]<1000&&genM>(massBin[i]-5)&&genM<(massBin[i]+5))||(massBin[i]>=1000&&genM>(massBin[i]*0.95)&&genM<(massBin[i]*1.05))))
//      if((zzsel>=100 && genM<pm && ZZ_pass_ID==1 && ZZ_pass_ISO==1 && ZZ_pass_SIP==1 && Z1Mass>40 && Z1Mass<120 && Z2Mass>4 && Z2Mass<120 && m4l>100 && x.getVal()>xMin[i] && x.getVal()<xMax[i])&&((massBin[i]<1000&&genM>(massBin[i]-5)&&genM<(massBin[i]+5))||(massBin[i]>=1000&&genM>(massBin[i]*0.95)&&genM<(massBin[i]*1.05)))) 
       dataset.add(ntupleVarSet, PUWeight*genHEPMCweight);
//       dataset.add(ntupleVarSet, weight*(wt->at(7)));
      //--------

    }
  }

  cout << "dataset n entries: " << dataset.sumEntries() << endl;

  RooDoubleCB* DCBall[100];
  RooRealVar* mean_ind[100];
  RooRealVar* sigma_ind[100];
  RooRealVar* a1_ind[100];
  RooRealVar* a2_ind[100];
  RooRealVar* n1_ind[100];
  RooRealVar* n2_ind[100];
  RooFitResult* fitres[100];
  RooDataSet* dataset_sub[100];

  for (int i=0; i<maxMassBin; i++) {
    char formulamass[200];

    mean_ind[i]= new RooRealVar("mean_CB","mean_CB",0.,-50., 50.) ;
    sigma_ind[i]= new RooRealVar("sigma_CB","sigma_CB",1, 0, 500);
    a1_ind[i]= new RooRealVar("a1_CB","a1_CB", 1., 0, 5.);
    n1_ind[i]= new RooRealVar("n1_CB","n1_CB", 1., 0, 5.);
    a2_ind[i]= new RooRealVar("a2_CB","a2_CB", 1., 0, 5.);
    n2_ind[i]= new RooRealVar("n2_CB","n2_CB", 1., 0, 20.);

    sprintf(tempmass,"mh%d",massBin[i]);
    sprintf(tempmass2,"massrc == massrc::mh%d",massBin[i]);
    DCBall[i] = new RooDoubleCB("DCBall","Double Crystal ball",x,*mean_ind[i],*sigma_ind[i],*a1_ind[i],*n1_ind[i],*a2_ind[i],*n2_ind[i]);
    dataset_sub[i]= (RooDataSet*)dataset.reduce(tempmass2);
    cout << "Individual fit, mass: "<<massBin[i]<<" , range: [" << xMin[i] << " , " << xMax[i] << "]." << endl<<endl;
    fitres[i] = (RooFitResult*)DCBall[i]->fitTo(*dataset_sub[i],Range(xMin[i],xMax[i]));
//    fitres[i] = (RooFitResult*)DCBall[i]->fitTo(*dataset_sub[i],SumW2Error(1),Range(xMin[i],xMax[i]),Strategy(2),NumCPU(8),Save(true));
    RooArgSet * params = DCBall[i]->getParameters(x);
    char paramfilename[100];
    sprintf(paramfilename,"SingleMassFit_ResoParam_MH%d_%s.txt",massBin[i],channel);
    params->writeToFile(paramfilename) ;

    TCanvas *c1 = new TCanvas("c1","c1",725,725);
    TPad *pad1 = new TPad("pad1","fit",0.05,0.35,0.95,0.97);
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","pull",0.05,0.02,0.95,0.35);
    pad2->Draw();

    int col;
    if(channel =="4mu") col=kOrange+7;
    if(channel =="4e") col=kAzure+2;
    if(channel =="2e2mu") col=kGreen+3;
    if(channel =="2l2q") col=kRed+3;

    char framename[100];
    sprintf(framename,"Resolution_MH%dGeV_%s",massBin[i],channel);

    RooPlot* xframe = x.frame(Range(xMin[i],xMax[i]),Bins(100),Title(framename)) ;
    xframe->GetYaxis()->SetTitleOffset(1.5);
    dataset_sub[i]->plotOn(xframe,DataError(RooAbsData::SumW2), MarkerStyle(kOpenCircle), MarkerSize(1.1));
    DCBall[i]->plotOn(xframe,LineColor(col),Slice(massrc,tempmass),ProjWData(massrc,dataset));

    RooHist* hpull = xframe->pullHist();
    RooPlot* frame2 = x.frame(Range(xMin[i],xMax[i]),Title("Pull Distribution")) ;
    frame2->addPlotable(hpull,"P");

    c1->cd();
    pad1->cd();
    xframe->Draw();
    pad2->cd() ;
    frame2->Draw() ;
    frame2->SetMinimum(-10);
    frame2->SetMaximum(10);
    char filename[100];
//    sprintf(filename,"/afs/cern.ch/user/c/cayou/www/HighMass/160217/resolution_test/Resolution_MH%d_%s_%s.png",massBin[i],fit,channel);
    sprintf(filename,"Resolution_MH%d_%s_singleMassFit.png",massBin[i],channel);

///////////////////////////THIS IS BREAKING???
    cout << endl << filename << endl;
    c1->SaveAs(filename);
    }
  
}

   void ResolutionFit_ind(){
//     int MassBin[] ={120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600};
//     int Inputfiles[]={115,120,124,125,126,130,135,140,145,150,155,160,165,170,175,180,190,200,210,230,250,270,300,350,400,450,500,550,600,700,800,900,1000,1500,2000,/*2500,*/3000};

//     int MassBin[] ={115,140,150,160,300,350,400,500,550,600,700,750,900};
//     int Inputfiles[]={115,140,150,160,300,350,400,500,550,600,700,750,900};

//     int MassBin[] ={750,800,1200,2000,3000,4000};
//     int Inputfiles[]={750,800,1200,2000,3000,4000};

     int MassBin[] ={1000};
     int Inputfiles[]={1000};
     int Nbins=sizeof(MassBin)/sizeof(*MassBin);
     int NFiles=sizeof(Inputfiles)/sizeof(*Inputfiles);

     IndiFit("2l2q",MassBin,Nbins,Inputfiles,NFiles);
//     IndiFit("4mu",MassBin,Nbins,Inputfiles,NFiles);
//     IndiFit("2e2mu",MassBin,Nbins,Inputfiles,NFiles);
//     IndiFit("4e",MassBin,Nbins,Inputfiles,NFiles);

   }
