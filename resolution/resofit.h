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
//#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
//#include "ZZAnalysis/AnalysisStep/interface/Category.h"
//#include "Math/GenVector/LorentzVector.h"
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
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>


//  float m4l,genM,pm,Z1Mass,Z2Mass;
  vector<float> *m4l=0,*Z1Mass=0,*Z2Mass=0;
  float genM;
  vector<float> *wt=0;
  vector<Short_t> *z1flav=0,*z2flav=0,*ZZsel=0,*ZZCandType=0,*Z1tau21=0,*Z1Pt=0,*Z2Pt=0;
//  Short_t z1flav, z2flav,zzsel;
  float weight = 0,PUWeight=0,genHEPMCweight=0;
  char cType[10];
  char tempmass[100];
  char tempmass2[100];
  double width[100];
  double xMinF=-2000;
  double xMaxF=1500;
  double xMin[100];
  double xMax[100];

  RooRealVar x("reso","m_{reco}-m_{true}",0.,xMinF,xMaxF,"GeV");
  RooRealVar w("myW","myW",1.0,-2000.,1500.);
  RooCategory massrc("massrc","massrc");
  RooDataSet* dataset;
  RooDataSet* dataset_sub[100];
  vector<TString> files;
  char inputfile[100];

  RooDoubleCB* DCBall[100];
  RooRealVar* mean_ind[100];
  RooRealVar* sigma_ind[100];
  RooRealVar* a1_ind[100];
  RooRealVar* a2_ind[100];
  RooRealVar* n1_ind[100];
  RooRealVar* n2_ind[100];
  RooFitResult* fitres[100];

  TString inputDir = "root://eoscms.cern.ch//eos/cms/store/user/covarell/2l2qTrees/160725/";
//  TString inputDir = "/afs/cern.ch/work/c/cayou/public/forWenzerRobert/2l2qsamples_new/";
//  TString inputDir = "/afs/cern.ch/work/c/cayou/public/forWenzerRobert/2l2qsamples_2bp/";
//  char dest[PATH_MAX];
//  sprintf(dest, "%s", gSystem->pwd());
//  char * pchar = strstr(dest, "prod");
//  strcpy(pchar, "prod/\0");
//  TString inputDir = dest;
  
  char channel[100]; 
//  int massBin[]={1000};
//  int massBin[]={400,450,500,550,600,700,750,800,900,1000,1200,1400,1600,1800,2000};
  int massBin[]={600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
//  int massBin[]={750,800,1000,1200,1600,2000,2500,3000,3500,4000};
//  int inputfiles[]={750,800,1200,2000,3000,4000};
//  int inputfiles[]={1000};
//  int inputfiles[]={400,450,500,550,600,700,750,800,900,1000,1500,2000};
  int inputfiles[]={600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500};
  short candType=1; //0 for merged jet (J), 1 for two resolved jets (jj)
  int exclude=0; //0 for exclude nothing; 1 exclude events with one jet type; 2 exclude events with both jet types
  int maxMassBin=sizeof(massBin)/sizeof(*massBin);;
  int Nfiles=sizeof(inputfiles)/sizeof(*inputfiles);

 float tau21bin[25] = {
    -0.00769231,
    0.0346154,  
    0.0769231,  
    0.119231,   
    0.161538,   
    0.203846,   
    0.246154,   
    0.288462,   
    0.330769,   
    0.373077,   
    0.415385,   
    0.457692,   
    0.5,        
    0.542308,   
    0.584615,   
    0.626923,   
    0.669231,   
    0.711538,   
    0.753846,   
    0.796154,   
    0.838461,   
    0.880769,   
    0.923077,   
    0.965385,   
    1.00769} ;
  float tau21corr[24] = {0,	      	   
    0      ,	   
    0,	   
    0.290173  , 
    -0.377686 , 
    0.0977722  ,
    0.412889   ,
    0.322422   ,
    -0.169658  ,
    -0.272581  ,
    -0.384607  ,
    0.109909   ,
    -0.123669  ,
    -0.158826  ,
    0.0764657  ,
    0.155703   ,
    -0.0750966 ,
    -0.0484963 ,
    0.239322   ,
    0.107814   ,
    -1.40918   ,
    0	   ,
    0	  , 
    0	   };
