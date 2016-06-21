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
  vector<Short_t> *z1flav=0,*z2flav=0,*zzsel=0,*candType=0;
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

  char dest[PATH_MAX];
  sprintf(dest, "%s", gSystem->pwd());
  char * pchar = strstr(dest, "prod");
  strcpy(pchar, "prod/\0");
  TString inputDir = dest;
  
  char channel[100]; 
  int massBin[]={350,400,450,500,600,700,750,850,1000,1200,2000};
  int inputfiles[]={350, 400, 450, 500, 600, 700, 750, 1000, 2000};
  short ZZCandType=1; //1 for merged jet (J), 2 for two resolved jets (jj)
  int exclude=2; //0 for exclude nothing; 1 exclude events with one jet type; 2 exclude events with both jet types
  int maxMassBin=sizeof(massBin)/sizeof(*massBin);;
  int Nfiles=sizeof(inputfiles)/sizeof(*inputfiles);

