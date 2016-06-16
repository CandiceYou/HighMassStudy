#include <TROOT.h>
#include <vector>
#include "TLorentzVector.h"
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/src/computeAngles.h>
#include <ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h>
#include "TLorentzRotation.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
//#include "loadMela.cc"

using namespace MEMNames;
using namespace std;

void testMEM(TString fin,TString fout){
  const double mu_mass(.105658);
  const double e_mass(.000511);
  const double tau_mass(1.77682);


  TChain *t = new TChain("candTree");
  t->Add(fin);

  TFile* file = new TFile(fout,"RECREATE");
  TTree* tree = (TTree*) t->CloneTree(0,"fast");

  MEMs test(13,125,"CTEQ6L");

//gen info

//  double Gen_p2bplus_VAJHU;
//  double Gen_p0plus_VAJHU;
//  double wt_2bp;

  Float_t GenLep1Pt,GenLep1Eta,GenLep1Phi;
  Float_t GenLep2Pt,GenLep2Eta,GenLep2Phi;
  Float_t GenLep3Pt,GenLep3Eta,GenLep3Phi;
  Float_t GenLep4Pt,GenLep4Eta,GenLep4Phi;
  float costheta1=0, costheta2=0, phi=0, costhetastar=0, phistar1=0;
  double Gen_mZ1=0,Gen_mZ2=0;

  Short_t GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;

//  tree->Branch("Gen_p2bplus_VAJHU",&Gen_p2bplus_VAJHU,"Gen_p2bplus_VAJHU/D");
//  tree->Branch("Gen_p0plus_VAJHU",&Gen_p0plus_VAJHU,"Gen_p0plus_VAJHU/D");
//  tree->Branch("wt_2bp",&wt_2bp,"wt_2bp/D");
  tree->Branch("Gen_costheta1",&costheta1,"Gen_costheta1/F");
  tree->Branch("Gen_costheta2",&costheta2,"Gen_costheta2/F");
  tree->Branch("Gen_phi",&phi,"Gen_phi/F");
  tree->Branch("Gen_costhetastar",&costhetastar,"Gen_costhetastar/F");
  tree->Branch("Gen_phistar1",&phistar1,"Gen_phistar1/F");
  tree->Branch("Gen_mZ1",&Gen_mZ1,"Gen_mZ1/D");
  tree->Branch("Gen_mZ2",&Gen_mZ2,"Gen_mZ2/D");

  t->SetBranchAddress("GenLep1Pt",&GenLep1Pt);
  t->SetBranchAddress("GenLep1Eta",&GenLep1Eta);
  t->SetBranchAddress("GenLep1Phi",&GenLep1Phi);
  t->SetBranchAddress("GenLep1Id",&GenLep1Id);

  t->SetBranchAddress("GenLep2Pt",&GenLep2Pt);
  t->SetBranchAddress("GenLep2Eta",&GenLep2Eta);
  t->SetBranchAddress("GenLep2Phi",&GenLep2Phi);
  t->SetBranchAddress("GenLep2Id",&GenLep2Id);

  t->SetBranchAddress("GenLep3Pt",&GenLep3Pt);
  t->SetBranchAddress("GenLep3Eta",&GenLep3Eta);
  t->SetBranchAddress("GenLep3Phi",&GenLep3Phi);
  t->SetBranchAddress("GenLep3Id",&GenLep3Id);

  t->SetBranchAddress("GenLep4Pt",&GenLep4Pt);
  t->SetBranchAddress("GenLep4Eta",&GenLep4Eta);
  t->SetBranchAddress("GenLep4Phi",&GenLep4Phi);
  t->SetBranchAddress("GenLep4Id",&GenLep4Id);

  for(int i=0 ; i<t->GetEntries(); i++){
//  for(int i=0 ; i<10; i++){
  vector<int> id;
  vector<TLorentzVector> p4;
  TLorentzVector lep1,lep2,lep3,lep4,z1,z2;

  t->GetEntry(i);
//  cout<<"entry:"<<i<<endl;
  double m1=0,m2=0,m3=0,m4=0; 
  if (abs(GenLep1Id)==13) m1=mu_mass;else if(abs(GenLep1Id)==11) m1=e_mass;else m1=tau_mass;
  if (abs(GenLep2Id)==13) m2=mu_mass;else if(abs(GenLep2Id)==11) m2=e_mass;else m2=tau_mass;
  if (abs(GenLep3Id)==13) m3=mu_mass;else if(abs(GenLep3Id)==11) m3=e_mass;else m3=tau_mass;
  if (abs(GenLep4Id)==13) m4=mu_mass;else if(abs(GenLep4Id)==11) m4=e_mass;else m4=tau_mass;

//  cout<<"Id="<<GenLep1Id<<","<<GenLep2Id<<","<<GenLep3Id<<","<<GenLep4Id<<endl;
//cout<<"mass="<<m1<<","<<m2<<","<<m3<<","<<m4<<endl;

  lep1.SetPtEtaPhiM(GenLep1Pt,GenLep1Eta,GenLep1Phi,m1);
  lep2.SetPtEtaPhiM(GenLep2Pt,GenLep2Eta,GenLep2Phi,m2);
  lep3.SetPtEtaPhiM(GenLep3Pt,GenLep3Eta,GenLep3Phi,m3);
  lep4.SetPtEtaPhiM(GenLep4Pt,GenLep4Eta,GenLep4Phi,m4);

  z1=lep1+lep2;
  z2=lep3+lep4;
  Gen_mZ1=z1.M();
  Gen_mZ2=z2.M();

  p4.push_back(lep1);
  p4.push_back(lep2);
  p4.push_back(lep3);
  p4.push_back(lep4);

  id.push_back(GenLep1Id);
  id.push_back(GenLep2Id);
  id.push_back(GenLep3Id);
  id.push_back(GenLep4Id);

//  test.computeME(MEMNames::k2bplus, MEMNames::kJHUGen ,p4,id,Gen_p2bplus_VAJHU);
//  test.computeME(MEMNames::kSMHiggs ,MEMNames::kJHUGen ,p4,id,Gen_p0plus_VAJHU); 

  mela::computeAngles(lep1,GenLep1Id,lep2,GenLep2Id,lep3,GenLep3Id,lep4,GenLep4Id,costhetastar,costheta1,costheta2,phi,phistar1);

//  wt_2bp=Gen_p2bplus_VAJHU/Gen_p0plus_VAJHU;

  tree->Fill();

  } 

  file->cd();
  tree->Write();
  file->Close();

}

void addGenMela(){

     testMEM("/afs/cern.ch/user/w/wqin/CMSSW_8_0_3/src/ZZAnalysis/AnalysisStep/test/prod/2016_gg2bp_merged/gg2bp_W0_M750/ZZ4lAnalysis.root","test.root");

/*
  TString  inputDir = "/afs/cern.ch/user/w/wqin/CMSSW_8_0_3/src/ZZAnalysis/AnalysisStep/test/prod/"; 

  int sampleMass_ggH[]={800};
  int Nfiles_ggH=sizeof(sampleMass_ggH)/sizeof(*sampleMass_ggH);
  char inputfile_ggH[100];
  char outputfile_ggH[100];

   for (int i=0; i<Nfiles_ggH; i++) {
     sprintf(inputfile_ggH,"2016_gg2bp_merged/gg2bp_W0_M%d/ZZ4lAnalysis.root",sampleMass_ggH[i]);
     sprintf(outputfile_ggH,"data/gg2bp_W0_M%d_ZZ4lAnalysis_addwt.root",sampleMass_ggH[i]);
     testMEM(inputDir+inputfile_ggH,outputfile_ggH);
   }
*/
}
