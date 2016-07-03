#include <string.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"


void selective_copy(char * fname, float pm){
  TString filename = fname;

  if (gSystem->AccessPathName(filename)) return;

  //open ROOT file and open candTree
  TFile *f = new TFile(filename,"update"); //option: NEW, CREATE, RECREATE, UPDATE, or READ
  TTree *T = (TTree*)f->Get("ZZTree/candTree");

  //create a new branch
  //F is for floating point: https://root.cern.ch/doc/master/classTBranch.html#ae30a4edb372811b4f02b485b7b0dfaca
  TBranch *polemass = T->Branch("PoleMass",&pm,"PoleMass/F");
  Long64_t nentries = T->GetEntries();

  //setup access to GenHMass
  float GenHMass;
  T->SetBranchAddress("GenHMass", &GenHMass);

  char * pchar = strstr(fname, ".root");
  strcpy(pchar, "_pruned.root\0");
  TString new_filename = fname;

  TFile *newfile = new TFile(new_filename,"recreate");
  TTree *newtree = T->CloneTree(0);

  for (Long64_t i=0;i<nentries;i++) {
      T->GetEntry(i);
      polemass->Fill();
      if(GenHMass <= pm) {
        newtree->Fill();
    }
    //else printf("Skipping branch w/ GenHMass = %d > %d\n", GenHMass, pm);
  }

  T->Print();
  T->Write();
  newtree->Write();
  delete f;
  delete newfile;
}

void addPMbranch(){
  char dest[PATH_MAX];
  sprintf(dest, "%s", gSystem->pwd());
  char * pchar = strstr(dest, "PT13TeV");
  strcpy(pchar, "PT13TeV/\0");
  //  /afs/cern.ch/user/r/rbarr/Analysis/CMSSW_8_0_9/src/ZZAnalysis/AnalysisStep/test/prod/RobertTest/mytree/PT13TeV/

  char inputfile_ggH[PATH_MAX];

  sprintf(inputfile_ggH,"%sggH%d/ZZ4lAnalysis.root",dest, 115);
  selective_copy(inputfile_ggH, 115);

}
