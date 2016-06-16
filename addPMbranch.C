void addbr(TString filename, double Pmass){
 if (!gSystem->AccessPathName(filename)){
   TFile *f = new TFile(filename,"update");
   TTree *T = (TTree*)f->Get("ZZTree/candTree");
   float pm;
   TBranch *polemass = T->Branch("PoleMass",&pm,"PoleMass/F");
   Long64_t nentries = T->GetEntries();
   for (Long64_t i=0;i<nentries;i++) {
      T->GetEntry(i);
      pm = Pmass;
      polemass->Fill();
   }
   T->Print();
   T->Write();
   delete f;
 }
}

void addPMbranch(){
  TString  inputDir = "/afs/cern.ch/work/c/cayou/HighMass_RunII/CMSSW_7_6_3_patch2/src/ZZAnalysis/AnalysisStep/test/prod/";
//  int sampleMass_ggH[]={1000,1500,2000,/*2500,*/3000};
  int sampleMass_ggH[]={1000};
  int sampleMass_VBF[]={1000,1500,2000,2500,3000};
//  int sampleMass_ggH[]={115,120,124,125,126,130,135,140,145,150,155,160,165,170,175,180,190,200,210,230,250,270,300,350,400,450,500,550,600,700,800,900,1000,1500,2000,/*2500,*/3000};
//  int sampleMass_VBF[]={115,120,124,125,126,/*130,*/135,140,145,150,155,160,165,170,175,180,190,200,210,230,250,270,300,350,400,450,500,550,600,700,800,900,1000,1500,2000,2500,3000};
  int Nfiles_ggH=sizeof(sampleMass_ggH)/sizeof(*sampleMass_ggH);
  int Nfiles_VBF=sizeof(sampleMass_VBF)/sizeof(*sampleMass_VBF);
  char inputfile_ggH[100];
  char inputfile_VBF[100];

   for (int i=0; i<Nfiles_ggH; i++) {
//     sprintf(inputfile_ggH,"highPtMuonId200/mytree/PT13TeV/ggH%d/ZZ4lAnalysis.root",sampleMass_ggH[i]);
     sprintf(inputfile_ggH,"highPtMuonId200/mytree/PT13TeV/ZZ4lAnalysis.root",sampleMass_ggH[i]);
     addbr(inputDir+inputfile_ggH,sampleMass_ggH[i]);
   }

/*
   for (int i=0; i<Nfiles_VBF; i++) {
     sprintf(inputfile_VBF,"highPtMuonId200/mytree/PT13TeV/VBFH%d/ZZ4lAnalysis.root",sampleMass_VBF[i]);
     addbr(inputDir+inputfile_VBF,sampleMass_VBF[i]);
   }
*/
}
