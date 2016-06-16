void copytree(){
//  TFile* file = TFile::Open("/afs/cern.ch/user/w/wqin/CMSSW_8_0_3/src/ZZAnalysis/AnalysisStep/test/prod/2016_gg2bp_merged/gg2bp_W0_M750/ZZ4lAnalysis.root");
 TFile* file = TFile::Open("test.root");
  TTree* originalTree = (TTree*)file->Get("candTree");
  TFile *newfile = new TFile("test_m750_w0.root","recreate");
  TTree* treeGrav = originalTree->CopyTree("ZZsel>=100&&Z1Flav*Z2Flav ==20449");
  treeGrav->Write();
}
