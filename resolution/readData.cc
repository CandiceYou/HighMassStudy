#include "resofit.h"
#include "style.cc"
void readData(char* channel="4e")
 {
  for (int i=0; i<maxMassBin; i++) {
    sprintf(tempmass,"mh%d",massBin[i]);
    massrc.defineType(tempmass,massBin[i]);

    if(strcmp(channel,"4e")==0) 
      width[i] = 1.9891+0.00554202*(massBin[i])+3.83558e-07*(massBin[i])*(massBin[i]);
    else if(strcmp(channel,"4mu")==0)
      width[i] = -4.58023+0.0191778*(massBin[i])+3.74327e-06*(massBin[i])*(massBin[i]);
    else if(strcmp(channel,"2e2mu")==0)
      width[i] = -3.28297+0.0153095*(massBin[i])+2.09897e-06*(massBin[i])*(massBin[i]);
    else if(strcmp(channel,"2l2q")==0){
      width[i] = 3.21246+0.0312538*(massBin[i])-7.29127e-07*(massBin[i])*(massBin[i]);
    }
    xMin[i] = width[i]*(-15);
    xMax[i] = width[i]*(10);
  }

  if (candType == 0)
    sprintf(cType,"J");
  else if (candType == 1)
    sprintf(cType,"jj");

  RooArgSet ntupleVarSet(x,w,massrc);
  dataset = new RooDataSet("resoM","resoM",ntupleVarSet,WeightVar("myW"));

  for (int i=0; i<Nfiles; i++) {
     sprintf(inputfile,"ggHiggs%d/ZZ2l2qAnalysis.root",inputfiles[i]);
    files.push_back(inputfile);
  }

  TChain *candTree = new TChain("ZZTree/candTree");
//  TChain *candTree = new TChain("candTree");

  for (vector<TString>::const_iterator file = files.begin(); file!=files.end(); ++file) {
   candTree->Add(inputDir+(*file));
   }

    int  nentries = candTree->GetEntries();

    //--- ggTree part
    candTree->SetBranchAddress("ZZMassRefit",&m4l);
    candTree->SetBranchAddress("GenHMass",&genM);
    candTree->SetBranchAddress("Z1Flav",&z1flav);
    candTree->SetBranchAddress("Z2Flav",&z2flav);
    candTree->SetBranchAddress("overallEventWeight",&weight);
    candTree->SetBranchAddress("PUWeight",&PUWeight);
    candTree->SetBranchAddress("genHEPMCweight",&genHEPMCweight);
//    candTree->SetBranchAddress("PoleMass",&pm); // not used
    candTree->SetBranchAddress("ZZsel",&ZZsel);
    candTree->SetBranchAddress("Z1Mass",&Z1Mass);
    candTree->SetBranchAddress("Z2Mass",&Z2Mass);
    candTree->SetBranchAddress("ZZCandType",&ZZCandType);
    candTree->SetBranchAddress("Z1tau21", &Z1tau21);

    for(int k=0; k<nentries; k++){
      candTree->GetEvent(k);
//      if (PUWeight*genHEPMCweight <= 0 ) cout << "Warning! Negative weight events" << endl;

//    Find candidate ID, prefer resolved
      int typ=-1 , candID=-1;
      for (int j = 0; j < ZZCandType->size(); j++) {
        for (int i=0; i<maxMassBin; i++) {
          if (((ZZCandType->at(j)==1 && Z1tau21->at(j)<=0.6)||ZZCandType->at(j)==2) && fabs(ZZsel->at(j))>=100 && x.getVal()>xMin[i] && x.getVal()<xMax[i] && Z1Mass->at(j)>=70 && Z1Mass->at(j)<=105 && Z2Mass->at(j)>=60 && ((massBin[i]<1000&&genM>(massBin[i]-25)&&genM<(massBin[i]+25))||(massBin[i]>=1000&&genM>(massBin[i]*0.95)&&genM<(massBin[i]*1.05)))){
            if (ZZCandType->at(j)==1) {typ=0; candID=j;}  //merged, SR
            else if (ZZCandType->at(j)==2) {typ=1; candID=j;break;}  //resolved, SR

            ntupleVarSet.setCatIndex("massrc",massBin[i]);
            ntupleVarSet.setRealValue("reso",(m4l->at(candID))-genM);
            ntupleVarSet.setRealValue("myW",PUWeight*genHEPMCweight);
          }
        }
      }
//    Sort entry into approrpriate dataset      
      if (typ==candType)
        dataset->add(ntupleVarSet, PUWeight*genHEPMCweight);
    } 

  cout << "dataset n entries: " << dataset->sumEntries() << endl;

  for (int i=0; i<maxMassBin; i++) {
    sprintf(tempmass2,"massrc == massrc::mh%d",massBin[i]);
    dataset_sub[i]= (RooDataSet*)dataset->reduce(tempmass2);   
  }
}
