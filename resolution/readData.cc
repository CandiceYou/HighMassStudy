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
    xMin[i] = width[i]*(-29);
    xMax[i] = width[i]*(24);
  }

  if (candType == 0)
    sprintf(cType,"J");
  else if (candType == 1)
    sprintf(cType,"jj");

  RooArgSet ntupleVarSet(x,w,massrc);
  dataset = new RooDataSet("resoM","resoM",ntupleVarSet,WeightVar("myW"));

  for (int i=0; i<Nfiles; i++) {
//    sprintf(inputfile,"ggHiggs%d/ZZ2l2qAnalysis.root",inputfiles[i]);
//    files.push_back(inputfile);
//    sprintf(inputfile,"VBFHiggs%d/ZZ2l2qAnalysis.root",inputfiles[i]);
//    files.push_back(inputfile);
    sprintf(inputfile,"BulkGrav%d/ZZ2l2qAnalysis.root",inputfiles[i]);
    files.push_back(inputfile);
  }

  TChain *candTree = new TChain("ZZTree/candTree");
//  TChain *candTree = new TChain("candTree");

  for (vector<TString>::const_iterator file = files.begin(); file!=files.end(); ++file) {
   candTree->Add(inputDir+(*file));
   }

    int  nentries = candTree->GetEntries();

    //--- ggTree part
    if(candType == 0)
      candTree->SetBranchAddress("ZZMass",&m4l);
    else if(candType == 1)
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
    candTree->SetBranchAddress("Z1Pt", &Z1Pt);
    candTree->SetBranchAddress("Z2Pt", &Z2Pt);

    int count = 0;
    for(int k=0; k<nentries; k++){
      candTree->GetEvent(k);
//      if (PUWeight*genHEPMCweight <= 0 ) cout << "Warning! Negative weight events" << endl;

//    Find candidate ID, prefer resolved
     int typ=-1 , lep=-1, tag = -1 , candID=-1 , candID_M=-1, candID_R=-1;

// Find candidate ID , prefer resolved
     for (int j = 0; j < ZZCandType->size(); j++) {
       if ( ((ZZCandType->at(j)==1 && Z1tau21->at(j)<=0.6)||ZZCandType->at(j)==2) && fabs(ZZsel->at(j))>=100 && Z1Mass->at(j)>=70 && Z1Mass->at(j)<=105 && Z2Mass->at(j)>=60){
         if (ZZCandType->at(j)==1) candID_M=j; //merged, SR
         else if (ZZCandType->at(j)==2) candID_R=j;  //resolved, SR
       }
     }

     if( (candID_M==-1) && (candID_R==-1)) {candID=-1 ; typ=-1;}
     else if ( (candID_M==-1) && (candID_R!=-1)) {candID=candID_R ; typ=1;}
     else if ( (candID_M!=-1) && (candID_R==-1)) {candID=candID_M ; typ=0;}
     else if ( (candID_M!=-1) && (candID_R!=-1)){
       if (Z1Pt->at(candID_M)>300 && Z2Pt->at(candID_M)>200) {typ=0; candID=candID_M;}
       else {typ=1; candID=candID_R;}
     }

     float t12weight = 1.;
     if (typ == 0) {
       for (int itau = 0; itau < 24; itau++) {
         if (Z1tau21->at(candID) > tau21bin[itau] && Z1tau21->at(candID) < tau21bin[itau+1]) t12weight = 1.+tau21corr[itau];
       }
     }

     for (int i=0; i<maxMassBin; i++) {
//       if (typ==candType && x.getVal()>xMin[i] && x.getVal()<xMax[i] && ((massBin[i]<1000&&genM>(massBin[i]-15)&&genM<(massBin[i]+15))||(massBin[i]>=1000&&genM>(massBin[i]*0.90)&&genM<(massBin[i]*1.10)))) {
       if (typ==candType && x.getVal()>xMin[i] && x.getVal()<xMax[i] && ((massBin[i]<1000&&genM>(massBin[i]-15)&&genM<(massBin[i]+15))||(massBin[i]>=1000&&genM>(massBin[i]*0.90)&&genM<(massBin[i]*1.10))) && ((massBin[i]<1200)||((massBin[i]>=1200)&&(((m4l->at(candID))-genM)>-500)))) {
          ntupleVarSet.setCatIndex("massrc",massBin[i]);
          ntupleVarSet.setRealValue("reso",(m4l->at(candID))-genM);
          ntupleVarSet.setRealValue("myW",weight*t12weight);
          dataset->add(ntupleVarSet,weight*t12weight);
        }
      }
    }

  cout << "dataset n entries: " << dataset->sumEntries() << endl;

  for (int i=0; i<maxMassBin; i++) {
    sprintf(tempmass2,"massrc == massrc::mh%d",massBin[i]);
    dataset_sub[i]= (RooDataSet*)dataset->reduce(tempmass2);   
  }
}
