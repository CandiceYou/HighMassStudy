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

  if (ZZCandType == 1)
    sprintf(cType,"J");
  else if (ZZCandType ==2)
    sprintf(cType,"jj");

  RooArgSet ntupleVarSet(x,w,massrc);
  dataset = new RooDataSet("resoM","resoM",ntupleVarSet,WeightVar("myW"));

  for (int i=0; i<Nfiles; i++) {
     sprintf(inputfile,"2016_2l2q/mytree/PT13TeV/ggH_zz2l2q_M%d/ZZ2l2qAnalysis.root",inputfiles[i]);
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
    candTree->SetBranchAddress("ZZsel",&zzsel);
    candTree->SetBranchAddress("Z1Mass",&Z1Mass);
    candTree->SetBranchAddress("Z2Mass",&Z2Mass);
    candTree->SetBranchAddress("ZZCandType",&candType);

    for(int k=0; k<nentries; k++){
      candTree->GetEvent(k);
//      if (PUWeight*genHEPMCweight <= 0 ) cout << "Warning! Negative weight events" << endl;

      switch (exclude){
      case 0:
        for (int w=0; w < (*candType).size(); w++) {
         if (candType->at(w)==ZZCandType) {
           for (int i=0; i<maxMassBin; i++) {
             ntupleVarSet.setCatIndex("massrc",massBin[i]);
             ntupleVarSet.setRealValue("reso",(m4l->at(w))-genM);
             ntupleVarSet.setRealValue("myW",PUWeight*genHEPMCweight);

             if((ZZCandType==1) && ((zzsel->at(w))>=100 && x.getVal()>xMin[i] && x.getVal()<xMax[i])&&((massBin[i]<1000&&genM>(massBin[i]-5)&&genM<(massBin[i]+5))||(massBin[i]>=1000&&genM>(massBin[i]*0.95)&&genM<(massBin[i]*1.05))))
               dataset->add(ntupleVarSet, PUWeight*genHEPMCweight);
             else if((ZZCandType==2) && (((Z1Mass->at(w))>70)&&((Z1Mass->at(w))<105)) && ((zzsel->at(w))>=100 && x.getVal()>xMin[i] && x.getVal()<xMax[i])&&((massBin[i]<1000&&genM>(massBin[i]-5)&&genM<(massBin[i]+5))||(massBin[i]>=1000&&genM>(massBin[i]*0.95)&&genM<(massBin[i]*1.05))) && ((massBin[i]<1200)||((massBin[i]>=1200)&&(((m4l->at(w))-genM)>-500))))
               dataset->add(ntupleVarSet, PUWeight*genHEPMCweight);
           }
         }
        }
        break;
      case 1:
        for (int w=0; w < (*candType).size(); w++) {
         if ((candType->at(w)==ZZCandType) && ((*candType).size()==2)) {
           for (int i=0; i<maxMassBin; i++) {
             ntupleVarSet.setCatIndex("massrc",massBin[i]);
             ntupleVarSet.setRealValue("reso",(m4l->at(w))-genM);
             ntupleVarSet.setRealValue("myW",PUWeight*genHEPMCweight);

             if((ZZCandType==1) && ((zzsel->at(w))>=100 && x.getVal()>xMin[i] && x.getVal()<xMax[i])&&((massBin[i]<1000&&genM>(massBin[i]-5)&&genM<(massBin[i]+5))||(massBin[i]>=1000&&genM>(massBin[i]*0.95)&&genM<(massBin[i]*1.05))))
               dataset->add(ntupleVarSet, PUWeight*genHEPMCweight);
             else if((ZZCandType==2) && (((Z1Mass->at(w))>70)&&((Z1Mass->at(w))<105)) && ((zzsel->at(w))>=100 && x.getVal()>xMin[i] && x.getVal()<xMax[i])&&((massBin[i]<1000&&genM>(massBin[i]-5)&&genM<(massBin[i]+5))||(massBin[i]>=1000&&genM>(massBin[i]*0.95)&&genM<(massBin[i]*1.05))) && ((massBin[i]<1200)||((massBin[i]>=1200)&&(((m4l->at(w))-genM)>-500))))
               dataset->add(ntupleVarSet, PUWeight*genHEPMCweight);
           }
         }
        }
        break;
      case 2:
       if (((*candType).size()==1) && (candType->at(0)==ZZCandType)) {
         for (int i=0; i<maxMassBin; i++) {
           ntupleVarSet.setCatIndex("massrc",massBin[i]);
           ntupleVarSet.setRealValue("reso",(m4l->at(0))-genM);
           ntupleVarSet.setRealValue("myW",PUWeight*genHEPMCweight);

           if((ZZCandType==1) && ((zzsel->at(0))>=100 && x.getVal()>xMin[i] && x.getVal()<xMax[i])&&((massBin[i]<1000&&genM>(massBin[i]-5)&&genM<(massBin[i]+5))||(massBin[i]>=1000&&genM>(massBin[i]*0.95)&&genM<(massBin[i]*1.05))))
             dataset->add(ntupleVarSet, PUWeight*genHEPMCweight);
           else if((ZZCandType==2) && (((Z1Mass->at(0))>70)&&((Z1Mass->at(0))<105)) && ((zzsel->at(0))>=100 && x.getVal()>xMin[i] && x.getVal()<xMax[i])&&((massBin[i]<1000&&genM>(massBin[i]-5)&&genM<(massBin[i]+5))||(massBin[i]>=1000&&genM>(massBin[i]*0.95)&&genM<(massBin[i]*1.05))) && ((massBin[i]<1200)||((massBin[i]>=1200)&&(((m4l->at(0))-genM)>-500))))
             dataset->add(ntupleVarSet, PUWeight*genHEPMCweight);
          }
        }
        break;
      default:
        cout << endl << "Bad value for \"exclude\"." << endl;
        break;
      }
    } 

  cout << "dataset n entries: " << dataset->sumEntries() << endl;

  for (int i=0; i<maxMassBin; i++) {
    sprintf(tempmass2,"massrc == massrc::mh%d",massBin[i]);
    dataset_sub[i]= (RooDataSet*)dataset->reduce(tempmass2);
  }
}
