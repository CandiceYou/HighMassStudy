#include "resofit.h"
#include "style.cc"
void  readData(char* channel="4e")                            
 {
  for (int i=0; i<maxMassBin; i++) {
    sprintf(tempmass,"mh%d",massBin[i]);
    massrc.defineType(tempmass,massBin[i]);
  if(channel=="4e") width[i] = 5.92172+0.00639677*(massBin[i]-600)-1.06763e-06*(massBin[i]-600)*(massBin[i]-600)+4.10566e-09*(massBin[i]-600)*(massBin[i]-600)*(massBin[i]-600);
  else if(channel=="4mu") width[i] = 8.8625+0.023313*(massBin[i]-600)+0.000014677*(massBin[i]-600)*(massBin[i]-600);
  else if(channel=="2e2mu") width[i] = 7.5026+0.0156385*(massBin[i]-600)+7.8748e-06*(massBin[i]-600)*(massBin[i]-600)+2.35478e-09*(massBin[i]-600)*(massBin[i]-600)*(massBin[i]-600);
  xMin[i] = width[i]*(-15);
  xMax[i] = width[i]*(10);
  }

  RooArgSet ntupleVarSet(x,w,massrc);
  dataset = new RooDataSet("resoM","resoM",ntupleVarSet,WeightVar("myW"));

  for (int i=0; i<Nfiles; i++) {
     sprintf(inputfile,"PT13TeV/ggH_zz2l2q_M%d/ZZ2l2qAnalysis.root",inputfiles[i]);
    files.push_back(inputfile);
  }

//TChain *ggTree = new TChain("ZZTree/candTree");
  TChain *candTree = new TChain("candTree");

  for (vector<TString>::const_iterator file = files.begin(); file!=files.end(); ++file) {
   candTree->Add(inputDir+(*file));
   }

    int  nentries = candTree->GetEntries();

    //--- ggTree part
    candTree->SetBranchAddress("ZZMass",&m4l);
    candTree->SetBranchAddress("GenHMass",&genM);
    candTree->SetBranchAddress("Z1Flav",&z1flav);
    candTree->SetBranchAddress("Z2Flav",&z2flav);
    candTree->SetBranchAddress("overallEventWeight",&weight);
    candTree->SetBranchAddress("PoleMass",&pm); // not used
    candTree->SetBranchAddress("ZZsel",&zzsel);
    candTree->SetBranchAddress("Z1Mass",&Z1Mass);
    candTree->SetBranchAddress("Z2Mass",&Z2Mass);

    for(int k=0; k<nentries; k++){
      candTree->GetEvent(k);

      if(channel=="4mu" && z1flav*z2flav != 28561) continue;
      if(channel=="4e" && z1flav*z2flav != 14641) continue;
      if(channel=="2e2mu" && z1flav*z2flav != 20449) continue;
      if (weight <= 0 ) cout << "Warning! Negative weight events" << endl;

     for (int i=0; i<maxMassBin; i++) {
      ntupleVarSet.setCatIndex("massrc",massBin[i]);
      ntupleVarSet.setRealValue("reso",m4l-genM);
      ntupleVarSet.setRealValue("myW",weight);
//      ntupleVarSet.setRealValue("myW",weight*(wt->at(7)));
for(auto && zz : zzsel)
      if((zz>=100 && x.getVal()>xMin[i] && x.getVal()<xMax[i])&&((massBin[i]<1000&&genM>(massBin[i]-5)&&genM<(massBin[i]+5))||(massBin[i]>=1000&&genM>(massBin[i]*0.95)&&genM<(massBin[i]*1.05))))

       dataset->add(ntupleVarSet, weight);
//       dataset.add(ntupleVarSet, weight*(wt->at(7)));

    }
  }

  cout << "dataset n entries: " << dataset->sumEntries() << endl;

  for (int i=0; i<maxMassBin; i++) {
    sprintf(tempmass2,"massrc == massrc::mh%d",massBin[i]);
    dataset_sub[i]= (RooDataSet*)dataset->reduce(tempmass2);
  }
}
