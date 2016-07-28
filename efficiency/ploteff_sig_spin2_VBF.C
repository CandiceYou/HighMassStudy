#include "TText.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveStats.h"

//spin : 0 spin0, 2 spin2.
//ch : 0 4mu, 1 4e, 2 2e2mu.

#define TAG_MACRO(j) {  float D_2jet = 1./(1.+0.06*(phjj_VAJHU_highestPTJets->at(j)/pvbf_VAJHU_highestPTJets->at(j)));     int tag=0;      if(nExtraJets >=2 && D_2jet > 0.5) tag = 1;     else {       if (local_ZZCandType==1) {          if (btag1stSubjet > 0.46 && btag2ndSubjet > 0.46)           tag = 2;       }       else {          if(btag1stJet > 0.46 && btag2ndJet > 0.46)           tag = 2;       }     }     if(tagged != tag) break;     hreco->Fill(GenHMass,(genHEPMCweight*PUWeight)); }

short local_ZZCandType; //1 for merged jet (J), 2 for two resolved jets (jj)
int exclude; //0 exclude nothing; 1 only consider ambiguous; 2 exclude ambiguous
int channel; //121 for eeqq, 169 for mumuqq 
int tagged; //0 for untagged, 1 for vbf-tagged, 2 for b-tagged

char * polytype="polyFunctot";
TString resolve;

TGraphErrors* makegr(int spin=0, int ch=0, Color_t color=2, int marker=20, int line=1){
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetOptStat(0000);
  gStyle->SetTitleFontSize(0.05);

  //get the right directory automatically
  TFile* f;
  char dest[PATH_MAX];
  sprintf(dest, "%s", gSystem->pwd());
  char * pchar = strstr(dest, "prod");
  strcpy(pchar, "prod/\0");
  TString inputDir = dest;
 
  int inputfiles_ggH[]={200,250,400,450,500,550,600,700,750,800,900,1000};
  
  int Nfiles_ggH=sizeof(inputfiles_ggH)/sizeof(*inputfiles_ggH);
  char inputfile[PATH_MAX];
  vector<TString> files_ggH;

   for (int i=0; i<Nfiles_ggH; i++) {
   sprintf(inputfile,"2016_VBF_2l2q/mytree/PT13TeV/VBFHiggs%d/ZZ2l2qAnalysis.root",inputfiles_ggH[i]);
   files_ggH.push_back(inputfile);
   }

  TChain *candTree = new TChain("ZZTree/candTree");
  for (vector<TString>::const_iterator file = files_ggH.begin(); file!=files_ggH.end(); ++file) 
  candTree->Add(inputDir+(*file));

  // draw raw histogram with 1GeV binning, 0GeV to 3050GeV.
  // const Int_t m=3050;
  const Int_t m=4500;
  int zzflav=28561;
  TH1F *hgen = new TH1F("hgen","hgen",m,0,m); //histogram
  TH1F *hreco = new TH1F("hreco","hreco",m,0,m);

  switch (ch){
  case 0: zzflav=28561;  break; //4mu
  case 1: zzflav=14641;  break; //4e
  case 2: zzflav=20449;  break; //2e2mu
  default:  cout<<"channel unknown\t"<<ch<<endl;
  }

  vector<short> *ZZsel=0,*ZZCandType=0;
  float GenHMass=0;
  vector<float> *ZZMass=0;
  short genFinalState;
  float PUWeight,genHEPMCweight;

  vector<Float_t> *JetPt = 0;
  vector<bool> *JetIsInZZCand = 0;
  vector<Float_t> *JetQGLikelihood = 0;
  vector<Float_t> *pvbf_VAJHU_highestPTJets = 0;
  vector<Float_t> *phjj_VAJHU_highestPTJets = 0;
  vector<Float_t> *JetBTagger = 0;
  vector<Float_t> *Z1Mass = 0;
  vector<Float_t> *Z2Mass = 0;
  vector<Float_t> *Z1tau21 = 0;
  vector<Short_t> *Z2Flav = 0;

  candTree->SetBranchAddress("ZZCandType",&ZZCandType);
  candTree->SetBranchAddress("ZZsel",&ZZsel);
  candTree->SetBranchAddress("genFinalState",&genFinalState);
  candTree->SetBranchAddress("PUWeight",&PUWeight);
  candTree->SetBranchAddress("genHEPMCweight",&genHEPMCweight);

  candTree->SetBranchAddress("ZZMassRefit",&ZZMass);
  candTree->SetBranchAddress("Z1Mass",&Z1Mass);
  candTree->SetBranchAddress("Z2Mass",&Z2Mass);
  candTree->SetBranchAddress("Z2Flav",&Z2Flav);
  candTree->SetBranchAddress("GenHMass",&GenHMass);
  candTree->SetBranchAddress("pvbf_VAJHU_highestPTJets", &pvbf_VAJHU_highestPTJets );
  candTree->SetBranchAddress("phjj_VAJHU_highestPTJets", &phjj_VAJHU_highestPTJets );
  candTree->SetBranchAddress("Z1tau21", &Z1tau21);
  candTree->SetBranchAddress("JetBTagger", &JetBTagger);
  candTree->SetBranchAddress("JetPt", &JetPt);
  candTree->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);
  candTree->SetBranchAddress("JetIsInZZCand", &JetIsInZZCand);

int nInJets = 0;
int nExtraJets = 0;

for (int i=0; i<candTree->GetEntries(); i++) {
  candTree->GetEntry(i);
  if(genFinalState!=ch) continue;
  hgen->Fill(GenHMass,(genHEPMCweight*PUWeight));

  float pt1stJet = 0.0001;
  float pt2ndJet = 0.0001;
  float btag1stJet = 0.;
  float btag2ndJet = 0.;
  float qglik1stJet = 0.;
  float qglik2ndJet = 0.;

  float pt1stSubjet = 0.0001;
  float pt2ndSubjet = 0.0001;
  float btag1stSubjet = 0.;
  float btag2ndSubjet = 0.;

  int nInJets = 0;
  int nExtraJets = 0;

  for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
    
    
    if (JetQGLikelihood->at(nJet) > -800.) {    // real jets
      if (JetIsInZZCand->at(nJet) ) {         
        if (pt1stJet < JetPt->at(nJet)) {
          pt2ndJet = pt1stJet;
          pt1stJet = JetPt->at(nJet);
          btag2ndJet = btag1stJet;
          btag1stJet = JetBTagger->at(nJet);
          qglik2ndJet = qglik1stJet;
          qglik1stJet = JetQGLikelihood->at(nJet);
        } else if (pt2ndJet < JetPt->at(nJet)) {
          pt2ndJet = JetPt->at(nJet);
          btag2ndJet = JetBTagger->at(nJet);
          qglik2ndJet = JetQGLikelihood->at(nJet);
        }
        nInJets++;
      } else nExtraJets++;  
    }
  } 

  for (unsigned int nJet=0; nJet<JetPt->size(); nJet++) {
    if (JetQGLikelihood->at(nJet) <= -800.) {            // subjets
      if (JetIsInZZCand->at(nJet) ) {         
        if (pt1stSubjet < JetPt->at(nJet)) {
          pt2ndSubjet = pt1stSubjet;
          pt1stSubjet = JetPt->at(nJet);
          btag2ndSubjet = btag1stSubjet;
          btag1stSubjet = JetBTagger->at(nJet);
        } else if (pt2ndSubjet < JetPt->at(nJet)) {
          pt2ndSubjet = JetPt->at(nJet);
                btag2ndSubjet = JetBTagger->at(nJet);
        }
      }
    }
  } 

//HERE THERE BE DRAGONS

  switch(exclude) {
case 0: //consider all
  if(ZZCandType->size() == 1) { // unambiguous
    if(local_ZZCandType==1) { //merged
      if( (ZZCandType->at(0) == local_ZZCandType) && (abs(ZZsel->at(0))>=100) && (Z2Mass->at(0) >= 60) && (Z1Mass->at(0)>=70) && (Z1Mass->at(0)<=105) && (abs(Z2Flav->at(0))==channel) && (Z1tau21->at(0) <= 0.6)) 
        TAG_MACRO(0);
      }

    else { //resolved
      if( (ZZCandType->at(0) == local_ZZCandType) && (abs(ZZsel->at(0))>=100) && (Z2Mass->at(0) >= 60) && (Z1Mass->at(0)>=70) && (Z1Mass->at(0)<=105) && (abs(Z2Flav->at(0))==channel))
        TAG_MACRO(0);
    }
  }

  if(ZZCandType->size() == 2) {  // ambiguous

    if(local_ZZCandType==1) { //merged
        if( (ZZCandType->at(0) == local_ZZCandType) && (abs(ZZsel->at(0))>=100) && (Z2Mass->at(0) >= 60) && (Z1Mass->at(0)>=70) && (Z1Mass->at(0)<=105) && (abs(Z2Flav->at(0))==channel) && (Z1tau21->at(0) <= 0.6)) 
        {
        TAG_MACRO(0);
      }
      
      else {
        
        if ((ZZCandType->at(1) == local_ZZCandType) && (ZZsel->at(1)>=100) && (Z2Mass->at(1) >= 60) && (Z1Mass->at(1)>=70) && (Z1Mass->at(1)<=105) && (abs(Z2Flav->at(1))==channel) && (Z1tau21->at(0) <= 0.6))
        TAG_MACRO(1);
      }
    }
    

    else { //resolved
      
      if( (ZZCandType->at(0) == local_ZZCandType) && (abs(ZZsel->at(0))>=100) && (Z2Mass->at(0) >= 60) && (Z1Mass->at(0)>=70) && (Z1Mass->at(0)<=105) && (abs(Z2Flav->at(0))==channel)) {
        TAG_MACRO(0);
      }

      else {
        
        if ((ZZCandType->at(1) == local_ZZCandType) && (ZZsel->at(1)>=100) && (Z2Mass->at(1) >= 60) && (Z1Mass->at(1)>=70) &&(Z1Mass->at(1)<=105) && (abs(Z2Flav->at(1))==channel))
        TAG_MACRO(1);
      }
    }
  }
  
  break;

case 2: // only consider unambiguous
  if(ZZCandType->size() != 1) break; 

    if(local_ZZCandType==1) { //merged
      if( (ZZCandType->at(0) == local_ZZCandType) && (abs(ZZsel->at(0))>=100) && (Z2Mass->at(0) >= 60) && (Z1Mass->at(0)>=70) && (Z1Mass->at(0)<=105) && (abs(Z2Flav->at(0))==channel) && (Z1tau21->at(0) <= 0.6)) 
        TAG_MACRO(0);
    }

    else { //resolved
      if( (ZZCandType->at(0) == local_ZZCandType) && (abs(ZZsel->at(0))>=100) && (Z2Mass->at(0) >= 60) && (Z1Mass->at(0)>=70) && (Z1Mass->at(0)<=105) && (abs(Z2Flav->at(0))==channel))
        TAG_MACRO(0);
    }
  } //switch
}//for



// bin contents of raw histograms
 double M_raw[m]={0};
 double gen_raw[m]={0};
 double reco_raw[m]={0};

// bin contents of merged histograms
 double M[]={300,350,400,450,500,550,600,650,700,750,800,850,900,1000};
 const Int_t n = sizeof(M)/sizeof(*M);

 double massE[n]={0};
 double gen[n]={0};
 double genE[n]={0};
 double reco[n]={0};
 double recoE[n]={0};
 double eff[n]={0};
 double effE[n]={0};

 double mass,width;

 for (int bin=1;bin<=m;bin++){
   M_raw[bin-1] = hreco->GetXaxis()->GetBinCenter(bin);
   gen_raw[bin-1] = hgen->GetBinContent(bin);
   reco_raw[bin-1] = hreco->GetBinContent(bin);
 }

 cout<<"spin="<<spin<<",ch="<<ch<<endl;

 for (int i=0;i<n;i++){
  mass = M[i];
  if(mass<=180) width=1;
  else if (mass>180 && mass<=300) width =5;
  else if (mass>300&&mass<600) width = 10;
  else if (mass>=600&&mass<1000) width = 25;
  else if (mass>=1000&&mass<2000) width = 100;
  else if (mass>=2000&&mass<2300) width=150;
  else width=250;

//merge
  for (int j=(mass-width);j<(mass+width);j++){
  //cout<<"mass="<<mass<<endl;
  //cout<<"gen_bin"<<i<<"="<<gen_raw[j]<<endl;
  //cout<<"reco_bin"<<i<<"="<<reco_raw[j]<<endl;
    gen[i]+=gen_raw[j];
    reco[i]+=reco_raw[j];
  }

  genE[i]=sqrt(gen[i]);
  recoE[i]=sqrt(reco[i]);

  // calculate efficiency and efficiency error
  double g,l,r;
  g=gen[i];
  r=reco[i];
  l=(g-r);
  if(g!=0){
  eff[i]=r/g;
  effE[i]=sqrt(l*l*r+r*r*l)/(g*g);
  }
  else {
  eff[i]=0;
  effE[i]=0;
  }
  cout<<eff[i]<<",";
}

  TF1 *polyFunctot= new TF1("polyFunctot","([0]+[1]*TMath::Erf( (x-[2])/[3] ))*([4]+[5]*x+[6]*x*x)+[7]*TMath::Gaus(x,[8],[9])", 300, 1000);
  polyFunctot->SetParameters(-4.42749e+00,4.61212e+0,-6.21611e+01,1.13168e+02,2.14321e+00,1.04083e-03,4.89570e-07, 0.1, 200, 30,0);
  polyFunctot->SetParLimits(1,0,.2);
  polyFunctot->SetParLimits(2,300,1000);
  polyFunctot->SetParLimits(8,500,800);
  polyFunctot->SetParLimits(9,10,170);
  polyFunctot->SetLineColor(color);
 
  if(spin==0)   polyFunctot->SetLineStyle(1);
  else if(spin==2) {
  polyFunctot->SetParLimits(7,0,0);
  polyFunctot->SetParLimits(8,0,0);
  polyFunctot->SetParLimits(9,0,0); 
  polyFunctot->SetLineStyle(2);
  }

  TGraphErrors *gr = new TGraphErrors (n,M,eff,massE,effE);

  gr->SetMarkerColor(color);
  gr->SetMarkerStyle(marker);
  gr->SetMarkerSize(1);
  cout<<"\n\nparameters:\n";
  gr->Fit(polytype);
  gr->SetLineColor(color);
  
  TString param_txt = resolve;
  param_txt += ".txt";

  FILE *fp = fopen(param_txt.Data(),"w");
  if (fp!=NULL && polytype =="polyFunctot" ) {
    for (int i=0;i<polyFunctot->GetNpar();i++) {
        Float_t value = polyFunctot->GetParameter(i);
        fprintf(fp,"p%d\t%f\n",i,value);
     }
  }
  fclose(fp);

  delete hreco; delete hgen; delete polyFunctot;
  return gr;
}


void ploteff_sig_spin2_80X_2(){
  TCanvas* c2 = new TCanvas("c2", "c2", 1000, 10, 1400, 800);
  //c2->SetLogx(); 
  c2->SetFillColor(0);
  c2->SetRightMargin(0.13);
  TMultiGraph *mg = new TMultiGraph();

  //auto resolve directory name
  char dest[PATH_MAX];  
  sprintf(dest, "%s", gSystem->pwd());

  char * pchar = strstr(dest, "/w/wqin");
  if(pchar != 0)
  strcpy(pchar, "/w/wqin/www/\0");

  pchar = strstr(dest, "/r/rbarr");
  if(pchar != 0)
  strcpy(pchar, "/r/rbarr/www/2l2q_VBF/June30/\0");

  pchar = strstr(dest, "/c/cayou");
  if(pchar != 0)
  strcpy(pchar, "/c/cayou/www/HighMass/\0");


  const char * string_ZZCandType = (local_ZZCandType==2)?"jj":"J";
  const char * string_channel = (channel==121)?"eeqq":"mumuqq";
  const char * string_exclude="";
    if(exclude==2) string_exclude = "reject_ambiguous";
    else if(exclude==1) string_exclude = "only_ambiguous";
    else if(exclude==0) string_exclude = "all";
  const char * string_tagged="";
    if(tagged==2) string_tagged = "b-tagged";
    else if(tagged==1) string_tagged = "vbf-tagged";
    else if(tagged==0) string_tagged = "untagged";

  sprintf(dest+strlen(dest), "vbf_efficiency_%s_%s_%s_%s", string_ZZCandType, string_channel, string_tagged, string_exclude);
  resolve = TString(dest);


  //TGraphErrors* ggH_4mu = makegr(0,0,kRed,20,1);
  //TGraphErrors* ggH_4e = makegr(0,1,kGreen,20,1);
  //TGraphErrors* ggH_2e2mu = makegr(0,2,kBlue,20,1);
  //TGraphErrors* spin2_4mu = makegr(2,0,kRed+2,24,2);
  //TGraphErrors* spin2_4e = makegr(2,1,kGreen+2,24,2);
  //TGraphErrors* spin2_2e2mu = makegr(2,2,kBlue+2,24,2);
  TGraphErrors* ggH_2l2q = makegr(0,99,kViolet,24,2);

  //mg->Add(ggH_4e);
  //mg->Add(ggH_4mu);
  //mg->Add(ggH_2e2mu);
  mg->Add(ggH_2l2q);
  //mg->Add(spin2_4e);
  //mg->Add(spin2_4mu);
  //mg->Add(spin2_2e2mu);
  mg->Draw("AP");
  mg->GetXaxis()->SetTitle("genHMass [GeV]");
  mg->GetYaxis()->SetTitle("efficiency*acceptance");
  mg->GetYaxis()->SetRangeUser(0,1);
  mg->GetXaxis()->SetRangeUser(100,4100);
  mg->GetYaxis()->SetRangeUser(0.,0.35);

  TLegend* leg3 = new TLegend(.9,0.3,0.99,.85);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  //leg3->AddEntry(ggH_4mu,"ggH->2b+ 4mu","pl");
  //leg3->AddEntry(ggH_2e2mu,"ggH->2b+ 2e2mu","pl");
  //leg3->AddEntry(ggH_4e,"ggH->2b+ 4e","pl");
  leg3->AddEntry(ggH_2l2q,"ggh->2b+ 2l2q","pl");

  //leg3->AddEntry(ggH_4mu,"ggH->0+ 4mu","pl");
  //leg3->AddEntry(ggH_2e2mu,"ggH->0+ 2e2mu","pl");
  //leg3->AddEntry(ggH_4e,"ggH->0+ 4e","pl");  
  //leg3->AddEntry(spin2_4mu,"ggH->2b+ 4mu","pl");
  //leg3->AddEntry(spin2_2e2mu,"ggH->2b+ 2e2mu","pl");
  //leg3->AddEntry(spin2_4e,"ggH->2b+ 4e","pl");
  //leg3->Draw();


  c2->Update();
  c2->SaveAs(resolve + ".png");
  c2->SaveAs(resolve + ".pdf");
}

void ploteff_sig_spin2_VBF() {
  local_ZZCandType=2;
  channel=169;
  tagged=1;
  exclude=(local_ZZCandType==2)? 0:2;
  ploteff_sig_spin2_80X_2(); 

  /*
  for(int c:{1,2}) {           //1 for merged jet (J), 2 for two resolved jets (jj)
    local_ZZCandType=c;
    exclude=(local_ZZCandType==2)? 0:2;   //give preference to resolved jets
  for (int ch:{121,169}) {               //121 for eeqq, 169 for mumuqq 
    channel=ch;
  for (int t:{0,1,2}) {                //0 for untagged, 1 for vbf-tagged, 2 for b-tagged
    tagged=t;
    printf("%d\t%d\t%d\t%d\n", local_ZZCandType, exclude, channel, tagged);
    ploteff_sig_spin2_80X_2();
  }}}
  */
}
