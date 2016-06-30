#include "TMath.h"
#include "TF1.h"
#include "TText.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "RooTFnBinding.h"
#include "fitFunction.c"
using namespace std;

TString resolve = "/afs/cern.ch/user/r/rbarr/www/bkg/bkg_EfficiencyFit";

void ploteff_bkg(){
 gStyle->SetPadLeftMargin(0.1);
 gStyle->SetPadRightMargin(0.15);
 gStyle->SetOptFit(0000);
 gStyle->SetOptStat(0000);
 gStyle->SetTitleFontSize(0.05);


TF1 *erf1 = new TF1("erf1", "0.5*[2]*(1+TMath::Erf( (x-[0]) / ([1]*sqrt(2)) ) )*([3]*x*x+[4]*x+[5])",100.0,1000.0);
erf1->SetParameters(118.90954335267955,27.544135827573065,9.960803718511713e-07,-0.25378550401112004,492.32846532876965,278086.23395067605);

TF1 *polyFunctot1 = new TF1("polyFunctot1","[0]+[1]*TMath::Sin((x-[2])/[3])*TMath::Gaus(x,[4],[5])+TMath::Erf((x-[6])/[7])*([8])", 100, 2000);

polyFunctot1->SetParameters(0,.1);
polyFunctot1->SetParLimits(0,0,0.25);
polyFunctot1->SetParameters(1,.04);
polyFunctot1->SetParLimits(1,0,0.05);
polyFunctot1->SetParameters(3,65);
polyFunctot1->SetParLimits(3,50,100);
polyFunctot1->SetParameters(4,800);
polyFunctot1->SetParLimits(4,600,1000);
polyFunctot1->SetParameters(5,200);
polyFunctot1->SetParLimits(5,100,300);
polyFunctot1->SetParameters(6,1);
polyFunctot1->SetParLimits(6,0,0.2);
polyFunctot1->SetParameters(7,50);
polyFunctot1->SetParLimits(7,0,200);
polyFunctot1->SetParameters(8,.1);
polyFunctot1->SetParLimits(8,0,.3);


TF1 *polyFunctot2 = new TF1("polyFunctot2","[0]+[1]*TMath::Sin((x-[2])/[3])*TMath::Gaus(x,[4],[5])+TMath::Erf((x-[6])/[7])*([8])", 100, 1000);

polyFunctot2->SetParameters(0,.1);
polyFunctot2->SetParLimits(0,0,0.25);
polyFunctot2->SetParameters(1,.04);
polyFunctot2->SetParLimits(1,0,0.07);
polyFunctot2->SetParameters(3,65);
polyFunctot2->SetParLimits(3,50,100);
polyFunctot2->SetParameters(4,800);
polyFunctot2->SetParLimits(4,500,1000);
polyFunctot2->SetParameters(5,200);
polyFunctot2->SetParLimits(5,100,300);
polyFunctot2->SetParameters(6,1);
polyFunctot2->SetParLimits(6,0,0.2);
polyFunctot2->SetParameters(7,50);
polyFunctot2->SetParLimits(7,0,200);
polyFunctot2->SetParameters(8,.1);
polyFunctot2->SetParLimits(8,0,.3);

//  TF1 *polyFunctot3= new TF1("polyFunctot3","TMath::Erf( (x-[0])/[1])*([2]+[3]*x+[4]*x*x)+[5]+[6]*x+[7]*x*x+[8]*x*x*x", 90., 2000);
//  polyFunctot3->SetParameters(1.79910e+02,1.38088e-01,3.49180e-01,-2.36380e-04,-5.12049e-06,5.39771e-06,-4.73109e-11,2.75448e-01,0);

  polyFunctot1->SetLineColor(3);
  polyFunctot2->SetLineColor(2);
//  polyFunctot3->SetLineColor(4);
 
 TFile* f1 = TFile::Open("root://lxcms03//data3/Higgs/160624/ggZZ4e/ZZ4lAnalysis.root"); 
 TFile* f2 = TFile::Open("root://lxcms03//data3/Higgs/160624/ggZZ4mu/ZZ4lAnalysis.root");
// TFile* f3 = TFile::Open("root://lxcms03//data3/Higgs/160624/ggZZ2e2mu/ZZ4lAnalysis.root"); 

 TTree * candTree1 = (TTree*) f1->Get("ZZTree/candTree");
 TTree * candTree2 = (TTree*) f2->Get("ZZTree/candTree");
// TTree * candTree3 = (TTree*) f3->Get("ZZTree/candTree");
 TH1F * hgen_4e = (TH1F*)f1->Get("PlotsZZ/hGenZZMass_4e");
 TH1F * hgen_4mu = (TH1F*)f2->Get("PlotsZZ/hGenZZMass_4mu");
// TH1F * hgen_2e2mu = (TH1F*)f3->Get("PlotsZZ/hGenZZMass_2e2mu");


 const Int_t m=2050;
 hgen_4e->SetBins(m,0,m);
 hgen_4mu->SetBins(m,0,m);
// hgen_2e2mu->SetBins(m,0,m);

 TH1F *hreco_4e = new TH1F("hreco_4e","hreco_4e",m,0,m);
 TH1F *hreco_4mu = new TH1F("hreco_4mu","hreco_4mu",m,0,m);
// TH1F *hreco_2e2mu = new TH1F("hreco_2e2mu","hreco_2e2mu",m,0,m);

 candTree1->Draw("GenHMass>>hreco_4e","(Z1Flav*Z2Flav == 14641)*(genHEPMCweight*PUWeight*dataMCWeight)");
 candTree2->Draw("GenHMass>>hreco_4mu","(Z1Flav*Z2Flav == 28561)*(genHEPMCweight*PUWeight*dataMCWeight)");
// candTree3->Draw("GenHMass>>hreco_2e2mu","(Z1Flav*Z2Flav == 20449)*(genHEPMCweight*PUWeight*dataMCWeight)");

cout << "Loaded root files\n";

double M[m]={0};
double gen_4e[m]={0};
double reco_4e[m]={0};
double gen_4mu[m]={0};
double reco_4mu[m]={0};
//double gen_2e2mu[m]={0};
//double reco_2e2mu[m]={0};

const Int_t n=65;
double M_new[n]={0};
for (int k=0;k<46;k++){
M_new[k]=(k+5)*20;
//cout<<"m["<<k<<"]="<<M_new[k]<<endl;
}
for (int q=46;q<n;q++){
M_new[q]=1000+(q-45)*50;
//cout<<"m["<<q<<"]="<<M_new[q]<<endl;
}

double massE[n]={0};

double gen_new_4e[n]={0};
double genE_new_4e[n]={0};
double reco_new_4e[n]={0};

double recoE_new_4e[n]={0};
double eff_4e[n]={0};
double effE_4e[n]={0};

double gen_new_4mu[n]={0};
double genE_new_4mu[n]={0};
double reco_new_4mu[n]={0};
double recoE_new_4mu[n]={0};
double eff_4mu[n]={0};
double effE_4mu[n]={0};

//double gen_new_2e2mu[n]={0};
//double genE_new_2e2mu[n]={0};
//double reco_new_2e2mu[n]={0};
//double recoE_new_2e2mu[n]={0};
//double eff_2e2mu[n]={0};
//double effE_2e2mu[n]={0};

double mass,width;

 for (int bin=1;bin<=m;bin++){
 M[bin-1] = hreco_4e->GetXaxis()->GetBinCenter(bin);
 gen_4e[bin-1] = hgen_4e->GetBinContent(bin);
 reco_4e[bin-1] = hreco_4e->GetBinContent(bin);

 gen_4mu[bin-1] = hgen_4mu->GetBinContent(bin);
 reco_4mu[bin-1] = hreco_4mu->GetBinContent(bin);

// gen_2e2mu[bin-1] = hgen_2e2mu->GetBinContent(bin);
// reco_2e2mu[bin-1] = hreco_2e2mu->GetBinContent(bin);
 }

 for (int i=0;i<n;i++){
  mass = M_new[i];
  if(mass<1000) width=10;
  else width=25;
  for (int j=(mass-width);j<(mass+width);j++){
       gen_new_4e[i]+=gen_4e[j];
       reco_new_4e[i]+=reco_4e[j];
       gen_new_4mu[i]+=gen_4mu[j];
       reco_new_4mu[i]+=reco_4mu[j];
//       gen_new_2e2mu[i]+=gen_2e2mu[j];
//       reco_new_2e2mu[i]+=reco_2e2mu[j];
      }
 genE_new_4e[i]=sqrt(gen_new_4e[i]);
 recoE_new_4e[i]=sqrt(reco_new_4e[i]);
 genE_new_4mu[i]=sqrt(gen_new_4mu[i]);
 recoE_new_4mu[i]=sqrt(reco_new_4mu[i]);
// genE_new_2e2mu[i]=sqrt(gen_new_2e2mu[i]);
// recoE_new_2e2mu[i]=sqrt(reco_new_2e2mu[i]);

 double g1,l1,r1,g2,l2,r2,g3,l3,r3;
 g1=gen_new_4e[i];
 r1=reco_new_4e[i];
 l1=(g1-r1);

 g2=gen_new_4mu[i];
 r2=reco_new_4mu[i];
 l2=(g2-r2);

// g3=gen_new_2e2mu[i];
// r3=reco_new_2e2mu[i];
 //l3=(g3-r3);

 
 if(g1!=0){
 eff_4e[i]=r1/g1;
 effE_4e[i]=sqrt(l1*l1*r1+r1*r1*l1)/(g1*g1);
 }
 else {
 eff_4e[i]=0;
 effE_4e[i]=0;
 }

 if(g2!=0){
 eff_4mu[i]=r2/g2;
 effE_4mu[i]=sqrt(l2*l2*r2+r2*r2*l2)/(g2*g2);
 }
 else {
 eff_4mu[i]=0;
 effE_4mu[i]=0;
 }

 if(g3!=0){
// eff_2e2mu[i]=r3/g3;
// effE_2e2mu[i]=sqrt(l3*l3*r3+r3*r3*l3)/(g3*g3);
 }
 else {
// eff_2e2mu[i]=0;
// effE_2e2mu[i]=0;
 }
 
 }

  TGraphErrors *gr31 = new TGraphErrors (n,M_new,eff_4e,massE,effE_4e);
  gr31->SetMarkerColor(3);
  gr31->GetXaxis()->SetTitle("genHMass");
  gr31->SetMarkerStyle(20);
  gr31->SetLineColor(3);
  gr31->SetMarkerSize(0.5);
  gr31->Fit("polyFunctot1");

  TGraphErrors *gr32 = new TGraphErrors (n,M_new,eff_4mu,massE,effE_4mu);
  gr32->SetMarkerColor(2);
  gr32->GetXaxis()->SetTitle("genHMass");
  gr32->SetLineColor(2);
  gr32->SetMarkerStyle(20);
  gr32->SetMarkerSize(0.5);
  gr32->Fit("polyFunctot2");

//  TGraphErrors *gr33 = new TGraphErrors (n,M_new,eff_2e2mu,massE,effE_2e2mu);
//  gr33->SetMarkerColor(4);
//  gr33->GetXaxis()->SetTitle("genHMass");
//  gr33->SetLineColor(4);
//  gr33->SetMarkerStyle(20);
//  gr33->SetMarkerSize(0.5);
//  gr33->Fit("polyFunctot3");


/*
TFile f ("bkg_eff.root","recreate");
f.cd();
gr31->SetName("bkgeff_4e"); gr31->Write();
gr32->SetName("bkgeff_4mu"); gr32->Write();
gr33->SetName("bkgeff_2e2mu"); gr33->Write();
f.Close();
*/


  TCanvas* c2 = new TCanvas("c2", "c2", 1000, 10, 1400, 800);
  c2->SetFillColor(0);

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr31);
  mg->Add(gr32);
//  mg->Add(gr33);
  //mg->SetTitle("backgroud efficiency*acceptance");
  mg->Draw("AP");
  mg->GetXaxis()->SetTitle("genHMass");
  mg->GetYaxis()->SetTitle("efficiency*acceptance");
  polyFunctot1->Draw("same");
  polyFunctot2->Draw("same");
//  polyFunctot3->Draw("same");
  TLegend* leg3 = new TLegend(.87,0.5,0.97,.85);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->AddEntry(gr31,"4e","p");
  leg3->AddEntry(gr32,"4mu","p");
//  leg3->AddEntry(gr33,"2e2mu","p");
  leg3->Draw();

  c2->Update();
  c2->SaveAs(resolve+".png");
  c2->SaveAs(resolve+".pdf");

  TString param_txt = resolve+"_4e.txt";

  FILE *fp = fopen(param_txt.Data(),"w");
  if (fp!=NULL) {
    for (int i=0;i<polyFunctot1->GetNpar();i++) {
        Float_t value = polyFunctot1->GetParameter(i);
        fprintf(fp,"p%d\t%f\n",i,value);
     }
  }
  fclose(fp);

  param_txt = resolve+"_4mu.txt";

  fp = fopen(param_txt.Data(),"w");
  if (fp!=NULL) {
    for (int i=0;i<polyFunctot2->GetNpar();i++) {
        Float_t value = polyFunctot2->GetParameter(i);
        fprintf(fp,"p%d\t%f\n",i,value);
     }
  }
  fclose(fp);


/*
 TH1F heff_4e = (*hreco_4e)/(*hgen_4e);
 TH1F heff_4mu = (*hreco_4mu)/(*hgen_4mu);
 TH1F heff_2e2mu = (*hreco_2e2mu)/(*hgen_2e2mu);
 //heff_4e.Rebin(20);
 //heff_4mu.Rebin(20);
 //heff_2e2mu.Rebin(20);
 TCanvas* c = new TCanvas("c", "c", 1000, 10, 800, 1000);
 TPad *pad1 = new TPad("pad1","gen",0.05,0.66,0.95,0.97);
 pad1->Draw();
 TPad *pad2 = new TPad("pad2","reco",0.05,0.33,0.95,0.66);
 pad2->Draw();
 TPad *pad3 = new TPad("pad3","eff",0.05,0.02,0.95,0.33);
 pad3->Draw();
 c->SetFillColor(0);
 //c->cd();
 pad1->cd();
 hgen_4mu->SetLineColor(kRed);
 hgen_4e->SetLineColor(kGreen);
 hgen_2e2mu->SetLineColor(kBlue);

 hgen_2e2mu->SetTitle("GenHMass (gen-level)");
 hgen_2e2mu->Draw();
 hgen_4mu->Draw("same");
 hgen_4e->Draw("same");

 pad2->cd();
 hreco_4mu->SetLineColor(kRed);
 hreco_4e->SetLineColor(kGreen);
 hreco_2e2mu->SetLineColor(kBlue);

 hreco_2e2mu->SetTitle("GenHMass (fullsim)");
 hreco_2e2mu->Draw("hist");
 hreco_4e->Draw("samehist");
 hreco_4mu->Draw("samehist");

 pad3->cd();
 heff_4mu.SetLineColor(kRed);
 heff_4e.SetLineColor(kGreen);
 heff_2e2mu.SetLineColor(kBlue);

 heff_4mu.SetTitle("efficiency*acceptance");
 heff_4mu.Draw("hist");
 heff_4e.Draw("samehist");
 heff_2e2mu.Draw("samehist");

 c->Update();
 char temp[150];
 sprintf(temp,"/afs/cern.ch/user/c/cayou/www/HighMass/160217/bkg_orignal.png");
 c->SaveAs(temp);
*/
}
