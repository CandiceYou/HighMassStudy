#include "TF1.h"
#include "TText.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveStats.h"
using namespace std;


void test(char* var1="Gen_costheta1", char* var2="helcosthetaZ1", char* title="test", int nbin=30, double min=-1, double max=1){
 gStyle->SetPadLeftMargin(0.1);
 gStyle->SetPadRightMargin(0.15);
 gStyle->SetOptFit(0000);
 gStyle->SetOptStat(0000);
 gStyle->SetTitleFontSize(0.05);

 TChain *candTree1 = new TChain("candTree");
//choose files with polemass>genmass
 candTree1->Add("../data/sample_rw2bp/ggH750_ZZ4lAnalysis_addwt.root"); //no width

 TChain *candTree2 = new TChain("SelectedTree");
 candTree2->Add("../data/ggH_2bp_zz4l_M750_Ga247_all.root");
// TChain *candTree2 = new TChain("newTree");
// candTree2->Add("../data/mH750_2bplus_gg_2e2mu.root");

 TH1F *hist1 = new TH1F("hist1","hist1",nbin,min,max);
 TH1F *hist2 = new TH1F("hist2","hist2",nbin,min,max);
 TH1F *hist3 = new TH1F("hist3","hist3",nbin,min,max);

 char hist1draw[500],hist2draw[500],hist3draw[500];
 sprintf(hist1draw,"%s>>hist1",var1);
 sprintf(hist2draw,"%s>>hist2",var2);
 sprintf(hist3draw,"%s>>hist3",var1);


 candTree1->Draw(hist1draw,"(GenHMass<=800&&GenHMass>=700&& genFinalState==2)*wt_2bp");
 candTree2->Draw(hist2draw,"ZZMass<=800&&ZZMass>=700&&interf==0");
 candTree1->Draw(hist3draw,"GenHMass<=800&&GenHMass>=700&& genFinalState==2");

 TCanvas* c = new TCanvas("c", "c", 800,800);
 c->SetFillColor(0);
 c->cd();

 hist1->SetLineWidth(2);
 hist2->SetLineWidth(2);
 hist3->SetLineWidth(2);

 hist1->SetLineColor(kRed);
 hist2->SetLineColor(kBlack);
 hist3->SetLineColor(kGreen+2);

 hist1->Scale(1/hist1->Integral("wt_2bp"));
 hist2->Scale(1/hist2->Integral());
 hist3->Scale(1/hist3->Integral());
// hist2->Scale(1/hist2->GetEntries());

 hist2->SetTitle(title);
 hist2->GetXaxis()->SetTitle(var2);
 if(var2!="Z1Mass"&&var2!="Z2Mass")
 hist2->GetYaxis()->SetRangeUser(0,0.04);
// hist2->GetYaxis()->SetTitle("Events");

 hist2->Draw("hist");
 hist1->Draw("samehist");
 hist3->Draw("samehist");

 c->Update();
 char temp[150];
 sprintf(temp,"%s.png",var2);
 c->SaveAs(temp);
 delete c;
 delete hist1;
 delete hist2;
 delete hist3;
}

void plotAngles(){
test("Gen_costheta1","helcosthetaZ1","",50,1, 1);
test("Gen_costheta2","helcosthetaZ2","",50,1, 1);
test("Gen_phi","helphi","",50,-3.14, 3.14);
test("Gen_costhetastar","costhetastar","",50,-1,1);
test("Gen_phistar1","phistarZ1","",50,-3.14,3.14);
test("Gen_mZ1","Z1Mass","",80,50,130);
test("Gen_mZ2","Z2Mass","",80,50,130);
}
