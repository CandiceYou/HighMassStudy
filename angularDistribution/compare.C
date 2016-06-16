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
 candTree1->Add("test_m750_w0.root");

 TChain *candTree2 = new TChain("SelectedTree");
 candTree2->Add("../data/ggH_2bp_zz4l_M750_Ga247_all.root");

 TChain *candTree3 = new TChain("newTree");
 candTree3->Add("../data/mH750_2bplus_gg_2e2mu.root");

 TChain *candTree4 = new TChain("SelectedTree");
 candTree4->Add("ggH_2bp_zz4l_M750_Ga247_all.root");

 TH1F *hist1 = new TH1F("hist1","hist1",nbin,min,max);
 TH1F *hist2 = new TH1F("hist2","hist2",nbin,min,max);
 TH1F *hist3 = new TH1F("hist3","hist3",nbin,min,max);
 TH1F *hist4 = new TH1F("hist4","hist4",nbin,min,max);

 char hist1draw[500],hist2draw[500],hist3draw[500],hist4draw[500];
 sprintf(hist1draw,"%s>>hist1",var1);
 sprintf(hist2draw,"%s>>hist2",var2);
 sprintf(hist3draw,"%s>>hist3",var2);
 sprintf(hist4draw,"%s>>hist4",var2);


 candTree1->Draw(hist1draw,"GenHMass<=770&&GenHMass>=730&& genFinalState==2 &&ZZsel>=100");
 candTree2->Draw(hist2draw,"ZZMass<=770&&ZZMass>=730&&flavortype==3");
 candTree3->Draw(hist3draw,"ZZMass<=770&&ZZMass>=730&&flavortype==3");
 candTree4->Draw(hist4draw,"ZZMass<=770&&ZZMass>=730&&flavortype==3");

 TCanvas* c = new TCanvas("c", "c", 800,800);
 c->SetFillColor(0);
 c->cd();

 hist1->SetLineWidth(2);
 hist2->SetLineWidth(2);
 hist3->SetLineWidth(2);
 hist4->SetLineWidth(2);

 hist1->SetLineColor(kBlack);
 hist2->SetLineColor(kRed);
 hist3->SetLineColor(kGreen+2);
 hist4->SetLineColor(kBlue);

 hist1->Scale(1/hist1->Integral());
 hist2->Scale(1/hist2->Integral());
 hist3->Scale(1/hist3->Integral());
 hist4->Scale(1/hist4->Integral());

 hist1->SetTitle(title);
 hist1->GetXaxis()->SetTitle(var2);
 if(var2!="Z1Mass"&&var2!="Z2Mass")
 hist1->GetYaxis()->SetRangeUser(0,0.04);
// hist2->GetYaxis()->SetTitle("Events");

 hist1->Draw("hist");
 hist2->Draw("samehist");
 hist3->Draw("samehist");
 //hist4->Draw("samehist");

  TLegend* leg = new TLegend(.65,0.65,0.85,.85);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hist1,"central","l");
  leg->AddEntry(hist2,"JHUGen 16June","l");
  leg->AddEntry(hist3,"JHUGen 15Winter","l");
 // leg->AddEntry(hist4,"0.3 width","l");
  leg->Draw();

 c->Update();
 char temp[150];
 sprintf(temp,"%s.png",var2);
 c->SaveAs(temp);
 delete c;
 delete hist1;
 delete hist2;
 delete hist3;
 delete hist4;
}

void compare(){
test("Gen_costheta1","helcosthetaZ1","",50,1, 1);
test("Gen_costheta2","helcosthetaZ2","",50,1, 1);
test("Gen_phi","helphi","",50,-3.14, 3.14);
test("Gen_costhetastar","costhetastar","",50,-1,1);
test("Gen_phistar1","phistarZ1","",50,-3.14,3.14);
test("Gen_mZ1","Z1Mass","",80,50,130);
test("Gen_mZ2","Z2Mass","",80,50,130);
}
