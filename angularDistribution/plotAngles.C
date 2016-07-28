#include "TF1.h"
#include "TText.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveStats.h"
using namespace std;


void test(const char* var1="Gen_costheta1", const char* var2="helcosthetaZ1", const char* title="test", int nbin=30, double min=-1, double max=1, int mass=750){
 gStyle->SetPadLeftMargin(0.1);
 gStyle->SetPadRightMargin(0.15);
 gStyle->SetOptFit(0000);
 gStyle->SetOptStat(0000);
 gStyle->SetTitleFontSize(0.05);

 TChain *candTree1 = new TChain("ZZTree/candTree");
//choose files with polemass>genmass

 candTree1->Add("../data/spin2.root"); //no width
// candTree1->Add("data_new/ggH800_ZZ4lAnalysis_addwt.root");
// candTree1->Add("data_new/ggH900_ZZ4lAnalysis_addwt.root");
// candTree1->Add("data_new/ggH1000_ZZ4lAnalysis_addwt.root");
// candTree1->Add("data_new/ggH1500_ZZ4lAnalysis_addwt.root");
// candTree1->Add("data_new/ggH2000_ZZ4lAnalysis_addwt.root");
// candTree1->Add("data_new/ggH3000_ZZ4lAnalysis_addwt.root");

 TChain *candTree2 = new TChain("ZZTree/candTree");
// candTree2->Add("data_new/ggH_2bp_zz4l_M750_Ga247_all.root");
 candTree2->Add("../data/2l2qsamples_2bp_madgraph/BulkGrav600/ZZ2l2qAnalysis.root");

 TH1F *hist1 = new TH1F("hist1","hist1",nbin,min,max);
 TH1F *hist2 = new TH1F("hist2","hist2",nbin,min,max);


 char hist1draw[500],hist2draw[500];
 sprintf(hist1draw,"%s>>hist1",var1);
 sprintf(hist2draw,"%s>>hist2",var2);


 candTree1->Draw(hist1draw,"(GenHMass<=650&&GenHMass>=550&& genFinalState==99)");
 candTree2->Draw(hist2draw,"(GenHMass<=650&&GenHMass>=550&& genFinalState==99)");

 TCanvas* c = new TCanvas("c", "c", 800,800);
 c->SetFillColor(0);
 c->cd();

 hist1->SetLineWidth(2);
 hist2->SetLineWidth(2);

 hist1->SetLineColor(kRed);
 hist2->SetLineColor(kBlack);


 hist1->Scale(1/hist1->Integral());
 hist2->Scale(1/hist2->Integral());

 hist1->SetTitle(var2);
 hist1->GetXaxis()->SetTitle(var2);
// if(var2!="Z1Mass"&&var2!="Z2Mass")
// hist2->GetYaxis()->SetRangeUser(0,0.04);
// hist2->GetYaxis()->SetTitle("Events");

 hist1->Draw("hist");
 hist2->Draw("samehist");


//Legend
 TLegend* leg = new TLegend(.65,0.8,0.85,.9);
 leg->SetFillColorAlpha(0,0);
 leg->SetBorderSize(0);
 leg->AddEntry(hist1,"spin2","l");
 leg->AddEntry(hist2,"spin2_madgraph","l");
 leg->Draw();
                           

 c->Update();
 char temp[150];
 sprintf(temp,"~/www/HighMassStudy/16.7.21/Mass%04d/%s.png",mass,var1);
 c->SaveAs(temp);
 delete c;
 delete hist1;
 delete hist2;
}


void plotAngles(){
test("helcosthetaZ1","helcosthetaZ1","",50,1, 1);
test("helcosthetaZ2","helcosthetaZ2","",50,1, 1);
test("helphi","helphi","",50,-3.14, 3.14);
test("costhetastar","costhetastar","",50,-1,1);
test("phistarZ1","phistarZ1","",50,-3.14,3.14);
test("Z1Mass","Z1Mass","",80,50,130);
test("Z2Mass","Z2Mass","",80,50,130);
}
