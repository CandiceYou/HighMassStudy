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
 candTree1->Add("~/Analysis/2l2q/CMSSW_8_0_7_wqin/src/ZZAnalysis/AnalysisStep/test/prod/2016_2l2q_spin2/output.root"); //no width

 TChain *candTree2 = new TChain("ZZTree/candTree");
 candTree2->Add("~/Analysis/2l2q/CMSSW_8_0_7_wqin/src/ZZAnalysis/AnalysisStep/test/prod/2l2qsamples_2bp_madgraph/output.root");

TChain *candTree3 = new TChain("ZZTree/candTree");
 candTree3->Add("/afs/cern.ch/user/r/rbarr/Analysis/2l2q/CMSSW_8_0_7_wqin/src/ZZAnalysis/AnalysisStep/test/prod/2016_2l2q/mytree/PT13TeV/ggHiggs600/ZZ2l2qAnalysis.root");


 TH1F *hist1 = new TH1F("hist1","hist1",nbin,min,max);
 TH1F *hist2 = new TH1F("hist2","hist2",nbin,min,max);
 TH1F *hist3 = new TH1F("hist3","hist3",nbin,min,max);

 char hist1draw[500],hist2draw[500],hist3draw[500], cut[500];
 sprintf(hist1draw,"%s>>hist1",var1);
 sprintf(hist2draw,"%s>>hist2",var2);
 sprintf(hist3draw,"%s>>hist3",var1);
 sprintf(cut,"(GenHMass<=%d&&GenHMass>=%d&& genFinalState==99)", mass+50, mass-50);


 candTree1->Draw(hist1draw,cut);
 candTree2->Draw(hist2draw,cut);
 candTree3->Draw(hist3draw,cut);

 TCanvas* c = new TCanvas("c", "c", 800,800);
 c->SetFillColor(0);
 c->cd();

 hist1->SetLineWidth(2);
 hist2->SetLineWidth(2);
 hist3->SetLineWidth(2);

 hist1->SetLineColor(kRed);
 hist2->SetLineColor(kBlack);
 hist3->SetLineColor(kMagenta+2);

 hist1->Scale(1/hist1->Integral());
 hist2->Scale(1/hist2->Integral());
 hist3->Scale(1); //integrals r hard
 cout << hist1->Integral() << endl <<hist2->Integral() << endl << hist3->Integral() << endl<< endl<< endl;

 hist1->SetTitle(var2);
 hist1->GetXaxis()->SetTitle(var2);
// if(var2!="Z1Mass"&&var2!="Z2Mass")
// hist2->GetYaxis()->SetRangeUser(0,0.04);
// hist2->GetYaxis()->SetTitle("Events");

 hist1->Draw("hist");
 hist2->Draw("samehist");
 hist3->Draw("whatdoiputhere");

 //Legend
 TLegend* leg = new TLegend(.65,0.8,0.85,.9);
 leg->SetFillColorAlpha(0,0);
 leg->SetBorderSize(0);
 leg->AddEntry(hist1,"spin2","l");
 leg->AddEntry(hist2,"spin2_madgraph","l");
 leg->AddEntry(hist3,"spin0","l");
 leg->Draw();

 c->Update();
 char temp[150];
 sprintf(temp,"~/www/HighMassStudy/16.7.21/Mass%04d/%s.png",mass,var1);
 c->SaveAs(temp);
 delete c;
 delete hist1;
 delete hist2;
 delete hist3;
}

void runner(int mass){
test("helcosthetaZ1","helcosthetaZ1","",50,1, 1,mass);
test("helcosthetaZ2","helcosthetaZ2","",50,1, 1,mass);
test("helphi","helphi","",50,-3.14, 3.14,mass);
test("costhetastar","costhetastar","",50,-1,1,mass);
test("phistarZ1","phistarZ1","",50,-3.14,3.14,mass);
test("Z1Mass","Z1Mass","",80,50,130,mass);
test("Z2Mass","Z2Mass","",80,50,130,mass);
}

void plotAngles() {
	int inputfiles[] = {600,800,1000,1400,1600,1800,2000,2500,3000,3500,4000,4500};
	const int n = sizeof(inputfiles)/sizeof(int);
	for(int i = 0; i < n;++i) runner(inputfiles[i]);
}
