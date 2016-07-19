#include <TROOT.h>
#include <vector>
#include "TText.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include <fstream>
#include <iostream>
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
using namespace std;

//mode : 0 ggH->0+, 1 ggH->2b+, 2 VBF->0+
//ch : 0 4mu, 1 4e, 2 2e2mu.

TGraphErrors* makegr(int mode=0, int eleSel=0, int ch=0, Color_t color=2, int marker=20, int line=1){
	gStyle->SetPadLeftMargin(0.1);
	gStyle->SetPadRightMargin(0.3);
	gStyle->SetOptStat(0000);
	gStyle->SetTitleFontSize(0.05);

	char* modeName[3] ={"ggH->0p","ggH->2bp","VBF->0p"};
	char* eleSelName[3]={"REG","TLE","RSE"};
	char* chan[3] ={"4mu","4e","2e2mu"};

	TH1F *hgen_pre;

	//Load Samples

	TFile* f;
	//  TString  inputDir = "root://lxcms03//data3/Higgs/160714/";
	//  TString  inputDir = "root://lxcms03//data3/Higgs/160718/";
	//  TString  inputDir = "/afs/cern.ch/work/c/cayou/HighMass_RunII/4l/160717/CMSSW_8_0_8/src/ZZAnalysis/AnalysisStep/test/prod/samples_2016mc_0718/";
	TString  inputDir = "/afs/cern.ch/work/c/cayou/HighMass_RunII/4l/160717/CMSSW_8_0_8/src/ZZAnalysis/AnalysisStep/test/prod/test/";
	//For samples below 300GeV, only events with GenHMass<PoleMass can be used.

	//  int inputfiles_ggHlow[]={115,120,124,125,/*126,*/130,135,140,145,150,155,160,165,170,175,180,200,210,230,250,270};
	//  int inputfiles_VBFlow[]={115,120,/*124,*/125,126,/*130,*/135,140,145,/*150,*/155,160,165,170,175,180,190,200,210,230,250,270};
	//  int inputfiles_ggH[]={300,350,400,450,500,550,600,700,750,800,900,1000,1500,2000,2500,3000};
	//  int inputfiles_VBF[]={300,350,400,450,500,550,600,700,/*750,*/800,900,1000,2000,2500,3000};
	//  int inputfiles_ggH_2bp[]={750,800,1200,2000,3000,4000};

	//debug
	int inputfiles_ggHlow[]={};
	int inputfiles_VBFlow[]={};
	//  int inputfiles_ggH[]={300,350,400,450,500,550,600,650,700,750,800,900,1500,2500,3000};
	int inputfiles_ggH[]={750};
	int inputfiles_VBF[]={};
	int inputfiles_ggH_2bp[]={};

	int Nfiles_ggH=sizeof(inputfiles_ggH)/sizeof(*inputfiles_ggH);
	int Nfiles_VBF=sizeof(inputfiles_VBF)/sizeof(*inputfiles_VBF);
	int Nfiles_ggHlow=sizeof(inputfiles_ggHlow)/sizeof(*inputfiles_ggHlow);
	int Nfiles_VBFlow=sizeof(inputfiles_VBFlow)/sizeof(*inputfiles_VBFlow);
	int Nfiles_ggH_2bp=sizeof(inputfiles_ggH_2bp)/sizeof(*inputfiles_ggH_2bp);


	char inputfile[100];
	vector<TString> files_ggH;
	vector<TString> files_VBF;

	for (int i=0; i<Nfiles_ggH; i++) {
		sprintf(inputfile,"ggH%d/ZZ4lAnalysis.root",inputfiles_ggH[i]);
		files_ggH.push_back(inputfile);
	}

	for (int i=0; i<Nfiles_ggHlow; i++) {
		//     sprintf(inputfile,"ggH%d/ZZ4lAnalysis_new.root",inputfiles_ggHlow[i]);
		sprintf(inputfile,"ggH%d/ZZ4lAnalysis.root",inputfiles_ggHlow[i]);
		files_ggH.push_back(inputfile);
	}

	for (int i=0; i<Nfiles_ggH_2bp; i++) {
		sprintf(inputfile,"ggGrav2PB_W0p3_M%d/ZZ4lAnalysis.root",inputfiles_ggH_2bp[i]);
		files_ggH.push_back(inputfile);
	}

	for (int i=0; i<Nfiles_VBF; i++) {
		sprintf(inputfile,"VBFH%d/ZZ4lAnalysis.root",inputfiles_VBF[i]);
		files_VBF.push_back(inputfile);
	}

	for (int i=0; i<Nfiles_VBFlow; i++) {
		//     sprintf(inputfile,"VBFH%d/ZZ4lAnalysis_new.root",inputfiles_VBFlow[i]);
		sprintf(inputfile,"VBFH%d/ZZ4lAnalysis.root",inputfiles_VBFlow[i]);
		files_VBF.push_back(inputfile);
	}

	TChain *candTree = new TChain;
	if (eleSel==0) candTree->SetName("ZZTree/candTree");
	else if (eleSel==1) candTree->SetName("ZZTreetle/candTree");
	else if (eleSel==2) candTree->SetName("ZZTreelooseEle/candTree");

	if (mode!=2){ //ggH
		for (vector<TString>::const_iterator file = files_ggH.begin(); file!=files_ggH.end(); ++file){ 
			candTree->Add(inputDir+(*file));
			//get gen hist from root files
			TFile* f = TFile::Open(inputDir+(*file));
			TH1F* h =(TH1F*)f->Get(Form("PlotsZZ/hGenZZMass_%s",chan[ch])); 
			if (file==files_ggH.begin() ) hgen_pre=h;
			else hgen_pre->Add(h);
		}
	}
	else{ //VBF
		for (vector<TString>::const_iterator file = files_VBF.begin(); file!=files_VBF.end(); ++file){
			candTree->Add(inputDir+(*file));
			//get gen hist from root files
			TFile* f = TFile::Open(inputDir+(*file));
			TH1F* h =(TH1F*)f->Get(Form("PlotsZZ/hGenZZMass_%s",chan[ch]));
			if (file==files_ggH.begin() ) hgen_pre=h;
			else hgen_pre->Add(h);
		}
	}

	// draw raw histogram with 1GeV binning, 0GeV to 3500GeV.
	const Int_t m=3500;
	TH1F *hreco = new TH1F("hreco","hreco",m,0,m);
	hgen_pre->SetBins(m,0,m);


	int zzflav,zzflav_tle;
	switch (ch){
		case 0: zzflav=28561;  break; //4mu
		case 1: zzflav=14641;zzflav_tle=29282;  break; //4e
		case 2: zzflav=20449;zzflav_tle=40898;  break; //2e2mu
		default:  cout<<"channel unknown.";
	}


	// Fill histograms

	short genFinalState=0,ZZsel=0,Z1Flav=0,Z2Flav=0;
	vector<short> *LepLepId=0,*LepisLoose=0;
	vector<float> *LepPt=0,*LepEta=0,*LepPhi=0,*reweightingweights=0;
	float GenHMass=0,genHEPMCweight=0,PUWeight=0,dataMCWeight=0;


	candTree->SetBranchAddress("GenHMass",&GenHMass);
	candTree->SetBranchAddress("genFinalState",&genFinalState);
	candTree->SetBranchAddress("ZZsel",&ZZsel);
	candTree->SetBranchAddress("Z1Flav",&Z1Flav);
	candTree->SetBranchAddress("Z2Flav",&Z2Flav);
	candTree->SetBranchAddress("genHEPMCweight",&genHEPMCweight);
	candTree->SetBranchAddress("dataMCWeight",&dataMCWeight);
	candTree->SetBranchAddress("PUWeight",&PUWeight);
	candTree->SetBranchAddress("LepLepId",&LepLepId);
	candTree->SetBranchAddress("LepPt",&LepPt);
	candTree->SetBranchAddress("LepEta",&LepEta);
	candTree->SetBranchAddress("LepPhi",&LepPhi);
	if(mode!=3)
		candTree->SetBranchAddress("reweightingweights",&reweightingweights);

	for (int i = 0; i < candTree->GetEntries(); i++) {
		candTree->GetEntry(i);

		if((eleSel==0||eleSel==2) && ZZsel==120){ //REG and RSE
			if (genFinalState==ch && fabs(Z1Flav*Z2Flav)==zzflav)
				hreco->Fill(GenHMass,(genHEPMCweight * PUWeight * dataMCWeight));
		}

		else if (eleSel==1 && ZZsel==120){ //TLE
			int TLE_index=-1, lep1_index=-1, lep2_index=-1;
			for (int i=0 ; i<LepLepId->size() ; i++){ 
				if (abs(LepLepId->at(i))==22) TLE_index=i;
			}
			if (TLE_index<=1) lep1_index=2, lep2_index=3;
			else lep1_index=0, lep2_index=1;
			TLorentzVector lep1,lep2,lep_TLE;
			lep1.SetPtEtaPhiM(LepPt->at(lep1_index),LepEta->at(lep1_index),LepPhi->at(lep1_index),0);
			lep2.SetPtEtaPhiM(LepPt->at(lep2_index),LepEta->at(lep2_index),LepPhi->at(lep2_index),0);
			lep_TLE.SetPtEtaPhiM(LepPt->at(TLE_index),LepEta->at(TLE_index),LepPhi->at(TLE_index),0);
			double TLE_dR_Z = lep_TLE.DeltaR(lep1+lep2);
			//  cout<<TLE_index<<","<<lep1_index<<","<<lep2_index<<",deltaR"<<TLE_dR_Z<<",pt"<<LepPt->at(TLE_index)<<endl;
			if (LepPt->at(TLE_index)>30 && TLE_dR_Z>1.6){
				if (genFinalState==ch && fabs(Z1Flav*Z2Flav)==zzflav_tle) 
					hreco->Fill(GenHMass,(genHEPMCweight * PUWeight * dataMCWeight));
			} 
		}
	}
	//Calculate Efficiency 

	// bin contents of raw histograms
	double M_raw[m]={0};
	double gen_raw[m]={0};
	double reco_raw[m]={0};

	// const Int_t n=35;
	// double M[n]={115,120,125,126,135,140,145,150,155,160,180,200,230,250,350,400,450,500,600,750,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400,3600};

	// const Int_t n=32;
	// double M[n]={115,120,125,126,135,140,145,150,155,160,180,200,230,250,350,400,450,500,600,750,900,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000};

	const Int_t n=1;
	double M[n]={750};

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
		gen_raw[bin-1] = hgen_pre->GetBinContent(bin);
		reco_raw[bin-1] = hreco->GetBinContent(bin);}

	for (int i=0;i<n;i++){
		mass = M[i];
		if(mass<=180) width=1;
		else if (mass>180 && mass<=300) width =5;
		else if (mass>300&&mass<600) width = 10;
		else if (mass>=600&&mass<1000) width = 25;
		else if (mass>=1000&&mass<2000) width = 150;
		else if (mass>=2000&&mass<2300) width=200;
		else if (mass>=2300&&mass<=3000) width=300;
		else width=400;
		//merge
		for (int j=(mass-width);j<(mass+width);j++){
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
		//cout<<i<<":"<<reco[i]<<","<<gen[i]<<","<<eff[i]<<endl;
	}

	TF1 *polyFunctot= new TF1("polyFunctot","([0]+[1]*TMath::Erf( (x-[2])/[3] ))*([4]+[5]*x+[6]*x*x+[10]*x*x*x)+[7]*TMath::Gaus(x,[8],[9])", 110., 3500);
	polyFunctot->SetParameters(-4.42749e+00,4.61212e+0,-6.21611e+01,1.13168e+02,2.14321e+00,1.04083e-03,4.89570e-07, 0.03, 200, 30,0);
	polyFunctot->SetParLimits(7,0,0.2);
	polyFunctot->SetParLimits(8,160,210);
	polyFunctot->SetParLimits(9,10,70);
	polyFunctot->SetLineColor(color);
	if(mode==0)   polyFunctot->SetLineStyle(1);
	else if(mode==1) {
		polyFunctot->SetParLimits(7,0,0);
		polyFunctot->SetParLimits(8,0,0);
		polyFunctot->SetParLimits(9,0,0); 
		polyFunctot->SetLineStyle(2);}

	TGraphErrors *gr = new TGraphErrors (n,M,eff,massE,effE);

	gr->SetMarkerColor(color);
	gr->SetMarkerStyle(marker);
	gr->SetMarkerSize(1);
	cout<<endl<<endl<<"parameters:"<<endl;
	TFitResultPtr r = gr->Fit("polyFunctot","S");
	gr->SetLineColor(color);

	ofstream fp;
	fp.open ("efficiency_param.txt",std::ofstream::out | std::ofstream::app);
	for (int i=0;i<11;i++) {
		Double_t value = r->Parameter(i);
		char par[200];
		sprintf(par,"%s_%s_param_%d: %E",modeName[mode],chan[ch],i,value);
		fp<<par<<endl;
	}
	fp<<endl;
	fp.close();

	delete hreco;
	return gr;
}

void eff_4l(){
	TCanvas* c2 = new TCanvas("c2", "c2", 1000, 10, 1400, 800);
	//c2->SetLogx(); 
	c2->SetFillColor(0);
	c2->SetRightMargin(0.18);
	TMultiGraph *mg = new TMultiGraph();
	//debug
	//TGraphErrors* ggH_4e = makegr(0,1,1,kGreen,20,1);
	TGraphErrors* ggH_2e2mu = makegr(0,2,2,kBlue,20,1);
	//mg->Add(ggH_4e);
	mg->Add(ggH_2e2mu);
	mg->Draw("AP");

	/*
	   cout<<"ggH_0+_4mu:"<<endl;
	   TGraphErrors* ggH_4mu = makegr(0,0,0,kRed,20,1);
	   cout<<"ggH_0+_4e:"<<endl;
	   TGraphErrors* ggH_4e = makegr(0,0,1,kGreen,20,1);
	   cout<<"ggH_0+_2e2mu:"<<endl;
	   TGraphErrors* ggH_2e2mu = makegr(0,0,2,kBlue,20,1);
	   cout<<"VBF_0+_4mu:"<<endl;
	   TGraphErrors* vbf_4mu = makegr(2,0,0,kRed+2,23,1);
	   cout<<"VBF_0+_4e:"<<endl;
	   TGraphErrors* vbf_4e = makegr(2,0,1,kGreen+2,23,1);
	   cout<<"VBF_0+_2e2mu:"<<endl;
	   TGraphErrors* vbf_2e2mu = makegr(2,0,2,kBlue+2,23,1);
	   cout<<"ggH_2b+_4mu:"<<endl;
	   TGraphErrors* spin2_4mu = makegr(1,0,0,kRed,24,2);
	   cout<<"ggH_2b+_4e:"<<endl;
	   TGraphErrors* spin2_4e = makegr(1,0,1,kGreen,24,2);
	   cout<<"ggH_2b+_2e2mu:"<<endl;
	   TGraphErrors* spin2_2e2mu = makegr(1,0,2,kBlue,24,2);

	   mg->Add(ggH_4e);
	   mg->Add(ggH_4mu);
	   mg->Add(ggH_2e2mu);
	   mg->Add(vbf_4e);
	   mg->Add(vbf_4mu);
	   mg->Add(vbf_2e2mu);
	   mg->Add(spin2_4e);
	   mg->Add(spin2_4mu);
	   mg->Add(spin2_2e2mu);
	   mg->Draw("AP");
	   mg->GetXaxis()->SetTitle("genHMass [GeV]");
	   mg->GetYaxis()->SetTitle("efficiency*acceptance");
	   mg->GetYaxis()->SetRangeUser(0,1);
	   mg->GetXaxis()->SetRangeUser(100,3000);

	   TLegend* leg3 = new TLegend(.83,0.3,0.99,.85);
	   leg3->SetFillColor(0);
	   leg3->SetBorderSize(0);
	   leg3->AddEntry(ggH_4mu,"ggH->0+ 4mu","pl");
	   leg3->AddEntry(ggH_2e2mu,"ggH->0+ 2e2mu","pl");
	   leg3->AddEntry(ggH_4e,"ggH->0+ 4e","pl");  
	   leg3->AddEntry(vbf_4mu,"VBF->0+ 4mu","pl");
	   leg3->AddEntry(vbf_2e2mu,"VBF->0+ 2e2mu","pl");
	   leg3->AddEntry(vbf_4e,"VBF->0+ 4e","pl");
	   leg3->AddEntry(spin2_4mu,"ggH->2b+ 4mu","pl");
	   leg3->AddEntry(spin2_2e2mu,"ggH->2b+ 2e2mu","pl");
	   leg3->AddEntry(spin2_4e,"ggH->2b+ 4e","pl");
	   leg3->Draw();
	   */
	c2->Update();
	c2->SaveAs("/afs/cern.ch/user/c/cayou/www/HighMass/test/sigEfficiencyFit_80X_4l_test.png");
	c2->SaveAs("/afs/cern.ch/user/c/cayou/www/HighMass/test/sigEfficiencyFit_80X_4l_test.pdf");
}
