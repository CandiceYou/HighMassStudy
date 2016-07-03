void 2l2q_eff(int prod=0){

 gStyle->SetPadLeftMargin(0.1);
 gStyle->SetPadRightMargin(0.15);
 gStyle->SetOptFit(0000);
 gStyle->SetOptStat(0000);
 gStyle->SetTitleFontSize(0.05);

// read sample files
  TChain *candTree = new TChain("ZZTree/candTree");

  TString  inputDir = "/afs/cern.ch/work/c/cayou/public/forWenzerRobert/2l2qsamples_new/";
  int inputfiles_ggH[]={200,250,350,400,450,500,600,750,1000,2000};
  int inputfiles_VBF[]={200,250,400,450,500,550,600,700,750,800,900,1000};
  int Nfiles_ggH=sizeof(inputfiles_ggH)/sizeof(*inputfiles_ggH);
  int Nfiles_VBF=sizeof(inputfiles_VBF)/sizeof(*inputfiles_VBF);
  char inputfile[100];
  vector<TString> files_ggH,files_VBF;

   for (int i=0; i<Nfiles_ggH; i++) {
     sprintf(inputfile,"ggHiggs%d/ZZ2l2qAnalysis.root",inputfiles_ggH[i]);
     files_ggH.push_back(inputfile);
   }

   for (int i=0; i<Nfiles_VBF; i++) {
     sprintf(inputfile,"VBFHiggs%d/ZZ2l2qAnalysis.root",inputfiles_VBF[i]);
     files_VBF.push_back(inputfile);
   }

 if (prod==0){ //ggH
  for (vector<TString>::const_iterator file = files_ggH.begin(); file!=files_ggH.end(); ++file)
  candTree->Add(inputDir+(*file));}
 else if (prod==1){ //VBF
  for (vector<TString>::const_iterator file = files_VBF.begin(); file!=files_VBF.end(); ++file)  
  candTree->Add(inputDir+(*file));}


// initialize histograms

 const int jetType = 2;
 const int LepFlav = 2;
 const int Tag = 3;
 char* jetName[2]={"Merged","Resolved"};
 char* lepName[2]={"eeqq","mumuqq"};
 char* legName[2]={"eeqq","#mu#muqq"};
 char* tagName[3]={"vbf-tagged","b-tagged","untagged"};

 TH1F* hreco[jetType][LepFlav][Tag];
 TH1F* heff[jetType][LepFlav][Tag];
 TGraphErrors *gr[jetType][LepFlav][Tag];

 for (int i=0;i<jetType;i++) {
   for (int j=0;j<LepFlav;j++) {
     for (int k=0;k<Tag;k++) {
        hreco[i][j][k] = new TH1F(Form("reco_jet%d_lep%d_tag%d",i,j,k),Form("reco_jet%d_lep%d_tag%d",i,j,k),3500,0,3500);
        heff[i][j][k] = new TH1F(Form("eff_jet%d_lep%d_tag%d",i,j,k),Form("eff_jet%d_lep%d_tag%d",i,j,k),3500,0,3500);
     }
   }
 }

 TH1F *hgen = new TH1F("hgen","hgen",3500,0,3500);
 
// read branches


  vector<short> *ZZsel=0,*ZZCandType=0,*Z2Flav=0;
  vector<float> *ZZMass=0,*Z1Mass=0,*Z2Mass=0,*Z1tau21=0;
  vector<float> *pvbf_VAJHU_highestPTJets=0,*phjj_VAJHU_highestPTJets=0,*pVAMCFM_qqZJJ_bkg=0,*pVAJHUGen_ggZZ_SM_sig=0; 
  vector<bool> *JetIsInZZCand=0;
  vector<float> *JetQGLikelihood=0,*JetBTagger=0,*JetPt=0;;
  float GenHMass=0,xsec=0,genHEPMCweight=0,PUWeight=0;


  candTree->SetBranchAddress("ZZCandType",&ZZCandType);
  candTree->SetBranchAddress("ZZsel",&ZZsel);
  candTree->SetBranchAddress("ZZMassRefit",&ZZMass);
  candTree->SetBranchAddress("Z1Mass",&Z1Mass);
  candTree->SetBranchAddress("Z2Mass",&Z2Mass);
  candTree->SetBranchAddress("pqqZJJ_VAMCFM",&pVAMCFM_qqZJJ_bkg);
  candTree->SetBranchAddress("p0plus_VAJHU",&pVAJHUGen_ggZZ_SM_sig);
  candTree->SetBranchAddress("phjj_VAJHU_highestPTJets",&phjj_VAJHU_highestPTJets);
  candTree->SetBranchAddress("pvbf_VAJHU_highestPTJets",&pvbf_VAJHU_highestPTJets);
  candTree->SetBranchAddress("GenHMass",&GenHMass);
  candTree->SetBranchAddress("xsec",&xsec);
  candTree->SetBranchAddress("genHEPMCweight",&genHEPMCweight);
  candTree->SetBranchAddress("PUWeight",&PUWeight);
  candTree->SetBranchAddress("Z1tau21", &Z1tau21);
  candTree->SetBranchAddress("Z2Flav", &Z2Flav);
  candTree->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood);
  candTree->SetBranchAddress("JetBTagger", &JetBTagger);
  candTree->SetBranchAddress("JetIsInZZCand", &JetIsInZZCand);
  candTree->SetBranchAddress("JetPt", &JetPt);



 for (int i = 0; i < candTree->GetEntries(); i++) {
     candTree->GetEntry(i);

// Fill gen histo
     hgen->Fill(GenHMass,(genHEPMCweight * PUWeight));

     int typ=-1 , lep=-1, tag = -1 , candID=-1;

// Find candidate ID , prefer resolved
     for (int j = 0; j < ZZCandType->size(); j++) {
         if ( ((ZZCandType->at(j)==1 && Z1tau21->at(j)<=0.6)||ZZCandType->at(j)==2) && fabs(ZZsel->at(j))>=100 && Z1Mass->at(j)>=70 && Z1Mass->at(j)<=105 && Z2Mass->at(j)>=60){
            if (ZZCandType->at(j)==1) {typ=0; candID=j;}  //merged, SR
            else if (ZZCandType->at(j)==2) {typ=1; candID=j;break;}  //resolved, SR
      }
     }

// Find leading jets
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
	if (JetQGLikelihood->at(nJet) > -800.) {            // real jets
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
	if (JetQGLikelihood->at(nJet) < -800.) {            // subjets
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

// VBF tagging 

     if (candID!=-1){
        float D_2jet=-1; 
        float phjj=phjj_VAJHU_highestPTJets->at(candID);
        float pvbf=pvbf_VAJHU_highestPTJets->at(candID); 
        if (pvbf>=0 && phjj>=0 && (pvbf+phjj)>0) D_2jet = 1./(1.+0.06*phjj/pvbf);
       
//lepton flavor
        if (abs(Z2Flav->at(candID))==121) lep=0;
        else if (abs(Z2Flav->at(candID))==169) lep=1;
        else cout<<"Error! undefined flavor!"<<endl;

//vbf-tagging && b-tagging
        if (nExtraJets>=2 && D_2jet>0.5) tag=0;  //vbf tagged
        else if ((typ==1 && btag1stJet>0.46 && btag2ndJet>0.46)||(typ==0 && btag1stSubjet>0.46 && btag2ndSubjet>0.46)) tag=1; //b-tagged
        else tag=2;  //untagged

//fill histo
        hreco[typ][lep][tag]->Fill(GenHMass,(genHEPMCweight * PUWeight));

     }
  }


//divide to get efficiency
/*
 for (int i=0;i<jetType;i++) {
   for (int j=0;j<LepFlav;j++) {
     for (int k=0;k<Tag;k++) {
        heff[i][j][k]=hreco[i][j][k]; 
        heff[i][j][k]->Divide(hgen); 
     }
   }
 }
*/

// get efficiency with error 


 for (int i=0;i<jetType;i++) {
   for (int j=0;j<LepFlav;j++) {
     for (int k=0;k<Tag;k++) {

   // bin contents of raw histograms
    const int m=3500;
    double M_raw[m]={0};
    double gen_raw[m]={0};
    double reco_raw[m]={0};
   
    double M[]={200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500,1600,1800,2000,2200,2400,2600,2800,3000};
   // int n=sizeof(M)/sizeof(*M);
    const int n = 30;
    double massE[n]={0}, gen[n]={0}, genE[n]={0}, reco[n]={0}, recoE[n]={0}, eff[n]={0}, effE[n]={0};
    double mass,width;
   
    for (int bin=1;bin<=m;bin++){
    M_raw[bin-1] = hreco[i][j][k]->GetXaxis()->GetBinCenter(bin);
    gen_raw[bin-1] = hgen->GetBinContent(bin);
    reco_raw[bin-1] = hreco[i][j][k]->GetBinContent(bin);}
   
    for (int i=0;i<n;i++){
     mass = M[i];
     if(mass<1000) width=25;
     else if (mass>=1000 && mass<1500) width =100;
     else if (mass>=1500 && mass<2000) width = 200;
     else if (mass>=2000) width = 400;
   
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
   //cout<<eff[i]<<",";
   }
   
    gr[i][j][k] = new TGraphErrors (n,M,eff,massE,effE);
   }
  }
 }


//Fitting Function

  TF1 *polyFunctot= new TF1("polyFunctot","([0]+[1]*TMath::Erf( (x-[2])/[3] ))*([4]+[5]*x+[6]*x*x+[10]*x*x*x)+[7]*TMath::Gaus(x,[8],[9])", 110., 1200);
  polyFunctot->SetParameters(0,0,1000,10,0,0,0,0,0,0,0);

//plotting

int col[2][3]={{2,3,4},{634,418,602}};

 TFile f ("2l2q_Efficiency.root","recreate");
 f.cd();
 for (j=0;j<2;j++){
  TCanvas* c = new TCanvas("c","c", 1600,800);
  c->SetFillColor(0);
  c->SetRightMargin(0.18);
  c->cd();
  TMultiGraph * mg = new TMultiGraph();
  TLegend* leg = new TLegend(.83,0.3,0.99,.85);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
 
  for (int i=0;i<jetType;i++) {
       for (int k=0;k<Tag;k++) {
         if (i==0) gr[i][j][k]->SetMarkerStyle(24);
         else if (i==1) gr[i][j][k]->SetMarkerStyle(20);
         gr[i][j][k]->SetMarkerSize(1);
         gr[i][j][k]->SetMarkerColor(col[i][k]);
         gr[i][j][k]->SetLineColor(col[i][k]);
         gr[i][j][k]->SetLineWidth(2);
         polyFunctot->SetLineColor(col[i][k]);
 //        gr[i][j][k]->Fit("polyFunctot");
         leg->AddEntry(gr[i][j][k],Form("%s_%s_%s",legName[j],jetName[i],tagName[k]),"lp");
         mg->Add(gr[i][j][k]);
 //write to file
         if (prod==0) gr[i][j][k]->SetName(Form("ggH_%s_%s_%s",lepName[j],jetName[i],tagName[k])); 
         else if (prod==1) gr[i][j][k]->SetName(Form("VBF_%s_%s_%s",lepName[j],jetName[i],tagName[k]));
         gr[i][j][k]->Write();
      }
   }
   mg->Draw("APC");
   mg->GetXaxis()->SetTitle("genHMass [GeV]");
   mg->GetYaxis()->SetTitle("efficiency*acceptance");
   mg->GetYaxis()->SetRangeUser(0,0.25);
   mg->GetXaxis()->SetRangeUser(200,3000);
   leg->Draw();
   c->Update();
   if (prod==0){ 
    c->SaveAs(Form("~/www/HighMass/160702/eff_ggH_all_%s.png",lepName[j]));
    c->SaveAs(Form("~/www/HighMass/160702/eff_ggH_all_%s.pdf",lepName[j]));}
   else if (prod==1){
    c->SaveAs(Form("~/www/HighMass/160702/eff_VBF_all_%s.png",lepName[j]));
    c->SaveAs(Form("~/www/HighMass/160702/eff_VBF_all_%s.pdf",lepName[j]));}
 }
 f.Close();
}
