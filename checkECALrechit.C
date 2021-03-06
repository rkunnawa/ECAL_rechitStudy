// macro to check the ECAL rechits with photons (selected by isolation cuts) near it.

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TRandom.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TEventList.h>
#include <TSystem.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"

#define pi 3.14159265359

float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = fabs(phi1 - phi2);
  if(dphi > pi)dphi -= 2*pi;
  float dr = sqrt(pow(deta,2) + pow(dphi,2));
  return dr;
}


#define NOBJECT_MAX 16384
#define MAXHITS 100000

using namespace std;

void checkECALrechit(int startfile = 0,
		     int endfile = 1,
		     std::string kFoname="output")

{


  TStopwatch timer;
  timer.Start();
  
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gStyle->SetOptStat(0);

  bool printDebug = false;
  bool doIsoPhoMatching = true;
  bool doJetMatching = false;
    
  TDatime date;

  std::string infile_Forest;

  infile_Forest = "pbpb_expressdata.txt";
  std::ifstream instr_Forest(infile_Forest.c_str(),std::ifstream::in);
  std::string filename_Forest;
  
  if(printDebug)cout<<"reading from "<<startfile<<" to "<<endfile<<endl;
  
  for(int ifile = 0;ifile<startfile;ifile++){
    instr_Forest>>filename_Forest;
  }

  const int N = 7; //6

  TChain * jetTree[N];

  string dir[N];
  dir[0] = "hltanalysis";
  dir[1] = "skimanalysis";
  dir[2] = "hiEvtAnalyzer";
  dir[3] = "akPu4CaloJetAnalyzer";
  dir[4] = "ggHiNtuplizer";
  dir[5] = "rechitanalyzer";
  dir[6] = "rechitanalyzer";

  string trees[N] = {
    "HltTree",
    "HltTree",
    "HiTree",
    "t",
    "EventTree",
    "eb",
    "ee"
  };

  for(int t = 0;t<N;++t){
    jetTree[t] = new TChain(string(dir[t]+"/"+trees[t]).data());
  }//tree loop ends
  
  for(int ifile = startfile; ifile<endfile; ++ifile){

    instr_Forest>>filename_Forest;
    cout<<"filename: "<<filename_Forest<<endl;

    for(int t = 0; t<N; ++t){
      jetTree[t]->Add(filename_Forest.c_str());    
      if(printDebug)cout << "Tree loaded  " << string(dir[t]+"/"+trees[t]).data() << endl;
      if(printDebug)cout << "Entries : " << jetTree[t]->GetEntries() << endl;
    }
    
    cout<<"Total number of events loaded in HiForest = "<<jetTree[2]->GetEntries()<<endl;

  }

  //hlt
  int minbiasBit;
  
  // jet
  int nref_F;
  float pt_F[1000];
  float rawpt_F[1000];
  float eta_F[1000];
  float phi_F[1000];

  // event
  float vz_F;
  int evt_F;
  int run_F;
  int lumi_F;
  int hiBin_F;

  // skim
  int pcollisionEventSelection_F;
  int pHBHENoiseFilter_F;

  // photon 
  std::vector<float> *ecalIso = 0;
  std::vector<float> *hcalIso = 0;
  std::vector<float> *trackIso = 0;
  std::vector<float> *hadronicOverEm = 0;
  std::vector<float> *phoSigmaIEtaIEta = 0;
  std::vector<float> *phoE = 0;
  std::vector<float> *phoEt = 0;
  std::vector<float> *phoeta = 0;
  std::vector<float> *phophi = 0;
  std::vector<float> *pho_swissCrx = 0;
  std::vector<float> *pho_seedTime = 0;
  
  // rechit ee,eb
  int eb_n;
  float eb_e[MAXHITS];
  float eb_et[MAXHITS];
  float eb_phi[MAXHITS];
  float eb_eta[MAXHITS];
  float eb_perp[MAXHITS];
  float eb_chi2[MAXHITS];
  float eb_err[MAXHITS];
  
  int ee_n;
  float ee_e[MAXHITS];
  float ee_et[MAXHITS];
  float ee_phi[MAXHITS];
  float ee_eta[MAXHITS];
  float ee_perp[MAXHITS];
  float ee_chi2[MAXHITS];
  float ee_err[MAXHITS];

  // set the branch address

  jetTree[0]->SetBranchAddress("HLT_HIL1MinimumBiasHF1ANDExpress_v1",&minbiasBit);
  
  jetTree[1]->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection_F);
  jetTree[1]->SetBranchAddress("pHBHENoiseFilterResultProducer",&pHBHENoiseFilter_F);

  jetTree[2]->SetBranchAddress("evt",&evt_F);
  jetTree[2]->SetBranchAddress("run",&run_F);
  jetTree[2]->SetBranchAddress("lumi",&lumi_F);
  jetTree[2]->SetBranchAddress("vz",&vz_F);
  jetTree[2]->SetBranchAddress("hiBin",&hiBin_F);

  jetTree[3]->SetBranchAddress("nref",&nref_F);
  jetTree[3]->SetBranchAddress("jtpt",pt_F);
  jetTree[3]->SetBranchAddress("jteta",eta_F);
  jetTree[3]->SetBranchAddress("rawpt",eta_F);
  jetTree[3]->SetBranchAddress("jtphi",phi_F);

  jetTree[4]->SetBranchAddress("pho_ecalClusterIsoR4",&ecalIso);
  jetTree[4]->SetBranchAddress("pho_hcalRechitIsoR4",&hcalIso);
  jetTree[4]->SetBranchAddress("pho_trackIsoR4PtCut20",&trackIso);
  jetTree[4]->SetBranchAddress("phoHoverE",&hadronicOverEm);
  jetTree[4]->SetBranchAddress("phoSigmaIEtaIEta",&phoSigmaIEtaIEta);
  jetTree[4]->SetBranchAddress("phoE",&phoE);
  jetTree[4]->SetBranchAddress("phoEt",&phoEt);
  jetTree[4]->SetBranchAddress("phoPhi",&phophi);
  jetTree[4]->SetBranchAddress("phoEta",&phoeta);
  jetTree[4]->SetBranchAddress("pho_swissCrx",&pho_swissCrx);
  jetTree[4]->SetBranchAddress("pho_seedTime",&pho_seedTime);

  jetTree[5]->SetBranchAddress("n",&eb_n);
  jetTree[5]->SetBranchAddress("e",eb_e);
  jetTree[5]->SetBranchAddress("et",eb_et);
  jetTree[5]->SetBranchAddress("phi",eb_phi);
  jetTree[5]->SetBranchAddress("eta",eb_eta);
  jetTree[5]->SetBranchAddress("perp",eb_perp);
  jetTree[5]->SetBranchAddress("chi2",eb_chi2);
  jetTree[5]->SetBranchAddress("eError",eb_err);

  jetTree[6]->SetBranchAddress("n",&ee_n);
  jetTree[6]->SetBranchAddress("e",ee_e);
  jetTree[6]->SetBranchAddress("et",ee_et);
  jetTree[6]->SetBranchAddress("phi",ee_phi);
  jetTree[6]->SetBranchAddress("eta",ee_eta);
  jetTree[6]->SetBranchAddress("perp",ee_perp);
  jetTree[6]->SetBranchAddress("chi2",ee_chi2);
  jetTree[6]->SetBranchAddress("eError",ee_err);


  // now define the output file and start the loops
  TFile *fout = new TFile(Form("%s_%d.root",kFoname.c_str(),endfile),"RECREATE");
  fout->cd();

  // define the histograms:
  TH1F * hChi2_rechit_isopho_ee = new TH1F("hChi2_rechit_isopho_ee","chi2 of rechits near isolated photons, endcap",70, 0, 70);
  TH1F * hChi2_rechit_e5_isopho_ee = new TH1F("hChi2_rechit_e5_isopho_ee","chi2 of rechits (e>5GeV) near isolated photons, endcap",70, 0, 70);
  TH1F * hChi2_rechit_isopho_eb = new TH1F("hChi2_rechit_isopho_eb","chi2 of rechits near isolated photons, barrel",70, 0, 70);
  TH1F * hChi2_rechit_e5_isopho_eb = new TH1F("hChi2_rechit_e5_isopho_eb","chi2 of rechits (e>5GeV) near isolated photons, barrel",70, 0, 70);

  TH1F * hErr_rechit_isopho_ee = new TH1F("hErr_rechit_isopho_ee","Err of rechits near isolated photons, endcap",100, 0, 0.5);
  TH1F * hErr_rechit_e5_isopho_ee = new TH1F("hErr_rechit_e5_isopho_ee","Err of rechits (e>5GeV) near isolated photons, endcap",100, 0, 0.5);
  TH1F * hErr_rechit_isopho_eb = new TH1F("hErr_rechit_isopho_eb","Err of rechits near isolated photons, barrel",100, 0, 0.5);
  TH1F * hErr_rechit_e5_isopho_eb = new TH1F("hErr_rechit_e5_isopho_eb","Err of rechits (e>5GeV) near isolated photons, barrel",100, 0, 0.5);

  TH2F * hChi2_vs_Err_isopho_ee = new TH2F("hChi2_vs_Err_isopho_ee","Endcap Rechits near isolated photons Endcap;Err;Chi2",100, 0, 0.5, 70, 0, 70);
  TH2F * hChi2_vs_Err_e5_isopho_ee = new TH2F("hChi2_vs_Err_e5_isopho_ee","Endcap Rechits (e>5GeV) near isolated photons Endcap;Err;Chi2",100, 0, 0.5, 70, 0, 70);
  TH2F * hChi2_vs_Err_isopho_eb = new TH2F("hChi2_vs_Err_isopho_eb","Endcap Rechits near isolated photons Barrel;Err;Chi2",100, 0, 0.5, 70, 0, 70);
  TH2F * hChi2_vs_Err_e5_isopho_eb = new TH2F("hChi2_vs_Err_e5_isopho_eb","Endcap Rechits (e>5GeV) near isolated photons Barrel;Err;Chi2",100, 0, 0.5, 70, 0, 70);


  TH1F * hChi2_rechit_Jet_ee = new TH1F("hChi2_rechit_Jet_ee","chi2 of rechits near jets , endcap",70, 0, 70);
  TH1F * hChi2_rechit_e5_Jet_ee = new TH1F("hChi2_rechit_e5_Jet_ee","chi2 of rechits (e>5GeV) near jets , endcap",70, 0, 70);
  TH1F * hChi2_rechit_Jet_eb = new TH1F("hChi2_rechit_Jet_eb","chi2 of rechits near jets , barrel",70, 0, 70);
  TH1F * hChi2_rechit_e5_Jet_eb = new TH1F("hChi2_rechit_e5_Jet_eb","chi2 of rechits (e>5GeV) near jets , barrel",70, 0, 70);

  TH1F * hErr_rechit_Jet_ee = new TH1F("hErr_rechit_Jet_ee","Err of rechits near jets , endcap",100, 0, 0.5);
  TH1F * hErr_rechit_e5_Jet_ee = new TH1F("hErr_rechit_e5_Jet_ee","Err of rechits (e>5GeV) near jets , endcap",100, 0, 0.5);
  TH1F * hErr_rechit_Jet_eb = new TH1F("hErr_rechit_Jet_eb","Err of rechits near jets , barrel",100, 0, 0.5);
  TH1F * hErr_rechit_e5_Jet_eb = new TH1F("hErr_rechit_e5_Jet_eb","Err of rechits (e>5GeV) near jets , barrel",100, 0, 0.5);

  TH2F * hChi2_vs_Err_Jet_ee = new TH2F("hChi2_vs_Err_Jet_ee","Endcap Rechits near jets  Endcap;Err;Chi2",100, 0, 0.5, 70, 0, 70);
  TH2F * hChi2_vs_Err_e5_Jet_ee = new TH2F("hChi2_vs_Err_e5_Jet_ee","Endcap Rechits (e>5GeV) near jets  Endcap;Err;Chi2",100, 0, 0.5, 70, 0, 70);
  TH2F * hChi2_vs_Err_Jet_eb = new TH2F("hChi2_vs_Err_Jet_eb","Endcap Rechits near jets  Barrel;Err;Chi2",100, 0, 0.5, 70, 0, 70);
  TH2F * hChi2_vs_Err_e5_Jet_eb = new TH2F("hChi2_vs_Err_e5_Jet_eb","Endcap Rechits (e>5GeV) near jets  Barrel;Err;Chi2",100, 0, 0.5, 70, 0, 70);

  
  TH1F * hEvtno = new TH1F("hEvtno","",2, 0, 2);

  TH1F * hCent = new TH1F("hcent","",200, 0, 200);

  TH1F * hPhotonSpectra_bad_010_barrel = new TH1F("hPhotonSpectra_bad_010_barrel","",100, 0, 300);
  TH1F * hPhotonSpectra_bad_1030_barrel = new TH1F("hPhotonSpectra_bad_1030_barrel","",100, 0, 300);
  TH1F * hPhotonSpectra_bad_3050_barrel = new TH1F("hPhotonSpectra_bad_3050_barrel","",100, 0, 300);
  TH1F * hPhotonSpectra_bad_50100_barrel = new TH1F("hPhotonSpectra_bad_50100_barrel","",100, 0, 300);
  TH1F * hPhotonSpectra_bad_barrel = new TH1F("hPhotonSpectra_bad_barrel","",100, 0, 300);

  TH1F * hPhotonSpectra_good_010_barrel = new TH1F("hPhotonSpectra_good_010_barrel","",100, 0, 300);
  TH1F * hPhotonSpectra_good_1030_barrel = new TH1F("hPhotonSpectra_good_1030_barrel","",100, 0, 300);
  TH1F * hPhotonSpectra_good_3050_barrel = new TH1F("hPhotonSpectra_good_3050_barrel","",100, 0, 300);
  TH1F * hPhotonSpectra_good_50100_barrel = new TH1F("hPhotonSpectra_good_50100_barrel","",100, 0, 300);
  TH1F * hPhotonSpectra_good_barrel = new TH1F("hPhotonSpectra_good_barrel","",100, 0, 300);

  TH1F * hPhotonSpectra_bad_010_endcap = new TH1F("hPhotonSpectra_bad_010_endcap","",100, 0, 300);
  TH1F * hPhotonSpectra_bad_1030_endcap = new TH1F("hPhotonSpectra_bad_1030_endcap","",100, 0, 300);
  TH1F * hPhotonSpectra_bad_3050_endcap = new TH1F("hPhotonSpectra_bad_3050_endcap","",100, 0, 300);
  TH1F * hPhotonSpectra_bad_50100_endcap = new TH1F("hPhotonSpectra_bad_50100_endcap","",100, 0, 300);
  TH1F * hPhotonSpectra_bad_endcap = new TH1F("hPhotonSpectra_bad_endcap","",100, 0, 300);

  TH1F * hPhotonSpectra_good_010_endcap = new TH1F("hPhotonSpectra_good_010_endcap","",100, 0, 300);
  TH1F * hPhotonSpectra_good_1030_endcap = new TH1F("hPhotonSpectra_good_1030_endcap","",100, 0, 300);
  TH1F * hPhotonSpectra_good_3050_endcap = new TH1F("hPhotonSpectra_good_3050_endcap","",100, 0, 300);
  TH1F * hPhotonSpectra_good_50100_endcap = new TH1F("hPhotonSpectra_good_50100_endcap","",100, 0, 300);
  TH1F * hPhotonSpectra_good_endcap = new TH1F("hPhotonSpectra_good_endcap","",100, 0, 300);

  
  if(printDebug) cout<<"Running through all the events now"<<endl;
  Long64_t nentries = jetTree[0]->GetEntries();
  //Long64_t nentries = 50000;
  if(printDebug) nentries = 10;
  TRandom rnd;

  for(int nEvt = 0; nEvt < nentries; ++ nEvt) {

    if(nEvt%1000 == 0)cout<<nEvt<<"/"<<nentries<<endl;
    if(printDebug)cout<<"nEvt = "<<nEvt<<endl;
    
    for(int t = 0;t<N;++t) jetTree[t]->GetEntry(nEvt);

    if(pcollisionEventSelection_F==0) continue;
    if(pHBHENoiseFilter_F == 0) continue;
    if(fabs(vz_F)>15) continue;
    if(!minbiasBit) continue; 

    hEvtno->Fill(1);
    hCent->Fill(hiBin_F);
    
    // minbias selection
    
    // photon isolation cuts

    std::vector<float> isopho_E;
    std::vector<float> isopho_Et;
    std::vector<float> isopho_eta;
    std::vector<float> isopho_phi;

    for(unsigned ipho = 0; ipho<phoE->size(); ++ipho){

      bool passedSpikeRejection = (   phoSigmaIEtaIEta->at(ipho) > 0.002
				   && pho_swissCrx->at(ipho) < 0.9
				   && TMath::Abs(pho_seedTime->at(ipho)) < 3);
      
      bool passedIsolation = (   ecalIso->at(ipho) < 4.2
			      && hcalIso->at(ipho)   < 2.2
			      && trackIso->at(ipho) < 2
			      && hadronicOverEm->at(ipho) < 0.1 );

      if(passedSpikeRejection && passedIsolation) {

	isopho_E.push_back(phoE->at(ipho));
	isopho_Et.push_back(phoEt->at(ipho));
	isopho_eta.push_back(phoeta->at(ipho));
	isopho_phi.push_back(phophi->at(ipho));
	
      }// isolated photon selection
      
    }// photon loop

    if(isopho_E.size()!=0 && doIsoPhoMatching){
    // now do the delta R matching with the ecal rechits and jets
      for(unsigned isop = 0; isop<isopho_E.size(); ++isop){

	bool Skip = false;
	// ecal rechit barrel
	for(int nrec = 0; nrec<eb_n; ++nrec){
	  if(Skip) break;	  
	  float delR = deltaR(isopho_eta[isop], isopho_phi[isop], eb_eta[nrec], eb_phi[nrec]);
	  if(delR<0.2){
	    // fill the histograms you want
	    hChi2_rechit_isopho_eb->Fill(eb_chi2[nrec]);
	    if(eb_e[nrec]>=5.0) hChi2_rechit_e5_isopho_eb->Fill(eb_chi2[nrec]);
	    hErr_rechit_isopho_eb->Fill(eb_err[nrec]);
	    if(eb_e[nrec]>=5.0) hErr_rechit_e5_isopho_eb->Fill(eb_err[nrec]);
	    hChi2_vs_Err_isopho_eb->Fill(eb_err[nrec], eb_chi2[nrec]);
	    if(eb_e[nrec]>=5.0) hChi2_vs_Err_e5_isopho_eb->Fill(eb_err[nrec], eb_chi2[nrec]);

	    if(eb_chi2[nrec]>=64.0){
	      if(hiBin_F <20) hPhotonSpectra_bad_010_barrel->Fill(isopho_Et[isop]);
	      if(hiBin_F >=20 && hiBin_F<60) hPhotonSpectra_bad_1030_barrel->Fill(isopho_Et[isop]);
	      if(hiBin_F >=60 && hiBin_F<100) hPhotonSpectra_bad_3050_barrel->Fill(isopho_Et[isop]);
	      if(hiBin_F >=100) hPhotonSpectra_bad_50100_barrel->Fill(isopho_Et[isop]);
	      Skip = true;
	    }
	    if(eb_chi2[nrec]<64.0){
	      if(hiBin_F <20) hPhotonSpectra_good_010_barrel->Fill(isopho_Et[isop]);
	      if(hiBin_F >=20 && hiBin_F<60) hPhotonSpectra_good_1030_barrel->Fill(isopho_Et[isop]);
	      if(hiBin_F >=60 && hiBin_F<100) hPhotonSpectra_good_3050_barrel->Fill(isopho_Et[isop]);
	      if(hiBin_F >=100) hPhotonSpectra_good_50100_barrel->Fill(isopho_Et[isop]);
	      Skip = true;
	    }
	    
	  }	  
	}// rechit barrel loop
	
	Skip = false;
	// ecal rechit endcap
	for(int nrec = 0; nrec<ee_n; ++nrec){  
	  if(Skip) break;	  
	  float delR = deltaR(isopho_eta[isop], isopho_phi[isop], ee_eta[nrec], ee_phi[nrec]);
	  if(delR<0.2){
	    // fill the histograms you want
	    hChi2_rechit_isopho_ee->Fill(ee_chi2[nrec]);
	    if(ee_e[nrec]>=5.0) hChi2_rechit_e5_isopho_ee->Fill(ee_chi2[nrec]);
	    hErr_rechit_isopho_ee->Fill(ee_err[nrec]);
	    if(ee_e[nrec]>=5.0) hErr_rechit_e5_isopho_ee->Fill(ee_err[nrec]);
	    hChi2_vs_Err_isopho_ee->Fill(ee_err[nrec], ee_chi2[nrec]);
	    if(ee_e[nrec]>=5.0) hChi2_vs_Err_e5_isopho_ee->Fill(ee_err[nrec], ee_chi2[nrec]);	    

	    if(ee_chi2[nrec]>=64.0){
	      if(hiBin_F <20) hPhotonSpectra_bad_010_endcap->Fill(isopho_Et[isop]);
	      if(hiBin_F >=20 && hiBin_F<60) hPhotonSpectra_bad_1030_endcap->Fill(isopho_Et[isop]);
	      if(hiBin_F >=60 && hiBin_F<100) hPhotonSpectra_bad_3050_endcap->Fill(isopho_Et[isop]);
	      if(hiBin_F >=100) hPhotonSpectra_bad_50100_endcap->Fill(isopho_Et[isop]);
	      Skip = true;
	    }
	    if(ee_chi2[nrec]<64.0){
	      if(hiBin_F <20) hPhotonSpectra_good_010_endcap->Fill(isopho_Et[isop]);
	      if(hiBin_F >=20 && hiBin_F<60) hPhotonSpectra_good_1030_endcap->Fill(isopho_Et[isop]);
	      if(hiBin_F >=60 && hiBin_F<100) hPhotonSpectra_good_3050_endcap->Fill(isopho_Et[isop]);
	      if(hiBin_F >=100) hPhotonSpectra_good_50100_endcap->Fill(isopho_Et[isop]);
	      Skip = true;
	    }
	  }
	}// rechit endcap loop

      }// isolated photon loop

    }// iso photon size

    isopho_E.clear();
    isopho_Et.clear();
    isopho_eta.clear();
    isopho_phi.clear();

    // now do the matching with jets.

    if(doJetMatching){
      for(int njet = 0; njet<nref_F; ++njet){

	if(pt_F[njet] <= 30.0) continue;
	
	// ecal rechit barrel
	for(int nrec = 0; nrec<eb_n; ++nrec){  
	  float delR = deltaR(eta_F[njet], phi_F[njet], eb_eta[nrec], eb_phi[nrec]);
	  if(delR<0.2){
	    // fill the histograms you want
	    hChi2_rechit_Jet_eb->Fill(eb_chi2[nrec]);
	    if(eb_e[nrec]>=5.0) hChi2_rechit_e5_Jet_eb->Fill(eb_chi2[nrec]);
	    hErr_rechit_Jet_eb->Fill(eb_err[nrec]);
	    if(eb_e[nrec]>=5.0) hErr_rechit_e5_Jet_eb->Fill(eb_err[nrec]);
	    hChi2_vs_Err_Jet_eb->Fill(eb_err[nrec], eb_chi2[nrec]);
	    if(eb_e[nrec]>=5.0) hChi2_vs_Err_e5_Jet_eb->Fill(eb_err[nrec], eb_chi2[nrec]);
	    
	  }
	}// rechit barrel loop

	// ecal rechit endcap
	for(int nrec = 0; nrec<ee_n; ++nrec){	  
	  float delR = deltaR(eta_F[njet], phi_F[njet], ee_eta[nrec], ee_phi[nrec]);
	  if(delR<0.2){
	    // fill the histograms you want
	    hChi2_rechit_Jet_ee->Fill(ee_chi2[nrec]);
	    if(ee_e[nrec]>=5.0) hChi2_rechit_e5_Jet_ee->Fill(ee_chi2[nrec]);
	    hErr_rechit_Jet_ee->Fill(ee_err[nrec]);
	    if(ee_e[nrec]>=5.0) hErr_rechit_e5_Jet_ee->Fill(ee_err[nrec]);
	    hChi2_vs_Err_Jet_ee->Fill(ee_err[nrec], ee_chi2[nrec]);
	    if(ee_e[nrec]>=5.0) hChi2_vs_Err_e5_Jet_ee->Fill(ee_err[nrec], ee_chi2[nrec]);
	    
	  }
	}// rechit endcap loop

      }
    }//jet matching
    
  }// event loop

  fout->Write();
  
  timer.Stop();
  cout<<"Macro finished: "<<endl;
  cout<<"CPU time (min)  = "<<(Float_t)timer.CpuTime()/60<<endl;
  cout<<"Real time (min) = "<<(Float_t)timer.RealTime()/60<<endl;
  
  
}
