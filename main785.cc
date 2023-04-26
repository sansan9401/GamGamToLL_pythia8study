// main78.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Ilkka Helenius <ilkka.m.helenius@jyu.fi>.

// Keywords: photon beam; UPC; photoproduction; exclusive;

// Main program to demonstrate how to generate different types of
// photon-initiated dilepton events in proton-proton collisions.

#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH2D.h"

using namespace Pythia8;

map<TString,TH1*> maphist;
TH1* GetHist(TString histname){
  TH1* h=NULL;
  auto mapit=maphist.find(histname);
  if(mapit!=maphist.end()) return mapit->second;
  return h;
}
void FillHist(TString histname,double value_x,double weight,int n_binx, const double *xbins){
  TH1D *this_hist = (TH1D*)GetHist(histname);
  if( !this_hist ){
    this_hist = new TH1D(histname, histname, n_binx, xbins);
    this_hist->SetDirectory(NULL);
    maphist[histname] = this_hist;
  }
  this_hist->Fill(value_x, weight);
}
void FillHist(TString histname,double value_x,double weight,int n_binx,double x_min,double x_max){
  TH1D *this_hist = (TH1D*)GetHist(histname);
  if( !this_hist ){
    this_hist = new TH1D(histname, histname, n_binx, x_min, x_max);
    this_hist->SetDirectory(NULL);
    maphist[histname] = this_hist;
  }
  this_hist->Fill(value_x, weight);
}
void FillHist(TString histname,double value_x, double value_y,double weight,int n_binx, const double *xbins,int n_biny, const double *ybins){
  TH2D *this_hist = (TH2D*)GetHist(histname);
  if( !this_hist ){
    this_hist = new TH2D(histname, histname, n_binx, xbins, n_biny, ybins);
    this_hist->SetDirectory(NULL);
    maphist[histname] = this_hist;
  }
  this_hist->Fill(value_x, value_y, weight);
}
void FillHist(TString histname,double value_x, double value_y,double weight,int n_binx, const double *xbins,int n_biny, double y_min, double y_max){
  TH2D *this_hist = (TH2D*)GetHist(histname);
  if( !this_hist ){
    this_hist = new TH2D(histname, histname, n_binx, xbins, n_biny, y_min, y_max);
    this_hist->SetDirectory(NULL);
    maphist[histname] = this_hist;
  }
  this_hist->Fill(value_x, value_y, weight);
}
double GetCosThetaCS(Vec4 l0,Vec4 l1){
  double l0pp=(l0.e()+l0.pz())/sqrt(2);
  double l0pm=(l0.e()-l0.pz())/sqrt(2);
  double l1pp=(l1.e()+l1.pz())/sqrt(2);
  double l1pm=(l1.e()-l1.pz())/sqrt(2);
  double dimass=m(l0,l1);
  double dipt=(l0+l1).pT();
  int direction=(l0+l1).pz()>0?1:-1;
  return direction*2*(l0pp*l1pm-l0pm*l1pp)/sqrt(dimass*dimass*(dimass*dimass+dipt*dipt));
}
// The main program.
int main(int argc,char **argv){
  if(argc!=3) return 1;
  int processType=atoi(argv[1]);
  int seed=atoi(argv[2]);
  // Generator.
  Pythia pythia;
  Info info;

  // Decrease the output.
  pythia.readString("Init:showChangedSettings = on");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Next:numberCount = 1000");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 1");
  pythia.readString("Next:numberShowEvent = 1");

  // Beam parameters.
  pythia.readString("Beams:frameType = 1");
  pythia.readString("Beams:eCM = 13000");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");

  // PDF
  if(processType==0) pythia.readString("PDF:pSet = LHAPDF6:MMHT2015qed_nlo");
  else pythia.readString("PDF:pSet = LHAPDF6:MMHT2015qed_nlo_inelastic");
  PDFPtr photon_elastic=make_shared<LHAPDF>(2212,"LHAPDF6:MMHT2015qed_nlo_elastic/0",&info);

  // Four possible contributions: double-dissociative = 1;
  // single-dissociative on side A (B) = 2 (3); elastic-elastic = 4.

  // Enable elastic photon beams according to processType above.
  // For dissociative beams photons from PDFs (NNPDF2.3 by default) are used.
  if ( processType == 4 || processType == 3 ){
    pythia.readString("PDF:beamA2gamma = on");
    pythia.readString("PDF:Proton2gammaSet = 2");
  }
  if ( processType == 4 || processType == 2 ){
    pythia.readString("PDF:beamB2gamma = on");
    pythia.readString("PDF:Proton2gammaSet = 2");
  }
  if ( processType == 5 ){
    pythia.readString("PDF:beamB2gamma = on");
    pythia.setPhotonFluxPtr(0,photon_elastic);
    pythia.readString("PDF:Proton2gammaSet = 0");
  }
  if ( processType == 6 ){
    pythia.readString("PDF:beamA2gamma = on");
    pythia.setPhotonFluxPtr(photon_elastic,0);
    pythia.readString("PDF:Proton2gammaSet = 0");
  }
  if ( processType == 7 ){
    pythia.readString("PDF:beamA2gamma = on");
    pythia.readString("PDF:beamB2gamma = on");
    pythia.setPhotonFluxPtr(photon_elastic,photon_elastic);
    pythia.readString("PDF:Proton2gammaSet = 0");
  }

  // Need to use virtuality dependent photon flux to obtain realistic
  // acoplanarity distribution also for elastic-elastic case.

  // Set outgoing lepton-pair id and switch on relevant process.
  int idLep = 13;
  if      (idLep == 11) pythia.readString("PhotonCollision:gmgm2ee = on");
  else if (idLep == 13) pythia.readString("PhotonCollision:gmgm2mumu = on");

  // Need to prepare for MPIs only for double-dissociative production.
  if ( processType > 1 ) pythia.readString("PartonLevel:MPI = off");

  // Use dipole shower as appropriate for the case where no colour connection
  // between remnants.
  pythia.readString("SpaceShower:dipoleRecoil = on");
  pythia.readString("SpaceShower:pTmaxMatch = 2");
  pythia.readString("SpaceShower:pTdampMatch = 1");

  // Kinematical cuts.
  pythia.readString("PhaseSpace:mHatMin = 50.0");
  pythia.readString("PhaseSpace:pTHatMinDiverge = 1.");

  //seed
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d",seed));

  // Initialize the histogram for acoplanarity.
  vector<double> mbins={50,60,80,100,150,200,400,1000,3000};

  // Initialize the generator.
  pythia.init();

  // Number of events.
  int nEvent = 10000;

  // Begin event loop. Skip if fails.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate next event.
    if (!pythia.next()) continue;

    // Event weight.
    map<TString,double> weights;
    weights[""]=pythia.info.weight();

    for(auto [suffix,weight]:weights){
      if(suffix=="") continue;
      if(!std::isnormal(weight)) weights[suffix]=1.;
    }
    // Find the final muons (last ones in the event record).
    int iPosPlus = 0, iPosMinus = 0;
    for (int i = 0; i < pythia.event.size();++i) {
      if (pythia.event[i].id() == idLep) iPosPlus = i;
      if (pythia.event[i].id() == -idLep) iPosMinus = i;
    }

    // Derive 4-momenta of leptons.
    Vec4 p1 = pythia.event[iPosPlus].p();
    Vec4 p2 = pythia.event[iPosMinus].p();

    // Fill histrogram with possible weights.
    if(m(p1,p2)<50) continue;
    for(auto [suffix,weight]:weights){
      FillHist("mass"+suffix,m(p1,p2),weight,mbins.size()-1,&mbins[0]);
      FillHist("rapidity"+suffix,(p1+p2).rap(),weight,100,-5,5);
      FillHist("pt"+suffix,(p1+p2).pT(),weight,200,0.,400.);
      FillHist("cost"+suffix,GetCosThetaCS(p1,p2),weight,20,-1,1.);
      if((p1+p2).rap()*(p1.eta()-p2.eta())>0) FillHist("forward"+suffix,m(p1,p2),weight,mbins.size()-1,&mbins[0]);
      else FillHist("backward"+suffix,m(p1,p2),weight,mbins.size()-1,&mbins[0]);

      if(TMath::Max(p1.pT(),p2.pT())>20&&TMath::Min(p1.pT(),p2.pT())>10&&fabs(p1.eta())<2.5&&fabs(p2.eta())<2.5){
	FillHist("mass_ac"+suffix,m(p1,p2),weight,mbins.size()-1,&mbins[0]);
	FillHist("rapidity_ac"+suffix,(p1+p2).rap(),weight,100,-5,5);
	FillHist("pt_ac"+suffix,(p1+p2).pT(),weight,200,0.,400.);
	FillHist("cost_ac"+suffix,GetCosThetaCS(p1,p2),weight,20,-1,1.);
	if((p1+p2).rap()*(p1.eta()-p2.eta())>0) FillHist("forward_ac"+suffix,m(p1,p2),weight,mbins.size()-1,&mbins[0]);
	else FillHist("backward_ac"+suffix,m(p1,p2),weight,mbins.size()-1,&mbins[0]);
      }
    }
  } // End of event loop.

  // Show statistics.
  pythia.stat();

  // Normalize to cross section [fb].
  double sigmaNorm = pythia.info.sigmaGen() / pythia.info.weightSum() * 1.e12;
  
  TFile f("hists.root","recreate");
  for(auto [name,hist]:maphist){
    hist->Scale(sigmaNorm);
    hist->Write();
  }
  f.Close();

  // Done.
  return 0;
}
