// main78.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Ilkka Helenius <ilkka.m.helenius@jyu.fi>.

// Keywords: photon beam; UPC; photoproduction; exclusive;

// Main program to demonstrate how to generate different types of
// photon-initiated dilepton events in proton-proton collisions.

#include "Pythia8/Pythia.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TFile.h"

using namespace Pythia8;
int seed=0;
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


void Run(double cut) {
  for(auto [name,hist]:maphist){
    if(hist) delete hist;
  }
  maphist.clear();

  // Generator.
  Pythia pythia;

  // Decrease the output.
  pythia.readString("Init:showChangedSettings = on");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Next:numberCount = 1000");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Beam parameters.
  pythia.readString("Beams:frameType = 1");
  pythia.readString("Beams:eCM = 13000");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");

  // PDF
  pythia.readString("PDF:pSet = LHAPDF6:NNPDF31_nlo_as_0118_luxqed");

  // Four possible contributions: double-dissociative = 1;
  // single-dissociative on side A (B) = 2 (3); elastic-elastic = 4.
  int processType = 1;

  // Enable elastic photon beams according to processType above.
  // For dissociative beams photons from PDFs (NNPDF2.3 by default) are used.
  if ( processType == 4 || processType == 3 )
    pythia.readString("PDF:beamA2gamma = on");
  if ( processType == 4 || processType == 2 )
    pythia.readString("PDF:beamB2gamma = on");

  // Need to use virtuality dependent photon flux to obtain realistic
  // acoplanarity distribution also for elastic-elastic case.
  pythia.readString("PDF:Proton2gammaSet = 2");

  // Set outgoing lepton-pair id and switch on relevant process.
  int idLep = 13;
  if      (idLep == 11) pythia.readString("PhotonCollision:gmgm2ee = on");
  else if (idLep == 13) pythia.readString("PhotonCollision:gmgm2mumu = on");

  // Need to prepare for MPIs only for double-dissociative production.
  if ( processType != 1 ) pythia.readString("PartonLevel:MPI = off");

  // Use dipole shower as appropriate for the case where no colour connection
  // between remnants.
  pythia.readString("SpaceShower:dipoleRecoil = on");
  pythia.readString("SpaceShower:pTmaxMatch = 2");
  pythia.readString("SpaceShower:pTdampMatch = 1");

  // Kinematical cuts.
  pythia.readString("PhaseSpace:mHatMin = 50.0");
  pythia.readString(Form("PhaseSpace:pTHatMinDiverge = %f",cut));

  //seed
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d",seed));

  // Initialize the generator.
  pythia.init();

  // Number of events.
  int nEvent = 3000;

  vector<double> mbins={50,60,80,100,150,200,400,1000,3000};
  TString suffix=Form("_cut%g",cut);

  // Begin event loop. Skip if fails.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate next event.
    if (!pythia.next()) continue;

    // Event weight.
    double weight = pythia.info.weight();

    // Find the final muons (last ones in the event record).
    int iPosPlus = 0, iPosMinus = 0;
    for (int i = 0; i < pythia.event.size();++i) {
      if (pythia.event[i].id() == idLep) iPosPlus = i;
      if (pythia.event[i].id() == -idLep) iPosMinus = i;
    }

    // Derive 4-momenta of leptons.
    Vec4 p1 = pythia.event[iPosPlus].p();
    Vec4 p2 = pythia.event[iPosMinus].p();
    
    if(m(p1,p2)<50) continue;
    FillHist("mass"+suffix,m(p1,p2),weight,mbins.size()-1,&mbins[0]);
    FillHist("rapidity"+suffix,(p1+p2).rap(),weight,100,-5,5);
    FillHist("pt"+suffix,(p1+p2).pT(),weight,200,0.,400.);
    if(TMath::Max(p1.pT(),p2.pT())>20&&TMath::Min(p1.pT(),p2.pT())>10&&fabs(p1.eta())<2.5&&fabs(p2.eta())<2.5){
      FillHist("mass_ac"+suffix,m(p1,p2),weight,mbins.size()-1,&mbins[0]);
      FillHist("rapidity_ac"+suffix,(p1+p2).rap(),weight,100,-5,5);
      FillHist("pt_ac"+suffix,(p1+p2).pT(),weight,200,0.,400.);
    }
  } // End of event loop.

  // Show statistics.
  pythia.stat();

  // Normalize to cross section [fb].
  double sigmaNorm = pythia.info.sigmaGen() / pythia.info.weightSum() * 1.e12;
  
  TFile f("hists.root","update");
  for(auto [name,hist]:maphist){
    hist->Scale(sigmaNorm);
    hist->Write();
  }
  f.Close();

  // Done.
  return;
}

// The main program.
int main(int argc,char** argv) {
  if(argc>1) seed=atoi(argv[1]);
  
  vector<double> cuts={0.5,1,2,3,5,7,10,15,20,30};
  for(double cut:cuts){
    Run(cut);
  }
  return 0;
}

