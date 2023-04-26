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

// The main program.
int main(int argc, char **argv){
  if(argc!=2) return 1;
  // Four possible contributions: double-dissociative = 1;
  // single-dissociative on side A (B) = 2 (3); elastic-elastic = 4.
  int processType = atoi(argv[1]);  

  // Generator.
  Pythia pythia;

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
  pythia.readString("PDF:pSet = 13");


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
  pythia.readString("PhaseSpace:pTHatMinDiverge = 1.");

  //seed
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");

  // Initialize the histogram for acoplanarity.
  vector<double> mbins={50,60,80,100,150,200,400,1000,3000};
  TH2D mpt("mpt","mpt",mbins.size()-1,&mbins[0],300,0.,600.);
  TH2D mpt_ac("mpt_ac","mpt_ac",mbins.size()-1,&mbins[0],300,0.,600.);

  // Initialize the generator.
  pythia.init();

  // Number of events.
  int nEvent = 10000;

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

    // Fill histrogram with possible weights.
    mpt.Fill(m(p1,p2),(p1+p2).pT(),weight);
    if(TMath::Max(p1.pT(),p2.pT())>20&&TMath::Min(p1.pT(),p2.pT())>10&&fabs(p1.eta())<2.5&&fabs(p2.eta())<2.5)
      mpt_ac.Fill(m(p1,p2),(p1+p2).pT(),weight);
  } // End of event loop.

  // Show statistics.
  pythia.stat();

  // Normalize to cross section [fb].
  double sigmaNorm = pythia.info.sigmaGen() / pythia.info.weightSum() * 1.e12;
  mpt.Scale(sigmaNorm);
  mpt_ac.Scale(sigmaNorm);
  
  TFile f("hists.root","recreate");
  mpt.Write();
  mpt_ac.Write();
  f.Close();

  // Done.
  return 0;
}
