////////////////////////////
// Author: James McIntyre //
// Pi0 Analysis STEP #5   //
////////////////////////////
The TSelector of STEP5 applies various cuts to the 2 gamma dataset (2578 ROOT files), using PROOF,
to reduce background and apply mass sideband subtraction weights when filling the main histograms of interest (plots of delta phi in various photon energy and |t| bins). These plots are saved to a ROOT file
with the name <TSelector name>Analysis.root. Additional histograms are also saved to this ROOT file,
as mentioned below, for use in future STEPs or for plots to be placed in my thesis. Subsequent STEPs
will further reduce the background using delta phi sideband subtraction and then determine the Pi0 cross-section.



The t| bin values are stored in STEP2 ROOT file as:
   TH1D hPiBinDW    // Lower bound
   TH1D hPiBinUP    // Upper bound

////////////////////////////////////
// Changing Photon Energy Binning //
////////////////////////////////////
Number of Energy Bins (currently 4):
   TH1D EChanBINS

Tagger Channels used in each Photon Energy Bin (Lower & Upper Channels are listed):
   TH1D EChanRange

Histogram with corresponding photon energies (GeV) for the range listed above("EChanRange"):
   TH1D EChanValue

For reference, the Hall B tagger channel energies (GeV) used are:
   // High energy end of bin in GeV
   // Tagchan# energy value is border between tagger channel # & (# - 1)
   double Tagchan0 = 5.385;
   double Tagchan1 = 5.335;
   double Tagchan2 = 5.280;
   double Tagchan3 = 5.220;
   double Tagchan4 = 5.165;
   double Tagchan5 = 5.115;
   double Tagchan6 = 5.065;
   double Tagchan7 = 5.015;
   double Tagchan8 = 4.965;
   double Tagchan9 = 4.915;
   double Tagchan10 = 4.860;
   double Tagchan11 = 4.800;
   double Tagchan12 = 4.740;
   double Tagchan13 = 4.675;
   double Tagchan14 = 4.615;
   double Tagchan15 = 4.570;
   double Tagchan16 = 4.530;
   double Tagchan17 = 4.485;
   double Tagchan18 = 4.435;
   double Tagchan19 = 4.385;

Mass values of the sideband subtraction regions:
   TH1D *Pi0_SB_Regions_allt_E0;
   TH1D *Pi0_SB_Regions_allt_E1;
   TH1D *Pi0_SB_Regions_allt_E2;
   TH1D *Pi0_SB_Regions_allt_E3;


Previous analysis procedure fitted TSplines to the various invariant mass sideband subtraction weights.
These TSplines come from STEP4's root file:
SIGNAL REGION
   TSpline *w0_spline_E0;
   TSpline *w0_spline_E1;
   TSpline *w0_spline_E2;
   TSpline *w0_spline_E3;

LOWER MASS SIDEBAND REGION
   TSpline *w1_spline_E0;
   TSpline *w1_spline_E1;
   TSpline *w1_spline_E2;
   TSpline *w1_spline_E3;

HIGHER MASS SIDEBAND REGION
   TSpline *w2_spline_E0;
   TSpline *w2_spline_E1;
   TSpline *w2_spline_E2;
   TSpline *w2_spline_E3;


This TSelector saves the following histograms to an output ROOT file:
   //////////// E0 BIN ////////////
   |t| from E_{#gamma} = 5.115 to 5.385 GeV for 0.1 <= |t| < 1.9 GeV^2
      TH1D *tPi0;
   #Delta#phi from E_{#gamma} = 5.115 to 5.385 GeV for 0.1 <= |t| < 1.9 GeV^2
      TH1D *dphi_Pi0;
   #Delta#phi from E_{#gamma} = 5.115 to 5.385 GeV for individual |t| bins
   TH1D *Pi_dphi_E0_t[20];

   //////////// E1 BIN ////////////
   |t| from E_{#gamma} = 4.860 to 5.115 GeV for 0.1 <= |t| < 1.9 GeV^2
      TH1D *tPi1;
   #Delta#phi from E_{#gamma} = 4.860 to 5.115 GeV for 0.1 <= |t| < 1.9 GeV^2
      TH1D *dphi_Pi1;
   #Delta#phi from E_{#gamma} = 4.860 to 5.115 GeV for individual |t| bins
   TH1D *Pi_dphi_E1_t[20];

   //////////// E2 BIN ////////////
   |t| from E_{#gamma} = 4.615 to 4.860 GeV for 0.1 <= |t| < 1.9 GeV^2
      TH1D *tPi2;
   #Delta#phi from E_{#gamma} = 4.615 to 4.860 GeV for 0.1 <= |t| < 1.9 GeV^2
      TH1D *dphi_Pi2;
   #Delta#phi from E_{#gamma} = 4.615 to 4.860 GeV for individual |t| bins
   TH1D *Pi_dphi_E2_t[20];

   //////////// E3 BIN ////////////
   |t| from E_{#gamma} = 4.358 to 4.615 GeV for 0.1 <= |t| < 1.9 GeV^2
      TH1D *tPi3;
   #Delta#phi from E_{#gamma} = 4.358 to 4.615 GeV for 0.1 <= |t| < 1.9 GeV^2
      TH1D *dphi_Pi3;
   #Delta#phi from E_{#gamma} = 4.358 to 4.615 GeV for individual |t| bins
   TH1D *Pi_dphi_E3_t[20];
