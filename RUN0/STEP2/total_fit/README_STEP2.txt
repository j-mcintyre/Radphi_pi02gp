////////////////////////////
// Author: James McIntyre //
// Pi0 Analysis STEP #2   //
////////////////////////////
The C++ macro in STEP2 fits histograms from STEP1, prints PNGs, and creates a ROOT file with parameters
used in the subsequent analysis steps.

The fitted histograms (listed below) cover mandelstam |t| from 0.1 - 1.9 GeV^2 and are binned in
four photon energy ranges (E0 - E3):
   TH1D Allt_E0
   TH1D Allt_E1
   TH1D Allt_E2
   TH1D Allt_E3

//////////////////////////
// Changing |t| binning //
//////////////////////////
The |t| range and bins are listed on line 201 (BinM2g[20][2]). To change |t| binning for the 
analysis you need to change "BinM2g" in this file and "M2g_bins" (used to define histo names) and
"BinM2g" (used in histogram filling) in STEP1 TSelector.
Current range:
   double BinM2g[20][2] = {{0.100, 0.130}, {0.130, 0.160}, {0.160, 0.190}, {0.190, 0.230}, {0.230, 0.270},
                           {0.270, 0.310}, {0.310, 0.360}, {0.360, 0.410}, {0.410, 0.460}, {0.460, 0.535},
                           {0.535, 0.610}, {0.610, 0.685}, {0.685, 0.785}, {0.785, 0.885}, {0.885, 0.985},
                           {0.985, 1.135}, {1.135, 1.285}, {1.285, 1.460}, {1.460, 1.660}, {1.660, 1.900}};

The above |t| bin values are stored in the ROOT file as:
   TH1D hBinDW    // Lower bound
   TH1D hBinUP    // Upper bound

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

/////////////////////
// Fitting Process //
/////////////////////
Invariant mass peaks for each signal change with photon energy and |t| 
(particularly at high and low |t|). The following fitting steps are really
useful for individual |t| bin fits (STEP3), but are useful here too.
1) Create and fit two TF1 single gaussian fits, one for the range of the Pi0 and 
   one for the Eta. Save the mean values as a double.
2) Create TF1 total fit (ftot) from a combination of double gaussians (Pi0, Eta, and Omega),
   single gaussian (Eta prime), and a third order polynomial with a scaling factor (pol3 * Scaler).
3) Fix fit parameters to values "somewhat" around what they should be, except pol3 & scaler set them to 1.
   Too close values can cause the fitter to fail.
4) Run "fixed parameter" fit with options "RB".
5) Release all parameters, except the Scaler. Keep that set to 1.
6) Set parameter limits (SetParLimits). Use the mean values +/- some value for the Pi0 and 
   Eta mean limits. Set widths to a fairly wide range. Keep the widths for a single double gaussian 
   the same range.
7) Fit the histogram of interest with options "SBRM".
8) Save pol3 parameters in a ROOT file, along with the other values mentioned at the beginning of this file.
9) Print fitting histograms, create directories, and save/organize PNG's in the created subdirectories.
10)Once the macro completes a visual inspection of the fit histograms saved as PNG’s is required. Look for inappropriate fit. 

/////////////////////
// Troubleshooting //
/////////////////////
If the fitting completes, but gives an odd fit:
1) Write an IF statement for that particular histogram (i.e. Echannel == 1) and alter the fit limits.

If the fits fails try the following (in the order listed):
1) Remove "M" from the options in step (7) above, then try to run the macro again.
   The fitter sometimes gets stuck minimizing the fit with this option.
2) Try changing the initial fixed fit parameter values that are set in step 3 above.
   Values too close to the correct value can cause the fitter to fail. Use an IF Statement
   to change the initial values for the failed histogram only. Sometimes there is only one 
   histogram that is the problem child. This step usually fixes those fittings that fail 
   even though the fit looks really good (e.g. initial parameters are too close and 
   the fitter gets stuck while trying to minimize the fit).
3) If the previous step fails, add an IF statement and implement different parameter limits
   in step (6) above. Base your changes off of the output text file created from your
   previous fit (2g_E_binned_Pi0_fit_allt_*.txt). The IF statement is added to change the
   parameter limits for the failed histogram and those there after. Limit the number of elements
   you change at once, several attempts may be required.
4) If all of the above fails to fix things, have a look at the failed TH1D. The single gaussian fit
   for Pi0 or Eta might give a mean that's way off (i.e. from a signal too weak or too broad).
   Also, check the text file "2g_E_binned_Pi0_fit_allt_*.txt" for the actual fitted values.
   If these means are way off then your SetParLimit values will be out of range of the signal.
