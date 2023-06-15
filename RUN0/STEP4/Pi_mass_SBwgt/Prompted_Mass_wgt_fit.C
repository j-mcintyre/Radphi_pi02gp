///////////////////////////////////////////////////
// Author: J. McIntyre
// 4 June 2023
///////////////////////////////////////////////////
// This file is used to find a TSpline fit of mass
// sideband subtraction weights and save to a TFile.
// ------------------------------------------------

/*
   //----------------------
   // Add line to histogram
   //----------------------
   TLine *zeroline3 = new TLine(0.74, 0, 1.18, 0);
   zeroline3->SetLineColor(1);
   zeroline3->SetLineStyle(10);
   zeroline3->SetLineWidth(2);
   zeroline3->Draw("hist");

   //----------------
   // Drawing Options
   //----------------
   gStyle->SetOptStat(0);                 // No Legend
   gStyle->SetOptTitle(0);                // No title
   histo->SetLineWidth(5);                // Line thickness
   histo->SetLineStyle(5);                // Line type (1=solid)
   histo->SetLineColor(kBlue);            // Line color
   histo->SetFillColor(kBlue);            // Filling color to histo (38 = light blue)
   histo->SetFillColorAlpha(kBlue,0.25);  // Suppose to be semi-transparent
   //----------------
*/

// Regarding #include, using <> tells the compiler to search for the .h file in the system include directories,
// while using "" has it search for the .h file in the current directory and then in the system include directories.
#include <iostream>
#include <TF1.h>
#include <cmath>
#include <TLegend.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TVirtualFitter.h>
#include <Math/MinimizerOptions.h>
#include <Math/Minimizer.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <fstream>
#include <math.h>
#include <TMath.h>
#include <string>
#include <TString.h>
#include <stdlib.h>
#include <Math/Functor.h>
#include <TSystem.h>
#include <TROOT.h>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <stdio.h>
#include <TPad.h>
#include <TText.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TColor.h>
#include <TAxis.h>
#include <TObjString.h>
#include <TList.h>
#include <Math/RootFinderAlgorithms.h>
#include <Math/RootFinder.h>
#include <Math/GSLMultiRootFinder.h>
#include <Math/WrappedMultiTF1.h>
#include <TKey.h>
#include <TGraphPainter.h>
#include <TAttFill.h>
#include <TSpline.h>

/*
#include <TLorentzVector.h>
#include <TVector3.h>
#include <vector>
#include <TH2.h>
#include <TF2.h>
#include <TGraphSmooth.h>

#include <TArrow.h>
#include <TPave.h>
#include <TPaveText.h>
#include <TApplication.h>
#include <TRint.h>
#include <TDirectory.h>
#include <TChain.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <TMarker.h>
#include "TRandom3.h"
*/

// using namespace std;

#define PI 3.14159265359
#define sigmaRANGE 2.25

TH1D *remass;
TF1  *fsig1;
TF1  *fsig2;
TF1  *fsig3;
TF1  *fsig4;
TF1  *fback;
TF1  *ftot;
TH1D *h_W0;
TH1D *h_W1;
TH1D *h_W2;

void Prompted_Mass_wgt_fit() {

   /*
	double tbins[20][2] = {{0.100, 0.130}, {0.130, 0.160}, {0.160, 0.190}, {0.190, 0.230}, {0.230, 0.270},
                            {0.270, 0.310}, {0.310, 0.360}, {0.360, 0.410}, {0.410, 0.460}, {0.460, 0.535},
                            {0.535, 0.610}, {0.610, 0.685}, {0.685, 0.785}, {0.785, 0.885}, {0.885, 0.985},
                            {0.985, 1.135}, {1.135, 1.285}, {1.285, 1.460}, {1.460, 1.660}, {1.660, 1.900}};
   */

   // STEP2 ROOT file
   TFile *fStep2 = new TFile("/nfs/direct/annex/mcintyre/GitRepo/Radphi_pi02gp/RUN0/STEP2/total_fit/TwoG_may_12_2023_STEP2.root");

   // Importing |t| bin values
   TH1D *hBinDW  = (TH1D*)fStep2->Get("hBinDW");
   TH1D *hBinUP  = (TH1D*)fStep2->Get("hBinUP");

   // Import number of |t| bins
   int tBinNum = hBinUP->GetNbinsX();

   double tbins[tBinNum][2];
   for (int b1 = 0; b1 < tBinNum; b1++) {
      int bs = b1 + 1;
      tbins[b1][0] = hBinDW->GetBinContent(bs);
      tbins[b1][1] = hBinUP->GetBinContent(bs);
   }

   // Import number of photon energy bins used
   TH1D *EChanBINS = (TH1D*)fStep2->Get("EChanBINS");
   int NumEChan = EChanBINS->GetBinContent(1);

   fStep2->Close();

   // Open STEP3 ROOT file to import sideband subtraction weights and errors
   TFile *fStep3 = new TFile("/nfs/direct/annex/mcintyre/GitRepo/Radphi_pi02gp/RUN0/STEP3/indiv_fits/TwoG_may_12_2023_STEP3.root");

   // Total number of bins (e.g. weights to import)
   int totalBins = NumEChan * tBinNum;
   double X_pi_sig_t[totalBins];
   double X_err_sig_t[totalBins];
   double X_pi_sb1_t[totalBins];
   double X_err_sb1_t[totalBins];
   double X_pi_sb2_t[totalBins];
   double X_err_sb2_t[totalBins];

   double pi_sig_t[tBinNum];
   double err_sig_t[tBinNum];
   double pi_sb1_t[tBinNum];
   double err_sb1_t[tBinNum];
   double pi_sb2_t[tBinNum];
   double err_sb2_t[tBinNum];

   int WeightCount = 0; 
   for (int i = 0; i < NumEChan; i++) {
      // Histogram from STEP3
      TString hist0Name = Form("h_W0_E%d", i);
      TH1D *h_W0 = (TH1D*)fStep3->Get(hist0Name);

      // Histogram from STEP3
      TString hist1Name = Form("h_W1_E%d", i);
      TH1D *h_W1 = (TH1D*)fStep3->Get(hist1Name);

      // Histogram from STEP3
      TString hist2Name = Form("h_W2_E%d", i);
      TH1D *h_W2 = (TH1D*)fStep3->Get(hist2Name);

      for (int j = 0; j < tBinNum; j++) {
         WeightCount = i * tBinNum + j;
         int bbq = 2 * j + 1;
	    	int acq = bbq + 1;
         X_pi_sig_t[WeightCount]  = h_W0->GetBinContent(bbq);
         X_err_sig_t[WeightCount] = h_W0->GetBinContent(acq);

         X_pi_sb1_t[WeightCount]  = h_W1->GetBinContent(bbq);
         X_err_sb1_t[WeightCount] = h_W1->GetBinContent(acq);

         X_pi_sb2_t[WeightCount]  = h_W2->GetBinContent(bbq);
         X_err_sb2_t[WeightCount] = h_W2->GetBinContent(acq);
      }
   }
   fStep3->Close();

   // Create a ROOT file for further analysis
  	TFile *MyFile = new TFile("TwoG_may_12_2023_STEP4.root", "NEW");
   MyFile->Close();

	// TSpline nodes
	int n = tBinNum;

   double PARTmeson_X[n];
   double w_sig_t[n];
   double w_err_sig_t[n];
   double w_sb1_t[n];
   double w_err_sb1_t[n];
   double w_sb2_t[n];
   double w_err_sb2_t[n];
   double meson_X[tBinNum];
   double ex[tBinNum];
	TString BkdgNAME;
	TString BkN;
   TString BkF;
   TString ENAME;
   TString EFILE;
	TString rangeNAME;
	TString ErrorNAME;
   TString hSig;
   TString hspline;
   TString hspline2;
   TString hxx;
   TString hSB;
   TString hSigpng;
   TString hSBpng;

   // TCanvas setup for capturing png files
   TCanvas *c101 = new TCanvas("c101", "", 800, 600);
   TCanvas *c102 = new TCanvas("c102", "", 800, 600);
   TCanvas *c103 = new TCanvas("c103", "", 800, 600);
   TCanvas *c105 = new TCanvas("c105", "", 800, 600);
   TCanvas *c106 = new TCanvas("c106", "", 800, 600);
   TCanvas *c107 = new TCanvas("c107", "", 800, 600);

   int HistBinCount = 0;
   for(int Echannel = 0; Echannel < NumEChan; Echannel++){

      int mi = 0;
      int NodeCount = 0;
		for(int i = 0; i < tBinNum; i++) {
         // Corresponding histogram bin entry
         HistBinCount = Echannel * tBinNum + i;

			double YY = tbins[i][1] - tbins[i][0];
			meson_X[i] = tbins[i][0] + (YY/2.0);
			ex[i] = YY/2.0;

         /* vvvvvvvv PART USED TO PLOT A TSPLINE vvvvvvvv */
         if (((Echannel == 0) && (i == 0 || i == 1 || i == 2 || i == 3 || i == 4 || i == 8 || i == 9 || i == 12 || i == 15 || i == 18 || i == 19)) ||
             ((Echannel == 1) && (i == 0 || i == 2 || i == 4 || i == 7 || i == 4 || i == 11 || i == 12 || i == 16 || i == 18 || i == 19)) ||
				 ((Echannel == 2) && (i == 0 || i == 1 || i == 2 || i == 4 || i == 9 || i == 10 || i == 11 || i == 14 || i == 16 || i == 17 || i == 19)) ||
				 ((Echannel == 3) && (i == 0 || i == 3 || i == 5 || i == 10 || i == 13 || i == 14 || i == 15 || i == 17 || i == 19))) {
				double ZZ = tbins[i][1] - tbins[i][0];
				PARTmeson_X[mi] = tbins[i][0] + (ZZ / 2.0);

            // Get a subset of bins to fit with a TSpline
				w_sig_t[mi] = X_pi_sig_t[HistBinCount];
				w_err_sig_t[mi] = X_err_sig_t[HistBinCount];

				w_sb1_t[mi]  = X_pi_sb1_t[HistBinCount];
				w_err_sb1_t[mi] = X_err_sb1_t[HistBinCount];

				w_sb2_t[mi]  = X_pi_sb2_t[HistBinCount];
				w_err_sb2_t[mi] = X_err_sb2_t[HistBinCount];

            cout << "w0     = " << w_sig_t[mi] << endl;
            cout << "w0 +/- = " << w_err_sig_t[mi] << endl;
            cout << "w1     = " << w_sb1_t[mi] << endl;
            cout << "w1 +/- = " << w_err_sb1_t[mi] << endl;
            cout << "w2     = " << w_sb2_t[mi] << endl;
            cout << "w2 +/- = " << w_err_sb2_t[mi] << endl;
            cout << "|t|    = " << PARTmeson_X[mi] << endl;
            cout << "----------------------" << endl;
            mi += 1;
            NodeCount += 1;
         }
         /* ^^^^^^^^ PART USED TO PLOT A TSPLINE ^^^^^^^^ */
		}

      /* vvvvvvvv PART USED TO PLOT A TSPLINE vvvvvvvv */ 
		// Modify location (y-axis) of the TSpline nodes.
      // This is a nudge in the right direction for a better TSpline fit.
      /* [] NUMBER IS NOT |t| BIN NUMBER, BUT THE "mi" NUMBER */
      if (Echannel == 0) {
   		//w_sig_t[0] -= 0.01;

   		w_sb2_t[0] = -0.94;
   		//w_sb1_t[6] += 0.002;

         // Sets the location of the first node's X to just below the low end of the bin.
         // Otherwise the first node would be in the center of the first bin
         // and create a sharp curve in the TSpline.
         PARTmeson_X[0] = 0.09;
      }
      else if (Echannel == 1) {
   		//w_sig_t[0] += 0.005;
   		//w_sig_t[18] -= 0.003;
   		w_sb1_t[0] -= 0.00025;
   		w_sb1_t[18] += 0.0003;
   		//w_sb1_t[6] += 0.005;

         //PARTmeson_X[0] = 0.09;
      }
      else if (Echannel == 2) {
   		w_sig_t[0] = 1.025;
   		w_sig_t[9] += 0.0025;
   		//w_sig_t[2] += 0.0035;
   		//w_sig_t[3] -= 0.008;
   		//w_sig_t[5] += 0.0015;
   		//w_sig_t[6] -= 0.0015;
   		//w_sig_t[7] += 0.0005;
   		//w_sig_t[8] += 0.0005;

   		w_sb1_t[0] = -0.1045;
   		//w_sb1_t[1] -= 0.0003;
   		//w_sb1_t[2] -= 0.007;
   		w_sb2_t[0] = -0.905;
   		//w_sb1_t[6] += 0.005;

         PARTmeson_X[0] = 0.09;
      }
      else if (Echannel == 3) {
   		w_sig_t[0] = 1.022;
   		w_sig_t[1] -= 0.001;
   		//w_sig_t[5] -= 0.006;
   		//w_sig_t[6] -= 0.006;
   		//w_sig_t[7] -= 0.004;
   		//w_sig_t[11] += 0.008;

   		//w_sb1_t[0] = -1.18;

   		//w_sb1_t[0] = -1.28;
   		w_sb1_t[7] += 0.0002;
   		//w_sb1_t[2] += 0.003;
   		//w_sb1_t[3] += 0.003;
   		//w_sb1_t[4] -= 0.002;
   		//w_sb1_t[5] += 0.006;
   		//w_sb2_t[6] += 0.002;
   		//w_sb1_t[7] -= 0.0035;

         PARTmeson_X[0] = 0.09;
      }

      // TSpline nodes
      // ALLOW DIFFERENT NUMBER OF TSPLINE NODES
      n = NodeCount;


		c101->cd();
		TGraph *grSig = new TGraph(n, PARTmeson_X, w_sig_t);
		grSig->SetTitle("#pi^{0} Signal Region Weights;|t| [GeV^2/c^{4}];Fill Weight");
		//grSig->GetXaxis()->CenterTitle(true);
		grSig->GetXaxis()->SetTitle("|t| [GeV^2/c^{4}]");
		grSig->GetYaxis()->CenterTitle(true);
		grSig->GetYaxis()->SetTitle("Fill Weight");
		grSig->Draw("ALP");

		TSpline3 *s = new TSpline3("",grSig);
		s->SetLineColor(kRed);
		s->Draw("same");

      // Update ROOT file for further analysis
      TFile *MyFile = new TFile("TwoG_may_12_2023_STEP4.root", "UPDATE");

		gFile = MyFile;
      hspline.Form("w0_spline_E%d", Echannel);
		gDirectory->WriteObject(s, hspline);


		c105->cd();
		TGraph *grSigw1 = new TGraph(n, PARTmeson_X, w_sb1_t);
		grSigw1->SetTitle("#pi^{0} SB Region 1 Weights (SB1);|t| [GeV^2/c^{4}];Fill Weight");
		//grSigw1->GetXaxis()->CenterTitle(true);
		grSigw1->GetXaxis()->SetTitle("|t| [GeV^2/c^{4}]");
		grSigw1->GetYaxis()->CenterTitle(true);
		grSigw1->GetYaxis()->SetTitle("Fill Weight");
		grSigw1->Draw("ALP");

		TSpline3 *st = new TSpline3("",grSigw1);
		st->SetLineColor(kRed);
		st->Draw("same");

		gFile = MyFile;
      hspline.Form("w1_spline_E%d", Echannel);
		gDirectory->WriteObject(st, hspline);

		c106->cd();
		TGraph *grSigw2 = new TGraph(n, PARTmeson_X, w_sb2_t);
		grSigw2->SetTitle("#pi^{0} SB Region 2 Weights (SB2);|t| [GeV^2/c^{4}];Fill Weight");
		//grSigw2->GetXaxis()->CenterTitle(true);
		grSigw2->GetXaxis()->SetTitle("|t| [GeV^2/c^{4}]");
		grSigw2->GetYaxis()->CenterTitle(true);
		grSigw2->GetYaxis()->SetTitle("Fill Weight");
		grSigw2->Draw("ALP");

		TSpline3 *st2 = new TSpline3("",grSigw2);
		st2->SetLineColor(kRed);
		st2->Draw("same");

		gFile = MyFile;
      hspline.Form("w2_spline_E%d", Echannel);
		gDirectory->WriteObject(st2, hspline);
      /* ^^^^^^^^ PART USED TO PLOT A TSPLINE ^^^^^^^^ */

      for (int wrt = 0; wrt < tBinNum; wrt++) {
         int placeholder = Echannel * tBinNum + wrt;
         pi_sig_t[wrt] = X_pi_sig_t[placeholder];
         err_sig_t[wrt] = X_err_sig_t[placeholder];
         pi_sb1_t[wrt] = X_pi_sb1_t[placeholder];
         err_sb1_t[wrt] = X_err_sb1_t[placeholder];
         pi_sb2_t[wrt] = X_pi_sb2_t[placeholder];
         err_sb2_t[wrt] = X_err_sb2_t[placeholder];
      }

      int fullN = tBinNum;
		c103->cd();
		TGraphErrors *gr3 = new TGraphErrors(fullN, meson_X, pi_sig_t, ex, err_sig_t);
		gr3->SetTitle("#pi^{0} Signal Region Weights;|t| [GeV^2/c^{4}];Fill Weight");
		gr3->GetYaxis()->CenterTitle(true);
		gr3->Draw("ALP");

		hSigpng.Form("Pi0_w0_E%d.png", Echannel);
      c103->Print(hSigpng);

      s->Draw("same");

		hSig.Form("Pi0_w0_spline_E%d.png", Echannel);
      c103->Print(hSig);

		gFile = MyFile;
		hSB.Form("Pi0_w0_E%d", Echannel);
		gDirectory->WriteObject(gr3, hSB);


		c102->cd();
		TGraphErrors *gr2 = new TGraphErrors(fullN,meson_X,pi_sb1_t,ex,err_sb1_t);
		gr2->SetTitle("#pi^{0} SB Region 1 Weights (SB1);|t| [GeV^2/c^{4}];Fill Weight");
		gr2->GetYaxis()->CenterTitle(true);
		gr2->Draw("ALP");

		hSigpng.Form("Pi0_w1_E%d.png", Echannel);
      c102->Print(hSigpng);

		st->Draw("same");

		hSig.Form("Pi0_w1_spline_E%d.png", Echannel);
      c102->Print(hSig);

		gFile = MyFile;
		hSB.Form("Pi0_w1_E%d", Echannel);
		gDirectory->WriteObject(gr2, hSB);


		c107->cd();
		TGraphErrors *gr4 = new TGraphErrors(fullN, meson_X, pi_sb2_t, ex, err_sb2_t);
		gr4->SetTitle("#pi^{0} SB Region 2 Weights (SB2);|t| [GeV^2/c^{4}];Fill Weight");
		gr4->GetYaxis()->CenterTitle(true);
		gr4->Draw("ALP");

		hSigpng.Form("Pi0_w2_E%d.png", Echannel);
      c107->Print(hSigpng);

		st2->Draw("same");

		hSig.Form("Pi0_w2_spline_E%d.png", Echannel);
      c107->Print(hSig);

		gFile = MyFile;
		hSB.Form("Pi0_w2_E%d", Echannel);
		gDirectory->WriteObject(gr4, hSB);

      MyFile->Close();
   }
   delete c101;
   delete c102;
   delete c103;
   delete c105;
   delete c106;
   delete c107;

   // Organize info into folders
   TString FolderN;
   for (int i = 0; i < NumEChan; i++) {
      FolderN.Form("mkdir -p histos/E%d", i);
      gSystem->Exec(FolderN);
   }
   gSystem->Exec("mkdir -p histos/Spline");
   gSystem->Exec("mv *_spline*.png histos/Spline/.");

   gSystem->Exec("mv *E0.png histos/E0/.");
   gSystem->Exec("mv *E1.png histos/E1/.");
   gSystem->Exec("mv *E2.png histos/E2/.");
   gSystem->Exec("mv *E3.png histos/E3/.");
}
