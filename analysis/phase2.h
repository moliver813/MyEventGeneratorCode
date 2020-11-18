#ifndef PHASE2_H
#define PHASE2_H

#ifdef RE_STYLE
#include <readableStyle.h>
#endif

#define PI TMath::Pi()

#include <vector>
#include <math.h>
#include <cmath>
#include <iostream>
#include <cstdlib>

//#include <fastjet/ClusterSequence.hh>

#include <TFile.h>
#include <TStyle.h>
#include <TColor.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TF2.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TText.h>
#include <TLatex.h>

//#include "analysis_params.h"

using namespace std;


// Needed for peak fitting
struct minMax
{
  minMax(double _min, double _max): min(_min), max(_max) {}
  double min;
  double max;
};


//  int iP2TriggerType = 0; // 0 for JetH, 1 for Pi0H

// Analysis variables 
  extern int nJetPtBins;
  extern int nParticlePtBins;

  extern int nPi0PtBin;
  extern int nZtBins;


  // Define the pt bins
  // Size of the bins is from STAR Jet-h paper
  // Desired bins are 0 - 2 based on the convention below.
//std::vector <double> jetPtBins = {10,15,20,40};
  extern std::vector <double> jetPtBins;
  // Desired bins are 0 - 10 based on the convention below.
  extern std::vector <double> particlePtBins;

  // Pi0-h bins:
  extern std::vector <double> pi0PtBins;
  extern std::vector <double> zTBins; 
  // may want another particle pt bins for pi0 analysis
  //extern std::vector <double> pi0ParticlePtBins;


	// Minimum pt for all particles
  extern	double particle_pt_min;
	// Minimum pt for particles included in Jet Reconstruction
	// Warning: this can be changed in run-time by phase1.cc 
  extern double jet_constituent_cut;
	// Minimum pt for at least one constituent of a jet with a hard core
  // used for specifc analyses.  For a global hard core cut on jets, see
  // globalHardCoreCut
  extern double jet_hard_core_cut;
  extern double eta_cut;

  // jet cone radius / resolution parameter
	extern double R;

	extern double eta_range;
	extern double delta_eta_range;

  extern bool requireHardCore;
  extern double globalHardCoreCut;

  extern double c_width;
  extern double c_height;

  extern bool use2DFit;
  extern bool useSecHFit;
  extern bool useModGausFit;
  extern bool useZYAM;

  extern int bkg2DMethod;
    const int aNoBkgSub = 0;
    const int a2DFitSub = 1;
    const int aDEtaFitSub = 2;
    const int aDEtaFarSub = 3;
		const int aDEtaGenGausFitSub = 4;
		const int aDEtaMixedGausFitSub = 5;
    const int aDEtaFarZYAMSub = 6;
  extern int dPhiFitMethod;
    const int a1G1G = 0;
    const int a1M1M = 1;
    const int a2G1M = 2;
    const int a2G1GG = 3;
		const int a1GG1GG = 4;
    const int a2G1SH = 5;
		const int a2G2G = 6;
  extern int dPhiBkgMethod;
/*    extern int aNoDPhiBkg; // B = 0
    extern int aFreeDPhiBkg; // B is a free parameter
    extern int aZYAM;  // B is fixed by ZYAM
*/

    const int aNoDPhiBkg = 0; // B = 0
    const int aFreeDPhiBkg = 1; // B is a free parameter
    const int aZYAM = 2;  // B is fixed by ZYAM



	extern bool plotSigmaFits;



TF1 * fitPeak(TH1D * hist, minMax peak, std::string name);

TF1 * phiCorrFit(TH1D * hist, std::string name);

TF1 * phiCorrFit_2Gaus(TH1D * hist, std::string name);


TF1 * phiCorrFit_1Gaus_1Gaus(TH1D * hist, std::string name);
TF1 * phiCorrFit_1Gaus_1Gaus_Fwhm(TH1D * hist, std::string name);

TF1 * phiCorrFit_1ModGaus_1ModGaus(TH1D * hist, std::string name);
TF1 * phiCorrFit_1ModGaus_1ModGaus_Fwhm(TH1D * hist, std::string name);

TF1 * phiCorrFit_2Gaus_1ModGaus(TH1D * hist, std::string name);
TF1 * phiCorrFit_2Gaus_1ModGaus_Fwhm(TH1D * hist, std::string name);

TF1 * phiCorrFit_2Gaus_1GenGaus(TH1D * hist, std::string name);
TF1 * phiCorrFit_2Gaus_1GenGaus_Fwhm(TH1D * hist, std::string name);

TF1 * phiCorrFit_2Gaus_2Gaus(TH1D * hist, std::string name);
TF1 * phiCorrFit_2Gaus_2Gaus_Fwhm(TH1D * hist, std::string name);

TF1 * phiCorrFit_1GenGaus_1GenGaus(TH1D * hist, std::string name);
TF1 * phiCorrFit_1GenGaus_1GenGaus_Fwhm(TH1D * hist, std::string name);

TF1 * phiCorrFit_2Gaus_SecH(TH1D * hist, std::string name);
TF1 * phiCorrFit_2Gaus_SecH_Fwhm(TH1D * hist, std::string name);

/** Newer method: 2 gaussians for nearside, and a Voigt profile for the 
 * awayside
 *  To be coded
 */

TF2 * etaPhiCorr_fit(TH2F * hist, std::string name);

/**
 * Method that takes a dPhi projection and returns the TF2 that matches the 
 * combinatorial background of the unprojected fit.
 */
TF1 * dEta_fit(TH1D * hist, std::string name, double deltaEtaRange, bool doFit);


TF1 * dEta_GenGaus_fit(TH1D * hist, std::string name, bool doFit);


TF1 * swiftJetFit(TH1F * hist, std::string name);

TF1 * swiftParticleFit(TH1F * hist, std::string name);

TF2 * createSwiftBackground(TF1 *jetFit, TF1 * particleFit, std::string name);








void parseInput(int argc, char * argv[], TString & inputFilename, TString & outputFilename);

void printHelp(std::string input = "");



#endif
