// phase2 file
// Take merged histograms, plot, fit 
// Create sigma, D_AA, R_AA plots
// 
// Can compile with:
//  make phase2
//
//  or
//
//  g++ -Wl,--no-as-needed -std=c++11 phase2.cc -o phase2 -I`root-config --incdir` `root-config --libs` `fastjet-config --cxxflags --libs --plugins`
//
//  or (if you have readableStyle.h)
//
//  make USER_DEFINED=-DRE_STYLE=1 basicAnalysis

#ifdef RE_STYLE
#include <readableStyle.h>
#endif

#define PI TMath::Pi()

#include <vector>
#include <math.h>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <iomanip>


#include <fastjet/ClusterSequence.hh>

#include <TFile.h>
#include <TStyle.h>
#include <TLine.h>

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
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TText.h>
#include <TLatex.h>
#include <TSystem.h>

#include "phase2.h"
#include "analysis_params.h"
//#include "fitAlgos.h"

using namespace std;

// If nonempty, a root file will be ready from this path to get the systematic 
// uncertainties on yields, rms, etc.
TString SystematicUncertInFilePath = "";

// If nonempty, a root file will be created at this path to store the systematic
// uncertainties, which will be calculated by comparing different event planes
// This should only be run on sims with no EP dependence
TString SystematicUncertOutFilePath = "";

TString outputDirPath = "output/";

void set_plot_style() {
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}



void normalizeTH1F(TH1F * hist) {
	double integral = hist->Integral("width");
	if (integral) hist->Scale(1./integral);
}


// General method to return a 2D-histogram whose columns have been normalized 
// to 1 (* bin size?)

TH2F * normalizeHistogramColumns(TH2F *inHist) {

  TH2F * normHist = (TH2F *) inHist->Clone(Form("%sNorm",inHist->GetName()));
  normHist->SetTitle(Form("%s (Normalized)",inHist->GetTitle()));
  TH2F * scaleHist = (TH2F *) inHist->Clone(Form("%sScale",inHist->GetName()));
  scaleHist->SetTitle(Form("%s (Scale)",inHist->GetTitle()));
  scaleHist->Sumw2(0);
  int nbinsx = inHist->GetNbinsX();
  int nbinsy = inHist->GetNbinsY();

  for (int i = 0; i < nbinsx+1; i++) {
    double integral = inHist->Integral(i,i,0,nbinsy+1,"");
  //  double integral = inHist->Integral(i,i,0,nbinsy+1,"width");
//    double integral = inHist->Integral(i,i+1,0,nbinsy+1,"width");
    integral = integral ? 1./integral : 0; // lol
    for (int j = 0; j < nbinsy+1; j++) {
      scaleHist->SetBinContent(i,j,integral);
      scaleHist->SetBinError(i,j,0);
    }
  }

  normHist->Multiply(scaleHist);
  delete scaleHist;
  return normHist;
}


TH1F * convolveHistos(TH1F * hist1, TH1F * hist2, TString name) {
	double hist1min = hist1->GetXaxis()->GetXmin();
	double hist1max = hist1->GetXaxis()->GetXmax();
	double hist2min = hist2->GetXaxis()->GetXmin();
	double hist2max = hist2->GetXaxis()->GetXmax();
	
	int hist1Bins = hist1->GetNbinsX();
	int hist2Bins = hist2->GetNbinsX();
	
	double hist1BinSize = hist1->GetXaxis()->GetBinWidth(1);
	double hist2BinSize = hist2->GetXaxis()->GetBinWidth(1);

	double convMax = ((hist1max - hist1min) + (hist2max - hist2min)) / 2.;
	double convMin = - convMax;

// guess for now FIXME
//	int convBins = hist1Bins + hist2Bins - 1;
//	int convBins = (hist1Bins + hist2Bins)/2 - 1;
	int convBins = (hist1Bins + hist2Bins - 1)*4;
//	int convBins = hist1Bins + hist2Bins;

//	TH1F * conv = (TH1F *) hist1->Clone(name);
	
	TH1F * conv = new TH1F(name,name,convBins,convMin,convMax);


	printf("CONVOLVE: h1 bins: %d   h2 bins: %d    Hc bins %d \n",hist1Bins,hist2Bins,convBins);
	printf("CONVOLVE: h1 Range: (%.2f  -   %.2f)\n",hist1min,hist1max);
	printf("CONVOLVE: h2 Range: (%.2f  -   %.2f)\n",hist2min,hist2max);
	printf("CONVOLVE: hC Range: (%.2f  -   %.2f)\n",convMin,convMax);


	// Method 1, fill from histograms (possible aliasing issues
/*
	// loop over bins in hist1 
		for (int j = 1; j <= hist2Bins; j++) {
	for (int i = 1; i <= hist1Bins; i++) { // remember that root histograms have 0 = underflow,  hist1Bins+1 = overflow
		double hist1BinContent = hist1->GetBinContent(i);
	//	double hist1BinCenter = hist1->GetBinCenter(i);
		double hist1BinCenter = hist1->GetXaxis()->GetBinCenter(i);

		// loop over bins in hist1
			double hist2BinContent = hist2->GetBinContent(j);
		//	double hist2BinCenter = hist2->GetBinCenter(j);
			double hist2BinCenter = hist2->GetXaxis()->GetBinCenter(j);
			
		printf("CONVOLVE: h2(%.3f,%.3f) - h1(%.3f,%.3f) = (%.3f,%.3e)\t\t",hist2BinCenter,hist2BinContent,hist1BinCenter,hist2BinContent,hist2BinCenter - hist1BinCenter,hist1BinContent * hist2BinContent);		
		printf("h2(%d) - h1(%d) = Hc(%d)\n",j,i,hist2->FindBin(hist2BinCenter - hist1BinCenter));
//		printf("h2(%d) - h1(%d) = Hc(%d)\n",j,i,hist2->FindFixBin(hist2BinCenter - hist1BinCenter));

	
			conv->Fill(hist2BinCenter - hist1BinCenter, hist1BinContent * hist2BinContent);
//			conv->Fill(hist2BinCenter - hist1BinCenter,1);
		}
  }
*/

	// Method 2: proper integral
	for (int i = 1; i <= convBins; i++) {
		double ConvBinCenter = conv->GetBinCenter(i); 
		// ConvBinCenter is delta
		double sum = 0;
		// sum of weights to get uncertainty right
		// double sumW = 0;

		for (int j = 1; j <= hist1Bins; j++) {
			int hist2Bin = hist2->FindBin(ConvBinCenter - hist1->GetBinCenter(j));  
			//check for overflow/underflow?

	//		sum += hist1->GetBinContent(j) * hist2->GetBinContent(hist2Bin);
			conv->Fill(ConvBinCenter,hist1->GetBinContent(j) * hist2->GetBinContent(hist2Bin));


		}
		// conv->Fill(ConvBinCenter,sum);
		// or 
		// conv->SetBinContent(i,sum); 
		// then do something with uncertainties

	}


//	conv->Rebin(3);

	return conv;
}

/** Calculate Fixed v2 from normalize histogram, via fitting.
  * The phase of the curve is fixed.
	* What did fixed mean? And who cares about the phase?
  */
//double calculateFixedV2AndError(TH1F * phiProj, double &v2Error) {
void calculateVnAndError(TH1F * phiProj, double &v2, double &v2Error, double &v4, double &v4Error, double &v6, double &v6Error) {

	printf("Doing a Vn Calculation...\n");
//	TF1 * fit = new TF1("fit","[0]*(1 + 2.*[1]*cos(2.*x)+2*[2]*cos(4.*x))",phiProj->GetXaxis()->GetXmin(),phiProj->GetXaxis()->GetXmax());
	TF1 * fit = new TF1("fit","[0]*(1+2.*[1]*cos(2.*x)+2*[2]*cos(4.*x)+2.*[3]*cos(6.*x))",phiProj->GetXaxis()->GetXmin(),phiProj->GetXaxis()->GetXmax());

	double fB_Guess = phiProj->Integral() / phiProj->GetNbinsX();

	fit->SetParameter(0,fB_Guess);
	//fit->SetParameter(0,0.5/PI);
	fit->SetParameter(1,0.1);
	fit->SetParameter(2,0.05);
	fit->SetParameter(3,0.01);
	phiProj->Fit(fit);

	v2 = fit->GetParameter(1);
	v2Error = fit->GetParError(1);
	v4 = fit->GetParameter(2);
	v4Error = fit->GetParError(2);
	v6 = fit->GetParameter(3);
	v6Error = fit->GetParError(3);
}




// FIXME make this work?
TF2 * dEtaToDPhiDEta(TH1F * hist, TString name) {
	TF2 * func = new TF2(name,"0*x + 0*y",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax(),-PI/2,3*PI/2);
	

	return func;
}	


/**
  * Loads systematic uncertainty tgraphs from input file applies them to the given matrix of tgraphs
  * Adds this systematic uncertainty to the present statistical uncertainty in quadrature. 
  * This systematic uncertainty is assumed to be uncorrelated 
  */
void	ApplySystematicUncertainty(vector <vector<TGraphErrors * >> tGraphs, TFile * f){ 
	// tGraphs expected to be indexed by k = 0,1, .., nEPBins, then j = 0,1, .. nJetPtBins

	for (int k = 0; k < nEPBins; k++) {
		for (int i = 0; i < nJetPtBins; i++) {
			TGraph * SysErr = (TGraph *) f->Get(Form("%s_SysUncert",tGraphs[k][i]->GetName()));		
			if (!SysErr) {
				fprintf(stderr,"%s_SysUncert Not Found in %s!\n",tGraphs[k][i]->GetName(),f->GetName());
				return;
			}
			int n = SysErr->GetN();
			// Cycle through points
			for (int j = 0; j < n; j++) {
				TGraphErrors * targetGraph = tGraphs[k][i];
				double x,y;
				targetGraph->GetPoint(j,x,y);
				double sysUn;
				SysErr->GetPoint(j,x,sysUn);
	
				double newError = sqrt(pow(targetGraph->GetErrorY(j),2) + pow(y*sysUn,2));
				targetGraph->SetPointError(j,targetGraph->GetErrorX(j),newError);
			}
		}
	}
}


/**
	* Similar to above, but intended for event plane ratios and differences (so no EP axis)
  * Loads systematic uncertainty tgraphs from input file applies them to the given matrix of tgraphs
  * Adds this systematic uncertainty to the present statistical uncertainty in quadrature. 
  */
void	ApplySystematicUncertainty(vector<TGraphErrors * > tGraphs, TFile * f){ 
	// tGraphs expected to be indexed by k = 0,1, .., nEPBins, then j = 0,1, .. nJetPtBins

	for (int i = 0; i < nJetPtBins; i++) {
		TGraph * SysErr = (TGraph *) f->Get(Form("%s_SysUncert",tGraphs[i]->GetName()));		
		if (!SysErr) {
			fprintf(stderr,"%s_SysUncert Not Found in %s!\n",tGraphs[i]->GetName(),f->GetName());
			return;
		}
		int n = SysErr->GetN();
		// Cycle through points
		for (int j = 0; j < n; j++) {
			TGraphErrors * targetGraph = tGraphs[i];
			double x,y;
			targetGraph->GetPoint(j,x,y);
			double sysUn;
			SysErr->GetPoint(j,x,sysUn);

			double newError = sqrt(pow(targetGraph->GetErrorY(j),2) + pow(y*sysUn,2));
			targetGraph->SetPointError(j,targetGraph->GetErrorX(j),newError);
		}
	}
}

/** For pp (where no EP dependence should exist)
  * Calculates the systematic fractional uncertainty in the tgraphs in the input matrices by using k =1,2, nEPBins as separate
  * trials.  Saves them to the parameter output file
  */
void	CalculateSystematicUncertainty(vector <vector<TGraphErrors * >> tGraphs, TFile * f){ 
	// tGraphs expected to be indexed by k = 0,1, .., nEPBins, then j = 0,1, .. nJetPtBins

	int Num = tGraphs[0].size(); 
//	for (int i = 0; i < nJetPtBins; i++) {
	for (int i = 0; i < Num; i++) {
		int n;
		TGraph * SysErr;
		n = tGraphs[0][i]->GetN();
		SysErr = new TGraph(n);
		SysErr->SetName(Form("%s_SysUncert",tGraphs[0][i]->GetName()));
		for (int j = 0; j < n; j++) {
			double sumOfSquares = 0;
			double sum = 0;
			double x = 0;
			double y;
			for (int k = 1; k <= nEPBins - 1; k++) {
				y=0;
				tGraphs[k][i]->GetPoint(j,x,y);
				sumOfSquares += y*y;
				sum += y;
			}
      if (y != 0) SysErr->SetPoint(j,x,(sqrt((sumOfSquares - sum*sum/(nEPBins-1))/(nEPBins-2.))/y));
		}
		f->Add(SysErr);
	}

}





void phase2(TString inputFilename = "root/pp.hist.root", TString
    outputFilename = "Final.root")
{
  // Open input file
  //	inputFilename = Form("%s.root", inputFilename.Data());
  TFile * fIn = TFile::Open(inputFilename.Data(),"READ");

  // Check if file is open
  if (fIn->IsOpen() == kFALSE)
  {
    Printf("inputFilename \"%s\" is not found!", inputFilename.Data());
    std::exit(1);
  }


  TString initialDirectory(gSystem->WorkingDirectory());
  gSystem->ChangeDirectory(outputDirPath.Data());

	gStyle->SetTitleSize(0.04,"xy");

  gStyle->SetOptStat(0);

	TString combinedClass = "";
	
  int nEvents = 0;
  int nEventsPerRun = 0;
  double nEventsPerRunDouble = 0;
  double totalWeight = 0.0;
  double totalWeightPerRun = 0.0;

  double averageCrossSection = 0;
  double crossSectionPerRun = 0;

	// Determine if Jet-Hadron or Pi0-Hadron 

	int iTriggerType = 0; // 0 for jets, 1 for pi0, 2 for pi0 with zt
//	int iParticleBinChoice = 0; // 0 for jetH pt bins, 1 for pi0H pt bins 2 for pi0H zt bins

	// really ad-hoc
	//TString sExampleHist = "pi0Hadron_pi0Pt_5_7_particlePt_0.15_0.25";
	TString sExampleHist = "pi0Hadron_pi0Pt_5_7_particlePt_0.20_0.40";
	TH2F * fExampleHist = 0;
	fExampleHist = (TH2F *) fIn->Get(sExampleHist);
	if (fExampleHist != 0) {
		iTriggerType = 1;
//		iParticleBinChoice = 1;
	} else {
		sExampleHist = "pi0Hadron_pi0Pt_5_7_zT_0.00_0.17";
		fExampleHist = (TH2F *) fIn->Get(sExampleHist);
		if (fExampleHist != 0) {
			iTriggerType = 2; 
//			iParticleBinChoice = 2;
		}
	}
	


  // Get runInfo TTree
  TTree * runInfo = (TTree *) fIn->Get("runInfo");
  if (!runInfo) {
    fprintf(stderr,"Error: TTree runInfo not found.  Exitting ... \n");
    exit(1);
  }
  // Get parameter TTree
  TTree * parameters = (TTree *) fIn->Get("parameters");
  if (!parameters) {
    printf("Parameters TTree not found!  Using parameters from analysis_params.h\n");
  } else {
    printf("Parameters TTree found!\n");
    if (parameters->GetEntries() ) {
      parameters->SetBranchAddress("jetR",&R);
      parameters->SetBranchAddress("constCut",&jet_constituent_cut);
      parameters->SetBranchAddress("etaCut",&eta_cut);
      parameters->SetBranchAddress("particlePtCut",&particle_pt_min);
      parameters->GetEntry(0);
			delta_eta_range = jet_eta_cut + eta_cut;
      if (R< 0.01) { // FIXME temporary measure for float bug
        R = 0.3;
        jet_constituent_cut = 3;
        eta_cut = 1;
        particle_pt_min = 0;
        printf("Using temporary correction\n");
      }
    } else {
      printf("Parameters TTree empty!  Using parameters from analysis_params.h");
    }
  }

	TString sTriggerName = "jet";
	TString sTriggerTitle = "Jet";

  TString jetClass = "jetPt_%.0f_%.0f";
  TString particleClass = "particlePt_%.2f_%.2f";

	TString sAssocPtTitle = "p_{T}^{associated} (GeV/c)";

	int nTriggerPtBins = nJetPtBins;
	std::vector<double> fTriggerPtBins = jetPtBins;

	int nParticleBins = nJetAssocPtBins;
	std::vector<double> assocParticlePtBins = {};

	switch (iTriggerType) {
		case 2:
			printf("Doing Pi0-hadron analysis with zT bins\n");
			sTriggerName = "pi0";
			sTriggerTitle = "#pi^{0}";
			jetClass = "pi0Pt_%.0f_%.0f";
			particleClass = "zT_%.2f_%.2f";
			sAssocPtTitle = "z_{T}";
			nTriggerPtBins = nPi0PtBins;
			fTriggerPtBins = pi0PtBins;
			nAssocParticleBins = nZtBins;
			assocParticlePtBins = zTBins;
			break;
		case 1:
			printf("Doing Pi0-hadron analysis with pT assoc bins\n");
			sTriggerName = "pi0";
			sTriggerTitle = "#pi^{0}";
			jetClass = "pi0Pt_%.0f_%.0f";
			nTriggerPtBins = nPi0PtBins;
			fTriggerPtBins = pi0PtBins;
			nAssocParticleBins = nPi0ParticlePtBins;
      assocParticlePtBins = pi0ParticlePtBins; 
			break;
		default:
		case 0:
			printf("Doing Jet-hadron analysis\n");
			assocParticlePtBins = particlePtBins; // the jet-h ones
			printf("Jet Radius Parameter:  %f\n",R);
			printf("Constituent Cut     :  %f\n",jet_constituent_cut); 
	}

  printf("Eta Cut             :  %f\n",eta_cut);
  printf("Particle Pt Cut     :  %f\n",particle_pt_min); 



  TBranch *bNEvents = runInfo->GetBranch("nEvents");
  TBranch *bTotalWeight = runInfo->GetBranch("totalWeight");
  TBranch *bCrossSection = runInfo->GetBranch("crossSection");

	// Check type of bNEvents for compatibility.  This is Kirill's fault.
	printf("DEBUG: bNEvents title = %s\n",bNEvents->GetTitle());
	bool useInt = !strcmp(bNEvents->GetTitle(),"nEvents/I");
	printf("useInt = %d\n",useInt);
	if (useInt)	bNEvents->SetAddress(&nEventsPerRun);
	else bNEvents->SetAddress(&nEventsPerRunDouble);
  bTotalWeight->SetAddress(&totalWeightPerRun);
  if ( bCrossSection) bCrossSection->SetAddress(&crossSectionPerRun);

  int nRuns = runInfo->GetEntries();
  for (int i = 0; i < nRuns; i++) {
    runInfo->GetEntry(i);
    nEvents += nEventsPerRun + (int) nEventsPerRunDouble;
    totalWeight += totalWeightPerRun;
    averageCrossSection += crossSectionPerRun;
  }
  // from picobarns to millibarns
  if (runInfo->GetEntries()) averageCrossSection = 1e-9 * averageCrossSection / runInfo->GetEntries();


  printf("According to the TTree, I have %d events from %d files, with total weight %f and average cross section %f millibarns.\n",nEvents,nRuns,totalWeight, averageCrossSection);
	// Checking the number of Event Plane Bins (either 1 for no event plane binning, or 4 for all,in,mid,out)
	if (fIn->Get("leadingJetPtEP_3")) {
		printf("It appears the input file used binning by event plane.\n");
		nEPBins = 4;
	}	

  // Get histograms
  TH2F * dEtaDPhi = (TH2F *) fIn->Get("dEtaDPhi");
  TH2F * etaPhi = (TH2F *) fIn->Get("etaPhi");
  TH1D * leadingJetPt = (TH1D *) fIn->Get("leadingJetPt");
  TH1D * jetPt = (TH1D *) fIn->Get("jetPt");
  TH1D * recJetPt = (TH1D *) fIn->Get("recJetPt");
  TH1D * trackPt = (TH1D *) fIn->Get("trackPt");
  TH1D * chargedTrackPt = (TH1D *) fIn->Get("chargedTrackPt");
	TH1D * leadingTrackPt = (TH1D *) fIn->Get("leadingTrackPt");
  TH1D * gammaEt = (TH1D *) fIn->Get("gammaEt");
  TH1D * jetPtUnweighted = (TH1D *) fIn->Get("jetPtUnweighted");
  TH1D * particleEnergy = (TH1D *) fIn->Get("particleEnergy");
  TH2F * jetEtaPhi = (TH2F *) fIn->Get("jetEtaPhi");
  TH2F * leadingJetEtaPhi = (TH2F *) fIn->Get("leadingJetEtaPhi");
  TH1F * dijetAj = (TH1F *) fIn->Get("dijetAj");
  TH2F * dijetAjJetPt = (TH2F *) fIn->Get("dijetAjJetPt");
  TH1F * dijetXj = (TH1F *) fIn->Get("dijetXj");
  TH2F * dijetXjJetPt = (TH2F *) fIn->Get("dijetXjJetPt");
  TH1I * hWeight = (TH1I *) fIn->Get("hWeight");


	TH2F * phiPt = (TH2F *) fIn->Get("phiPt");
	TH2F * jetPhiPt = (TH2F *) fIn->Get("jetPhiPt");


	// Event plane 
  TH1D * leadingJetPtEP[nEPBins];
  for (int k = 0; k < nEPBins; k++) {
		combinedClass = Form("leadingJetPtEP_%d", k);
		leadingJetPtEP[k] = (TH1D *) fIn->Get(combinedClass);
	}

  TH1F * vertexR = (TH1F *) fIn->Get("vertexR");
  TH1F * vertexRHardCore = (TH1F *) fIn->Get("vertexRHardCore");
  TH1F * vertexRNoHardCore = (TH1F *) fIn->Get("vertexRNoHardCore");
  TH1F * vertexRDijetCut = (TH1F *) fIn->Get("vertexRDijetCut");
  TH2F * vertexXY = (TH2F *) fIn->Get("vertexXY");
  TH2F * vertexXYHardCore = (TH2F *) fIn->Get("vertexXYHardCore");
  TH2F * vertexXYNoHardCore = (TH2F *) fIn->Get("vertexXYNoHardCore");
  TH2F * vertexXYRot = (TH2F *) fIn->Get("vertexXYRot");
  TH2F * vertexXYRotHardCore = (TH2F *) fIn->Get("vertexXYRotHardCore");
  TH2F * vertexXYRotNoHardCore = (TH2F *) fIn->Get("vertexXYRotNoHardCore");
  TH2F * vertexXYRotDijet = (TH2F *) fIn->Get("vertexXYRotDijet");
  TH2F * vertexHighestPtXRot = (TH2F *) fIn->Get("vertexHighestPtXRot");
  TH2F * vertexXRotJetPt = (TH2F *) fIn->Get("vertexXRotJetPt");
  TH2F * vertexXRotJetPtHardCore = (TH2F *) fIn->Get("vertexXRotJetPtHardCore");
  TH2F * vertexXRotTriggerPt = (TH2F *) fIn->Get("vertexXRotTriggerPt");


  TH2F * vertexJetPtR = (TH2F *) fIn->Get("vertexJetPtR");
  TH2F * vertexJetPtRHardCore = (TH2F *) fIn->Get("vertexJetPtRHardCore");
  TH2F * vertexJetPtRDijet = (TH2F *) fIn->Get("vertexJetPtRDijet");
  TH2F * vertexDijetAjR = (TH2F *) fIn->Get("vertexDijetAjR");
  TH2F * vertexHighestPtR = (TH2F *) fIn->Get("vertexHighestPtR");
  TH2F * vertexHighestJetZR = (TH2F *) fIn->Get("vertexHighestJetZR");
  TH2F * vertexHighestJetZXRot = (TH2F *) fIn->Get("vertexHighestJetZXRot");
  TH3F * vertexXYRotAjDijet = (TH3F *) fIn->Get("vertexXYRotAjDijet");

  TH2F * ptLossVertexRLeadingJet = (TH2F *) fIn->Get("ptLossVertexRLeadingJet");
  TH2F * ptLossVertexRSubleadingJet = (TH2F *) fIn->Get("ptLossVertexRSubleadingJet");

  TH1F * qqbar_pt = (TH1F *) fIn->Get("qqbar_pt");
  TH1F * qq_pt = (TH1F *) fIn->Get("qq_pt");
  TH1F * gg_pt = (TH1F *) fIn->Get("gg_pt");
  TH1F * qg_pt = (TH1F *) fIn->Get("qg_pt");

  TH1F * qqbar_PartonPtByJetBin[nJetPtBins];
  TH1F * qq_PartonPtByJetBin[nJetPtBins];
  TH1F * gq_PartonPtByJetBin[nJetPtBins];
  TH1F * qg_PartonPtByJetBin[nJetPtBins];
  TH1F * gg_PartonPtByJetBin[nJetPtBins];

  TH1F * PartonPtByJetBinEP[nEPBins][nJetPtBins];
  TH1F * qqbar_PartonPtByJetBinEP[nEPBins][nJetPtBins];
  TH1F * qq_PartonPtByJetBinEP[nEPBins][nJetPtBins];
  TH1F * gq_PartonPtByJetBinEP[nEPBins][nJetPtBins];
  TH1F * qg_PartonPtByJetBinEP[nEPBins][nJetPtBins];
  TH1F * gg_PartonPtByJetBinEP[nEPBins][nJetPtBins];

	TH1F * ptLossLeadingJetPtBinEP[nEPBins][nJetPtBins];
	TH1F * energyLossLeadingJetPtBinEP[nEPBins][nJetPtBins];
	TH1F * ptLossLeadingPartonPtBinEP[nEPBins][nJetPtBins];
	TH1F * energyLossLeadingPartonPtBinEP[nEPBins][nJetPtBins];

	TH1F * ptLossSubleadingJetPtBinEP[nEPBins][nJetPtBins];
	TH1F * energyLossSubleadingJetPtBinEP[nEPBins][nJetPtBins];
	TH1F * ptLossSubleadingPartonPtBinEP[nEPBins][nJetPtBins];
	TH1F * energyLossSubleadingPartonPtBinEP[nEPBins][nJetPtBins];




  TH2F * hsPtLeadingJetPt = (TH2F *) fIn->Get("hsPtLeadingJetPt");
  TH2F * hsPtSubLeadingJetPt = (TH2F *) fIn->Get("hsPtSubLeadingJetPt");
  TH2F * hsDEtaDPhiLeadingJet = (TH2F *) fIn->Get("hsDEtaDPhiLeadingJet");
  TH2F * hsDEtaDPhiSubLeadingJet = (TH2F *) fIn->Get("hsDEtaDPhiSubLeadingJet");


	TH2F * hs1hs2Pt = (TH2F *) fIn->Get("hs1hs2Pt");
	TH1F * hs1hs2dPhi = (TH1F *) fIn->Get("hs1hs2dPhi");


  // Integral? GetEffectiveEntries? GetEntries?
  // StatOverflows?


  if (! (dEtaDPhi && etaPhi && jetPt && particleEnergy && jetEtaPhi && hWeight) ) {
    fprintf(stderr,"Error: a histogram was not found!\n");
    return;
  }

  double nevents = hWeight->GetEntries();
  double njets = jetPt->GetEntries();
  printf("This file has %.f events with %.f %ss.\n",nevents,njets,sTriggerName.Data());

  // Jet-hadron is done more differentially

  TH2F * jetHadron[nTriggerPtBins][nAssocParticleBins];
//  TH2F * jetHadron[nJetPtBins][nAssocParticleBins];
  TH2F * jetHadronScatCent[nTriggerPtBins][nAssocParticleBins];
//  TH2F * jetHadronRecoil[nJetPtBins][nAssocParticleBins];
//  TH2F * jetHadronNotRecoil[nJetPtBins][nAssocParticleBins];
  TH2F * jetHadronOnlyJet[nTriggerPtBins][nAssocParticleBins];
  TH2F * jetHadronExcludeJet[nTriggerPtBins][nAssocParticleBins];
  TH1F * ptBinHadronPt[nTriggerPtBins];

  TH1F * ptBinHadronPtEP[nEPBins][nTriggerPtBins];
	TH2F * jetHadronEP[nEPBins][nTriggerPtBins][nAssocParticleBins];
	TH2F * jetHadronScatCentEP[nEPBins][nTriggerPtBins][nAssocParticleBins];


  // for Swift background
  TH1F *swiftJetYield[nTriggerPtBins];
  TH1F *swiftParticleYield[nTriggerPtBins][nAssocParticleBins];

  combinedClass = "";
  double njets_bin = 0;


  njets_bin = jetPtUnweighted->Integral(jetPt->FindBin(0),jetPt->FindBin(100));
  printf("I have %f binned %ss.",njets_bin,sTriggerName.Data());
  njets_bin = jetPt->Integral(jetPt->FindBin(0),jetPt->FindBin(100));
  printf("I have %f binned %ss after weighting..\n",njets_bin,sTriggerName.Data());



	for (unsigned int i = 0; i < nJetPtBins; i++)
	{
    // Parton info
    combinedClass = Form("qqbar_PartonPtByJetBin_" + jetClass,jetPtBins.at(i),jetPtBins.at(i+1));
    qqbar_PartonPtByJetBin[i] = (TH1F *) fIn->Get(combinedClass);
    combinedClass = Form("qq_PartonPtByJetBin_" + jetClass,jetPtBins.at(i),jetPtBins.at(i+1));
    qq_PartonPtByJetBin[i] = (TH1F *) fIn->Get(combinedClass);
    combinedClass = Form("gq_PartonPtByJetBin_" + jetClass,jetPtBins.at(i),jetPtBins.at(i+1));
    gq_PartonPtByJetBin[i] = (TH1F *) fIn->Get(combinedClass);
    combinedClass = Form("qg_PartonPtByJetBin_" + jetClass,jetPtBins.at(i),jetPtBins.at(i+1));
    qg_PartonPtByJetBin[i] = (TH1F *) fIn->Get(combinedClass);
    combinedClass = Form("gg_PartonPtByJetBin_" + jetClass,jetPtBins.at(i),jetPtBins.at(i+1));
    gg_PartonPtByJetBin[i] = (TH1F *) fIn->Get(combinedClass);
	}
	printf("Finished loading some parton histograms\n");	

  for (unsigned int i = 0; i < nTriggerPtBins; i++)
  {
    njets_bin = jetPt->Integral(jetPt->FindBin(fTriggerPtBins[i]),jetPt->FindBin(fTriggerPtBins[i+1]));

    combinedClass = Form("ptBinHadronPt_%sPt_%.0f_%.0f",sTriggerName.Data(),fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
    //combinedClass = Form("ptBinHadronPt_" + jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
		printf("Searching for object %s\n",combinedClass.Data());
    ptBinHadronPt[i] = (TH1F *) fIn->Get(combinedClass);

    combinedClass = Form("swiftJetYield_%sPt_%.0f_%.0f",sTriggerName.Data(),fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
    //combinedClass = Form("swiftJetYield_" + jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
		printf("Searching for object %s\n",combinedClass.Data());
    swiftJetYield[i] = (TH1F *) fIn->Get(combinedClass);
    for (unsigned int j = 0; j < nAssocParticleBins; j++)
    {
      combinedClass = Form("%sHadron_" + jetClass + "_" + particleClass, sTriggerName.Data(), fTriggerPtBins.at(i), fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
			printf("Searching for object %s\n",combinedClass.Data());
      jetHadron[i][j] = (TH2F *) fIn->Get(combinedClass);
      jetHadron[i][j]->Scale(1./jetHadron[i][j]->GetXaxis()->GetBinWidth(1)); // normalized to dN/dEta

      combinedClass = Form("%sHadronScatCent_" + jetClass + "_" + particleClass, sTriggerName.Data(), fTriggerPtBins.at(i), fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
			printf("Searching for object %s\n",combinedClass.Data());
      jetHadronScatCent[i][j] = (TH2F *) fIn->Get(combinedClass);
      jetHadronScatCent[i][j]->Scale(1./jetHadronScatCent[i][j]->GetXaxis()->GetBinWidth(1));

      // FIXME temporarily making the recoil, not recoil things duplicates of the
      // regular one, to have backwards compatibility for testing

//      combinedClass = Form("jetHadronRecoil_" + jetClass + "_" + particleClass, jetPtBins.at(i), jetPtBins.at(i+1), particlePtBins.at(j), particlePtBins.at(j+1));
//      jetHadronRecoil[i][j] = (TH2F *) jetHadron[i][j]->Clone();
//      jetHadronRecoil[i][j]->SetName(combinedClass);

      //			jetHadronRecoil[i][j] = (TH2F *) fIn->Get(combinedClass);
      //	jetHadronRecoil[i][j]->Scale(1./jetHadronRecoil[i][j]->GetXaxis()->GetBinWidth(1));

  //    combinedClass = Form("jetHadronNotRecoil_" + jetClass + "_" + particleClass, jetPtBins.at(i), jetPtBins.at(i+1), particlePtBins.at(j), particlePtBins.at(j+1));
//      jetHadronNotRecoil[i][j] = (TH2F *) jetHadron[i][j]->Clone();
//      jetHadronNotRecoil[i][j]->SetName(combinedClass);
      //		jetHadronNotRecoil[i][j] = (TH2F *) fIn->Get(combinedClass);
      //			jetHadronNotRecoil[i][j]->Scale(1./jetHadronNotRecoil[i][j]->GetXaxis()->GetBinWidth(1));

      combinedClass = Form("%sHadronOnlyJet_" + jetClass + "_" + particleClass, sTriggerName.Data(), fTriggerPtBins.at(i), fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
			printf("Searching for object %s\n",combinedClass.Data());
      jetHadronOnlyJet[i][j] = (TH2F *) fIn->Get(combinedClass);
      jetHadronOnlyJet[i][j]->Scale(1./jetHadronOnlyJet[i][j]->GetXaxis()->GetBinWidth(1));

      combinedClass = Form("%sHadronExcludeJet_" + jetClass + "_" + particleClass, sTriggerName.Data(), fTriggerPtBins.at(i), fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
			printf("Searching for object %s\n",combinedClass.Data());
      jetHadronExcludeJet[i][j] = (TH2F *) fIn->Get(combinedClass);
      jetHadronExcludeJet[i][j]->Scale(1./jetHadronExcludeJet[i][j]->GetXaxis()->GetBinWidth(1));

      combinedClass = Form("swiftParticleYield_" + jetClass + "_" + particleClass, sTriggerName.Data(), fTriggerPtBins.at(i),fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
			printf("Searching for object %s\n",combinedClass.Data());
      swiftParticleYield[i][j] = (TH1F *) fIn->Get(combinedClass);
    }
  }
 	printf("Finished loading the 2D histograms without EP\n");

	for (unsigned int k = 0; k < nEPBins; k++)
  {
		for (unsigned int i = 0; i < nJetPtBins; i++)
		{

			// Parton info
			combinedClass = Form("PartonPtByJetBin_EP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			PartonPtByJetBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
			combinedClass = Form("qqbar_PartonPtByJetBin_EP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			qqbar_PartonPtByJetBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
			combinedClass = Form("qq_PartonPtByJetBin_EP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			qq_PartonPtByJetBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
			combinedClass = Form("gq_PartonPtByJetBin_EP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			gq_PartonPtByJetBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
			combinedClass = Form("qg_PartonPtByJetBin_EP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			qg_PartonPtByJetBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
			combinedClass = Form("gg_PartonPtByJetBin_EP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			gg_PartonPtByJetBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);

			// pt and energy loss
			combinedClass = Form("ptLossLeadingJetPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			ptLossLeadingJetPtBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
			combinedClass = Form("energyLossLeadingJetPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			energyLossLeadingJetPtBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
			combinedClass = Form("ptLossLeadingPartonPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			ptLossLeadingPartonPtBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
			//	combinedClass = Form("energyLeadingLossPartonPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			combinedClass = Form("energyLossLeadingPartonPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			energyLossLeadingPartonPtBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
			// Temporary bug fix
			if (!energyLossLeadingPartonPtBinEP[k][i]) {
				combinedClass = Form("energyLeadingLossPartonPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
				energyLossLeadingPartonPtBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
				energyLossLeadingPartonPtBinEP[k][i]->SetName(Form("energyLossLeadingPartonPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1)));
			}

			combinedClass = Form("ptLossSubleadingJetPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			ptLossSubleadingJetPtBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
			combinedClass = Form("energyLossSubleadingJetPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			energyLossSubleadingJetPtBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
			combinedClass = Form("ptLossSubleadingPartonPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			ptLossSubleadingPartonPtBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
			combinedClass = Form("energyLossSubleadingPartonPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
			energyLossSubleadingPartonPtBinEP[k][i] = (TH1F *) fIn->Get(combinedClass);
		}

		for (unsigned int i = 0; i < nTriggerPtBins; i++)
		{

			combinedClass = Form("ptBinHadronPt_EP_%d_" + jetClass,k,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinHadronPtEP[k][i] = (TH1F *) fIn->Get(combinedClass);
			if (!ptBinHadronPtEP[k][i]) {
				fprintf(stderr,"Could not find %s in input file!\n",combinedClass.Data());
				exit(1);
			}
			for (unsigned int j = 0; j < nAssocParticleBins; j++)
			{
				combinedClass = Form("%sHadron_EP_%d_" + jetClass + "_" + particleClass, sTriggerName.Data(), k, fTriggerPtBins.at(i), fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
				jetHadronEP[k][i][j] = (TH2F *) fIn->Get(combinedClass);
				if (jetHadronEP[k][i][j]) { 
					jetHadronEP[k][i][j]->Scale(1./jetHadronEP[k][i][j]->GetYaxis()->GetBinWidth(1)); // scaling by delta phi bin width
					TString tempName = jetHadronEP[k][i][j]->GetName();
					tempName.ReplaceAll(".","");
					jetHadronEP[k][i][j]->SetName(tempName);
				}

				combinedClass = Form("%sHadronScatCent_EP_%d_" + jetClass + "_" + particleClass, sTriggerName.Data(), k, fTriggerPtBins.at(i), fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
				jetHadronScatCentEP[k][i][j] = (TH2F *) fIn->Get(combinedClass);
				if (jetHadronScatCentEP[k][i][j]) { 
					jetHadronScatCentEP[k][i][j]->Scale(1./jetHadronScatCentEP[k][i][j]->GetYaxis()->GetBinWidth(1));
					TString tempName = jetHadronScatCentEP[k][i][j]->GetName();
					tempName.ReplaceAll(".","");
					jetHadronScatCentEP[k][i][j]->SetName(tempName);
				}

			}
		}
	}

  // Soft Drop histograms 
  TH2F * hJetPtZ = (TH2F *) fIn->Get("hJetPtZ");
  TH2F * hJetPtMass = (TH2F *) fIn->Get("hJetPtMass");
  TH1D * hSDJetPt = (TH1D *) fIn->Get("hSDJetPt");
  TH1D * hSDJetPtBkgSub = (TH1D *) fIn->Get("hSDJetPtBkgSub");
	TH2F * hJetPtDR = (TH2F *) fIn->Get("hJetPtDR");  
	TH2F * hJetPtZ_DRCut = (TH2F *) fIn->Get("hJetPtZ_DRCut");
	TH2F * hJetPtSD2Mass = (TH2F *) fIn->Get("hJetPtSD2Mass");



  // -----------------------------------------------------------------------------
  // NORMALIZATION
  // -----------------------------------------------------------------------------


  //now make these averages:
  double nEventsPerRun_f = nEvents/nRuns;
  totalWeightPerRun = totalWeight/nRuns;


  printf("nruns 			 	= %d\n",nRuns);
  printf("total events 	= %d\n",nEvents);
  printf("total weight 	= %f\n",totalWeight);
  printf("events/run 	 	= %f\n",nEventsPerRun_f);
  printf("weight/run   	= %f\n",totalWeightPerRun);


  //Normalizing all histograms by nfiles, so that we are effectively averaging 
  // between them all

  bool weightingOn = true; //?
  if (weightingOn) {
    double rescale = 1;
    if (totalWeight == nEvents) {
      printf("I think this is not weighted\n.");
      if (averageCrossSection != 0.) {
        rescale = averageCrossSection * (1.0/nRuns) * (1.0 / totalWeightPerRun); // =  1.0/totalWeight
      } else {
        rescale = (1.0/nRuns) * (1.0 / totalWeightPerRun); // =  1.0/totalWeight
      }

    } else {
      if (averageCrossSection != 0.) {
        rescale = averageCrossSection * (1.0/nRuns) * (1.0/totalWeightPerRun); // something about weighting?
      } else {
				printf("Ignoring Cross Section\n");
        rescale = (1.0/nRuns) * (1.0/totalWeightPerRun); // something about weighting?
      }	

      //note: it is actually the same for both cases
    }

    TObjLink *lnk = gDirectory->GetList()->FirstLink();
    TObject *obj;
    TH1 *hobj;
    while (lnk) {
      obj = gDirectory->FindObject(lnk->GetObject()->GetName());    
      lnk = lnk->Next();
      if (obj->InheritsFrom("TH1") && strcmp(obj->GetName(),"hWeight") && strcmp(obj->GetName(),"jetPtUnweighted")) {
        hobj = (TH1 *) obj;
        hobj->Scale(rescale);
      }    
    }
  } 
	printf("Finished applying general weighting\n");

	// Creating leadingJetPt histogram for normalization
  if (!leadingJetPt) leadingJetPt = (TH1D *) hsPtLeadingJetPt->ProjectionY("leadingJetPt");

  for (unsigned int i = 0; i < nTriggerPtBins; i++) {
    njets_bin = leadingJetPt->Integral(jetPt->FindBin(fTriggerPtBins[i]),jetPt->FindBin(fTriggerPtBins[i+1]));
    if (njets_bin) {
      for (unsigned int j = 0; j < nAssocParticleBins; j++) {
        jetHadron[i][j]->Scale(1./njets_bin);
        jetHadronScatCent[i][j]->Scale(1./njets_bin);
//        jetHadronRecoil[i][j]->Scale(1./njets_bin);
 //       jetHadronNotRecoil[i][j]->Scale(1./njets_bin);
        swiftParticleYield[i][j]->Scale(1./njets_bin);
      }
    }
  }

	bool bUseSeparateNormalizations = true; // whether to normalize by N_{jet} or N_{jet within EP bin}
	// Correct is to use Separate normalizations. This is also what was done for 2.76
	// This is what you want for the final normalization. If RPF is done, don't do in-bin normalization yet

	if (bUseSeparateNormalizations) printf("Normalizing by number of jets in each EP bin separately.\n");

	// Event Plane analysis
	if (leadingJetPtEP[0]) {
		for (unsigned int k = 0; k < nEPBins; k++) {
			for (unsigned int i = 0; i < nTriggerPtBins; i++) {
				printf("k = %d i = %d\n",k,i);
				// Normalizing by 1 / (N_jets within EP Bin)
				printf("leadingJetPtEP [%d] has name %s\n",k,leadingJetPtEP[k]->GetName());
				printf("jetPt has name %s\n",jetPt->GetName());
				if (bUseSeparateNormalizations) {
					njets_bin = leadingJetPtEP[k]->Integral(jetPt->FindBin(fTriggerPtBins[i]),jetPt->FindBin(fTriggerPtBins[i+1]));
				} else { 
					njets_bin = leadingJetPtEP[0]->Integral(jetPt->FindBin(fTriggerPtBins[i]),jetPt->FindBin(fTriggerPtBins[i+1]));
				}
				printf("Found integral %f\n",njets_bin);
				if (njets_bin) {
					for (unsigned int j = 0; j < nAssocParticleBins; j++) {
						printf("jetHadron EP has name %s\n",jetHadronEP[k][i][j]->GetName());
						jetHadronEP[k][i][j]->Scale(1./njets_bin);
            if (jetHadronScatCentEP[k][i][j]) jetHadronScatCentEP[k][i][j]->Scale(1./njets_bin);
					}
				}
			}
		}
	}
	printf("what\n");
  // Normalizing histograms:
  // jetPt should be normalized by nevents
  // jet-h correlations should be normalized by njets
  //  jetPt->Scale(1./nevents);    
  //  recJetPt->Scale(1./nevents);    

  TH2F * hJetPtZNorm = 0;
  if(hJetPtZ) {
		// Setting 0 bin in z to 0 for 
		for (int i = 0; i < hJetPtZ->GetNbinsX(); i++) {
	//		double x_cent = hJetPtZ->GetXaxis()->GetBinCenter(i+1);
			hJetPtZ->SetBinContent(i+1,1,0.0);
		}
    hJetPtZ->SetTitle("Soft Drop Groomed z");
    hJetPtZNorm = normalizeHistogramColumns(hJetPtZ);
    hJetPtZNorm->GetYaxis()->SetTitle("z_{g}");
    hJetPtZNorm->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
  }
	TH2F * hJetPtZ_DRCutNorm = 0;
	TH2F * hJetPtDRNorm = 0;
	if (hJetPtZ_DRCut) {
		for (int i = 0; i < hJetPtZ_DRCut->GetNbinsX(); i++) {
			hJetPtZ_DRCut->SetBinContent(i+1,1,0.0);
		}
		hJetPtZ_DRCut->SetTitle("Soft Drop Groomed z");
		hJetPtZ_DRCutNorm = normalizeHistogramColumns(hJetPtZ_DRCut);
    hJetPtZ_DRCutNorm->GetYaxis()->SetTitle("z_{g}");
    hJetPtZ_DRCutNorm->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");

		hJetPtDR->SetTitle("Soft Drop #Delta R");
		hJetPtDR->GetYaxis()->SetRangeUser(0.1,hJetPtDR->GetYaxis()->GetXmax());
		hJetPtDRNorm = normalizeHistogramColumns(hJetPtDR);
		hJetPtDRNorm->GetYaxis()->SetTitle("#Delta R");
		hJetPtDRNorm->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");

	}
	TH2F * hJetPtSD2MassNorm = 0;
	if ( hJetPtSD2Mass) {
		hJetPtSD2MassNorm = normalizeHistogramColumns(hJetPtSD2Mass);
	}




  TH2F * hJetPtMassNorm = 0;
  if(hJetPtMass) {
    hJetPtMass->GetXaxis()->SetTitle("Mass (GeV/c^2)");
    hJetPtMassNorm = normalizeHistogramColumns(hJetPtMass);
    hJetPtMassNorm->SetName(hJetPtMassNorm->GetName());
    hJetPtMassNorm->GetYaxis()->SetTitle("Mass");
    hJetPtMassNorm->GetXaxis()->SetTitle("m^{jet} (GeV/c^2)");
  }

	printf("Finished Normalization\n");

  // -----------------------------------------------------------------------------
  // END OF NORMALIZATION (or is it?)
  // -----------------------------------------------------------------------------

	gSystem->ChangeDirectory(initialDirectory.Data());
	gSystem->mkdir(outputDirPath.Data(),1);
  gSystem->ChangeDirectory(outputDirPath.Data());
	// Create output file within the output directory
  TFile * fOut = TFile::Open(outputFilename.Data(),"RECREATE");
  if (fOut->IsOpen() == kFALSE)
  {
    Printf("outputFilename \"%s\" failed to open!", outputFilename.Data());
    std::exit(1);
  }
  fOut->Add(runInfo);
  if (parameters) fOut->Add(parameters);


  // Draw and save histograms
  TCanvas *canvas = new TCanvas("canvas","canvas",c_width,c_height);
  // Plot dEtaDPhi for multiple types of drawings
  dEtaDPhi->GetXaxis()->SetTitle("#Delta#eta");
  dEtaDPhi->GetYaxis()->SetTitle("#Delta#phi");
  dEtaDPhi->Draw("surf1");
  canvas->Print("dEtaDPhi.surf.pdf");
  dEtaDPhi->Draw("colz");
  canvas->Print("dEtaDPhi.pdf");

  fOut->Add(dEtaDPhi);

  // etaPhi plot
  etaPhi->GetXaxis()->SetRangeUser(-eta_cut,eta_cut);
  etaPhi->GetXaxis()->SetTitle("#eta");
  etaPhi->GetYaxis()->SetTitle("#phi");
  etaPhi->Draw("colz");
  canvas->Print("etaPhi.pdf");
  fOut->Add(etaPhi);
  // Overall energy of particle
  particleEnergy->GetXaxis()->SetTitle("E (GeV)");
  particleEnergy->Draw();
  canvas->Print("particleEnergy.pdf");
  fOut->Add(particleEnergy);
  // Jet spectra
  TH1F *jetPtCrossSection = (TH1F *) jetPt->Clone();
  jetPtCrossSection->SetName("jetPtCrossSection");
  jetPtCrossSection->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  jetPtCrossSection->GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}d#eta (mb c/GeV)");
  jetPtCrossSection->Scale(1.0/(eta_range * jetPtCrossSection->GetXaxis()->GetBinWidth(1)));
  jetPtCrossSection->Draw();
  canvas->SetLogy(1);
  canvas->Print("jetPtCrossSection.pdf");
  canvas->SetLogy(0);
	fOut->Add(leadingJetPt);
  fOut->Add(jetPt);
  fOut->Add(jetPtUnweighted);

  // Reconstructed (Background Subtracted) Jet spectrum
  TH1F *recJetPtCrossSection = (TH1F *) recJetPt->Clone();
  recJetPtCrossSection->SetName("recJetPtCrossSection");
  recJetPtCrossSection->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  recJetPtCrossSection->GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}d#eta (mb c/GeV)");
  recJetPtCrossSection->Scale(1.0/(eta_range * recJetPtCrossSection->GetXaxis()->GetBinWidth(1)));
  recJetPtCrossSection->Draw();
  canvas->SetLogy(1);
  canvas->Print("recJetPtCrossSection.pdf");
  canvas->SetLogy(0);

	int JetBinsMarkerStyle[5] = {21,20,22,23,34};

	// Draw the Jet Eta-Phi, maybe find v2 of the triggers
	// maybe calculate v2 in the trigger pt bins
	//jetPhiPt->Draw("COLZ");
	TF1 * FitJetEtaPtBins[nTriggerPtBins];
	TH1D * JetEtaPtBins[nTriggerPtBins];
	TLegend * legJetEtaPtBins = new TLegend(0.6,0.7,0.85,0.9);
	for (int i = 0; i < nTriggerPtBins; i++) {
		int iTriggerBinLow = jetPhiPt->GetYaxis()->FindFixBin(fTriggerPtBins[i]);
		int iTriggerBinHigh = jetPhiPt->GetYaxis()->FindFixBin(fTriggerPtBins[i+1]) - 1;
		JetEtaPtBins[i] = (TH1D *) jetPhiPt->ProjectionX(Form("%sPhiPtBin%d",sTriggerName.Data(),i),iTriggerBinLow,iTriggerBinHigh);
		JetEtaPtBins[i]->SetTitle(Form("%s #phi Distribution",sTriggerTitle.Data())); 
		JetEtaPtBins[i]->Rebin(4);
		JetEtaPtBins[i]->GetYaxis()->SetTitle("dN/d#phi (a.u.)");
		JetEtaPtBins[i]->SetLineColor(cTriggerColorList[i]);
		JetEtaPtBins[i]->SetMarkerColor(cTriggerColorList[i]);
		JetEtaPtBins[i]->SetMarkerStyle(JetBinsMarkerStyle[i]);
		double scale = JetEtaPtBins[i]->Integral("width");
		if (scale != 0) JetEtaPtBins[i]->Scale(1./scale);
		//FitJetEtaPtBins[i] = new TF1(Form("TriggerV2FitBin%d",i),"(1./(2*TMath::Pi())) * (1 + 2 * [0] * TMath::Cos(2*x))",0,2*PI);
		FitJetEtaPtBins[i] = new TF1(Form("TriggerV2FitBin%d",i),"(1./(2*TMath::Pi())) * (1 + 2 * [0] * TMath::Cos(2.*x) + 2*[1]*TMath::Cos(4.*x))",0,2*PI);
		JetEtaPtBins[i]->Fit(FitJetEtaPtBins[i],"0");
		FitJetEtaPtBins[i]->SetLineColor(cTriggerColorList[i]);
	}
	for (int i = 0; i < nTriggerPtBins; i++) {
		if (i == 0) JetEtaPtBins[i]->Draw(); 
		else JetEtaPtBins[i]->Draw("SAME");
		FitJetEtaPtBins[i]->Draw("SAME");
		legJetEtaPtBins->AddEntry(JetEtaPtBins[i],Form("%.0f #leq p_{T} < %.0f",fTriggerPtBins[i],fTriggerPtBins[i+1]),"lp");
	}
	legJetEtaPtBins->Draw("SAME");
	canvas->Print("TriggerPhiDistPtBins.pdf");
	canvas->Print("TriggerPhiDistPtBins.C");
	TGraphErrors * TriggerV2Graph = new TGraphErrors(nTriggerPtBins);
	TGraphErrors * TriggerV4Graph = new TGraphErrors(nTriggerPtBins);
	for (int i = 0; i < nTriggerPtBins; i++) {
		TriggerV2Graph->SetPoint(i,(fTriggerPtBins[i] + fTriggerPtBins[i+1])/2.,FitJetEtaPtBins[i]->GetParameter(0));
		TriggerV2Graph->SetPointError(i,(fTriggerPtBins[i+1] - fTriggerPtBins[i])/2.,FitJetEtaPtBins[i]->GetParError(0));


		TriggerV4Graph->SetPoint(i,(fTriggerPtBins[i] + fTriggerPtBins[i+1])/2.,FitJetEtaPtBins[i]->GetParameter(1));
		TriggerV4Graph->SetPointError(i,(fTriggerPtBins[i+1] - fTriggerPtBins[i])/2.,FitJetEtaPtBins[i]->GetParError(1));

	}
	TriggerV2Graph->SetMarkerStyle(kFullSquare);
	TriggerV2Graph->SetName("Trigger_V2");
	TriggerV2Graph->SetTitle("Trigger v_{2};p_{T} GeV/c;v_{2}");

	TriggerV4Graph->SetMarkerStyle(kFullSquare);
	TriggerV4Graph->SetName("Trigger_V4");
	TriggerV4Graph->SetTitle("Trigger v_{4};p_{T} GeV/c;v_{4}");

	TriggerV2Graph->Draw("ALP");
	canvas->Print("TriggerV2.pdf");
	canvas->Print("TriggerV2.C");
	fOut->Add(TriggerV2Graph);

	TriggerV4Graph->Draw("ALP");
	canvas->Print("TriggerV4.pdf");
	canvas->Print("TriggerV4.C");
	fOut->Add(TriggerV4Graph);

  if (dijetAj) {
    dijetAj->Rebin(4);
    dijetAj->Scale(1.0/dijetAj->Integral("width"));
    canvas->Clear();
    dijetAj->Draw();
    canvas->Print("dijetAj.pdf");
    fOut->Add(dijetAj);
  }
  if (dijetAjJetPt) {
    fOut->Add(dijetAjJetPt);
  }

  if (dijetXj) {
    dijetXj->Rebin(4);
    dijetXj->Scale(1.0/dijetXj->Integral("width"));
    canvas->Clear();
    dijetXj->Draw();
    canvas->Print("dijetXj.pdf");
    fOut->Add(dijetXj);
  }
  if (dijetXjJetPt) {

    fOut->Add(dijetXjJetPt);
  }

  if (vertexR) {
    fOut->Add(vertexR);
    if (vertexXY) fOut->Add(vertexXY);
    if (vertexXYHardCore) fOut->Add(vertexXYHardCore);
    if (vertexXYNoHardCore) fOut->Add(vertexXYNoHardCore);
    if (vertexRDijetCut) fOut->Add(vertexRDijetCut);
    if (vertexXYRot) fOut->Add(vertexXYRot);
    if (vertexXYRotHardCore) fOut->Add(vertexXYRotHardCore);
    if (vertexXYRotNoHardCore) fOut->Add(vertexXYRotNoHardCore);
    if (vertexXYRotDijet) fOut->Add(vertexXYRotDijet);
    if (vertexJetPtRHardCore) fOut->Add(vertexJetPtRHardCore);
    if (vertexJetPtRDijet) fOut->Add(vertexJetPtRDijet);
    if (vertexDijetAjR) fOut->Add(vertexDijetAjR);
    if (vertexJetPtR) { 
			fOut->Add(vertexJetPtR);
			TH2F * vertexJetPtRNorm = normalizeHistogramColumns(vertexJetPtR);
		}
    if (vertexXRotJetPt) {
				fOut->Add(vertexXRotJetPt);
				TH2F * vertexXRotJetPtNorm = normalizeHistogramColumns(vertexXRotJetPt);
				TProfile * vertexXRotJetPt_pfx = vertexXRotJetPt->ProfileX();
				vertexXRotJetPt_pfx->GetYaxis()->SetTitle("<x_{rot}> (fm)");
//				fOut->Add(vertexXRotJetPt_pfx);
		}
    if (vertexXRotJetPtHardCore) fOut->Add(vertexXRotJetPtHardCore);
    if (vertexXRotTriggerPt) {
      fOut->Add(vertexXRotTriggerPt);
      canvas->Clear();
      vertexXRotTriggerPt->Draw("COLZ");
      canvas->Print("vertexXRotTriggerPt.pdf");

      TH2F * vertexXRotTriggerPtNorm = normalizeHistogramColumns(vertexXRotTriggerPt);
      vertexXRotTriggerPtNorm->Draw("COLZ");
      canvas->SetLogz();
      canvas->Print("vertexXRotTriggerPtNorm.pdf");
      canvas->Clear();

    }


    TCanvas *cVertexCmp = new TCanvas("cVertexCmp","cVertexCmp",c_width,c_height);
    cVertexCmp->Divide(2,1);

    TH2F *vertexXYHCOverAll = (TH2F *) vertexXYHardCore->Clone();
    TH2F *vertexXYNoHCOverAll = (TH2F *) vertexXYNoHardCore->Clone();

    TH2F *vertexXY_Clone = (TH2F *) vertexXY->Clone();

    int rebinVertex = 2;	
    vertexXY_Clone->Rebin2D(rebinVertex,rebinVertex);
    vertexXYHCOverAll->Rebin2D(rebinVertex,rebinVertex);
    vertexXYNoHCOverAll->Rebin2D(rebinVertex,rebinVertex);

    vertexXY_Clone->Scale(1./vertexXY_Clone->Integral("width"));
    vertexXYHCOverAll->Scale(1./vertexXYHCOverAll->Integral("width"));
    vertexXYNoHCOverAll->Scale(1./vertexXYNoHCOverAll->Integral("width"));

    cVertexCmp->cd(1);
    vertexXYHCOverAll->Draw("COLZ");
    cVertexCmp->cd(2);
    vertexXYNoHCOverAll->Draw("COLZ");
    cVertexCmp->Print("vertexCMP.pdf");

    cVertexCmp->Clear();
    cVertexCmp->Divide(2,2);


    if(vertexXYRotDijet) {

      TH2F *vertexXYRot_Clone = (TH2F *) vertexXYRot->Clone();
      TH2F *vertexXYRotHC_Clone = (TH2F *) vertexXYRotHardCore->Clone();
      TH2F *vertexXYRotNoHC_Clone = (TH2F *) vertexXYRotNoHardCore->Clone();
      TH2F *vertexXYRotDijet_Clone = (TH2F *) vertexXYRotDijet->Clone();
      vertexXYRot_Clone->Rebin2D(rebinVertex,rebinVertex);
      vertexXYRotHC_Clone->Rebin2D(rebinVertex,rebinVertex);
      vertexXYRotNoHC_Clone->Rebin2D(rebinVertex,rebinVertex);
      vertexXYRotDijet_Clone->Rebin2D(rebinVertex,rebinVertex);
      vertexXYRot_Clone->Scale(1./vertexXYRot_Clone->Integral("width"));
      vertexXYRotHC_Clone->Scale(1./vertexXYRotHC_Clone->Integral("width"));
      vertexXYRotNoHC_Clone->Scale(1./vertexXYRotNoHC_Clone->Integral("width"));
      vertexXYRotDijet_Clone->Scale(1./vertexXYRotDijet_Clone->Integral("width"));

      cVertexCmp->cd(1);
      vertexXYRotHC_Clone->SetTitle("Hard Core Cut");
      vertexXYRotHC_Clone->Draw("COLZ");
      cVertexCmp->cd(2);
      vertexXYRotDijet_Clone->SetTitle("Dijet Cut");
      vertexXYRotDijet_Clone->Draw("COLZ");
      cVertexCmp->cd(3);
      //		vertexXYRotNoHC_Clone->Draw("COLZ");
      vertexXYRot_Clone->SetTitle("Inclusive");
      vertexXYRot_Clone->Draw("COLZ");
      cVertexCmp->Print("vertexRotCMP.pdf");


      //projecting:
      double N_near,N_away;
      TText *text = new TText();
      TLatex *latex = new TLatex();
      // Inclusive
      TH1F * vertexRot_Clone_px = (TH1F *) vertexXYRot_Clone->ProjectionX("_px");
      vertexRot_Clone_px->SetTitle("Jet Axis Projection");
      TH1F * vertexRot_Clone_py = (TH1F *) vertexXYRot_Clone->ProjectionY("_py");
      vertexRot_Clone_py->SetTitle("Perpendicular Axis Projection");

      //calculating s = N_{near} / N_{away}
      // including overflow, underflow bins
      N_near = vertexRot_Clone_px->Integral(0,vertexRot_Clone_px->GetXaxis()->FindBin(0.));
      N_away = vertexRot_Clone_px->Integral(vertexRot_Clone_px->GetXaxis()->FindBin(0.),vertexRot_Clone_px->GetNbinsX() + 1);

      cVertexCmp->Clear();
      cVertexCmp->Divide(2,2);
      cVertexCmp->cd(2);
      vertexXYRot_Clone->Draw("COLZ");
      cVertexCmp->cd(1);
      vertexRot_Clone_py->Draw("COLZ");
      cVertexCmp->cd(4);
      vertexRot_Clone_px->Draw("COLZ");
      cVertexCmp->cd(3);
      latex->DrawLatex(0.5,0.6,"s = N_{near}/N_{away}");
      text->DrawText(0.5,0.5,Form("s = %f",N_near / N_away));
      cVertexCmp->Print("vertexProjInc.pdf");

      // Hard Core Cut
      TH1F * vertexRotHC_Clone_px = (TH1F *) vertexXYRotHC_Clone->ProjectionX("_px");
      vertexRotHC_Clone_px->SetTitle("Jet Axis Projection");
      TH1F * vertexRotHC_Clone_py = (TH1F *) vertexXYRotHC_Clone->ProjectionY("_py");
      vertexRotHC_Clone_py->SetTitle("Perpendicular Axis Projection");

      N_near = vertexRotHC_Clone_px->Integral(0,vertexRotHC_Clone_px->GetXaxis()->FindBin(0.));
      N_away = vertexRotHC_Clone_px->Integral(vertexRotHC_Clone_px->GetXaxis()->FindBin(0.),vertexRotHC_Clone_px->GetNbinsX() + 1);

      cVertexCmp->Clear();
      cVertexCmp->Divide(2,2);
      cVertexCmp->cd(2);
      vertexXYRotHC_Clone->Draw("COLZ");
      cVertexCmp->cd(1);
      vertexRotHC_Clone_py->Draw("COLZ");
      cVertexCmp->cd(4);
      vertexRotHC_Clone_px->Draw("COLZ");
      cVertexCmp->cd(3);
      latex->DrawLatex(0.5,0.6,"s = N_{near}/N_{away}");
      text->DrawText(0.5,0.5,Form("s = %f",N_near / N_away));
      cVertexCmp->Print("vertexProjHC.pdf");

      //Dijet Cut
      TH1F * vertexRotDijet_Clone_px = (TH1F *) vertexXYRotDijet_Clone->ProjectionX("_px");
      vertexRotDijet_Clone_px->SetTitle("Jet Axis Projection");
      TH1F * vertexRotDijet_Clone_py = (TH1F *) vertexXYRotDijet_Clone->ProjectionY("_py");
      vertexRotDijet_Clone_py->SetTitle("Perpendicular Axis Projection");

      N_near = vertexRotDijet_Clone_px->Integral(0,vertexRotDijet_Clone_px->GetXaxis()->FindBin(0.));
      N_away = vertexRotDijet_Clone_px->Integral(vertexRotDijet_Clone_px->GetXaxis()->FindBin(0.),vertexRotDijet_Clone_px->GetNbinsX() + 1);

      cVertexCmp->Clear();
      cVertexCmp->Divide(2,2);
      cVertexCmp->cd(2);
      vertexXYRotDijet_Clone->Draw("COLZ");
      cVertexCmp->cd(1);
      vertexRotDijet_Clone_py->Draw("COLZ");
      cVertexCmp->cd(4);
      vertexRotDijet_Clone_px->Draw("COLZ");
      cVertexCmp->cd(3);
      latex->DrawLatex(0.5,0.6,"s = N_{near}/N_{away}");
      text->DrawText(0.5,0.5,Form("s = %f",N_near / N_away));
      cVertexCmp->Print("vertexProjDijet.pdf");
    }



    //now using R
    cVertexCmp->Clear();
    int rebinVertexR = 1;	

    vertexR->Rebin(rebinVertexR);
    vertexRHardCore->Rebin(rebinVertexR);
    //	vertexRNoHardCore->Rebin(rebinVertexR);
    vertexRDijetCut->Rebin(rebinVertexR);

    vertexR->Scale(1.0/vertexR->Integral("width"));
    vertexRHardCore->Scale(1.0/vertexRHardCore->Integral("width"));
    //		vertexRNoHardCore->Scale(1.0/vertexRNoHardCore->Integral("width"));
    vertexRDijetCut->Scale(1.0/vertexRDijetCut->Integral("width"));

    vertexR->SetLineColor(kBlack);
    vertexR->SetMarkerColor(kBlack);
    vertexR->SetMarkerStyle(kFullSquare);
    vertexRDijetCut->SetLineColor(kRed);
    vertexRDijetCut->SetMarkerColor(kRed);
    vertexRDijetCut->SetMarkerStyle(kFullStar);
    vertexRDijetCut->SetMarkerSize(2);
    //		vertexRNoHardCore->SetLineColor(kBlue);
    vertexRHardCore->SetLineColor(kBlue);
    vertexRHardCore->SetMarkerColor(kBlue);
    vertexRHardCore->SetMarkerStyle(8);
    vertexRHardCore->SetMarkerSize(1);

    vertexR->Draw("LP");
    vertexRHardCore->Draw("LP SAME");
    //	vertexRNoHardCore->Draw("LP SAME");
    vertexRDijetCut->Draw("LP SAME");


    TLegend * legVertex = new TLegend(0.50, 0.65, 0.90, 0.85);
    legVertex->AddEntry(vertexR, "All", "lp");
    legVertex->AddEntry(vertexRHardCore, "Hard Core (>6 GeV/c) Jets", "lp");
    //	legVertex->AddEntry(vertexRNoHardCore, "Removed Jets", "l");
    legVertex->AddEntry(vertexRDijetCut, "Dijet Cut","lp");
    legVertex->Draw("same");

    cVertexCmp->Print("vertexRCMP.pdf");


    // the vertexJetPtR is stored in a TH2, so I have to project it into different
    // histograms to easily scale it, and make useful plots.

    if (vertexJetPtR ) {
			
			// FIXME some of these might be using the pi0 binnings
      int rebinVertexJetPtR = 1;	
      TH1F * vertexJetPtRBin[nJetPtBins];
      TH1F * vertexJetPtRHardCoreBin[nJetPtBins];
      TH1F * vertexJetPtRDijetBin[nJetPtBins];
      TLegend * legVertexJetPtRBin[nJetPtBins];
      for (int i = 0; i < nJetPtBins; i++) {
        combinedClass = Form("vertexR_" + jetClass, jetPtBins.at(i), jetPtBins.at(i+1));
        vertexJetPtRBin[i] = (TH1F *) vertexJetPtR->ProjectionY(combinedClass,i,i+1,"e");
        vertexJetPtRBin[i]->Rebin(rebinVertexJetPtR);
        double scale = vertexJetPtRBin[i]->Integral("width");
        scale = scale != 0 ? 1./scale : 0;
        vertexJetPtRBin[i]->Scale(scale);
        vertexJetPtRBin[i]->SetTitle(Form("%0.f < p_{T}^{jet} < %0.f",jetPtBins.at(i),jetPtBins.at(i+1)));			


        combinedClass = Form("vertexRHardCore_" + jetClass, jetPtBins.at(i), jetPtBins.at(i+1));
        vertexJetPtRHardCoreBin[i] = (TH1F *) vertexJetPtRHardCore->ProjectionY(combinedClass,i,i+1,"e");
        vertexJetPtRHardCoreBin[i]->Rebin(rebinVertexJetPtR);
        scale = vertexJetPtRHardCoreBin[i]->Integral("width");
        scale = scale != 0 ? 1./scale : 0;
        vertexJetPtRHardCoreBin[i]->Scale(scale);

        combinedClass = Form("vertexRDijet_" + jetClass, jetPtBins.at(i), jetPtBins.at(i+1));
        vertexJetPtRDijetBin[i] = (TH1F *) vertexJetPtRDijet->ProjectionY(combinedClass,i,i+1,"e");
        vertexJetPtRDijetBin[i]->Rebin(rebinVertexJetPtR);
        scale = vertexJetPtRDijetBin[i]->Integral("width");
        scale = scale != 0 ? 1./scale : 0;
        vertexJetPtRDijetBin[i]->Scale(scale);

        legVertexJetPtRBin[i] = new TLegend(0.50,0.65,0.90,0.85); 

      }
      TCanvas *cVertexJetPtR = new TCanvas("cVertexJetPtR","cVertexJetPtR",c_width,c_height);
      cVertexJetPtR->Divide(TMath::CeilNint(TMath::Sqrt(nJetPtBins)), 
          TMath::FloorNint(TMath::Sqrt(nJetPtBins)));


      //actually, should make it so it compares hard core to none?

      //		TLegend * legVertex2 = new TLegend(0.50, 0.65, 0.90, 0.85);
      //		TString style = "LP";
      for (int i = 0; i < nJetPtBins; i++) {
        cVertexJetPtR->cd(i+1);

        vertexJetPtRBin[i]->SetLineColor(1);
        vertexJetPtRBin[i]->SetMarkerStyle(33);
        vertexJetPtRBin[i]->SetMarkerColor(1);
        vertexJetPtRHardCoreBin[i]->SetLineColor(2);
        vertexJetPtRHardCoreBin[i]->SetMarkerStyle(33);
        vertexJetPtRHardCoreBin[i]->SetMarkerColor(2);
        vertexJetPtRDijetBin[i]->SetLineColor(3);
        vertexJetPtRDijetBin[i]->SetMarkerStyle(33);
        vertexJetPtRDijetBin[i]->SetMarkerColor(3);

        //			if (i == 1) style = "SAME LP";
        vertexJetPtRBin[i]->Draw("LP");
        vertexJetPtRHardCoreBin[i]->Draw("LP SAME");

        if (i>1) vertexJetPtRDijetBin[i]->Draw("LP SAME");
        //			vertexJetPtRBin[i]->GetYaxis()->SetRangeUser(0,1.1*max(vertexJetPtRBin[i]->GetBinContent(vertexJetPtRBin[i]->GetMaximumBin()),max(vertexJetPtRHardCoreBin[i]->GetBinContent(vertexJetPtRHardCoreBin[i]->GetMaximumBin()),vertexJetPtRDijetBin[i]->GetBinContent(vertexJetPtRDijetBin[i]->GetMaximumBin()))));			


        legVertexJetPtRBin[i]->AddEntry(vertexJetPtRBin[i],"All","lp");			
        legVertexJetPtRBin[i]->AddEntry(vertexJetPtRHardCoreBin[i],"Hard Core (>6 GeV/c)","lp");			
        if (i>1) legVertexJetPtRBin[i]->AddEntry(vertexJetPtRDijetBin[i],"Dijet Cut","lp");			

        legVertexJetPtRBin[i]->Draw("SAME");

      }
      cVertexJetPtR->Print("vertexJetPtRCmp.pdf");


      // Comparing R vs. jet Pt
      cVertexJetPtR->Clear();
      TLegend * legVertexJetPtR = new TLegend(0.50,0.65,0.90,0.85); 

      TString storeTitle = vertexR->GetTitle();
      vertexR->SetTitle("Vertex R by p^{jet}_{T}");
      vertexR->SetMarkerStyle(kOpenSquare);
      vertexR->SetMarkerColor(kBlack);
      vertexR->SetLineColor(kBlack);    
      vertexR->SetLineStyle(1);    


      vertexR->Draw("LP");
      legVertexJetPtR->AddEntry(vertexR,"All Vertices","lp");


      int colors[6] = {1,2,3,7,4,6};

      for (int i = 0; i < nJetPtBins; i++) {
        vertexJetPtRBin[i]->SetMarkerColor(colors[cTriggerColorList[i]]);
        vertexJetPtRBin[i]->SetLineColor(colors[cTriggerColorList[i]]);
        //    if (!i)  vertexJetPtRBin[i]->Draw("LP");
        vertexJetPtRBin[i]->Draw("LP SAME");
        legVertexJetPtR->AddEntry(vertexJetPtRBin[i],vertexJetPtRBin[i]->GetTitle(),"lp"); 
      }
      legVertexJetPtR->Draw("SAME");
      cVertexJetPtR->Print("vertexJetPtR.pdf");
      vertexR->SetTitle(storeTitle);



      // vertexR highest Pt stuff
      // FLEP
      //vertexHighestPtR
      TH2F *vertexHCCutR = (TH2F *) vertexHighestPtR->Clone();
      vertexHCCutR->SetName("vertexHCCutR");		
      vertexHCCutR->SetTitle("Vertex R vs Hard Core Cut");		

      TCanvas *cVertexHighestPtR = new TCanvas("cVertexHighestPtR","cVertexHighestPtR",c_width,c_height);
      cVertexHighestPtR->cd();
      vertexHighestPtR->GetXaxis()->SetRangeUser(2.,30.);
      vertexHighestPtR->Draw("COLZ");
      cVertexHighestPtR->SetLogz();
      cVertexHighestPtR->Print("vertexHighestPtR.pdf");
      vertexHighestPtR->GetXaxis()->SetRangeUser(0.,100.);

      //each bin (pt,R) in vertexHCCutR will be filled with the integral 
      // \int_{pt}^{\infty} d(pt') vertexHighestPtR(pt',R)

      // this is inefficient with respect to memory reading, but saves
      // time and memory.
      int nbinsx = vertexHighestPtR->GetNbinsX();
      int nbinsy = vertexHighestPtR->GetNbinsY();

      for (int j = 0; j < nbinsy; j++) {
        //manually do the integral starting from infinity
        // filling in vertexHCCutR along the way 
        double integral;
        // add in overflow
        integral = vertexHighestPtR->GetBinContent(vertexHighestPtR->GetBin(nbinsx+1,j));
        for (int i = nbinsx ; i >= 0; i--) {
          // do integral
          integral += vertexHighestPtR->GetBinContent(vertexHighestPtR->GetBin(i,j));
          // set value
          vertexHCCutR->SetBinContent(vertexHCCutR->GetBin(i,j),integral); 				
        }
      }

      cVertexHighestPtR->cd();
      vertexHCCutR->GetXaxis()->SetRangeUser(2.,30.);
      vertexHCCutR->Draw("COLZ");
      cVertexHighestPtR->Print("vertexHCCutR.pdf");
      vertexHCCutR->GetXaxis()->SetRangeUser(0.,100.);

      // Normalized version
      TH2F *vertexHCCutRNorm = (TH2F *) vertexHCCutR->Clone();
      vertexHCCutRNorm->SetName("vertexHCCutRNorm");
      vertexHCCutRNorm->SetTitle("Vertex R vs Hard Core Cut (normalized)");
      TH2F *vertexHCCutRScale = (TH2F *) vertexHCCutR->Clone();
      vertexHCCutRScale->SetName("vertexHCCutRScale");
      vertexHCCutRScale->SetTitle("vertexHCCutRScale");

      for (int i = 0; i < nbinsx+1; i++) {
        double integral = vertexHCCutR->Integral(i,i+1,0,nbinsy+1,"width");
        integral = integral ? 1./integral : 0; // lol
        for (int j = 0; j < nbinsy+1; j++) {
          vertexHCCutRScale->SetBinContent(i,j,integral);
        }
      }
      //scaling so each column is normalized
      vertexHCCutRNorm->Multiply(vertexHCCutRScale);
      //drawing the norms as a test
      cVertexHighestPtR->cd(); cVertexHighestPtR->Clear();
      // showing profile as well
      vertexHCCutRNorm->GetXaxis()->SetRangeUser(0.,20.);
      vertexHCCutRNorm->GetYaxis()->SetRangeUser(0.,10.);
      vertexHCCutRNorm->Draw("COLZ");
      cVertexHighestPtR->SetLogz(1);
      TProfile *vertexHCCutRNorm_px = vertexHCCutRNorm->ProfileX("_pfx",1,-1);
     // TProfile *vertexHCCutRNorm_px = vertexHCCutRNorm->ProfileX("_pfx",1,-1,"s");
      //vertexHCCutRNorm_px->SetMarkerStyle(kOpenSquare);
      vertexHCCutRNorm_px->SetLineWidth(2);
      vertexHCCutRNorm_px->SetLineColor(1);
      vertexHCCutRNorm_px->SetMarkerColor(1);
      vertexHCCutRNorm_px->Draw("SAME");

      cVertexHighestPtR->Print("vertexHCCutRNorm.pdf");
      vertexHCCutRNorm->GetXaxis()->SetRangeUser(0.,100.);


      // vertexHighestPtXRot stuff
      if (vertexHighestPtXRot) {

        TH2F *vertexHCCutXRot = (TH2F *) vertexHighestPtXRot->Clone();
        vertexHCCutXRot->SetName("vertexHCCutXRot");		
        vertexHCCutXRot->SetTitle("Vertex x_{rot} vs Hard Core Cut");		

        TCanvas *cVertexHighestPtXRot = new TCanvas("cVertexHighestPtXRot","cVertexHighestPtXRot",c_width,c_height);
        cVertexHighestPtXRot->cd();
        vertexHighestPtXRot->GetXaxis()->SetRangeUser(2.,30.);
        vertexHighestPtXRot->Draw("COLZ");
        cVertexHighestPtXRot->SetLogz();
        cVertexHighestPtXRot->Print("vertexHighestPtXRot.pdf");
        vertexHighestPtXRot->GetXaxis()->SetRangeUser(0.,100.);

        //each bin (pt,R) in vertexHCCutR will be filled with the integral 
        // \int_{pt}^{\infty} d(pt') vertexHighestPtR(pt',R)

        // this is inefficient with respect to memory reading, but saves
        // time and memory.
        nbinsx = vertexHighestPtXRot->GetNbinsX();
        nbinsy = vertexHighestPtXRot->GetNbinsY();

        for (int j = 0; j < nbinsy; j++) {
          //manually do the integral starting from infinity
          // filling in vertexHCCutR along the way 
          double integral;
          // add in overflow
          integral = vertexHighestPtXRot->GetBinContent(vertexHighestPtXRot->GetBin(nbinsx+1,j));
          for (int i = nbinsx ; i >= 0; i--) {
            // do integral
            integral += vertexHighestPtXRot->GetBinContent(vertexHighestPtXRot->GetBin(i,j));
            // set value
            vertexHCCutXRot->SetBinContent(vertexHCCutR->GetBin(i,j),integral); 				
          }
        }

        cVertexHighestPtXRot->cd();
        vertexHCCutXRot->GetXaxis()->SetRangeUser(2.,30.);
        vertexHCCutXRot->Draw("COLZ");
        cVertexHighestPtXRot->Print("vertexHCCutXRot.pdf");
        vertexHCCutXRot->GetXaxis()->SetRangeUser(0.,100.);

        // Normalized version
        TH2F *vertexHCCutXRotNorm = (TH2F *) vertexHCCutXRot->Clone();
        vertexHCCutXRotNorm->SetName("vertexHCCutXRotNorm");
        vertexHCCutXRotNorm->SetTitle("Vertex x_{rot} vs Hard Core Cut (normalized)");
        vertexHCCutXRotNorm->GetXaxis()->SetTitle("Const. Cut (GeV/c)");
        TH2F *vertexHCCutXRotScale = (TH2F *) vertexHCCutXRot->Clone();
        vertexHCCutXRotScale->SetName("vertexHCCutXRotScale");
        vertexHCCutXRotScale->SetTitle("vertexHCCutXRotScale");

        for (int i = 0; i < nbinsx+1; i++) {
          double integral = vertexHCCutXRot->Integral(i,i+1,0,nbinsy+1,"width");
          integral = integral ? 1./integral : 0; // lol
          for (int j = 0; j < nbinsy+1; j++) {
            vertexHCCutXRotScale->SetBinContent(i,j,integral);
          }
        }
        //scaling so each column is normalized
        vertexHCCutXRotNorm->Multiply(vertexHCCutXRotScale);
        //drawing the norms as a test
        cVertexHighestPtXRot->cd(); cVertexHighestPtXRot->Clear();
        // showing profile as well
        vertexHCCutXRotNorm->GetXaxis()->SetRangeUser(0.,20.);
        vertexHCCutXRotNorm->GetYaxis()->SetRangeUser(-10.,10.);
        vertexHCCutXRotNorm->Draw("COLZ");
        cVertexHighestPtXRot->SetLogz(1);
        TProfile *vertexHCCutXRotNorm_px = vertexHCCutXRotNorm->ProfileX("_pfx",1,-1);
     //   TProfile *vertexHCCutXRotNorm_px = vertexHCCutXRotNorm->ProfileX("_pfx",1,-1,"s");
        //vertexHCCutRNorm_px->SetMarkerStyle(kOpenSquare);
        vertexHCCutXRotNorm_px->SetLineWidth(2);
        vertexHCCutXRotNorm_px->SetLineColor(1);
        vertexHCCutXRotNorm_px->SetMarkerColor(1);
        vertexHCCutXRotNorm_px->Draw("SAME");

        cVertexHighestPtXRot->Print("vertexHCCutXRotNorm.pdf");
        vertexHCCutXRotNorm->GetXaxis()->SetRangeUser(0.,100.);
      }



      if(vertexHighestJetZR) {
        fOut->Add(vertexHighestJetZR);
        cVertexHighestPtR->Clear();
        cVertexHighestPtR->SetLogz();
        vertexHighestJetZR->Draw("COLZ");
        cVertexHighestPtR->Print("vertexHighestJetZR.pdf");
        cVertexHighestPtR->Clear();
        TH2F * vertexHighestJetZRNorm = normalizeHistogramColumns(vertexHighestJetZR);
        TProfile * vertexHighestJetZRNorm_pfx = vertexHighestJetZRNorm->ProfileX("_pfx",1,-1);
    //    TProfile * vertexHighestJetZRNorm_pfx = vertexHighestJetZRNorm->ProfileX("_pfx",1,-1,"s");
        vertexHighestJetZRNorm_pfx->SetLineColor(kBlack);
        vertexHighestJetZRNorm->Draw("COLZ");
        vertexHighestJetZRNorm_pfx->Draw("SAME");
        cVertexHighestPtR->Print("vertexHighestJetZRNorm.pdf");
      }

      //		legVertex2->AddEntry(vertexJetPtRBin[i],Form("%.0f < p_{T}^{jet} < %.0f",jetPtBins.at(i),jetPtBins.at(i+1)),"lp");
    }
    //	TString jetClass = "jetPt_%.0f_%.0f";
    //		legVertex2->Draw("SAME");

  }


  if (trackPt) fOut->Add(trackPt);
	if (chargedTrackPt) fOut->Add(chargedTrackPt);
	if (leadingTrackPt) fOut->Add(leadingTrackPt);
  if(gammaEt) fOut->Add(gammaEt);
  if(recJetPt) fOut->Add(recJetPt);


  canvas->Clear();
  jetEtaPhi->GetXaxis()->SetRangeUser(-eta_cut,eta_cut);
  jetEtaPhi->GetXaxis()->SetTitle("#eta");
  jetEtaPhi->GetYaxis()->SetTitle("#phi");
  jetEtaPhi->Draw("colz");
  canvas->Print("jetEtaPhi.pdf");
  fOut->Add(jetEtaPhi);

	// Calculating Particle Event plane v2
	if (phiPt) {
		int nPhiBins = phiPt->GetNbinsY(); 

		// FIXME use same pt bins as associated??

		TGraphErrors * particleV2 = new TGraphErrors(nPhiBins);
		particleV2->SetName("Track_V2");
		particleV2->SetTitle("Track V2");
		TGraphErrors * particleV4 = new TGraphErrors(nPhiBins);
		particleV4->SetName("Track_V4");
		particleV4->SetTitle("Track V4");
		TGraphErrors * particleV6 = new TGraphErrors(nPhiBins);
		particleV6->SetName("Track_V6");
		particleV6->SetTitle("Track V6");
		

		for (int i = 0; i < nPhiBins; i++) {
			TH1F * phiProj = (TH1F *) phiPt->ProjectionX("_proj_phi",i+1,i+1);
			double v2,v2Error,v4,v4Error,v6,v6Error;
			normalizeTH1F(phiProj);
//			v2 = calculateFixedV2AndError(phiProj,v2Error);
//			calculateVnAndError(phiProj,v2,v2Error,v4,v4Error);
			calculateVnAndError(phiProj,v2,v2Error,v4,v4Error,v6,v6Error);
			particleV2->SetPoint(i,phiPt->GetYaxis()->GetBinCenter(i+1),v2);
			particleV2->SetPointError(i,phiPt->GetYaxis()->GetBinWidth(i+1)/2.,v2Error);

			particleV4->SetPoint(i,phiPt->GetYaxis()->GetBinCenter(i+1),v4);
			particleV4->SetPointError(i,phiPt->GetYaxis()->GetBinWidth(i+1)/2.,v4Error);

			particleV6->SetPoint(i,phiPt->GetYaxis()->GetBinCenter(i+1),v6);
			particleV6->SetPointError(i,phiPt->GetYaxis()->GetBinWidth(i+1)/2.,v6Error);

			delete phiProj;
		}
		particleV2->SetMarkerStyle(kFullSquare);
		particleV2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		particleV2->GetYaxis()->SetTitle("V_{2}");
		particleV4->SetMarkerStyle(kFullSquare);
		particleV4->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		particleV4->GetYaxis()->SetTitle("V_{4}");
		particleV6->SetMarkerStyle(kFullSquare);
		particleV6->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		particleV6->GetYaxis()->SetTitle("V_{6}");



		particleV2->Draw("ALP");
		canvas->Print("TrackV2.pdf");
		canvas->Print("TrackV2.C");
		particleV4->Draw("ALP");
		canvas->Print("TrackV4.pdf");
		canvas->Print("TrackV4.C");



		fOut->Add(phiPt);
		fOut->Add(particleV2);
		fOut->Add(particleV4);
		fOut->Add(particleV6);
	}
	


  if (leadingJetEtaPhi) {
    canvas->Clear();
    leadingJetEtaPhi->GetXaxis()->SetRangeUser(-eta_cut,eta_cut);
    leadingJetEtaPhi->GetXaxis()->SetTitle("#eta");
    leadingJetEtaPhi->GetYaxis()->SetTitle("#phi");
    leadingJetEtaPhi->Draw("colz");
    canvas->Print("leadingJetEtaPhi.pdf");
    fOut->Add(leadingJetEtaPhi);
  }

  // JetHadron
  canvas->Clear();

  Int_t index = 1;

  for (unsigned int i = 0; i < nTriggerPtBins; i++)
  {
    index = 1;
    canvas->Clear();
    canvas->Divide(TMath::CeilNint(TMath::Sqrt(nAssocParticleBins)), 
        TMath::FloorNint(TMath::Sqrt(nAssocParticleBins)));
    combinedClass = Form("jetHadron_" + jetClass + ".pdf",
        fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));

    for (unsigned int j = 0; j < nAssocParticleBins; j++)
    {
      canvas->cd(index);
      jetHadron[i][j]->SetContour(20);     
      jetHadron[i][j]->Draw("COLZ");
      //		canvas->GetPad(index)->SetLogz();
      //	jetHadron[i][j]->Draw("colz");
      index++;
    }
    canvas->Print(combinedClass);
  }
  //	canvas->Print("jetHadron.pdf");
  // Reset the divide canvas;
  canvas->cd();

  // Fit projections and draw
  // 2 particle correlations
  // It has to be a TH1D according to root documentation
  TH1D * phiProjection = dEtaDPhi->ProjectionY("phiProjection");
  minMax nearSideLimits(-PI/4, PI/4);
  minMax awaySideLimits(3.0*PI/4, 5.0*PI/4);

  // Perform fit
  TF1 * fullFit = phiCorrFit(phiProjection, "fullFit");
  //	TF1 * nearSideFit = fitPeak(phiProjection, nearSideLimits, "nearSide");
  //	TF1 * awaySideFit = fitPeak(phiProjection, awaySideLimits, "awaySide");


  // NEW METHOD
  //	phiCorrFit_2Gaus

  TF1 * fullFit2 = phiCorrFit_2Gaus(phiProjection,"fullFit2");



  // Set fit parameters
  //	nearSideFit->SetLineColor(kRed);
  //	awaySideFit->SetLineColor(kBlue);

  fullFit->SetLineColor(kViolet);
  fullFit->SetLineStyle(2);
  fullFit2->SetLineColor(kViolet);

  // Draw projection and fit
  phiProjection->Draw();
  //	nearSideFit->Draw("same");
  //	awaySideFit->Draw("same");
  //  fullFit->Draw("SAME");
  fullFit2->Draw("SAME");


  // Save projection
  canvas->Print("phiProjection.pdf");

  // For eta within +/- 0.5
  TH1D * phiProjectionEtaCut = dEtaDPhi->ProjectionY("phiProjectionEtaCut", 25, 75);
  TF1 * etaCutFullFit = phiCorrFit(phiProjectionEtaCut, "etaCutFullFit");
  //	TF1 * nearSideEtaCutFit = fitPeak(phiProjectionEtaCut, nearSideLimits, "nearSideEtaCut");
  //	TF1 * awaySideEtaCutFit = fitPeak(phiProjectionEtaCut, awaySideLimits, "awaySideEtaCut");

  // Set fit parameters
  //	nearSideEtaCutFit->SetLineColor(kRed);
  //	awaySideEtaCutFit->SetLineColor(kBlue);

  etaCutFullFit->SetLineColor(kViolet);

  // Draw projection and fit
  phiProjectionEtaCut->Draw();
  //	nearSideEtaCutFit->Draw("same");
  //	awaySideEtaCutFit->Draw("same");
  etaCutFullFit->Draw("SAME");

  // Save projection
  canvas->Print("phiProjectionEtaCut.pdf");

  // Jet Hadron
  minMax nearSideJetHLimits(-PI/6, PI/6);
  minMax awaySideJetHLimits(3.0*PI/4, 5.0*PI/4);

  // Setup Canvas
  canvas->Clear();

  

  // Define Histograms
  TH1D * jetHProjection[nEPBins][nTriggerPtBins][nAssocParticleBins];
  TH1D * jetHdEtaProjection[nEPBins][nTriggerPtBins][nAssocParticleBins];

	// Projections especially for export to GH code
	TH1D * jetHFullDEtaProjection[nEPBins][nTriggerPtBins][nAssocParticleBins];
  TH1D * jetHNearDEtaProjection[nEPBins][nTriggerPtBins][nAssocParticleBins];
  TH1D * jetHFarDEtaProjection[nEPBins][nTriggerPtBins][nAssocParticleBins];





  TH1F * jetHSCdEtaProjection[nEPBins][nTriggerPtBins][nAssocParticleBins];

  TH1D * jetHRecoilProjection[nEPBins][nTriggerPtBins][nAssocParticleBins];
  TH1D * jetHNotRecoilProjection[nEPBins][nTriggerPtBins][nAssocParticleBins];

  // Fit Functions
  TF1 * jetHdEtaProjectionFit[nEPBins][nTriggerPtBins][nAssocParticleBins];
  TF2 * jetHdEtadPhiFit[nEPBins][nTriggerPtBins][nAssocParticleBins];
  TF2 * jetHFit[nEPBins][nTriggerPtBins][nAssocParticleBins]; //one used to subtract
  TF1 * jetHProjectionFit[nEPBins][nTriggerPtBins][nAssocParticleBins];
  TF1 * jetHProjectionFit2[nEPBins][nTriggerPtBins][nAssocParticleBins];
	TH1F * fastMixedHist[nEPBins][nTriggerPtBins][nAssocParticleBins];
  TH2F * jetHadronBkgSub[nEPBins][nTriggerPtBins][nAssocParticleBins];


  TCanvas *cJetHBkgSub = new TCanvas("cJetBkgSub","cJetBkgSub",c_width,c_height);
  TCanvas *cJetHdEtaFit = new TCanvas("cJetHdEtaFit","cJetHdEtaFit",c_width,c_height);
  TCanvas *cJetHdEtadPhiFit = new TCanvas("cJetHdEtadPhiFit","cJetHdEtadPhiFit",c_width,c_height);

  TString combinedNames;




	TF2 * fAccCorr;

	// Acceptance Correction
	if (acceptanceCorrection) {
		TString functionString = "1 * (TMath::Abs(x) <= [1])";
		functionString += "+(TMath::Abs(x)>[1])*(TMath::Abs(x)<=(2*[0]-[1]))*0.5*(2*[0]-[1]-TMath::Abs(x))/([0]-[1])";
		functionString += "+0*y";
		fAccCorr = new TF2("fAccCorr",functionString,jetHadronEP[0][0][0]->GetXaxis()->GetXmin(),jetHadronEP[0][0][0]->GetXaxis()->GetXmax(),jetHadronEP[0][0][0]->GetYaxis()->GetXmin(),jetHadronEP[0][0][0]->GetYaxis()->GetXmax());	
//		fAccCorr = new TF2("fAccCorr",functionString,jetHadron[0][0]->GetXaxis()->GetXmin(),jetHadron[0][0]->GetXaxis()->GetXmax(),jetHadron[0][0]->GetYaxis()->GetXmin(),jetHadron[0][0]->GetYaxis()->GetXmax());	
		fAccCorr->SetParameter(0,eta_cut);							
		fAccCorr->SetParameter(1,R);							

		canvas->Clear();
		canvas->SetLogz(0);
		fAccCorr->Draw("surf");
		canvas->Print("acceptanceFunction.pdf");

	}


	vector <vector<TGraphErrors *> > ptBinASWidthsGraphsEP;
	vector <vector<TGraphErrors *> > ptBinASRmsGraphsEP;
	vector <vector<TGraphErrors *> > ptBinASIntegralsGraphsEP;
	vector <vector<TGraphErrors *> > ptBinASFwhmGraphsEP;
	vector <vector<TGraphErrors *> > ptBinASFwhmRescaledGraphsEP;
	vector <vector<TGraphErrors *> > ptBinASYieldsGraphsEP;
	vector <vector<TGraphErrors *> > ptBinASMeansGraphsEP;
	vector <vector<TGraphErrors *> > ptBinASBetasGraphsEP;
	vector <vector<TGraphErrors *> > ptBinBsGraphsEP;
	vector <vector<TGraph *> > 			ptBinChiSqGraphsEP;
	vector <vector<TGraph *> >			  ptBinChiSqOverNDFGraphsEP;


	vector <vector<TGraphErrors *> > ptBinNSRmsGraphsEP;
	vector <vector<TGraphErrors *> > ptBinNSIntegralsGraphsEP;


	// Loop over event plane bins
	for (int k = 0; k < nEPBins; k++) {



		// Vectors for widths


		// Awayside Variables
		vector < vector <double> > ptBinASWidths;
		vector < vector <double> > ptBinASWidthsErr;
		vector < vector <double> > ptBinASRms;
		vector < vector <double> > ptBinASRmsErr;
		vector < vector <double> > ptBinASFwhm;
		vector < vector <double> > ptBinASFwhmErr;
		vector < vector <double> > ptBinASFwhmRescaled;
		vector < vector <double> > ptBinASFwhmRescaledErr;
		vector < vector <double> > ptBinASYields;
		vector < vector <double> > ptBinASYieldsErr;
		vector < vector <double> > ptBinASIntegrals;
		vector < vector <double> > ptBinASIntegralsErr;
		vector < vector <double> > ptBinASMeans;
		vector < vector <double> > ptBinASMeansErr;
		vector < vector <double> > ptBinASBetas;
		vector < vector <double> > ptBinASBetasErr;

		vector < vector <double> > ptBinNSRms;
		vector < vector <double> > ptBinNSRmsErr;
		vector < vector <double> > ptBinNSIntegrals;
		vector < vector <double> > ptBinNSIntegralsErr;
		vector < vector <double> > ptBinNSYields;
		vector < vector <double> > ptBinNSYieldsErr;


		// Background Variables
		vector < vector <double> > ptBinBs;
		vector < vector <double> > ptBinBsErr;

		// Goodness of Fit
		vector < vector <double> > ptBinChiSq;
		vector < vector <double> > ptBinChiSqOverNDF;


		for (int i = 0; i < nTriggerPtBins; i++ ) {
			vector <double> ASWidthsRow;
			vector <double> ASWidthsErrRow;
			vector <double> ASRmsRow;
			vector <double> ASRmsErrRow;
			vector <double> ASFwhmRow;
			vector <double> ASFwhmErrRow;
			vector <double> ASFwhmRescaledRow;
			vector <double> ASFwhmRescaledErrRow;
			vector <double> ASYieldsRow;
			vector <double> ASYieldsErrRow;
			vector <double> ASIntegralsRow;
			vector <double> ASIntegralsErrRow;
			vector <double> ASMeansRow;
			vector <double> ASMeansErrRow;
			vector <double> ASBetasRow;
			vector <double> ASBetasErrRow;
			vector <double> BsRow;
			vector <double> BsErrRow;
			vector <double> ChiSqRow;
			vector <double> ChiSqOverNDFRow;

			vector <double> NSRmsRow;
			vector <double> NSRmsErrRow;
			vector <double> NSIntegralsRow;
			vector <double> NSIntegralsErrRow;
			vector <double> NSYieldsRow;
			vector <double> NSYieldsErrRow;


			for(int j = 0; j < nAssocParticleBins; j++) {
				ASWidthsRow.push_back(0.0);
				ASWidthsErrRow.push_back(0.0);
				ASRmsRow.push_back(0.0);
				ASRmsErrRow.push_back(0.0);
				ASFwhmRow.push_back(0.0);
				ASFwhmErrRow.push_back(0.0);
				ASFwhmRescaledRow.push_back(0.0);
				ASFwhmRescaledErrRow.push_back(0.0);
				ASMeansRow.push_back(0.0);
				ASMeansErrRow.push_back(0.0);
				ASYieldsRow.push_back(0.0);
				ASYieldsErrRow.push_back(0.0);
				ASIntegralsRow.push_back(0.0);
				ASIntegralsErrRow.push_back(0.0);
				ASBetasRow.push_back(0.0);
				ASBetasErrRow.push_back(0.0);
				BsRow.push_back(0.0);
				BsErrRow.push_back(0.0);
				ChiSqRow.push_back(0.0);
				ChiSqOverNDFRow.push_back(0.0);

				NSIntegralsRow.push_back(0.0);
				NSIntegralsErrRow.push_back(0.0);
				NSYieldsRow.push_back(0.0);
				NSYieldsErrRow.push_back(0.0);
				NSRmsRow.push_back(0.0);
				NSRmsErrRow.push_back(0.0);


			}

			ptBinASWidths.push_back(ASWidthsRow);
			ptBinASWidthsErr.push_back(ASWidthsErrRow);
			ptBinASRms.push_back(ASRmsRow);
			ptBinASRmsErr.push_back(ASRmsErrRow);
			ptBinASFwhm.push_back(ASFwhmRow);
			ptBinASFwhmErr.push_back(ASFwhmErrRow);
			ptBinASFwhmRescaled.push_back(ASFwhmRescaledRow);
			ptBinASFwhmRescaledErr.push_back(ASFwhmRescaledErrRow);
			ptBinASYields.push_back(ASYieldsRow);
			ptBinASYieldsErr.push_back(ASYieldsErrRow);
			ptBinASIntegrals.push_back(ASIntegralsRow);
			ptBinASIntegralsErr.push_back(ASIntegralsErrRow);
			ptBinASMeans.push_back(ASMeansRow);
			ptBinASMeansErr.push_back(ASMeansErrRow);
			ptBinASBetas.push_back(ASBetasRow);
			ptBinASBetasErr.push_back(ASBetasErrRow);
			ptBinBs.push_back(BsRow);
			ptBinBsErr.push_back(BsErrRow);
			ptBinChiSq.push_back(ChiSqRow);
			ptBinChiSqOverNDF.push_back(ChiSqOverNDFRow);

			ptBinNSRms.push_back(NSRmsRow);
			ptBinNSRmsErr.push_back(NSRmsErrRow);
			ptBinNSIntegrals.push_back(NSIntegralsRow);
			ptBinNSIntegralsErr.push_back(NSIntegralsErrRow);
			ptBinNSYields.push_back(NSYieldsRow);
			ptBinNSYieldsErr.push_back(NSYieldsErrRow);

		}

		TF1 * fitSwiftJet[nTriggerPtBins];
		TF1 * fitSwiftParticle[nTriggerPtBins][nAssocParticleBins];
		TF2 * swiftBackground[nTriggerPtBins][nAssocParticleBins];

		TCanvas *cBkgSub = new TCanvas("cBkgSub","cBkgSub",c_width,c_height); // canvas for AS with bkgsub
		TCanvas *cProjCmp = new TCanvas("cProjCmp","cProjCmp",c_width,c_height); // canvas for comparing 
		// all, recoil, notRecoil dPhi projections
		TCanvas *cSwiftJet = new TCanvas("cSwiftJet","cSwiftJet",c_width,c_height);
		TCanvas *cSwiftParticle = new TCanvas("cSwiftParticle","cSwiftParticle",c_width,c_height);
		TCanvas *cSwiftBackground = new TCanvas("cSwiftBackground","cSwiftBackground",c_width,c_height);

		TString sEPLabel = Form("EP_%d",k);

		TString sJetHProjection2 = "";
		TString sBkgSub = "";
		TString sSwiftJet = Form("jetYield_%s.pdf",sEPLabel.Data());
//		TString sSwiftJet = "output/jetYield.pdf";
		TString sSwiftParticle = "";
		TString sSwiftBackground = "";
		TString sJetHBkgSub = "";	
		TString sJetHdEtaFit = "";
		TString sJetHdEtadPhiFit = "";

		cSwiftJet->Divide(TMath::CeilNint(TMath::Sqrt(nTriggerPtBins)),
				TMath::FloorNint(TMath::Sqrt(nTriggerPtBins)));


		TH1F * constantAssocYield;
		if (fastMixed && fastMixingMode == 1) {
			 constantAssocYield = (TH1F *) swiftParticleYield[0][0]->Clone("constantAssocYield");
			for (int n = 1; n <= constantAssocYield->GetNbinsX(); n++) {
				constantAssocYield->SetBinContent(n,1);
			}
		}


		int n_y_SubDiv = TMath::FloorNint(TMath::Sqrt(nAssocParticleBins));
//		int n_x_SubDiv = TMath::CeilNint(TMath::Sqrt(nAssocParticleBins));
		int n_x_SubDiv = TMath::CeilNint(nAssocParticleBins*1.0/n_y_SubDiv);
		printf("DEBUG: n_x = %d, n_y = %d (from nAssocPtBins = %d)\n",n_x_SubDiv,n_y_SubDiv,nAssocParticleBins);
		for (unsigned int i = 0; i < nTriggerPtBins; i++)
		{
			index = 0;
			canvas->Clear();
			canvas->Divide(n_x_SubDiv,n_y_SubDiv);
			cBkgSub->Clear();
			cBkgSub->Divide(n_x_SubDiv,n_y_SubDiv);
			cProjCmp->Clear();
			cProjCmp->Divide(n_x_SubDiv,n_y_SubDiv);
			cSwiftParticle->Clear();
			cSwiftParticle->Divide(n_x_SubDiv,n_y_SubDiv);
			cSwiftBackground->Clear();
			cSwiftBackground->Divide(n_x_SubDiv,n_y_SubDiv);
			cJetHBkgSub->Clear();
			cJetHBkgSub->Divide(n_x_SubDiv,n_y_SubDiv);
			cJetHdEtaFit->Clear();   
			cJetHdEtaFit->Divide(n_x_SubDiv,n_y_SubDiv);
			cJetHdEtadPhiFit->Clear();   
			cJetHdEtadPhiFit->Divide(n_x_SubDiv,n_y_SubDiv);

			combinedClass = Form("jetHProjection_" + jetClass + sEPLabel + ".pdf",
					fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			sBkgSub = Form("jetHProjectionBkgSub_" + jetClass + sEPLabel + ".pdf",
					fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			sSwiftParticle = Form("particleYield_" + jetClass + sEPLabel + ".pdf",
					fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));		
			sSwiftBackground = Form("swiftBackground_" + jetClass + sEPLabel + ".pdf",
					fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));		
			sJetHBkgSub = Form("jetHBkgSub_" + jetClass + sEPLabel,
					fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));		
			sJetHdEtaFit = Form("jetHdEtaFit_" + jetClass + sEPLabel,
					fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));		
			sJetHdEtadPhiFit = Form("jetHdEtadPhiFit_" + jetClass + sEPLabel + ".pdf",
					fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));		

	/*		combinedClass = Form("jetHProjection_" + jetClass + ".pdf",
					jetPtBins.at(i),jetPtBins.at(i+1));
			sBkgSub = Form("jetHProjectionBkgSub_" + jetClass + ".pdf",
					jetPtBins.at(i),jetPtBins.at(i+1));
			sSwiftParticle = Form("particleYield_" + jetClass + ".pdf",
					jetPtBins.at(i),jetPtBins.at(i+1));		
			sSwiftBackground = Form("swiftBackground_" + jetClass + ".pdf",
					jetPtBins.at(i),jetPtBins.at(i+1));		
			sJetHBkgSub = Form("jetHBkgSub_" + jetClass + ".pdf",
					jetPtBins.at(i),jetPtBins.at(i+1));		
			sJetHdEtaFit = Form("jetHdEtaFit_" + jetClass + ".pdf",
					jetPtBins.at(i),jetPtBins.at(i+1));		
			sJetHdEtadPhiFit = Form("jetHdEtadPhiFit_" + jetClass + ".pdf",
					jetPtBins.at(i),jetPtBins.at(i+1));		
*/


			cSwiftJet->cd(i+1);
			swiftJetYield[i]->SetMarkerStyle(1);
			swiftJetYield[i]->SetMarkerColor(6);
			swiftJetYield[i]->SetMarkerSize(1);
			swiftJetYield[i]->SetLineColor(6);

			swiftJetYield[i]->SetTitle("");	
			swiftJetYield[i]->Draw();

			combinedNames = swiftJetYield[i]->GetName();
			combinedNames += "_fit";
			fitSwiftJet[i] = swiftJetFit(swiftJetYield[i],combinedNames.Data());
			fitSwiftJet[i]->SetLineColor(9);
			fitSwiftJet[i]->SetLineStyle(9);
			fitSwiftJet[i]->Draw("SAME");


			for (unsigned int j = 0; j < nAssocParticleBins; j++)
			{
				// Set appropraite part of canvas
				index++;

				printf("======================================================================\n");
				printf("|| Fitting JetHadron Jet %d (%.1f <= p_{T}^{jet} < %.1f GeV/c) Hadron %d (%.1f <= p_{T}^{assoc} < %.1f GeV/c) \n",i,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1),j,assocParticlePtBins.at(j),assocParticlePtBins.at(j+1));
				printf("======================================================================\n");

				// SWIFT
				// -------------------------------------------------------
				cSwiftParticle->cd(index);
				swiftParticleYield[i][j]->SetMarkerStyle(1);
				swiftParticleYield[i][j]->SetMarkerColor(kBlue);
				swiftParticleYield[i][j]->SetMarkerSize(1);
				swiftParticleYield[i][j]->SetLineColor(kBlue);
				swiftParticleYield[i][j]->SetTitle(Form("Particle Yield for %.1f < p_{T}^{jet} < %.f GeV/c, %.2f < p_{T}^{part} < %.2f",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1),assocParticlePtBins.at(j),assocParticlePtBins.at(j+1)));
				swiftParticleYield[i][j]->Draw();


				combinedNames = swiftParticleYield[i][j]->GetName();
				combinedNames += "_fit";
							fitSwiftParticle[i][j] = swiftParticleFit(swiftParticleYield[i][j],combinedNames.Data());
							fitSwiftParticle[i][j]->SetLineColor(kCyan);
							fitSwiftParticle[i][j]->SetLineStyle(kCyan);
	//      			fitSwiftParticle[i][j]->Draw("SAME");





        if (jetHadronScatCentEP[k][i][j]) {
          combinedNames = jetHadronScatCentEP[k][i][j]->GetName();
          combinedNames += "_dEtaProjection";
         // jetHSCdEtaProjection[k][i][j] = jetHadronScatCentEP[k][i][j]->ProjectionX(combinedNames.Data(),dPhi_Bkg_min_bin,dPhi_Bkg_max_bin,"e");
          jetHSCdEtaProjection[k][i][j] = (TH1F *) jetHadronScatCentEP[k][i][j]->ProjectionX(combinedNames.Data());
          jetHSCdEtaProjection[k][i][j]->GetXaxis()->SetTitle("#Delta#eta");
          jetHSCdEtaProjection[k][i][j]->GetYaxis()->SetTitle("1/N_{jet} dN/d#eta");
          jetHSCdEtaProjection[k][i][j]->SetTitle(combinedNames.Data());
          // averaging over number of bins in phi.
 //         jetHSCdEtaProjection[k][i][j]->Scale(1./(1+dPhi_Bkg_max_bin - dPhi_Bkg_min_bin));
        }

				// Create the background functions.

				combinedNames = swiftParticleYield[i][j]->GetName();
				combinedNames += "_swiftBackground";
			//  			swiftBackground[i][j] = createSwiftBackground(fitSwiftJet[i],fitSwiftParticle[i][j],jetHadron[i][j],combinedNames.Data());
		 //   		swiftBackground[i][j]->SetLineColor(6);
				//

				TH2F * fMainHist = jetHadronEP[k][i][j]; // for convenience

				// Apply Acceptance Correction
				if (acceptanceCorrection) {
		//			TF2 * fAccCorr = new TF2("fAccCorr",functionString,jetHadron[0][0]->GetXaxis()->GetXmin(),jetHadron[0][0]->GetXaxis()->GetXmax(),jetHadron[0][0]->GetYaxis()->GetXmin(),jetHadron[0][0]->GetYaxis()->GetXmax());	
//					jetHadron[i][j]->Divide(fAccCorr);
						jetHadronEP[k][i][j]->Divide(fAccCorr);
				}	 				
				if (fastMixed)  { // use fast mixed method
					combinedNames = jetHadronEP[k][i][j]->GetName();
					combinedNames += "_2DFit";
					if (fastMixingMode == 1) {
						fastMixedHist[k][i][j] = convolveHistos(swiftJetYield[i],constantAssocYield,Form("%s_fastMixedDEta",jetHadronEP[k][i][j]->GetName()));
					} else if (fastMixingMode == 2) {
						fastMixedHist[k][i][j] = convolveHistos(swiftJetYield[i],swiftParticleYield[i][j],Form("%s_fastMixedDEta",jetHadronEP[k][i][j]->GetName()));
					} else if (fastMixingMode == 3 && jetHadronScatCentEP[k][i][j]) {
						fastMixedHist[k][i][j] = convolveHistos(swiftJetYield[i],jetHSCdEtaProjection[k][i][j],Form("%s_fastMixedDEta",jetHadronEP[k][i][j]->GetName()));
					}
					// scaling to unity at deta = 0
					fastMixedHist[k][i][j]->Sumw2(false);
					fastMixedHist[k][i][j]->Scale(1./fastMixedHist[k][i][j]->GetBinContent(fastMixedHist[k][i][j]->FindBin(0.0)));


					// FIXME make a tf2 from this to divide the th2
	 //       jetHFit[k][i][j] = dEtaToDPhiDEta(fastMixedHist[k][i][j],combinedNames);

					// FIXME or make a th2, though this is memory intensive
//					TH2F * tempClone = (TH2F *) jetHadronEP[k][i][j]->Clone("temp");
					for (int m = 1; m <= jetHadronEP[k][i][j]->GetNbinsX(); m++) {
						for (int n = 1; n <= jetHadronEP[k][i][j]->GetNbinsY(); n++) {
							TH2F * hist = jetHadronEP[k][i][j];
							double hX = hist->GetXaxis()->GetBinCenter(m);
							double hY = hist->GetYaxis()->GetBinCenter(n);
							int binN = hist->FindFixBin(hX,hY);
		
							// how to check for outside the mixedHistogram
							double mixedVal = fastMixedHist[k][i][j]->GetBinContent(fastMixedHist[k][i][j]->FindFixBin(hX,hY));
							if (mixedVal != 0) {							
								hist->SetBinContent(binN,hist->GetBinContent(binN)/mixedVal);
								hist->SetBinError(binN,hist->GetBinError(binN)/mixedVal);
							}
						}
					}
				}



				


	//			double dPhi_Bkg_min = -PI/2;
//				double dPhi_Bkg_max = PI/2;
				double dPhi_Bkg_min = -PI/3;
				double dPhi_Bkg_max = PI/3;
				int dPhi_Bkg_min_bin = jetHadronEP[k][i][j]->GetYaxis()->FindFixBin(dPhi_Bkg_min);
				int dPhi_Bkg_max_bin = jetHadronEP[k][i][j]->GetYaxis()->FindFixBin(dPhi_Bkg_max);      
				combinedNames = jetHadronEP[k][i][j]->GetName();
				combinedNames += "_dEtaProjection";
				jetHdEtaProjection[k][i][j] = jetHadronEP[k][i][j]->ProjectionX(combinedNames.Data(),dPhi_Bkg_min_bin,dPhi_Bkg_max_bin,"e");
				jetHdEtaProjection[k][i][j]->GetXaxis()->SetTitle("#Delta#eta");
				jetHdEtaProjection[k][i][j]->GetYaxis()->SetTitle("1/N_{jet} dN/d#Delta#eta");
				jetHdEtaProjection[k][i][j]->SetTitle(combinedNames.Data());
				// averaging over number of bins in phi.
				jetHdEtaProjection[k][i][j]->Scale(1./(1+dPhi_Bkg_max_bin - dPhi_Bkg_min_bin));

        

        

        // FIXME experiment
//        jetHdEtaProjection[k][i][j]->Rebin(2);
//        jetHdEtaProjection[k][i][j]->Scale(1./2);

				// FIXME should the 2d histograms be scaled by DEta bin width early?

				// could vary this as a systematic
				// i've been using 1 for a while
				int nRebinFarEta = 1; // 8 // FIXME 

				// Save a projection of the far delta eta region.
//				double fDeltaEtaPeakCut = 0.6; // the cut around the NS peak
				//double fDeltaEtaMax = jetHadronEP[k][i][j]->GetXaxis()->GetXmax();
				int iNegativeDeltaEtaCutBin = fMainHist->GetXaxis()->FindFixBin(-fDeltaEtaPeakCut);
				int iPositiveDeltaEtaCutBin =  fMainHist->GetXaxis()->FindFixBin(fDeltaEtaPeakCut);

				int iMinDeltaEtaBin = fMainHist->GetXaxis()->FindFixBin(-delta_eta_cut);
				int iMaxDeltaEtaBin = fMainHist->GetXaxis()->FindFixBin(delta_eta_cut);

				// should the edge bin be included in the near eta?

				// Rescale so that the histograms are normalized to the number of total triggers, not the
				// number in the EP bin
				double ExportScale = 1;
				double fNJetsTotal = leadingJetPtEP[0]->Integral(jetPt->FindBin(fTriggerPtBins[i]),jetPt->FindBin(fTriggerPtBins[i+1]));
				double fNJetsEPBin = leadingJetPtEP[k]->Integral(jetPt->FindBin(fTriggerPtBins[i]),jetPt->FindBin(fTriggerPtBins[i+1]));
				if (fNJetsTotal > 0 && fNJetsEPBin > 0) {
					ExportScale = fNJetsEPBin / fNJetsTotal;
				}

			//		njets_bin = leadingJetPtEP[k]->Integral(jetPt->FindBin(fTriggerPtBins[i]),jetPt->FindBin(fTriggerPtBins[i+1]));
				//} else { 
			//		njets_bin = leadingJetPtEP[0]->Integral(jetPt->FindBin(fTriggerPtBins[i]),jetPt->FindBin(fTriggerPtBins[i+1]));

				// Produce a full projection especially for export to GH Phase4

				combinedNames = Form("Proj_PtBin%d_EP%d_FullDPhi_ObsBin%d",i+1,k-1,j);
				jetHFullDEtaProjection[k][i][j] = fMainHist->ProjectionY(combinedNames.Data(),iMinDeltaEtaBin,iMaxDeltaEtaBin,"e");
				jetHFullDEtaProjection[k][i][j]->SetTitle("Full #Delta#eta Region (#Delta#phi Projection)");
				double fFullEtaRangeNorm = (1 + iMaxDeltaEtaBin - iMinDeltaEtaBin);
				printf("debug: will scale Full histo by %f\n",fFullEtaRangeNorm);
				jetHFullDEtaProjection[k][i][j]->Scale(1./fFullEtaRangeNorm);
				jetHFullDEtaProjection[k][i][j]->Scale(ExportScale);


				// Producing the Near Eta region Projection
//				combinedNames = jetHadronEP[k][i][j]->GetName();
				combinedNames = Form("Proj_PtBin%d_EP%d_NearEtaDPhi_ObsBin%d",i+1,k-1,j);
//				combinedNames += "_NearEtaProj";
				jetHNearDEtaProjection[k][i][j] = fMainHist->ProjectionY(combinedNames.Data(),iNegativeDeltaEtaCutBin,iPositiveDeltaEtaCutBin,"e");
				jetHNearDEtaProjection[k][i][j]->SetTitle("Near #Delta#eta Region (#Delta#phi Projection)");
				double fNearEtaRangeNorm = (1 + iPositiveDeltaEtaCutBin - iNegativeDeltaEtaCutBin);
				printf("debug: will scale neareta histo by %f\n",fNearEtaRangeNorm);
				jetHNearDEtaProjection[k][i][j]->Scale(1./fNearEtaRangeNorm);
				jetHNearDEtaProjection[k][i][j]->Scale(ExportScale);
			

				// Producing the Far Eta Region Projection
//				combinedNames = jetHadronEP[k][i][j]->GetName();
				combinedNames = Form("Proj_PtBin%d_EP%d_FarEtaDPhi_ObsBin%d",i+1,k-1,j);
//				combinedNames += "_FarEtaProj";
				// 2020Mar16 switch to iNegativeDeltaEtaCutBin - 1
				jetHFarDEtaProjection[k][i][j] = fMainHist->ProjectionY(combinedNames.Data(),iMinDeltaEtaBin,iNegativeDeltaEtaCutBin-1,"e");
				jetHFarDEtaProjection[k][i][j]->SetTitle("Far #Delta#eta Region (#Delta#phi Projection)");
				
				TH1D * hTempFar = fMainHist->ProjectionY("hTempFar",iPositiveDeltaEtaCutBin+1,iMaxDeltaEtaBin,"e"); // I don't want to worry about resetting axes
				jetHFarDEtaProjection[k][i][j]->Add(hTempFar);
				delete hTempFar;

				// normalize by eta range?
//				double fEtaRangeNorm = fMainHist->GetXaxis()->GetBinUpEdge(fMainHist->GetXaxis()->GetNbins()) - fMainHist->GetXaxis()->GetBinLowEdge(1) - (fMainHist->GetXaxis()->GetBinLowEdge(iPositiveDeltaEtaCutBin)-fMainHist->GetXaxis()->GetBinUpEdge(iNegativeDeltaEtaCutBin));
				// Normalize by number of bins?
				//double fEtaRangeNorm = (1 + iMaxDeltaEtaBin - iMinDeltaEtaBin) - (1 + iPositiveDeltaEtaCutBin - iNegativeDeltaEtaCutBin);
				double fEtaRangeNorm = (iMaxDeltaEtaBin - iMinDeltaEtaBin) - (iPositiveDeltaEtaCutBin - iNegativeDeltaEtaCutBin);
				// will be 0 if PosDeltaEtaCutBin is too large
				printf("debug: will scale fareta histo by %f\n",fEtaRangeNorm);
				jetHFarDEtaProjection[k][i][j]->Scale(1./fEtaRangeNorm);
				// Now scaled 1 / NSelectedBins. This is what we want for getting what to subtract from the 2D histogram.
				jetHFarDEtaProjection[k][i][j]->Scale(ExportScale);

				// jetHadronEP already normalized by DeltaPhi bin width

				jetHFarDEtaProjection[k][i][j]->Rebin(nRebinFarEta); // rebinning in delta phi
				jetHFarDEtaProjection[k][i][j]->Scale(1./nRebinFarEta); // now back in dNdPhi

				double fastMixedIntegral = 0;
				double dEtaIntegral = 0;

				combinedNames = jetHadronEP[k][i][j]->GetName();
				combinedNames += "_2DFit";
				switch (bkg2DMethod) {
					case aNoBkgSub:
						jetHFit[k][i][j] = new TF2(combinedNames,"0*x+0*y",-2*eta_cut,2*eta_cut,-PI/2,3*PI/2);

						break;
					case a2DFitSub:
						cJetHdEtadPhiFit->cd(index);
						jetHdEtadPhiFit[k][i][j] = etaPhiCorr_fit(jetHadronEP[k][i][j], Form("%s_1",combinedNames.Data()));
						jetHdEtadPhiFit[k][i][j]->SetContour(20);
						jetHdEtadPhiFit[k][i][j]->Draw("SURF2");
			//			fOut->Add(jetHdEtadPhiFit[k][i][j]);
				

						jetHFit[k][i][j] = (TF2 *) jetHdEtadPhiFit[k][i][j]->Clone(combinedNames);
						jetHFit[k][i][j]->SetParameter("Y_NS_1",0); 
						jetHFit[k][i][j]->SetParameter("Y_NS_2",0); 
						jetHFit[k][i][j]->SetParameter("Y_AS_2",0); 

						break;
					case aDEtaFitSub:
						combinedNames = jetHadronEP[k][i][j]->GetName();
						combinedNames += "_dEtaFit";
						jetHdEtaProjectionFit[k][i][j] = dEta_fit(jetHdEtaProjection[k][i][j], combinedNames.Data(),delta_eta_cut,true);

						jetHFit[k][i][j] = new TF2(combinedNames,"[0]+0*x+0*y",-delta_eta_range,delta_eta_range,-PI/2,3*PI/2);
						jetHFit[k][i][j]->SetParameter(0,jetHdEtaProjectionFit[k][i][j]->GetParameter(0));

						break;
					case aDEtaFarSub:
						combinedNames = jetHadronEP[k][i][j]->GetName();
						combinedNames += "_dEtaFit";
						// FIXME update dEta range? could do that for all cases
						//jetHdEtaProjection[k][i][j]->GetXaxis()->SetRangeUser(-delta_eta_cut,delta_eta_cut);
						jetHdEtaProjectionFit[k][i][j] = dEta_fit(jetHdEtaProjection[k][i][j], combinedNames.Data(),delta_eta_cut,false);

						jetHFit[k][i][j] = new TF2(combinedNames,"[0]+0*x+0*y",-delta_eta_range,delta_eta_range,-PI/2,3*PI/2);
						jetHFit[k][i][j]->SetParameter(0,jetHdEtaProjectionFit[k][i][j]->GetParameter(0));

						break;
					case aDEtaGenGausFitSub:
						combinedNames = jetHadronEP[k][i][j]->GetName();
						combinedNames += "_dEtaFit";
						jetHdEtaProjectionFit[k][i][j] = dEta_GenGaus_fit(jetHdEtaProjection[k][i][j], combinedNames.Data(),true);

						jetHFit[k][i][j] = new TF2(combinedNames,"[0]+0*x+0*y",-delta_eta_range,delta_eta_range,-PI/2,3*PI/2);
						jetHFit[k][i][j]->SetParameter(0,jetHdEtaProjectionFit[k][i][j]->GetParameter(0));
									
						break;
          case aDEtaMixedGausFitSub:
						combinedNames = jetHadronEP[k][i][j]->GetName();
						combinedNames += "_dEtaFit";

            if (particlePtBins.at(j+1) <= jet_constituent_cut) {
              jetHdEtaProjectionFit[k][i][j] = dEta_GenGaus_fit(jetHdEtaProjection[k][i][j], combinedNames.Data(),true);
            } else {
              jetHdEtaProjectionFit[k][i][j] = dEta_fit(jetHdEtaProjection[k][i][j], combinedNames.Data(),delta_eta_cut,true);
            }
  
						jetHFit[k][i][j] = new TF2(combinedNames,"[0]+0*x+0*y",-delta_eta_range,delta_eta_range,-PI/2,3*PI/2);
						jetHFit[k][i][j]->SetParameter(0,jetHdEtaProjectionFit[k][i][j]->GetParameter(0));
      
            break;
					case aDEtaFarZYAMSub: // take the Far Eta region (rebinned and rescaled) and find ZYAM for d^2 N / dDeltaEta dDeltaPhi
						// FIXME do I need to do something with jetHdEtaProjectionFit
						printf("DEtaFar ZYAM found Value %f\n",jetHFarDEtaProjection[k][i][j]->GetBinContent(jetHFarDEtaProjection[k][i][j]->GetMinimumBin()));
						jetHFit[k][i][j] = new TF2(combinedNames,"[0]+0*x+0*y",-delta_eta_range,delta_eta_range,-PI/2,3*PI/2);
						jetHFit[k][i][j]->SetParameter(0,jetHFarDEtaProjection[k][i][j]->GetBinContent(jetHFarDEtaProjection[k][i][j]->GetMinimumBin()));			
						break;
					default:
						fprintf(stderr,"Error: Somehow, invalid choice for dPhi-dEta 2d Background");
						exit(1);
						break;
				}

				cJetHdEtaFit->cd(index);
				jetHdEtaProjection[k][i][j]->GetYaxis()->SetRangeUser(0,jetHdEtaProjection[k][i][j]->GetMaximum());
				jetHdEtaProjection[k][i][j]->Draw();

          printf("HEllo, w0rld\n");
       // if (jetHSCdEtaProjection[k][i][j]) {
        if (jetHadronScatCentEP[k][i][j]) {
          jetHSCdEtaProjection[k][i][j]->SetLineColor(16);
          jetHSCdEtaProjection[k][i][j]->SetMarkerColor(12);

 //         jetHdEtaProjection[k][i][j]->Add(jetHSCdEtaProjection[k][i][j],-1); // FIXME
//				jetHdEtaProjection[k][i][j]->GetYaxis()->SetRangeUser(0,jetHdEtaProjection[k][i][j]->GetMaximum());
//          jetHdEtaProjection[k][i][j]->Draw();

          jetHSCdEtaProjection[k][i][j]->Draw("SAME");
        }
          printf("HEllo, w0rld\n");

//				if (bkg2DMethod == aFastMixed) {

//					fastMixedHist[k][i][j]->SetLineColor(kGreen);
//					fastMixedHist[k][i][j]->Draw("SAME");
//					jetHdEtaProjection[k][i][j]->Divide(fastMixedHist[k][i][j]);
//					jetHdEtaProjection[k][i][j]->Draw("SAME");

	//				fastMixedHist[k][i][j]->Draw("");
				if (fastMixed) fOut->Add(fastMixedHist[k][i][j]);

				if (bkg2DMethod == aDEtaFitSub || bkg2DMethod == aDEtaFarSub || ((bkg2DMethod == aDEtaMixedGausFitSub) && (particlePtBins.at(j+1) > jet_constituent_cut))) {
					jetHdEtaProjectionFit[k][i][j]->SetLineColor(kViolet);
					jetHdEtaProjectionFit[k][i][j]->Draw("SAME");

					// Drawing subparts
					TF1 *doubleGauss_1 = (TF1 *) jetHdEtaProjectionFit[k][i][j]->Clone();
					doubleGauss_1->SetLineColor(kGreen);
					doubleGauss_1->SetLineStyle(2);
					doubleGauss_1->SetParameter(0,0);
					doubleGauss_1->SetParameter(3,0);
					doubleGauss_1->Draw("SAME");
					
					TF1 *doubleGauss_2 = (TF1 *) jetHdEtaProjectionFit[k][i][j]->Clone();
					doubleGauss_2->SetLineColor(kCyan);
					doubleGauss_2->SetLineStyle(2);
					doubleGauss_2->SetParameter(0,0);
					doubleGauss_2->SetParameter(1,0);
					doubleGauss_2->Draw("SAME");

					TF1 *tentFunc = (TF1 *) jetHdEtaProjectionFit[k][i][j]->Clone();
					tentFunc->SetLineColor(kRed);
					tentFunc->SetParameter(1,0);
					tentFunc->SetParameter(3,0);
					tentFunc->Draw("SAME");
				} else if (bkg2DMethod == aDEtaGenGausFitSub || ((bkg2DMethod == aDEtaMixedGausFitSub) && (particlePtBins.at(j+1) <= jet_constituent_cut))) {
					jetHdEtaProjectionFit[k][i][j]->SetLineColor(kViolet);
					jetHdEtaProjectionFit[k][i][j]->Draw("SAME");
					
					// Drawing subparts
					TF1 *genGauss_1 = (TF1 *) jetHdEtaProjectionFit[k][i][j]->Clone();
					genGauss_1->SetLineColor(kGreen);
					genGauss_1->SetLineStyle(2);
					genGauss_1->SetParameter(0,0);
					genGauss_1->Draw("SAME");
					
					TF1 *constFunc = (TF1 *) jetHdEtaProjectionFit[k][i][j]->Clone();
					constFunc->SetLineColor(kRed);
					constFunc->SetParameter(1,0);
					constFunc->Draw("SAME");

				}


				


				jetHFit[k][i][j]->SetTitle(combinedNames);


				jetHadronBkgSub[k][i][j] = (TH2F *) jetHadronEP[k][i][j]->Clone();
				jetHadronBkgSub[k][i][j]->SetName(Form("%s_2DBkgSub",jetHadronBkgSub[k][i][j]->GetName()));

	//      if (!useZYAM) {
	//      }
				fOut->Add(jetHFit[k][i][j]);
				fOut->Add(jetHadronEP[k][i][j]);
//				fOut->Add(jetHdEtaProjectionFit[k][i][j]);


				// Subtracting 2D Fit.
				jetHadronBkgSub[k][i][j]->Add(jetHFit[k][i][j],-1);
				cJetHBkgSub->cd(index);

				// Drawing it a bit nicer
				if (jetHadronBkgSub[k][i][j]->GetXaxis()->GetXmax() > 4) jetHadronBkgSub[k][i][j]->GetXaxis()->SetRangeUser(-4,4);

				jetHadronBkgSub[k][i][j]->Draw("COLZ");
				// FIXME draw the delta_eta_cut and near eta cut lines on top
				TLine * lLeftDeltaEtaCut = new TLine(-delta_eta_cut,-PI/2,-delta_eta_cut,3*PI/2);
				lLeftDeltaEtaCut->Draw("SAME");
				TLine * lRightDeltaEtaCut = new TLine(delta_eta_cut,-PI/2,delta_eta_cut,3*PI/2);
				lRightDeltaEtaCut->Draw("SAME");
				TLine * lLeftDeltaEtaPeakCut = new TLine(-fDeltaEtaPeakCut,-PI/2,-fDeltaEtaPeakCut,3*PI/2);
				lLeftDeltaEtaPeakCut->Draw("SAME");
				TLine * lRightDeltaEtaPeakCut = new TLine(fDeltaEtaPeakCut,-PI/2,fDeltaEtaPeakCut,3*PI/2);
				lRightDeltaEtaPeakCut->Draw("SAME");

				// Box around Far Eta NS Region
				


				// Perform projection in dPhi
				combinedNames = jetHadronEP[k][i][j]->GetName();
				combinedNames += "_Projection";
			//	jetHProjection[k][i][j] = jetHadronBkgSub[k][i][j]->ProjectionY(combinedNames.Data());
				// Now with Delta Phi Cut 
				// FIXME should this use delta_eta_cut
				//int xBinLow = jetHadronBkgSub[k][i][j]->GetXaxis()->FindBin(-fDeltaEtaPeakCut); //  projecting the peak region
				//int xBinHigh = jetHadronBkgSub[k][i][j]->GetXaxis()->FindBin(fDeltaEtaPeakCut);
				int xBinLow = jetHadronBkgSub[k][i][j]->GetXaxis()->FindBin(-delta_eta_cut); // projecting the whole analysis region
				int xBinHigh = jetHadronBkgSub[k][i][j]->GetXaxis()->FindBin(delta_eta_cut);
				jetHProjection[k][i][j] = jetHadronBkgSub[k][i][j]->ProjectionY(combinedNames.Data(),xBinLow,xBinHigh);

				// New: Normalize by nBinsDeltaEta. This is bad for the normalization of the nearside peak
				int nBinsDeltaEta = 1 + xBinHigh - xBinLow;
				printf("Debug: normalizing full projection by %d DeltaEta bins\n",nBinsDeltaEta);
				jetHProjection[k][i][j]->Scale(1./nBinsDeltaEta);
				// Now with Delta Phi Cut 
				




				//	jetHProjection[i][j] = jetHadron[i][j]->ProjectionY(combinedNames.Data());
				jetHProjection[k][i][j]->SetStats(kFALSE);
				jetHProjection[k][i][j]->GetXaxis()->SetTitle("#Delta#phi");
//				jetHProjection[k][i][j]->GetYaxis()->SetTitle("1/N_{jet} dN/d#Delta#phi");
				jetHProjection[k][i][j]->GetYaxis()->SetTitle("1/N_{jet} d^{2}N/d#Delta#phid#Delta#eta"); // FIXME
				jetHProjection[k][i][j]->SetTitle(Form("%s-h, %.0f #leq p^{%s}_{T} < %.0f, %.2f #leq %s < %.2f, EP Bin %d",sTriggerTitle.Data(),fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1),assocParticlePtBins.at(j),sAssocPtTitle.Data(),assocParticlePtBins.at(j+1),k));

				int rebinProj = 2;
				jetHProjection[k][i][j]->Rebin(rebinProj);
				jetHProjection[k][i][j]->Scale(1./rebinProj);



				combinedNames = jetHadronEP[k][i][j]->GetName();
				combinedNames += "_fit";
				printf("Executing sigma fit ...\n");
				switch (dPhiFitMethod) {
					case a1G1G:  // Simplest: 1 Gaussian for each peak
						jetHProjectionFit2[k][i][j] = phiCorrFit_1Gaus_1Gaus(jetHProjection[k][i][j],combinedNames.Data());
						break;
					case a1M1M: // 1 Mod Gaus for NS, 1 Mod Gaus for AS
						// FIXME
						jetHProjectionFit2[k][i][j] = phiCorrFit_1ModGaus_1ModGaus(jetHProjection[k][i][j],combinedNames.Data());
						break;
					case a2G1M: // 2 Gaus for NS, 1 Mod Gaus for AS
						// FIXME
						jetHProjectionFit2[k][i][j] = phiCorrFit_2Gaus_1ModGaus(jetHProjection[k][i][j],combinedNames.Data());
						break;
					case a2G1GG:
						jetHProjectionFit2[k][i][j] = phiCorrFit_2Gaus_1GenGaus(jetHProjection[k][i][j],combinedNames.Data());
						break;
					case a1GG1GG:
						jetHProjectionFit2[k][i][j] = phiCorrFit_1GenGaus_1GenGaus(jetHProjection[k][i][j],combinedNames.Data());
						break;
					case a2G1SH:
						jetHProjectionFit2[k][i][j] = phiCorrFit_2Gaus_SecH(jetHProjection[k][i][j],combinedNames.Data());
						break;
					case a2G2G:
						jetHProjectionFit2[k][i][j] = phiCorrFit_2Gaus_2Gaus(jetHProjection[k][i][j],combinedNames.Data());
						break;
					default:
						fprintf(stderr,"Error: Somehow, invalid choice for dPhi Fit Method");
						exit(1);
						break;
				}

				jetHProjectionFit2[k][i][j]->SetLineColor(kRed);
				jetHProjectionFit2[k][i][j]->SetLineWidth(1);
				jetHProjectionFit2[k][i][j]->SetNpx(200);
				fOut->Add(jetHProjectionFit2[k][i][j]);

				TH1F * clone_proj = (TH1F *) jetHProjection[k][i][j]->Clone(Form("%s_BkgSub",jetHProjection[k][i][j]->GetName()));
	 //     clone_proj->SetName(Form("%s_BkgSub",clone_proj->GetName()));

				clone_proj->GetXaxis()->SetRangeUser(PI/2, 3.*PI/2);

				//TF1 * bkg_near = (TF1 *) jetHProjectionFit[i][j]->Clone();
				TF1 * bkg_near_2 = (TF1 *) jetHProjectionFit2[k][i][j]->Clone();

				// old method: subtracting nearside
			//	bkg_near_2->SetParameter("Y_AS",0);

				// now, just subtracting background term
		//		bkg_near_2->SetParameter("Y_AS",0);
				if (bkg_near_2->GetParNumber("Y_AS") != -1) {
					bkg_near_2->SetParameter("Y_AS",0);
				} else {
					bkg_near_2->SetParameter("Y_AS_1",0);
					bkg_near_2->SetParameter("Y_AS_2",0);
				}
				if (bkg_near_2->GetParNumber("Y_NS") != -1 ) {
					bkg_near_2->SetParameter("Y_NS",0); 
				} else {
					bkg_near_2->SetParameter("Y_NS_1",0);
					bkg_near_2->SetParameter("Y_NS_2",0); 
				}


				clone_proj->Add(bkg_near_2,-1);

				canvas->cd(index);
				// Save values to array
				double sigma_as = jetHProjectionFit2[k][i][j]->GetParameter("Sigma_AS");
				double sigma_as_err = jetHProjectionFit2[k][i][j]->GetParError(jetHProjectionFit2[k][i][j]->GetParNumber("Sigma_AS"));
				ptBinASWidths[i][j] = sigma_as;
				ptBinASWidthsErr[i][j] = sigma_as_err;

				printf("Found: ChiSq = %f \t\t NDF = %d\n",jetHProjectionFit2[k][i][j]->GetChisquare(),jetHProjectionFit2[k][i][j]->GetNDF());

				ptBinChiSq[i][j] = jetHProjectionFit2[k][i][j]->GetChisquare();
				ptBinChiSqOverNDF[i][j] = ptBinChiSq[i][j] / jetHProjectionFit2[k][i][j]->GetNDF();

				clone_proj->GetXaxis()->SetRangeUser(PI - rmsRange, PI + rmsRange);
				ptBinASRms[i][j] = clone_proj->GetRMS();
				ptBinASRmsErr[i][j] = clone_proj->GetRMSError();
				clone_proj->GetXaxis()->SetRangeUser(PI/2, 3.*PI/2);

				Double_t integral_err = 0;
				ptBinASIntegrals[i][j] = clone_proj->IntegralAndError(clone_proj->GetXaxis()->FindBin(PI - rmsRange),clone_proj->GetXaxis()->FindBin(PI + rmsRange),integral_err,"width");
			//	ptBinASIntegrals[i][j] = clone_proj->IntegralAndError(clone_proj->GetXaxis()->FindBin(PI/2),clone_proj->GetXaxis()->FindBin(3.*PI/2),integral_err,"width");
			 // ptBinASIntegrals[i][j] = clone_proj->Integral(clone_proj->GetXaxis()->FindBin(PI/2),clone_proj->GetXaxis()->FindBin(3.*PI/2),&integral_err,"width");
	//MARKER
				ptBinASIntegralsErr[i][j]=integral_err;
	
				// Scaling by bin size
				ptBinASIntegrals[i][j] = ptBinASIntegrals[i][j] / (particlePtBins.at(j+1) - particlePtBins.at(j));
				ptBinASIntegralsErr[i][j] = ptBinASIntegralsErr[i][j] / (particlePtBins.at(j+1) - particlePtBins.at(j));



	//			clone_proj->GetXaxis()->SetRangeUser(-PI/2, PI/2);
				ptBinNSIntegrals[i][j] = clone_proj->IntegralAndError(clone_proj->GetXaxis()->FindBin(0 - rmsRange),clone_proj->GetXaxis()->FindBin(0+rmsRange),integral_err,"width");
	//MARKER
				ptBinNSIntegralsErr[i][j]=integral_err;
	
				// Scaling by bin size
				ptBinNSIntegrals[i][j] = ptBinNSIntegrals[i][j] / (particlePtBins.at(j+1) - particlePtBins.at(j));
				ptBinNSIntegralsErr[i][j] = ptBinNSIntegralsErr[i][j] / (particlePtBins.at(j+1) - particlePtBins.at(j));


			
				clone_proj->GetXaxis()->SetRangeUser(0 - rmsRange, 0 + rmsRange);
				ptBinNSRms[i][j] = clone_proj->GetRMS();
				ptBinNSRmsErr[i][j] = clone_proj->GetRMSError();


				clone_proj->GetXaxis()->SetRangeUser(-PI/2, 3*PI/2);



				combinedNames += "_reparam";

				
				// Choose dPhi Fit Method for the main, Omega/FWHM analysis
				printf("Executing omega fit ... \n");
				switch (dPhiFitMethod) {
					case a1G1G: // 1 Gaus for each peak
						jetHProjectionFit[k][i][j] = phiCorrFit_1Gaus_1Gaus_Fwhm(jetHProjection[k][i][j],combinedNames.Data());
						break;
					case a1M1M: // 1 Mod Gaus for NS, 1 Mod Gaus for AS
						jetHProjectionFit[k][i][j] = phiCorrFit_1ModGaus_1ModGaus_Fwhm(jetHProjection[k][i][j],combinedNames.Data());
						break;
					case a2G1M: // 2 Gaus for NS, 1 Mod Gaus for AS
						jetHProjectionFit[k][i][j] = phiCorrFit_2Gaus_1ModGaus_Fwhm(jetHProjection[k][i][j],combinedNames.Data());
						break;
					case a2G1GG:
						jetHProjectionFit[k][i][j] = phiCorrFit_2Gaus_1GenGaus_Fwhm(jetHProjection[k][i][j],combinedNames.Data());
						break;
					case a1GG1GG:
						jetHProjectionFit[k][i][j] = phiCorrFit_1GenGaus_1GenGaus_Fwhm(jetHProjection[k][i][j],combinedNames.Data());
						break;
					case a2G1SH:
						jetHProjectionFit[k][i][j] = phiCorrFit_2Gaus_SecH_Fwhm(jetHProjection[k][i][j],combinedNames.Data());
						break;
					case a2G2G:
						jetHProjectionFit[k][i][j] = phiCorrFit_2Gaus_2Gaus_Fwhm(jetHProjection[k][i][j],combinedNames.Data());
						break;
					default:
						fprintf(stderr,"Error: Somehow, invalid choice for dPhi Fit Method");
						exit(1);
						break;
				}
				jetHProjectionFit[k][i][j]->SetName(Form("fit1_%s_%d_%d",sEPLabel.Data(),i,j));
		

				double beta_as = 0;
				double beta_as_err = 0;	
				if (jetHProjectionFit[k][i][j]->GetParNumber("Beta_AS") >= 0) { // checking if this formula has a beta variable
					beta_as = jetHProjectionFit[k][i][j]->GetParameter("Beta_AS");
					beta_as_err = jetHProjectionFit[k][i][j]->GetParError(jetHProjectionFit[k][i][j]->GetParNumber("Beta_AS"));
				}
				double B_param = jetHProjectionFit[k][i][j]->GetParameter("B");
				double B_param_err = jetHProjectionFit[k][i][j]->GetParError(jetHProjectionFit[k][i][j]->GetParNumber("B"));


				jetHProjectionFit[k][i][j]->SetLineColor(kViolet);
				jetHProjectionFit[k][i][j]->SetLineWidth(1);
				jetHProjectionFit[k][i][j]->SetNpx(200);
				fOut->Add(jetHProjectionFit[k][i][j]); 
				ptBinASFwhm[i][j] = jetHProjectionFit[k][i][j]->GetParameter("Omega_AS");
				ptBinASFwhmErr[i][j] = jetHProjectionFit[k][i][j]->GetParError(jetHProjectionFit[k][i][j]->GetParNumber("Omega_AS"));

/*
				printf("Found: ChiSq = %f \t\t NDF = %d\n",jetHProjectionFit[k][i][j]->GetChisquare(),jetHProjectionFit[k][i][j]->GetNDF());

				ptBinChiSq[i][j] = jetHProjectionFit[k][i][j]->GetChisquare();
				ptBinChiSqOverNDF[i][j] = ptBinChiSq[i][j] / jetHProjectionFit[k][i][j]->GetNDF();
	*/

				ptBinASFwhmRescaled[i][j] = ptBinASFwhm[i][j] / (2.* TMath::Sqrt(2.* TMath::Log(2.)));
				ptBinASFwhmRescaledErr[i][j] = ptBinASFwhmErr[i][j] / (2.* TMath::Sqrt(2.* TMath::Log(2.)));

				ptBinASYields[i][j] = jetHProjectionFit[k][i][j]->GetParameter("Y_AS"); 
				ptBinASYieldsErr[i][j] = jetHProjectionFit[k][i][j]->GetParError(jetHProjectionFit[k][i][j]->GetParNumber("Y_AS")); 

        // Scaling Yields by bin width;
				ptBinASYields[i][j] = ptBinASYields[i][j]/ (particlePtBins.at(j+1) - particlePtBins.at(j));
				ptBinASYieldsErr[i][j] = ptBinASYieldsErr[i][j] / (particlePtBins.at(j+1) - particlePtBins.at(j));


				ptBinASBetas[i][j] = beta_as;
				ptBinASBetasErr[i][j] = beta_as_err;

				ptBinBs[i][j] = B_param;
				ptBinBsErr[i][j] = B_param_err;


//				ptBinHadronPt[i]->GetXaxis()->SetRangeUser(particlePtBins.at(j),particlePtBins.at(j+1));
				ptBinHadronPtEP[k][i]->GetXaxis()->SetRangeUser(particlePtBins.at(j),particlePtBins.at(j+1));

				// FIXME do this by event plane bin

				ptBinASMeans[i][j] = ptBinHadronPtEP[k][i]->GetMean();
				ptBinASMeansErr[i][j] = ptBinHadronPtEP[k][i]->GetMeanError();

			//	ptBinHadronPt[i]->GetXaxis()->UnZoom();  
				ptBinHadronPtEP[k][i]->GetXaxis()->UnZoom();  

				// Draw projection and fit
			 // jetHProjection[i][j]->GetYaxis()->SetRangeUser(0,1.1*jetHProjection[i][j]->GetBinContent(jetHProjection[i][j]->GetMaximumBin()));
				jetHProjection[k][i][j]->GetYaxis()->SetRangeUser(-0.1*jetHProjection[k][i][j]->GetBinContent(jetHProjection[k][i][j]->GetMaximumBin()),1.1*jetHProjection[k][i][j]->GetBinContent(jetHProjection[k][i][j]->GetMaximumBin()));
				jetHProjection[k][i][j]->Draw();

				// Different fit methods have different parameters:
				// 1M1M has only one NS yield'

				//Sigma Constituents
	//      TF1 *away_side_sigma = (TF1 *) jetHProjectionFit2[i][j]->Clone();
	//      away_side_sigma->SetParameter("B",0);
	//      TF1 *const_bg_sigma = (TF1 *) jetHProjectionFit2[i][j]->Clone();
	 //     TF1 *near_side_1_sigma = (TF1 *) jetHProjectionFit2[i][j]->Clone();
	//      TF1 *near_side_2_sigma = (TF1 *) jetHProjectionFit2[i][j]->Clone();
				
				TF1 *away_side, *away_side_2, *const_bg, *near_side_1, *near_side_2;

				if (plotSigmaFits) {// use Sigma Constituents for plot
					away_side = (TF1 *) jetHProjectionFit2[k][i][j]->Clone();
					away_side_2 = (TF1 *) jetHProjectionFit2[k][i][j]->Clone();
					const_bg = (TF1 *) jetHProjectionFit2[k][i][j]->Clone();
					near_side_1 = (TF1 *) jetHProjectionFit2[k][i][j]->Clone();
					near_side_2 = (TF1 *) jetHProjectionFit2[k][i][j]->Clone();
				} else { // Omega Constituents
					away_side = (TF1 *) jetHProjectionFit[k][i][j]->Clone();
					away_side_2 = (TF1 *) jetHProjectionFit2[k][i][j]->Clone();
					const_bg = (TF1 *) jetHProjectionFit[k][i][j]->Clone();
					near_side_1 = (TF1 *) jetHProjectionFit[k][i][j]->Clone();
					near_side_2 = (TF1 *) jetHProjectionFit[k][i][j]->Clone();
				}
				away_side->SetParameter("B",0);

				


				switch (dPhiFitMethod) {
					case a1G1G:
					case a1GG1GG:
					case a1M1M:
						away_side->SetParameter("Y_NS",0);
						away_side->SetLineColor(kGreen);
						away_side->SetLineStyle(9);
						near_side_1->SetParameter("B",0);
						near_side_1->SetParameter("Y_AS",0);
						near_side_1->SetLineColor(kCyan);
						near_side_1->SetLineStyle(3);
						near_side_1->Draw("SAME");
						const_bg->SetParameter("Y_NS",0);
						const_bg->SetParameter("Y_AS",0);

						break;
					case a2G2G:
						away_side->SetParameter("Y_NS_1",0);
						away_side->SetParameter("Y_NS_2",0);
						away_side->SetParameter("Y_AS_2",0);
						away_side->SetLineColor(kGreen);
						away_side->SetLineStyle(2);
	
						away_side_2->SetParameter("B",0);
						away_side_2->SetParameter("Y_AS_1",0);
						away_side_2->SetParameter("Y_NS_1",0);
						away_side_2->SetParameter("Y_NS_2",0);
						away_side_2->SetLineColor(kSpring);
						away_side_2->SetLineStyle(3);					
						away_side_2->Draw("SAME");	

						near_side_1->SetParameter("B",0);
						near_side_1->SetParameter("Y_NS_2",0);
						near_side_1->SetParameter("Y_AS_1",0);
						near_side_1->SetParameter("Y_AS_2",0);
						near_side_1->SetLineColor(kCyan);
						near_side_1->SetLineStyle(2);
						near_side_1->Draw("SAME");

						near_side_2->SetParameter("B",0);
						near_side_2->SetParameter("Y_NS_1",0);
						near_side_2->SetParameter("Y_AS_1",0);
						near_side_2->SetParameter("Y_AS_2",0);
						near_side_2->SetLineColor(kRed);
						near_side_2->SetLineStyle(3);
						near_side_2->Draw("SAME");

						const_bg->SetParameter("Y_NS_1",0);
						const_bg->SetParameter("Y_NS_2",0);
						const_bg->SetParameter("Y_AS_1",0); 
						const_bg->SetParameter("Y_AS_2",0); 

						break;
					default:
						away_side->SetParameter("Y_NS_1",0);
						away_side->SetParameter("Y_NS_2",0);
						away_side->SetLineColor(kGreen);
						away_side->SetLineStyle(2);

						near_side_1->SetParameter("B",0);
						near_side_1->SetParameter("Y_NS_2",0);
						near_side_1->SetParameter("Y_AS",0);
						near_side_1->SetLineColor(kCyan);
						near_side_1->SetLineStyle(3);
						near_side_1->Draw("SAME");

						near_side_2->SetParameter("B",0);
						near_side_2->SetParameter("Y_NS_1",0);
						near_side_2->SetParameter("Y_AS",0);
						near_side_2->SetLineColor(kRed);
						near_side_2->SetLineStyle(3);
						near_side_2->Draw("SAME");

						const_bg->SetParameter("Y_NS_1",0);
						const_bg->SetParameter("Y_NS_2",0);
						const_bg->SetParameter("Y_AS",0);
					break;
				}
				
			if (plotSigmaFits)  jetHProjectionFit2[k][i][j]->Draw("SAME"); else jetHProjectionFit[k][i][j]->Draw("SAME");
				
				const_bg->SetLineColor(15);
				const_bg->SetLineStyle(3);
				away_side->Draw("SAME");

	//      const_bg->SetParameter("B",0.15);

				const_bg->Update();
				const_bg->Draw("SAME");
	//      jetHProjectionFit2[i][j]->Draw("SAME");

				cBkgSub->cd(index);
				TH1F * clone_proj_2 = (TH1F *) jetHProjection[k][i][j]->Clone(Form("%s_BkgSub_2",jetHProjection[k][i][j]->GetName()));
	 //     clone_proj_2->SetName(Form("%s_BkgSub_2",clone_proj_2->GetName()));
				TF1 *bkg_near = (TF1 *) jetHProjectionFit[k][i][j]->Clone();
				bkg_near->SetParameter("Y_AS",0);
				clone_proj_2->Add(bkg_near,-1);
				clone_proj_2->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);

		//		if (plotSigmaFits) { // use Sigma Fits
				clone_proj->Draw("");
				clone_proj->GetYaxis()->SetRangeUser(0,clone_proj->GetYaxis()->GetXmax());
	//			} else {// use Omega fits
	//				clone_proj_2->Draw("");
	//			}
				away_side->Draw("SAME");

				//			delete bkg_near;

				// FIXME is this right?
			/*	if (index == 
						TMath::CeilNint(TMath::Sqrt(nAssocParticleBins))*
						TMath::FloorNint(TMath::Sqrt(nAssocParticleBins))) 
					index++;*/

			}
			canvas->Print(combinedClass);
			cBkgSub->Print(sBkgSub);
			cSwiftParticle->Print(sSwiftParticle);
			cJetHBkgSub->Print(sJetHBkgSub + ".pdf");
			cJetHdEtaFit->Print(sJetHdEtaFit + ".pdf");
			if (bkg2DMethod == a2DFitSub) cJetHdEtadPhiFit->Print(sJetHdEtadPhiFit);
		}
		cSwiftJet->Print(sSwiftJet);
		// Save projections
		//canvas->Print("jetHProjection.pdf");


		// Create vector for ptBins for TGraph
		std::vector <double> ptBinsForTGraph;	
		// Take the average of the two values for each bin location
		for (unsigned int i = 0; i < nAssocParticleBins; i++) { ptBinsForTGraph.push_back((assocParticlePtBins.at(i) + assocParticlePtBins.at(i+1))/2); }

		// Graph for widths
		vector <TGraphErrors *> ptBinASWidthsGraphs;
		vector <TGraphErrors *> ptBinASRmsGraphs;
		vector <TGraphErrors *> ptBinASIntegralsGraphs;
		vector <TGraphErrors *> ptBinASFwhmGraphs;
		vector <TGraphErrors *> ptBinASFwhmRescaledGraphs;
		vector <TGraphErrors *> ptBinASYieldsGraphs;
		vector <TGraphErrors *> ptBinASMeansGraphs;
		vector <TGraphErrors *> ptBinASBetasGraphs;
		vector <TGraphErrors *> ptBinBsGraphs;
		vector <TGraph *> 			ptBinChiSqGraphs;
		vector <TGraph *> 			ptBinChiSqOverNDFGraphs;

		vector <TGraphErrors *> ptBinNSRmsGraphs;
		vector <TGraphErrors *> ptBinNSIntegralsGraphs;


		for (int i = 0; i < nTriggerPtBins; i++) {
			TGraphErrors * graph = new TGraphErrors(nAssocParticleBins,&ptBinsForTGraph[0], &ptBinASWidths[i][0],0,&ptBinASWidthsErr[i][0]);
			ptBinASWidthsGraphs.push_back(graph);
			graph = new TGraphErrors(nAssocParticleBins,&ptBinsForTGraph[0], &ptBinASRms[i][0],0,&ptBinASRmsErr[i][0]);
			ptBinASRmsGraphs.push_back(graph);
			graph = new TGraphErrors(nAssocParticleBins,&ptBinsForTGraph[0], &ptBinASIntegrals[i][0],0,&ptBinASIntegralsErr[i][0]);
			ptBinASIntegralsGraphs.push_back(graph);
			graph = new TGraphErrors(nAssocParticleBins,&ptBinsForTGraph[0], &ptBinASFwhm[i][0],0,&ptBinASFwhmErr[i][0]);
			ptBinASFwhmGraphs.push_back(graph);
			graph = new TGraphErrors(nAssocParticleBins,&ptBinsForTGraph[0], &ptBinASFwhmRescaled[i][0],0,&ptBinASFwhmRescaledErr[i][0]);
			ptBinASFwhmRescaledGraphs.push_back(graph);
			graph = new TGraphErrors(nAssocParticleBins,&ptBinsForTGraph[0], &ptBinASYields[i][0],0,&ptBinASYieldsErr[i][0]);
			ptBinASYieldsGraphs.push_back(graph);
			graph = new TGraphErrors(nAssocParticleBins,&ptBinsForTGraph[0], &ptBinASMeans[i][0],0,&ptBinASMeansErr[i][0]);
			ptBinASMeansGraphs.push_back(graph);
			graph = new TGraphErrors(nAssocParticleBins,&ptBinsForTGraph[0], &ptBinASBetas[i][0],0,&ptBinASBetasErr[i][0]);
			ptBinASBetasGraphs.push_back(graph);
			graph = new TGraphErrors(nAssocParticleBins,&ptBinsForTGraph[0], &ptBinBs[i][0],0,&ptBinBsErr[i][0]);
			ptBinBsGraphs.push_back(graph);

			graph = new TGraphErrors(nAssocParticleBins,&ptBinsForTGraph[0], &ptBinNSIntegrals[i][0],0,&ptBinNSIntegralsErr[i][0]);
			ptBinNSIntegralsGraphs.push_back(graph);
			graph = new TGraphErrors(nAssocParticleBins,&ptBinsForTGraph[0], &ptBinNSRms[i][0],0,&ptBinNSRmsErr[i][0]);
			ptBinNSRmsGraphs.push_back(graph);



			TGraph * graph2 = new TGraph(nAssocParticleBins,&ptBinsForTGraph[0], &ptBinChiSq[i][0]);
			ptBinChiSqGraphs.push_back(graph2);
			graph2 = new TGraph(nAssocParticleBins,&ptBinsForTGraph[0], &ptBinChiSqOverNDF[i][0]);
			ptBinChiSqOverNDFGraphs.push_back(graph2);
			



		}


		for (int i = 0; i < nTriggerPtBins; i++) {
			ptBinASWidthsGraphs[i]->SetMarkerSize(1);
// FIXME
			//Width Graphs
			combinedClass = Form(jetClass +"_"+ sEPLabel + "_AS_Widths",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinASWidthsGraphs[i]->SetName(combinedClass.Data());
			ptBinASWidthsGraphs[i]->SetTitle("Awayside Widths");
			ptBinASWidthsGraphs[i]->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			ptBinASWidthsGraphs[i]->GetYaxis()->SetTitle("#sigma_{as}");

			//RMS Graphs
			ptBinASRmsGraphs[i]->SetMarkerSize(1);
			combinedClass = Form(jetClass +"_"+ sEPLabel + "_AS_Rms",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinASRmsGraphs[i]->SetName(combinedClass.Data());
			ptBinASRmsGraphs[i]->SetTitle("Away-side RMS");
			ptBinASRmsGraphs[i]->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			ptBinASRmsGraphs[i]->GetYaxis()->SetTitle("RMS");

			ptBinNSRmsGraphs[i]->SetMarkerSize(1);
			combinedClass = Form(jetClass +"_"+ sEPLabel + "_NS_Rms",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinNSRmsGraphs[i]->SetName(combinedClass.Data());
			ptBinNSRmsGraphs[i]->SetTitle("Near-side RMS");
			ptBinNSRmsGraphs[i]->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			ptBinNSRmsGraphs[i]->GetYaxis()->SetTitle("RMS");

			//Integral Graphs
			ptBinASIntegralsGraphs[i]->SetMarkerSize(1);
			combinedClass = Form(jetClass +"_"+ sEPLabel + "_AS_Integrals",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinASIntegralsGraphs[i]->SetName(combinedClass.Data());
			ptBinASIntegralsGraphs[i]->SetTitle("Away-side Yields (Integral)");
			ptBinASIntegralsGraphs[i]->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			ptBinASIntegralsGraphs[i]->GetYaxis()->SetTitle("dN/dp_{T} (GeV/c)^{-1}");

			ptBinNSIntegralsGraphs[i]->SetMarkerSize(1);
			combinedClass = Form(jetClass +"_"+ sEPLabel + "_NS_Integrals",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinNSIntegralsGraphs[i]->SetName(combinedClass.Data());
			ptBinNSIntegralsGraphs[i]->SetTitle("Near-side Yields (Integral)");
			ptBinNSIntegralsGraphs[i]->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			ptBinNSIntegralsGraphs[i]->GetYaxis()->SetTitle("dN/dp_{T} (GeV/c)^{-1}");


			//FWHM Graphs
			ptBinASFwhmGraphs[i]->SetMarkerSize(1);
			combinedClass = Form(jetClass +"_"+ sEPLabel + "_AS_Fwhm",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinASFwhmGraphs[i]->SetName(combinedClass.Data());
			ptBinASFwhmGraphs[i]->SetTitle("Awayside FWHM");
			ptBinASFwhmGraphs[i]->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			ptBinASFwhmGraphs[i]->GetYaxis()->SetTitle("FWHM");
			// Rescaled Version
			ptBinASFwhmRescaledGraphs[i]->SetMarkerSize(1);
			combinedClass = Form(jetClass +"_"+ sEPLabel + "_AS_Fwhm_Rescaled",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinASFwhmRescaledGraphs[i]->SetName(combinedClass.Data());
			ptBinASFwhmRescaledGraphs[i]->SetTitle("Awayside FWHM");
			ptBinASFwhmRescaledGraphs[i]->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			ptBinASFwhmRescaledGraphs[i]->GetYaxis()->SetTitle("FWHM");


			//Yield Graphs
			combinedClass = Form(jetClass +"_"+ sEPLabel + "_AS_Yields",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinASYieldsGraphs[i]->SetName(combinedClass.Data());
			ptBinASYieldsGraphs[i]->SetTitle("Awayside Yields");
			ptBinASYieldsGraphs[i]->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			ptBinASYieldsGraphs[i]->GetYaxis()->SetTitle("Y_{as}");

			//Mean Graphs
			combinedClass = Form(jetClass +"_"+ sEPLabel + "_Means",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinASMeansGraphs[i]->SetName(combinedClass.Data());
			ptBinASMeansGraphs[i]->SetTitle("Awayside Means");
			ptBinASMeansGraphs[i]->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			ptBinASMeansGraphs[i]->GetYaxis()->SetTitle("<p_{T}> GeV/c");

			//Beta Graphs
			combinedClass = Form(jetClass +"_"+ sEPLabel + "_AS_Betas",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinASBetasGraphs[i]->SetName(combinedClass.Data());
			ptBinASBetasGraphs[i]->SetTitle("Awayside Betas");
			ptBinASBetasGraphs[i]->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			ptBinASBetasGraphs[i]->GetYaxis()->SetTitle("#beta");

			//Background Graphs
			combinedClass = Form(jetClass +"_"+ sEPLabel + "_Bs",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinBsGraphs[i]->SetName(combinedClass.Data());
			ptBinBsGraphs[i]->SetTitle("Background Parameter");
			ptBinBsGraphs[i]->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			ptBinBsGraphs[i]->GetYaxis()->SetTitle("B");


			// Chi Squared Graphs
			combinedClass = Form(jetClass +"_"+ sEPLabel + "_ChiSq",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinChiSqGraphs[i]->SetName(combinedClass.Data());
			ptBinChiSqGraphs[i]->SetTitle("Fit ChiSquared");
			ptBinChiSqGraphs[i]->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			ptBinChiSqGraphs[i]->GetYaxis()->SetTitle("#Chi^{2}");

			combinedClass = Form(jetClass +"_"+ sEPLabel + "_ChiSqOverNDF",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			ptBinChiSqOverNDFGraphs[i]->SetName(combinedClass.Data());
			ptBinChiSqOverNDFGraphs[i]->SetTitle("Fit ChiSquared Over NDF");
			ptBinChiSqOverNDFGraphs[i]->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			ptBinChiSqOverNDFGraphs[i]->GetYaxis()->SetTitle("#Chi^{2}/NDF");


		}

		canvas->Clear();

		// Adding some plots


		// Widths Method comparison plot
		// Sigma, RMS, FWHM_rescaled
		TLegend *legAS = new TLegend(0.50,0.65,0.90,0.85); 
		double max_y = 0;
		double min_y = 0;
		TString style = "ALP";

		for(int i = 0; i < nTriggerPtBins; i++) {
			max_y = (min_y = 0);

			ptBinASWidthsGraphs[i]->SetLineColor(kBlack);
			ptBinASWidthsGraphs[i]->SetMarkerColor(kBlack);
			ptBinASWidthsGraphs[i]->SetMarkerStyle(33);

			ptBinASRmsGraphs[i]->SetLineColor(kRed);
			ptBinASRmsGraphs[i]->SetMarkerColor(kRed);
			ptBinASRmsGraphs[i]->SetMarkerStyle(33);

			ptBinASFwhmRescaledGraphs[i]->SetLineColor(kBlue);
			ptBinASFwhmRescaledGraphs[i]->SetMarkerColor(kBlue);
			ptBinASFwhmRescaledGraphs[i]->SetMarkerStyle(33);

			min_y = min(TMath::MinElement(ptBinASWidthsGraphs[i]->GetN(),ptBinASWidthsGraphs[i]->GetY()),TMath::MinElement(ptBinASRmsGraphs[i]->GetN(),ptBinASRmsGraphs[i]->GetY()));
			min_y = min(min_y,TMath::MinElement(ptBinASFwhmGraphs[i]->GetN(),ptBinASFwhmGraphs[i]->GetY()));
			max_y = max(TMath::MaxElement(ptBinASWidthsGraphs[i]->GetN(),ptBinASWidthsGraphs[i]->GetY()),TMath::MaxElement(ptBinASRmsGraphs[i]->GetN(),ptBinASRmsGraphs[i]->GetY()));
			max_y = max(max_y,TMath::MaxElement(ptBinASFwhmGraphs[i]->GetN(),ptBinASFwhmGraphs[i]->GetY()));

			legAS->AddEntry(ptBinASWidthsGraphs[i],"sigma", "lp");
			legAS->AddEntry(ptBinASRmsGraphs[i],"RMS after NS Subtraction", "lp");
			legAS->AddEntry(ptBinASFwhmRescaledGraphs[i],"FWHM (Rescaled by 1/2.35 )", "lp");

			//ptBinASWidthsGraphs[i]->GetYaxis()->SetRangeUser(min_y,max_y);
			ptBinASWidthsGraphs[i]->GetYaxis()->SetRangeUser(0.,max_y);
			ptBinASWidthsGraphs[i]->Draw("ALP");		
			ptBinASRmsGraphs[i]->Draw("SAME LP");		
			ptBinASFwhmRescaledGraphs[i]->Draw("SAME LP");		

			combinedClass = Form("as_Width_cmp_" + jetClass + sEPLabel + ".pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			legAS->Draw("SAME");
			canvas->Print(combinedClass);
			legAS->Clear();
			canvas->Clear();
		}
		canvas->SetLogy(0);





		// Widths Plot
		canvas->Clear();
		max_y = (min_y = 0);
		style = "ALP";
		for(int i = 0; i < nTriggerPtBins; i++) {
			ptBinASWidthsGraphs[i]->SetLineColor(cTriggerColorList[i]);
			ptBinASWidthsGraphs[i]->SetMarkerColor(cTriggerColorList[i]);
			ptBinASWidthsGraphs[i]->SetMarkerStyle(33);
			min_y = min(min_y,TMath::MinElement(ptBinASWidthsGraphs[i]->GetN(),ptBinASWidthsGraphs[i]->GetY()));
			max_y = max(max_y,TMath::MaxElement(ptBinASWidthsGraphs[i]->GetN(),ptBinASWidthsGraphs[i]->GetY()));
			if (i == 1) style = "SAME LP";
			ptBinASWidthsGraphs[i]->Draw(style);		
			combinedClass = Form("%0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1));
			legAS->AddEntry(ptBinASWidthsGraphs[i],combinedClass, "lp");
		}
		ptBinASWidthsGraphs[0]->GetYaxis()->SetRangeUser(min_y,max_y);
		//  canvas->SetLogy(1);
		legAS->Draw("SAME");
		canvas->Print(Form("as%s_widths.pdf",sEPLabel.Data()));
		canvas->Clear();
		legAS->Clear();
		canvas->SetLogy(0);

		// RMS plots
		max_y = (min_y = 0);
		style = "ALP";
		for(int i = 0; i < nTriggerPtBins; i++) {
			ptBinASRmsGraphs[i]->SetLineColor(cTriggerColorList[i]);
			ptBinASRmsGraphs[i]->SetMarkerColor(cTriggerColorList[i]);
			ptBinASRmsGraphs[i]->SetMarkerStyle(33);
			min_y = min(min_y,TMath::MinElement(ptBinASRmsGraphs[i]->GetN(),ptBinASRmsGraphs[i]->GetY()));
			max_y = max(max_y,TMath::MaxElement(ptBinASRmsGraphs[i]->GetN(),ptBinASRmsGraphs[i]->GetY()));
			if (i == 1) style = "SAME LP";
			ptBinASRmsGraphs[i]->Draw(style);		
			combinedClass = Form("%0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1));
			legAS->AddEntry(ptBinASRmsGraphs[i],combinedClass, "lp");
		}
		ptBinASRmsGraphs[0]->GetYaxis()->SetRangeUser(min_y,max_y);
		//  canvas->SetLogy(1);
		legAS->Draw("SAME");
		canvas->Print(Form("as%s_RMS.pdf",sEPLabel.Data()));
		canvas->Clear();
		legAS->Clear();
		canvas->SetLogy(0);


		// Yields Method comparison plot
		// Yield parameter, Integral of bkg-sub plot
		legAS = new TLegend(0.50,0.65,0.90,0.85); 
		max_y = 0;
		min_y = 0;
		style = "ALP";

		for(int i = 0; i < nTriggerPtBins; i++) {
			max_y = (min_y = 0);

			ptBinASYieldsGraphs[i]->SetLineColor(kBlack);
			ptBinASYieldsGraphs[i]->SetMarkerColor(kBlack);
			ptBinASYieldsGraphs[i]->SetMarkerStyle(33);

			ptBinASIntegralsGraphs[i]->SetLineColor(kRed);
			ptBinASIntegralsGraphs[i]->SetMarkerColor(kRed);
			ptBinASIntegralsGraphs[i]->SetMarkerStyle(33);

			min_y = min(TMath::MinElement(ptBinASYieldsGraphs[i]->GetN(),ptBinASYieldsGraphs[i]->GetY()),TMath::MinElement(ptBinASIntegralsGraphs[i]->GetN(),ptBinASIntegralsGraphs[i]->GetY()));
			max_y = max(TMath::MaxElement(ptBinASYieldsGraphs[i]->GetN(),ptBinASYieldsGraphs[i]->GetY()),TMath::MaxElement(ptBinASIntegralsGraphs[i]->GetN(),ptBinASIntegralsGraphs[i]->GetY()));

			legAS->AddEntry(ptBinASYieldsGraphs[i],"Yield Parameter", "lp");
			legAS->AddEntry(ptBinASIntegralsGraphs[i],"Integral after NS Subtraction", "lp");

			ptBinASYieldsGraphs[i]->GetYaxis()->SetRangeUser(0.,max_y);
			ptBinASYieldsGraphs[i]->Draw("ALP");		
			ptBinASIntegralsGraphs[i]->Draw("SAME LP");		

			legAS->Draw("SAME");
			combinedClass = Form("as_Yield_cmp_" + jetClass + sEPLabel + ".pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			canvas->Print(combinedClass);
			legAS->Clear();
			canvas->Clear();
		}
		canvas->SetLogy(0);



		// Integral Yield plots
		max_y = (min_y = 0);
		style = "ALP";
		for(int i = 0; i < nTriggerPtBins; i++) {
			ptBinASIntegralsGraphs[i]->SetLineColor(cTriggerColorList[i]);
			ptBinASIntegralsGraphs[i]->SetMarkerColor(cTriggerColorList[i]);
			ptBinASIntegralsGraphs[i]->SetMarkerStyle(33);
			min_y = min(min_y,TMath::MinElement(ptBinASIntegralsGraphs[i]->GetN(),ptBinASIntegralsGraphs[i]->GetY()));
			max_y = max(max_y,TMath::MaxElement(ptBinASIntegralsGraphs[i]->GetN(),ptBinASIntegralsGraphs[i]->GetY()));
			if (i == 1) style = "SAME LP";
			ptBinASIntegralsGraphs[i]->Draw(style);		
			combinedClass = Form("%0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1));
			legAS->AddEntry(ptBinASIntegralsGraphs[i],combinedClass, "lp");
		}
		ptBinASIntegralsGraphs[0]->GetYaxis()->SetRangeUser(min_y,max_y);
		//  canvas->SetLogy(1);
		legAS->Draw("SAME");
		canvas->Print(Form("as%s_Integrals.pdf",sEPLabel.Data()));
		canvas->Print(Form("as%s_Integrals.C",sEPLabel.Data()));
		canvas->Clear();
		legAS->Clear();
		canvas->SetLogy(0);



		// FWHM plots
		max_y = (min_y = 0);
		style = "ALP";
		for(int i = 0; i < nTriggerPtBins; i++) {
			ptBinASFwhmGraphs[i]->SetLineColor(cTriggerColorList[i]);
			ptBinASFwhmGraphs[i]->SetMarkerColor(cTriggerColorList[i]);
			ptBinASFwhmGraphs[i]->SetMarkerStyle(33);
			min_y = min(min_y,TMath::MinElement(ptBinASFwhmGraphs[i]->GetN(),ptBinASFwhmGraphs[i]->GetY()));
			max_y = max(max_y,TMath::MaxElement(ptBinASFwhmGraphs[i]->GetN(),ptBinASFwhmGraphs[i]->GetY()));
			if (i == 1) style = "SAME LP";
			ptBinASFwhmGraphs[i]->Draw(style);		
			combinedClass = Form("%0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1));
			legAS->AddEntry(ptBinASFwhmGraphs[i],combinedClass, "lp");
		}
		ptBinASFwhmGraphs[0]->GetYaxis()->SetRangeUser(min_y,max_y);
		//  canvas->SetLogy(1);
		legAS->Draw("SAME");
		canvas->Print(Form("as%s_FWHM.pdf",sEPLabel.Data()));
		canvas->Clear();
		legAS->Clear();
		canvas->SetLogy(0);


		// Yields Plot
		max_y = (min_y = 0);
		style = "ALP";
		for(int i = 0; i < nTriggerPtBins; i++) {
			ptBinASYieldsGraphs[i]->SetLineColor(cTriggerColorList[i]);
			ptBinASYieldsGraphs[i]->SetMarkerColor(cTriggerColorList[i]);
			ptBinASYieldsGraphs[i]->SetMarkerStyle(33);
			min_y = min(min_y,TMath::MinElement(ptBinASYieldsGraphs[i]->GetN(),ptBinASYieldsGraphs[i]->GetY()));
			max_y = max(max_y,TMath::MaxElement(ptBinASYieldsGraphs[i]->GetN(),ptBinASYieldsGraphs[i]->GetY()));
			if (i == 1) style = "SAME LP";
			ptBinASYieldsGraphs[i]->Draw(style);		
			combinedClass = Form("%0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1));
			legAS->AddEntry(ptBinASYieldsGraphs[i],combinedClass, "lp");
		}
		ptBinASYieldsGraphs[0]->GetYaxis()->SetRangeUser(min_y,max_y);
		legAS->Draw("SAME");
		//	canvas->SetLogy(1);
		canvas->Print(Form("as%s_yields.pdf",sEPLabel.Data()));
		//	canvas->SetLogy(0);
		canvas->Clear();
		legAS->Clear();


		// Print Beta Parameters for the awayside
		max_y = 0;
		min_y = 1;
		style = "ALP";
		for(int i = 0; i < nTriggerPtBins; i++) {
			ptBinASBetasGraphs[i]->SetLineColor(cTriggerColorList[i]);
			ptBinASBetasGraphs[i]->SetMarkerColor(cTriggerColorList[i]);
			ptBinASBetasGraphs[i]->SetMarkerStyle(33);
			min_y = min(min_y,TMath::MinElement(ptBinASBetasGraphs[i]->GetN(),ptBinASBetasGraphs[i]->GetY()));
			max_y = max(max_y,TMath::MaxElement(ptBinASBetasGraphs[i]->GetN(),ptBinASBetasGraphs[i]->GetY()));
			if (i == 1) style = "SAME LP";
			ptBinASBetasGraphs[i]->Draw(style);		
			combinedClass = Form("%0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1));
			legAS->AddEntry(ptBinASBetasGraphs[i],combinedClass, "lp");
		}
		ptBinASBetasGraphs[0]->GetYaxis()->SetRangeUser(min_y,max_y);
		legAS->Draw("SAME");
			canvas->SetLogy(1);
		canvas->Print(Form("as%s_Betas.pdf",sEPLabel.Data()));
			canvas->SetLogy(0);
		canvas->Clear();

		// Print Background Parameters
		legAS->Clear();
		max_y = 0;
		min_y = 1;
		style = "ALP";
		for(int i = 0; i < nTriggerPtBins; i++) {
			ptBinBsGraphs[i]->SetLineColor(cTriggerColorList[i]);
			ptBinBsGraphs[i]->SetMarkerColor(cTriggerColorList[i]);
			ptBinBsGraphs[i]->SetMarkerStyle(33);
	//    min_y = min(min_y,TMath::MinElement(ptBinBsGraphs[i]->GetN(),ptBinBsGraphs[i]->GetY()));
			max_y = max(max_y,TMath::MaxElement(ptBinBsGraphs[i]->GetN(),ptBinBsGraphs[i]->GetY()));
			if (i == 1) style = "SAME LP";
			ptBinBsGraphs[i]->Draw(style);		
			combinedClass = Form("%0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1));
			legAS->AddEntry(ptBinBsGraphs[i],combinedClass, "lp");
		}
		min_y = 1e-4;
		ptBinBsGraphs[0]->GetYaxis()->SetRangeUser(min_y,max_y);
		legAS->Draw("SAME");
			canvas->SetLogy(1);
		canvas->Print(Form("as%s_Bs.pdf",sEPLabel.Data()));
			canvas->SetLogy(0);
		canvas->Clear();


		for (int i = 0; i < nTriggerPtBins; i++ ) {
			swiftJetYield[i]->SetName(Form("jetYield%s_jetPt_%.1f_%.1f",sEPLabel.Data(),fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
			fOut->Add(swiftJetYield[i]);
		}

		cSwiftParticle->Clear();
		cSwiftParticle->cd();
		for (int j = 0; j < nAssocParticleBins; j++) {
			for (int i = 1; i < nTriggerPtBins; i++ ) {
				swiftParticleYield[0][j]->Add(swiftParticleYield[i][j]);
			}
			swiftParticleYield[0][j]->SetName(Form("particleYield%s_particlePt_%.2f_%.2f",sEPLabel.Data(),assocParticlePtBins.at(j),assocParticlePtBins.at(j+1)));
			swiftParticleYield[0][j]->Draw();
			cSwiftParticle->Print(Form("%s.pdf",swiftParticleYield[0][j]->GetName()));
			fOut->Add(swiftParticleYield[0][j]);
		}


			ptBinASWidthsGraphsEP.push_back(ptBinASWidthsGraphs);
			ptBinASRmsGraphsEP.push_back(ptBinASRmsGraphs);
			ptBinASIntegralsGraphsEP.push_back(ptBinASIntegralsGraphs);
			ptBinASFwhmGraphsEP.push_back(ptBinASFwhmGraphs);
			ptBinASYieldsGraphsEP.push_back(ptBinASYieldsGraphs);
			ptBinASMeansGraphsEP.push_back(ptBinASMeansGraphs);
			ptBinASBetasGraphsEP.push_back(ptBinASBetasGraphs);
			ptBinBsGraphsEP.push_back(ptBinBsGraphs);
			ptBinChiSqGraphsEP.push_back(ptBinChiSqGraphs);
			ptBinChiSqOverNDFGraphsEP.push_back(ptBinChiSqOverNDFGraphs);

			ptBinNSRmsGraphsEP.push_back(ptBinNSRmsGraphs);
			ptBinNSIntegralsGraphsEP.push_back(ptBinNSIntegralsGraphs);


	}

	if (leadingJetPt) fOut->Add(leadingJetPt);
  // Drawing, Saving Hard Scatter Things
  if (qqbar_pt) {
    fOut->Add(qqbar_pt);
    fOut->Add(qq_pt);
    fOut->Add(gg_pt);
    fOut->Add(qg_pt);
  }
  if (qqbar_PartonPtByJetBin[0]) {
    TH1F * partonPtByJetBin[nJetPtBins];
    for (int i = 0; i < nTriggerPtBins; i++ ) {

      partonPtByJetBin[i] = (TH1F *) qqbar_PartonPtByJetBin[i]->Clone();
      partonPtByJetBin[i]->SetName(Form("partonPtByJetBin_%.0f_%.0f",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
      partonPtByJetBin[i]->SetTitle(Form("Parton p_{T} for %.0f #leq p_{T}^{jet} < %.0f (GeV/c)",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
      partonPtByJetBin[i]->GetYaxis()->SetTitle("1/N_{jet} dN/dp_{T}");
      partonPtByJetBin[i]->SetLineColor(kBlue);
      partonPtByJetBin[i]->Add(qq_PartonPtByJetBin[i]);
      partonPtByJetBin[i]->Add(gq_PartonPtByJetBin[i]);
      partonPtByJetBin[i]->Add(qg_PartonPtByJetBin[i]);
      partonPtByJetBin[i]->Add(gg_PartonPtByJetBin[i]);

      // normalizing to unity
      if (double pIntegral = partonPtByJetBin[i]->Integral("width")) partonPtByJetBin[i]->Scale(1./pIntegral);
      // normalizing to 1/N_{jet}  ?
     // if (int nJetEntries = partonPtByJetBin[i]->GetEntries()) partonPtByJetBin[i]->Scale(1./nJetEntries);

    
      canvas->Clear();  
      partonPtByJetBin[i]->Draw();
      canvas->Print(Form("partonPtByJetBin_jetPt_%.0f_%.0f.pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));


    }



  } // End of Event Plane Loop

	if (ApplySystUncertBeforeRatio && SystematicUncertInFilePath.Length() && nEPBins > 1) {
		printf("Leading Systematic Uncertainties from Event Plane analysis.\n");	

		gSystem->ChangeDirectory(initialDirectory.Data());
		TFile * fSystematicUncertInFile = TFile::Open(SystematicUncertInFilePath.Data(),"READ");
		if (!fSystematicUncertInFile) {
			fprintf(stderr,"Error: Unable to open systematic error input file!\n");
			exit(1);
		}
		gSystem->ChangeDirectory(outputDirPath.Data());
	

		ApplySystematicUncertainty(ptBinASYieldsGraphsEP,fSystematicUncertInFile);
		ApplySystematicUncertainty(ptBinNSIntegralsGraphsEP,fSystematicUncertInFile);
		ApplySystematicUncertainty(ptBinASIntegralsGraphsEP,fSystematicUncertInFile);

		ApplySystematicUncertainty(ptBinASWidthsGraphsEP,fSystematicUncertInFile);
		ApplySystematicUncertainty(ptBinNSRmsGraphsEP,fSystematicUncertInFile);
		ApplySystematicUncertainty(ptBinASRmsGraphsEP,fSystematicUncertInFile);

	}



	for (int k = 0; k < nEPBins; k++) {
		for (int i = 0; i < nTriggerPtBins; i++) {

			fOut->Add(ptBinASWidthsGraphsEP[k][i]);
			fOut->Add(ptBinASRmsGraphsEP[k][i]);
//			fOut->Add(ptBinASIntegralsGraphsEP[k][i]);
			fOut->Add(ptBinASFwhmGraphsEP[k][i]);
			fOut->Add(ptBinASYieldsGraphsEP[k][i]);
			fOut->Add(ptBinASMeansGraphsEP[k][i]);
			fOut->Add(ptBinASBetasGraphsEP[k][i]);
			fOut->Add(ptBinBsGraphsEP[k][i]);
			fOut->Add(ptBinChiSqGraphsEP[k][i]);
			fOut->Add(ptBinChiSqOverNDFGraphsEP[k][i]);

			fOut->Add(ptBinNSRmsGraphsEP[k][i]);
//			fOut->Add(ptBinNSIntegralsGraphsEP[k][i]);


		}
	}


	// Making Event Plane Comparison Plots

	// Yield Ratios
	vector<TGraphErrors *> OutOverIn_AS_JetPtBin;
	vector<TGraphErrors *> MidOverIn_AS_JetPtBin;
	vector<TGraphErrors *> OutOverIn_NS_JetPtBin;
	vector<TGraphErrors *> MidOverIn_NS_JetPtBin;

	// Parameter Method
	vector<TGraphErrors *> OutOverInPar_AS_JetPtBin;
	vector<TGraphErrors *> MidOverInPar_AS_JetPtBin;

	// Yield Differences
	vector<TGraphErrors *> OutMinusIn_AS_JetPtBin;
	vector<TGraphErrors *> MidMinusIn_AS_JetPtBin;
	vector<TGraphErrors *> OutMinusIn_NS_JetPtBin;
	vector<TGraphErrors *> MidMinusIn_NS_JetPtBin;

	// RMS
	vector<TGraphErrors *> RMSOutOverIn_AS_JetPtBin;
	vector<TGraphErrors *> RMSMidOverIn_AS_JetPtBin;
	vector<TGraphErrors *> RMSOutOverIn_NS_JetPtBin;
	vector<TGraphErrors *> RMSMidOverIn_NS_JetPtBin;

	// Widths (Fit Parameter Sigma)
	vector<TGraphErrors *> WidthsOutOverIn_AS_JetPtBin;
	vector<TGraphErrors *> WidthsMidOverIn_AS_JetPtBin;
//	vector<TGraphErrors *> WidthsOutOverIn_NS_JetPtBin;
//	vector<TGraphErrors *> WidthsMidOverIn_NS_JetPtBin;


	if (nEPBins == 4) {
		// Could move this up
		Int_t colorList[4] = {kBlack, kBlue-4, kGreen-3, kRed+1};
		Int_t markerList[4] = {kOpenSquare,kFullSquare,kFullTriangleUp,kFullCircle};
		vector <TString> titleList = {"All","In-Plane","Mid-Plane","Out-of-Plane"};		

		canvas->Clear();
		TLegend *legAS = new TLegend(0.50,0.65,0.90,0.85); 

		TString style = "AP";


		// Parton Pt By Jet Bin Plots
		for (int i = 0; i < nJetPtBins; i++) {
			
		}


		// AwaySide

		// Integral Yield plots vs EP Bin
		for(int i = 0; i < nTriggerPtBins; i++) {
			TMultiGraph * mg = new TMultiGraph();
			for(int k = 0; k < nEPBins; k++) {
				ptBinASIntegralsGraphsEP[k][i]->SetLineColor(colorList[k]);
				ptBinASIntegralsGraphsEP[k][i]->SetMarkerColor(colorList[k]);
				ptBinASIntegralsGraphsEP[k][i]->SetMarkerStyle(markerList[k]);
				mg->Add(ptBinASIntegralsGraphsEP[k][i]);
				legAS->AddEntry(ptBinASIntegralsGraphsEP[k][i],titleList[k], "lp");
			}
			canvas->SetLogy(1);
			mg->Draw(style);
			mg->SetTitle(Form("%0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			mg->GetXaxis()->SetTitle(ptBinASIntegralsGraphsEP[0][i]->GetXaxis()->GetTitle());
			mg->GetYaxis()->SetTitle(ptBinASIntegralsGraphsEP[0][i]->GetYaxis()->GetTitle());
	//		ptBinASIntegralsGraphs[0]->GetYaxis()->SetRangeUser(min_y,max_y);
			legAS->Draw("SAME");
			canvas->Print(Form("EP_AS_Integrals_JetPt_%0.f_%0.f.pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
			canvas->Print(Form("EP_AS_Integrals_JetPt_%0.f_%0.f.C",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
			canvas->Clear();
			legAS->Clear();
			canvas->SetLogy(0);
		}

		// Parameter Yield plots vs EP Bin
		for(int i = 0; i < nTriggerPtBins; i++) {
			TMultiGraph * mg = new TMultiGraph();
			for(int k = 0; k < nEPBins; k++) {
				ptBinASYieldsGraphsEP[k][i]->SetLineColor(colorList[k]);
				ptBinASYieldsGraphsEP[k][i]->SetMarkerColor(colorList[k]);
				ptBinASYieldsGraphsEP[k][i]->SetMarkerStyle(markerList[k]);
				mg->Add(ptBinASYieldsGraphsEP[k][i]);
				legAS->AddEntry(ptBinASYieldsGraphsEP[k][i],titleList[k], "lp");
			}
			mg->Draw(style);
			mg->SetTitle(Form("%0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			mg->GetXaxis()->SetTitle(ptBinASYieldsGraphsEP[0][i]->GetXaxis()->GetTitle());
			mg->GetYaxis()->SetTitle(ptBinASYieldsGraphsEP[0][i]->GetYaxis()->GetTitle());
	//		ptBinASIntegralsGraphs[0]->GetYaxis()->SetRangeUser(min_y,max_y);
			//  canvas->SetLogy(1);
			legAS->Draw("SAME");
			canvas->Print(Form("EP_AS_Yields_" + jetClass + ".pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
			canvas->Print(Form("EP_AS_Yields_" + jetClass + ".C",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
			canvas->Clear();
			legAS->Clear();
			canvas->SetLogy(0);
		}



		for(int i = 0; i < nTriggerPtBins; i++) {
			TMultiGraph * mg = new TMultiGraph();
			for(int k = 0; k < nEPBins; k++) {
				ptBinASRmsGraphsEP[k][i]->SetLineColor(colorList[k]);
				ptBinASRmsGraphsEP[k][i]->SetMarkerColor(colorList[k]);
				ptBinASRmsGraphsEP[k][i]->SetMarkerStyle(markerList[k]);
				mg->Add(ptBinASRmsGraphsEP[k][i]);
				legAS->AddEntry(ptBinASRmsGraphsEP[k][i],titleList[k], "lp");
			}
			mg->Draw(style);
			mg->SetTitle(Form("%0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			mg->GetXaxis()->SetTitle(ptBinASRmsGraphsEP[0][i]->GetXaxis()->GetTitle());
			mg->GetYaxis()->SetTitle(ptBinASRmsGraphsEP[0][i]->GetYaxis()->GetTitle());
			legAS->Draw("SAME");
			canvas->Print(Form("EP_AS_Rms_JetPt_%0.f_%0.f.pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
			canvas->Print(Form("EP_AS_Rms_JetPt_%0.f_%0.f.C",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
			canvas->Clear();
			legAS->Clear();
			canvas->SetLogy(0);
//			delete mg;
		}


			//NearSide
		for(int i = 0; i < nTriggerPtBins; i++) {
			TMultiGraph * mg = new TMultiGraph();
			for(int k = 0; k < nEPBins; k++) {
				ptBinNSIntegralsGraphsEP[k][i]->SetLineColor(colorList[k]);
				ptBinNSIntegralsGraphsEP[k][i]->SetMarkerColor(colorList[k]);
				ptBinNSIntegralsGraphsEP[k][i]->SetMarkerStyle(markerList[k]);
				mg->Add(ptBinNSIntegralsGraphsEP[k][i]);
				legAS->AddEntry(ptBinNSIntegralsGraphsEP[k][i],titleList[k], "lp");
			}
			canvas->SetLogy(1);
			mg->Draw(style);
			mg->SetTitle(Form("%0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			mg->GetXaxis()->SetTitle(ptBinNSIntegralsGraphsEP[0][i]->GetXaxis()->GetTitle());
			mg->GetYaxis()->SetTitle(ptBinNSIntegralsGraphsEP[0][i]->GetYaxis()->GetTitle());
	//		ptBinASIntegralsGraphs[0]->GetYaxis()->SetRangeUser(min_y,max_y);
			legAS->Draw("SAME");
			canvas->Print(Form("EP_NS_Integrals_JetPt_%0.f_%0.f.pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
			canvas->Print(Form("EP_NS_Integrals_JetPt_%0.f_%0.f.C",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
			canvas->Clear();
			legAS->Clear();
			canvas->SetLogy(0);
//			delete mg;
		}


		for(int i = 0; i < nTriggerPtBins; i++) {
			TMultiGraph * mg = new TMultiGraph();
			for(int k = 0; k < nEPBins; k++) {
				ptBinNSRmsGraphsEP[k][i]->SetLineColor(colorList[k]);
				ptBinNSRmsGraphsEP[k][i]->SetMarkerColor(colorList[k]);
				ptBinNSRmsGraphsEP[k][i]->SetMarkerStyle(markerList[k]);
				mg->Add(ptBinNSRmsGraphsEP[k][i]);
				legAS->AddEntry(ptBinNSRmsGraphsEP[k][i],titleList[k], "lp");
			}
			mg->Draw(style);
			mg->SetTitle(Form("%0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			mg->GetXaxis()->SetTitle(ptBinNSRmsGraphsEP[0][i]->GetXaxis()->GetTitle());
			mg->GetYaxis()->SetTitle(ptBinNSRmsGraphsEP[0][i]->GetYaxis()->GetTitle());
			legAS->Draw("SAME");
			canvas->Print(Form("EP_NS_Rms_JetPt_%0.f_%0.f.pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
			canvas->Print(Form("EP_NS_Rms_JetPt_%0.f_%0.f.C",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
			canvas->Clear();
			legAS->Clear();
			canvas->SetLogy(0);

		}

	for (int k = 0; k < nEPBins; k++) {
		for (int i = 0; i < nTriggerPtBins; i++) {

			fOut->Add(ptBinASIntegralsGraphsEP[k][i]);
			fOut->Add(ptBinNSIntegralsGraphsEP[k][i]);
		}
	}

	// FIXME Calculate systematic uncertainties here
	// Note that this is for pp, when the event plane division is useful for estimating syst. uncert.
	if (SystematicUncertOutFilePath.Length() && nEPBins > 1) {
		printf("Creating Systematic Uncertainties from Event Plane analysis.\n");	

		gSystem->ChangeDirectory(initialDirectory.Data());
		TFile * fSystematicUncertOutFile = TFile::Open(SystematicUncertOutFilePath.Data(),"RECREATE");
		if (!fSystematicUncertOutFile) {
			fprintf(stderr,"Error: Unable to open systematic error output file!\n");
			exit(1);
		}
		gSystem->ChangeDirectory(outputDirPath.Data());

		CalculateSystematicUncertainty(ptBinASYieldsGraphsEP,fSystematicUncertOutFile);
		CalculateSystematicUncertainty(ptBinNSIntegralsGraphsEP,fSystematicUncertOutFile);
		CalculateSystematicUncertainty(ptBinASIntegralsGraphsEP,fSystematicUncertOutFile);

		CalculateSystematicUncertainty(ptBinASWidthsGraphsEP,fSystematicUncertOutFile);
		CalculateSystematicUncertainty(ptBinNSRmsGraphsEP,fSystematicUncertOutFile);
		CalculateSystematicUncertainty(ptBinASRmsGraphsEP,fSystematicUncertOutFile);


		fSystematicUncertOutFile->Write();
	}


	// Removing first point from TGraphs
	int nSkipPoints = 0; // How many points are we dropping.  2 for EP
	for (int l = 0; l < nSkipPoints; l++) { 
		for (int k = 0; k < nEPBins; k++) {
			for (int i = 0; i < nTriggerPtBins; i++) {
				ptBinASIntegralsGraphsEP[k][i]->RemovePoint(0);
				ptBinNSIntegralsGraphsEP[k][i]->RemovePoint(0);
				ptBinASYieldsGraphsEP[k][i]->RemovePoint(0);
				ptBinASRmsGraphsEP[k][i]->RemovePoint(0);
				ptBinNSRmsGraphsEP[k][i]->RemovePoint(0);
			}
		}
	}
		// Making out/in mid/in ratio plots
	//	vector< * TGraphErrors> OutOverInJetPtBin;
	//	vector< * TGraphErrors> MidOverInJetPtBin;

	TF1 * constOne = new TF1("constOne","1",0,15);
	constOne->SetLineColor(kGray);
	constOne->SetLineStyle(2);
	constOne->SetLineWidth(3);

	TF1 * constZero = new TF1("constZero","0",0,15);
	constZero->SetLineColor(kGray);
	constZero->SetLineStyle(2);
	constZero->SetLineWidth(3);
		
	float RatioMin = 0.5;
	float RatioMax = 1.5;

	float DiffMin = -0.8;
	float DiffMax = 1.45;

		//Away Side
		// Yields
		for(int i = 0; i < nTriggerPtBins; i++) {
			printf("Doing something for trigger pt bin %d\n",i);
			combinedNames = Form("OutOverIn_AS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			TGraphErrors * OutOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			OutOverIn->SetName(combinedNames);
			OutOverIn->SetTitle(Form("Away-Side Ratio, %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			combinedNames = Form("MidOverIn_AS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			TGraphErrors * MidOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			MidOverIn->SetName(combinedNames);
			MidOverIn->SetTitle(Form("Away-Side Ratio, %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));

			combinedNames = Form("OutMinusIn_AS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			TGraphErrors * OutMinusIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			OutMinusIn->SetName(combinedNames);
			OutMinusIn->SetTitle(Form("Away-Side D_{RP}(Out), %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			combinedNames = Form("MidMinusIn_AS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			TGraphErrors * MidMinusIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			MidMinusIn->SetName(combinedNames);
			MidMinusIn->SetTitle(Form("Away-Side D_{RP}(Mid), %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));

			double * InX = ptBinASIntegralsGraphsEP[1][i]->GetX();
			double * InY = ptBinASIntegralsGraphsEP[1][i]->GetY();
			double * InYErr = ptBinASIntegralsGraphsEP[1][i]->GetEY();

			double * MidY = ptBinASIntegralsGraphsEP[2][i]->GetY();
			double * MidYErr = ptBinASIntegralsGraphsEP[2][i]->GetEY();

			double * OutY = ptBinASIntegralsGraphsEP[3][i]->GetY();
			double * OutYErr = ptBinASIntegralsGraphsEP[3][i]->GetEY();
	
			for (int j = 0; j < nAssocParticleBins-nSkipPoints; j++) {
			//	double OutOverInPoint = ptBinASIntegralsGraphsEP[3][i]->Get
				if (InY != 0) {
					double OutOverInPoint = OutY[j] / InY[j];
					double OutOverInErr = OutOverInPoint * TMath::Sqrt(TMath::Power(OutYErr[j]/OutY[j],2) + TMath::Power(InYErr[j]/InY[j],2));
				
					OutOverIn->SetPoint(j,InX[j],OutOverInPoint);
					OutOverIn->SetPointError(j,0,OutOverInErr);

					double MidOverInPoint = MidY[j] / InY[j];
					double MidOverInErr = MidOverInPoint * TMath::Sqrt(TMath::Power(MidYErr[j]/MidY[j],2) + TMath::Power(InYErr[j]/InY[j],2));

					MidOverIn->SetPoint(j,InX[j],MidOverInPoint);
					MidOverIn->SetPointError(j,0,MidOverInErr);


					double OutMinusInPoint = OutY[j] - InY[j];
					double OutMinusInErr = TMath::Sqrt(TMath::Power(OutYErr[j],2) + TMath::Power(InYErr[j],2));
				
					OutMinusIn->SetPoint(j,InX[j],OutMinusInPoint);
					OutMinusIn->SetPointError(j,0,OutMinusInErr);

					double MidMinusInPoint = MidY[j] - InY[j];
					double MidMinusInErr = TMath::Sqrt(TMath::Power(MidYErr[j],2) + TMath::Power(InYErr[j],2));

					MidMinusIn->SetPoint(j,InX[j],MidMinusInPoint);
					MidMinusIn->SetPointError(j,0,MidMinusInErr);

				}
			}

//			canvas->Clear();
			OutOverIn->SetMarkerStyle(kFullSquare);
			OutOverIn->SetMarkerColor(kBlack);
			OutOverIn->SetFillColorAlpha(kBlack,0.33);
			OutOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			OutOverIn->GetYaxis()->SetTitle("Y_{Out}/Y_{In}");

			// FIXME Draw these after producing the whole set and applying the systematic uncertainties
/*
			OutOverIn->Draw("AP 3");
			constOne->Draw("SAME");
			OutOverIn->Draw("P 3 SAME");
			OutOverIn->GetYaxis()->SetRangeUser(RatioMin,RatioMax);
			canvas->Print(Form("EP_AS_OutOverIn_"+jetClass+".pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
//			canvas->Print(Form("EP_AS_OutOverIn_JetPt_%0.f_%0.f.eps",jetPtBins.at(i),jetPtBins.at(i+1)));
			canvas->Print(Form("EP_AS_OutOverIn_"+jetClass+".C",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
	*/		
//			canvas->Clear();
			MidOverIn->SetMarkerStyle(kFullSquare);
			MidOverIn->SetMarkerColor(kBlack);
			MidOverIn->SetFillColorAlpha(kBlack,0.33);
			MidOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			MidOverIn->GetYaxis()->SetTitle("Y_{Mid}/Y_{In}");

/*			MidOverIn->Draw("AP 3");
			constOne->Draw("SAME");
			MidOverIn->Draw("P 3 SAME");
			MidOverIn->GetYaxis()->SetRangeUser(RatioMin,RatioMax);
			canvas->Print(Form("EP_AS_MidOverIn_"+jetClass+".pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
			canvas->Print(Form("EP_AS_MidOverIn_"+jetClass+".C",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
*/
			OutOverIn_AS_JetPtBin.push_back(OutOverIn);
			MidOverIn_AS_JetPtBin.push_back(MidOverIn);


//			canvas->Clear();
			OutMinusIn->SetMarkerStyle(kFullSquare);
			OutMinusIn->SetMarkerColor(kBlack);
			OutMinusIn->SetFillColorAlpha(kBlack,0.33);
			OutMinusIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			OutMinusIn->GetYaxis()->SetTitle("Y_{Out} - Y_{In}");
/*
			OutMinusIn->Draw("AP 3");
			constZero->Draw("SAME");
			OutMinusIn->Draw("P 3 SAME");
			OutMinusIn->GetYaxis()->SetRangeUser(DiffMin,DiffMax);
			canvas->Print(Form("EP_AS_OutMinusIn_"+jetClass+".pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
//			canvas->Print(Form("EP_AS_OutMinusIn_JetPt_%0.f_%0.f.eps",jetPtBins.at(i),jetPtBins.at(i+1)));
			canvas->Print(Form("EP_AS_OutMinusIn_"+jetClass+".C",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
*/
			
			canvas->Clear();
			MidMinusIn->SetMarkerStyle(kFullSquare);
			MidMinusIn->SetMarkerColor(kBlack);
			MidMinusIn->SetFillColorAlpha(kBlack,0.33);
			MidMinusIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			MidMinusIn->GetYaxis()->SetTitle("Y_{Mid} - Y_{In}");
/*
			MidMinusIn->Draw("AP 3");
			constZero->Draw("SAME");
			MidMinusIn->Draw("P 3 SAME");
			MidMinusIn->GetYaxis()->SetRangeUser(DiffMin,DiffMax);
			canvas->Print(Form("EP_AS_MidMinusIn_"+jetClass+".pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
//			canvas->Print(Form("EP_AS_MidMinusIn_JetPt_%0.f_%0.f.eps",jetPtBins.at(i),jetPtBins.at(i+1)));
			canvas->Print(Form("EP_AS_MidMinusIn_"+jetClass+".C",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
*/

			OutMinusIn_AS_JetPtBin.push_back(OutMinusIn);
			MidMinusIn_AS_JetPtBin.push_back(MidMinusIn);


		}

		// AwaySide
		// Parameter Yield
		for(int i = 0; i < nTriggerPtBins; i++) {
			TGraphErrors * OutOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			combinedNames = Form("OutOverInPar_AS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			OutOverIn->SetName(combinedNames);
			OutOverIn->SetTitle("Out/In Yield Ratio (Away-Side)");
			TGraphErrors * MidOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			combinedNames = Form("MidOverInPar_AS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			MidOverIn->SetName(combinedNames);
			MidOverIn->SetTitle("Mid/In Yield Ratio (Away-Side)");
			
			double * InX = ptBinASYieldsGraphsEP[1][i]->GetX();
			double * InY = ptBinASYieldsGraphsEP[1][i]->GetY();
			double * InYErr = ptBinASYieldsGraphsEP[1][i]->GetEY();

			double * MidY = ptBinASYieldsGraphsEP[2][i]->GetY();
			double * MidYErr = ptBinASYieldsGraphsEP[2][i]->GetEY();

			double * OutY = ptBinASYieldsGraphsEP[3][i]->GetY();
			double * OutYErr = ptBinASYieldsGraphsEP[3][i]->GetEY();
	
			for (int j = 0; j < nAssocParticleBins-nSkipPoints; j++) {
			//	double OutOverInPoint = ptBinASIntegralsGraphsEP[3][i]->Get
				if (InY != 0) {
					double OutOverInPoint = OutY[j] / InY[j];
					double OutOverInErr = OutOverInPoint * TMath::Sqrt(TMath::Power(OutYErr[j]/OutY[j],2) + TMath::Power(InYErr[j]/InY[j],2));
				
					OutOverIn->SetPoint(j,InX[j],OutOverInPoint);
					OutOverIn->SetPointError(j,0,OutOverInErr);

					double MidOverInPoint = MidY[j] / InY[j];
					double MidOverInErr = MidOverInPoint * TMath::Sqrt(TMath::Power(MidYErr[j]/MidY[j],2) + TMath::Power(InYErr[j]/InY[j],2));

					MidOverIn->SetPoint(j,InX[j],MidOverInPoint);
					MidOverIn->SetPointError(j,0,MidOverInErr);


				}
			}

			canvas->Clear();
			OutOverIn->SetMarkerStyle(kFullSquare);
			OutOverIn->SetMarkerColor(kBlack);
			OutOverIn->SetFillColorAlpha(kBlack,0.33);
			OutOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			OutOverIn->GetYaxis()->SetTitle("Y_{Out}/Y_{In}");
/*			OutOverIn->Draw("AP 3");
			constOne->Draw("SAME");
			OutOverIn->Draw("P 3 SAME");
			OutOverIn->GetYaxis()->SetRangeUser(RatioMin,RatioMax);
			canvas->Print(Form("EP_AS_OutOverInPar_"+jetClass+".pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
			canvas->Print(Form("EP_AS_OutOverInPar_"+jetClass+".C",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
	*/		
			canvas->Clear();
			MidOverIn->SetMarkerStyle(kFullSquare);
			MidOverIn->SetMarkerColor(kBlack);
			MidOverIn->SetFillColorAlpha(kBlack,0.33);
			MidOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			MidOverIn->GetYaxis()->SetTitle("Y_{Mid}/Y_{In}");
/*			MidOverIn->Draw("AP 3");
			constOne->Draw("SAME");
			MidOverIn->Draw("P 3 SAME");
			MidOverIn->GetYaxis()->SetRangeUser(RatioMin,RatioMax);
			canvas->Print(Form("EP_AS_MidOverInPar_"+jetClass+".pdf",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			canvas->Print(Form("EP_AS_MidOverInPar_"+jetClass+".C",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
*/
			OutOverInPar_AS_JetPtBin.push_back(OutOverIn);
			MidOverInPar_AS_JetPtBin.push_back(MidOverIn);
		}

		//Away Side
		// RMS
		for(int i = 0; i < nTriggerPtBins; i++) {
			combinedNames = Form("RMSOutOverIn_AS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			TGraphErrors * RMSOutOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			RMSOutOverIn->SetName(combinedNames);
	//		OutOverIn->SetTitle("Out/In Yield Ratio (Away-Side)");
			RMSOutOverIn->SetTitle(Form("Away-Side RMS Ratio, %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			combinedNames = Form("RMSMidOverIn_AS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			TGraphErrors * RMSMidOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			RMSMidOverIn->SetName(combinedNames);
		//	MidOverIn->SetTitle("Out/In Yield Ratio (Away-Side)");
			RMSMidOverIn->SetTitle(Form("Away-Side RMS Ratio, %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			
			double * InX = ptBinASRmsGraphsEP[1][i]->GetX();
			double * InY = ptBinASRmsGraphsEP[1][i]->GetY();
			double * InYErr = ptBinASRmsGraphsEP[1][i]->GetEY();

			double * MidY = ptBinASRmsGraphsEP[2][i]->GetY();
			double * MidYErr = ptBinASRmsGraphsEP[2][i]->GetEY();

			double * OutY = ptBinASRmsGraphsEP[3][i]->GetY();
			double * OutYErr = ptBinASRmsGraphsEP[3][i]->GetEY();
	
			for (int j = 0; j < nAssocParticleBins-nSkipPoints; j++) {
			//	double OutOverInPoint = ptBinASIntegralsGraphsEP[3][i]->Get
				if (InY != 0) {
					double OutOverInPoint = OutY[j] / InY[j];
					double OutOverInErr = OutOverInPoint * TMath::Sqrt(TMath::Power(OutYErr[j]/OutY[j],2) + TMath::Power(InYErr[j]/InY[j],2));
				
					RMSOutOverIn->SetPoint(j,InX[j],OutOverInPoint);
					RMSOutOverIn->SetPointError(j,0,OutOverInErr);

					double MidOverInPoint = MidY[j] / InY[j];
					double MidOverInErr = MidOverInPoint * TMath::Sqrt(TMath::Power(MidYErr[j]/MidY[j],2) + TMath::Power(InYErr[j]/InY[j],2));

					RMSMidOverIn->SetPoint(j,InX[j],MidOverInPoint);
					RMSMidOverIn->SetPointError(j,0,MidOverInErr);
				}
			}

		//	canvas->Clear();
			RMSOutOverIn->SetMarkerStyle(kFullSquare);
			RMSOutOverIn->SetMarkerColor(kBlack);
			RMSOutOverIn->SetFillColorAlpha(kBlack,0.33);
			RMSOutOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			RMSOutOverIn->GetYaxis()->SetTitle("RMS_{Out}/RMS_{In}");

	//		canvas->Clear();
			RMSMidOverIn->SetMarkerStyle(kFullSquare);
			RMSMidOverIn->SetMarkerColor(kBlack);
			RMSMidOverIn->SetFillColorAlpha(kBlack,0.33);
			RMSMidOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			RMSMidOverIn->GetYaxis()->SetTitle("RMS_{Mid}/RMS_{In}");

			RMSOutOverIn_AS_JetPtBin.push_back(RMSOutOverIn);
			RMSMidOverIn_AS_JetPtBin.push_back(RMSMidOverIn);
		}

		//Away Side
		// Widths (Fit Parameter)
		for(int i = 0; i < nTriggerPtBins; i++) {
			combinedNames = Form("WidthsOutOverIn_AS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			TGraphErrors * WidthsOutOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			WidthsOutOverIn->SetName(combinedNames);
	//		OutOverIn->SetTitle("Out/In Yield Ratio (Away-Side)");
			WidthsOutOverIn->SetTitle(Form("Away-Side Widths Ratio, %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			combinedNames = Form("WidthsMidOverIn_AS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			TGraphErrors * WidthsMidOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			WidthsMidOverIn->SetName(combinedNames);
		//	MidOverIn->SetTitle("Out/In Yield Ratio (Away-Side)");
			WidthsMidOverIn->SetTitle(Form("Away-Side Widths Ratio, %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			
			double * InX = ptBinASWidthsGraphsEP[1][i]->GetX();
			double * InY = ptBinASWidthsGraphsEP[1][i]->GetY();
			double * InYErr = ptBinASWidthsGraphsEP[1][i]->GetEY();

			double * MidY = ptBinASWidthsGraphsEP[2][i]->GetY();
			double * MidYErr = ptBinASWidthsGraphsEP[2][i]->GetEY();

			double * OutY = ptBinASWidthsGraphsEP[3][i]->GetY();
			double * OutYErr = ptBinASWidthsGraphsEP[3][i]->GetEY();
	
			for (int j = 0; j < nAssocParticleBins-nSkipPoints; j++) {
			//	double OutOverInPoint = ptBinASIntegralsGraphsEP[3][i]->Get
				if (InY != 0) {
					double OutOverInPoint = OutY[j] / InY[j];
					double OutOverInErr = OutOverInPoint * TMath::Sqrt(TMath::Power(OutYErr[j]/OutY[j],2) + TMath::Power(InYErr[j]/InY[j],2));
				
					WidthsOutOverIn->SetPoint(j,InX[j],OutOverInPoint);
					WidthsOutOverIn->SetPointError(j,0,OutOverInErr);

					double MidOverInPoint = MidY[j] / InY[j];
					double MidOverInErr = MidOverInPoint * TMath::Sqrt(TMath::Power(MidYErr[j]/MidY[j],2) + TMath::Power(InYErr[j]/InY[j],2));

					WidthsMidOverIn->SetPoint(j,InX[j],MidOverInPoint);
					WidthsMidOverIn->SetPointError(j,0,MidOverInErr);
				}
			}
		//	canvas->Clear();
			WidthsOutOverIn->SetMarkerStyle(kFullSquare);
			WidthsOutOverIn->SetMarkerColor(kBlack);
			WidthsOutOverIn->SetFillColorAlpha(kBlack,0.33);
			WidthsOutOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			WidthsOutOverIn->GetYaxis()->SetTitle("#sigma_{Out}/#sigma_{In}");

		//	canvas->Clear();
			WidthsMidOverIn->SetMarkerStyle(kFullSquare);
			WidthsMidOverIn->SetMarkerColor(kBlack);
			WidthsMidOverIn->SetFillColorAlpha(kBlack,0.33);
			WidthsMidOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			WidthsMidOverIn->GetYaxis()->SetTitle("#sigma_{Mid}/#sigma_{In}");

			WidthsOutOverIn_AS_JetPtBin.push_back(WidthsOutOverIn);
			WidthsMidOverIn_AS_JetPtBin.push_back(WidthsMidOverIn);
		}



		//Near Side
		for(int i = 0; i < nTriggerPtBins; i++) {
			TGraphErrors * OutOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			combinedNames = Form("OutOverIn_NS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			OutOverIn->SetName(combinedNames);
			OutOverIn->SetTitle(Form("Near-Side Ratio, %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			TGraphErrors * MidOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			combinedNames = Form("MidOverIn_NS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			MidOverIn->SetName(combinedNames);
			MidOverIn->SetTitle(Form("Near-Side Ratio, %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			
			combinedNames = Form("OutMinusIn_NS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			TGraphErrors * OutMinusIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			OutMinusIn->SetName(combinedNames);
			OutMinusIn->SetTitle(Form("Near-Side D_{RP}(Out), %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			combinedNames = Form("MidMinusIn_NS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			TGraphErrors * MidMinusIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			MidMinusIn->SetName(combinedNames);
			MidMinusIn->SetTitle(Form("Near-Side D_{RP}(Mid), %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));



			double * InX = ptBinNSIntegralsGraphsEP[1][i]->GetX();
			double * InY = ptBinNSIntegralsGraphsEP[1][i]->GetY();
			double * InYErr = ptBinNSIntegralsGraphsEP[1][i]->GetEY();

			double * MidY = ptBinNSIntegralsGraphsEP[2][i]->GetY();
			double * MidYErr = ptBinNSIntegralsGraphsEP[2][i]->GetEY();

			double * OutY = ptBinNSIntegralsGraphsEP[3][i]->GetY();
			double * OutYErr = ptBinNSIntegralsGraphsEP[3][i]->GetEY();
	
			for (int j = 0; j < nAssocParticleBins-nSkipPoints; j++) {
			//	double OutOverInPoint = ptBinASIntegralsGraphsEP[3][i]->Get
				if (InY != 0) {
					double OutOverInPoint = OutY[j] / InY[j];
					double OutOverInErr = OutOverInPoint * TMath::Sqrt(TMath::Power(OutYErr[j]/OutY[j],2) + TMath::Power(InYErr[j]/InY[j],2));
				
					OutOverIn->SetPoint(j,InX[j],OutOverInPoint);
					OutOverIn->SetPointError(j,0,OutOverInErr);

					double MidOverInPoint = MidY[j] / InY[j];
					double MidOverInErr = MidOverInPoint * TMath::Sqrt(TMath::Power(MidYErr[j]/MidY[j],2) + TMath::Power(InYErr[j]/InY[j],2));

					MidOverIn->SetPoint(j,InX[j],MidOverInPoint);
					MidOverIn->SetPointError(j,0,MidOverInErr);

					double OutMinusInPoint = OutY[j] - InY[j];
					double OutMinusInErr = TMath::Sqrt(TMath::Power(OutYErr[j],2) + TMath::Power(InYErr[j],2));
				
					OutMinusIn->SetPoint(j,InX[j],OutMinusInPoint);
					OutMinusIn->SetPointError(j,0,OutMinusInErr);

					double MidMinusInPoint = MidY[j] - InY[j];
					double MidMinusInErr = TMath::Sqrt(TMath::Power(MidYErr[j],2) + TMath::Power(InYErr[j],2));

					MidMinusIn->SetPoint(j,InX[j],MidMinusInPoint);
					MidMinusIn->SetPointError(j,0,MidMinusInErr);

				}
			}

//			canvas->Clear();
			OutOverIn->SetMarkerStyle(kFullSquare);
			OutOverIn->SetMarkerColor(kBlack);
			OutOverIn->SetFillColorAlpha(kBlack,0.33);
			OutOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			OutOverIn->GetYaxis()->SetTitle("Y_{Out}/Y_{In}");

//			canvas->Clear();
			MidOverIn->SetMarkerStyle(kFullSquare);
			MidOverIn->SetMarkerColor(kBlack);
			MidOverIn->SetFillColorAlpha(kBlack,0.33);
			MidOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			MidOverIn->GetYaxis()->SetTitle("Y_{Mid}/Y_{In}");

			OutOverIn_NS_JetPtBin.push_back(OutOverIn);
			MidOverIn_NS_JetPtBin.push_back(MidOverIn);


			canvas->Clear();
			OutMinusIn->SetMarkerStyle(kFullSquare);
			OutMinusIn->SetMarkerColor(kBlack);
			OutMinusIn->SetFillColorAlpha(kBlack,0.33);
			OutMinusIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			OutMinusIn->GetYaxis()->SetTitle("Y_{Out} - Y_{In}");
/*			OutMinusIn->Draw("AP 3");
			constZero->Draw("SAME");
			OutMinusIn->Draw("P 3 SAME");
			OutMinusIn->GetYaxis()->SetRangeUser(DiffMin,DiffMax);
			canvas->Print(Form("EP_NS_OutMinusIn_"+jetClass+".pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
//			canvas->Print(Form("EP_NS_OutMinusIn_JetPt_%0.f_%0.f.eps",jetPtBins.at(i),jetPtBins.at(i+1)));
			canvas->Print(Form("EP_NS_OutMinusIn_"+jetClass+".C",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
	*/		
			canvas->Clear();
			MidMinusIn->SetMarkerStyle(kFullSquare);
			MidMinusIn->SetMarkerColor(kBlack);
			MidMinusIn->SetFillColorAlpha(kBlack,0.33);
			MidMinusIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			MidMinusIn->GetYaxis()->SetTitle("Y_{Mid} - Y_{In}");
/*			MidMinusIn->Draw("AP 3");
			constZero->Draw("SAME");
			MidMinusIn->Draw("P 3 SAME");
			MidMinusIn->GetYaxis()->SetRangeUser(DiffMin,DiffMax);
			canvas->Print(Form("EP_NS_MidMinusIn_"+jetClass+".pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
	//		canvas->Print(Form("EP_NS_MidMinusIn_JetPt_%0.f_%0.f.eps",jetPtBins.at(i),jetPtBins.at(i+1)));
			canvas->Print(Form("EP_NS_MidMinusIn_"+jetClass+".C",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
*/

			OutMinusIn_NS_JetPtBin.push_back(OutMinusIn);
			MidMinusIn_NS_JetPtBin.push_back(MidMinusIn);




		}


		// Near Side
		// RMS
		for(int i = 0; i < nTriggerPtBins; i++) {
			combinedNames = Form("RMSOutOverIn_NS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			TGraphErrors * RMSOutOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			RMSOutOverIn->SetName(combinedNames);
			RMSOutOverIn->SetTitle(Form("Near-Side RMS Ratio, %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			combinedNames = Form("RMSMidOverIn_NS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
			TGraphErrors * RMSMidOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
			RMSMidOverIn->SetName(combinedNames);
			RMSMidOverIn->SetTitle(Form("Near-Side RMS Ratio, %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
			
			double * InX = ptBinNSRmsGraphsEP[1][i]->GetX();
			double * InY = ptBinNSRmsGraphsEP[1][i]->GetY();
			double * InYErr = ptBinNSRmsGraphsEP[1][i]->GetEY();

			double * MidY = ptBinNSRmsGraphsEP[2][i]->GetY();
			double * MidYErr = ptBinNSRmsGraphsEP[2][i]->GetEY();

			double * OutY = ptBinNSRmsGraphsEP[3][i]->GetY();
			double * OutYErr = ptBinNSRmsGraphsEP[3][i]->GetEY();
			for (int j = 0; j < nAssocParticleBins-nSkipPoints; j++) {
				if (InY != 0) {
					double OutOverInPoint = OutY[j] / InY[j];
					double OutOverInErr = OutOverInPoint * TMath::Sqrt(TMath::Power(OutYErr[j]/OutY[j],2) + TMath::Power(InYErr[j]/InY[j],2));
				
					RMSOutOverIn->SetPoint(j,InX[j],OutOverInPoint);
					RMSOutOverIn->SetPointError(j,0,OutOverInErr);

					double MidOverInPoint = MidY[j] / InY[j];
					double MidOverInErr = MidOverInPoint * TMath::Sqrt(TMath::Power(MidYErr[j]/MidY[j],2) + TMath::Power(InYErr[j]/InY[j],2));

					RMSMidOverIn->SetPoint(j,InX[j],MidOverInPoint);
					RMSMidOverIn->SetPointError(j,0,MidOverInErr);
				}
			}

			RMSOutOverIn->SetMarkerStyle(kFullSquare);
			RMSOutOverIn->SetMarkerColor(kBlack);
			RMSOutOverIn->SetFillColorAlpha(kBlack,0.33);
			RMSOutOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			RMSOutOverIn->GetYaxis()->SetTitle("RMS_{Out}/RMS_{In}");

			RMSMidOverIn->SetMarkerStyle(kFullSquare);
			RMSMidOverIn->SetMarkerColor(kBlack);
			RMSMidOverIn->SetFillColorAlpha(kBlack,0.33);
			RMSMidOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
			RMSMidOverIn->GetYaxis()->SetTitle("RMS_{Mid}/RMS_{In}");

			RMSOutOverIn_NS_JetPtBin.push_back(RMSOutOverIn);
			RMSMidOverIn_NS_JetPtBin.push_back(RMSMidOverIn);
		}

	// Near Side
	// Widths (Fit Parameter)
/* To be implemented when NS Width is defined
	for(int i = 0; i < nTriggerPtBins; i++) {
		combinedNames = Form("WidthsOutOverIn_NS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
		TGraphErrors * WidthsOutOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
		WidthsOutOverIn->SetName(combinedNames);
		WidthsOutOverIn->SetTitle(Form("Near-Side Widths Ratio, %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
		combinedNames = Form("WidthsMidOverIn_NS_"+jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
		TGraphErrors * WidthsMidOverIn = new TGraphErrors(nAssocParticleBins-nSkipPoints);
		WidthsMidOverIn->SetName(combinedNames);
		WidthsMidOverIn->SetTitle(Form("Near-Side Widths Ratio, %0.f #leq p_{T}^{%s} #leq %0.f",fTriggerPtBins.at(i),sTriggerTitle.Data(),fTriggerPtBins.at(i+1)));
		
		double * InX = ptBinNSWidthsGraphsEP[1][i]->GetX();
		double * InY = ptBinNSWidthsGraphsEP[1][i]->GetY();
		double * InYErr = ptBinNSWidthsGraphsEP[1][i]->GetEY();

		double * MidY = ptBinNSWidthsGraphsEP[2][i]->GetY();
		double * MidYErr = ptBinNSWidthsGraphsEP[2][i]->GetEY();

		double * OutY = ptBinNSWidthsGraphsEP[3][i]->GetY();
		double * OutYErr = ptBinNSWidthsGraphsEP[3][i]->GetEY();
		for (int j = 0; j < nAssocParticleBins-nSkipPoints; j++) {
			if (InY != 0) {
				double OutOverInPoint = OutY[j] / InY[j];
				double OutOverInErr = OutOverInPoint * TMath::Sqrt(TMath::Power(OutYErr[j]/OutY[j],2) + TMath::Power(InYErr[j]/InY[j],2));
			
				WidthsOutOverIn->SetPoint(j,InX[j],OutOverInPoint);
				WidthsOutOverIn->SetPointError(j,0,OutOverInErr);

				double MidOverInPoint = MidY[j] / InY[j];
				double MidOverInErr = MidOverInPoint * TMath::Sqrt(TMath::Power(MidYErr[j]/MidY[j],2) + TMath::Power(InYErr[j]/InY[j],2));

				WidthsMidOverIn->SetPoint(j,InX[j],MidOverInPoint);
				WidthsMidOverIn->SetPointError(j,0,MidOverInErr);
			}
		}

		WidthsOutOverIn->SetMarkerStyle(kFullSquare);
		WidthsOutOverIn->SetMarkerColor(kBlack);
		WidthsOutOverIn->SetFillColorAlpha(kBlack,0.33);
		WidthsOutOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
		WidthsOutOverIn->GetYaxis()->SetTitle("#sigma_{Out}/#sigma_{In}");

		WidthsMidOverIn->SetMarkerStyle(kFullSquare);
		WidthsMidOverIn->SetMarkerColor(kBlack);
		WidthsMidOverIn->SetFillColorAlpha(kBlack,0.33);
		WidthsMidOverIn->GetXaxis()->SetTitle(sAssocPtTitle.Data());
		WidthsMidOverIn->GetYaxis()->SetTitle("#sigma_{Mid}/#sigma_{In}");

		WidthsOutOverIn_NS_JetPtBin.push_back(WidthsOutOverIn);
		WidthsMidOverIn_NS_JetPtBin.push_back(WidthsMidOverIn);
	}
  */



	// Apply Systematic Uncertainties
	if (!ApplySystUncertBeforeRatio && SystematicUncertInFilePath.Length() && nEPBins > 1) {
		printf("Loading Systematic Uncertainties from Event Plane analysis.\n");	

		gSystem->ChangeDirectory(initialDirectory.Data());
		TFile * fSystematicUncertInFile = TFile::Open(SystematicUncertInFilePath.Data(),"READ");
		if (!fSystematicUncertInFile) {
			fprintf(stderr,"Error: Unable to open systematic error input file!\n");
			exit(1);
		}
		gSystem->ChangeDirectory(outputDirPath.Data());

		ApplySystematicUncertainty(OutOverIn_AS_JetPtBin,fSystematicUncertInFile);
		ApplySystematicUncertainty(MidOverIn_AS_JetPtBin,fSystematicUncertInFile);
		ApplySystematicUncertainty(OutOverInPar_AS_JetPtBin,fSystematicUncertInFile);
		ApplySystematicUncertainty(MidOverInPar_AS_JetPtBin,fSystematicUncertInFile);
		ApplySystematicUncertainty(OutOverIn_NS_JetPtBin,fSystematicUncertInFile);
		ApplySystematicUncertainty(MidOverIn_NS_JetPtBin,fSystematicUncertInFile);

		ApplySystematicUncertainty(RMSOutOverIn_AS_JetPtBin,fSystematicUncertInFile);
		ApplySystematicUncertainty(RMSMidOverIn_AS_JetPtBin,fSystematicUncertInFile);
		ApplySystematicUncertainty(RMSOutOverIn_NS_JetPtBin,fSystematicUncertInFile);
		ApplySystematicUncertainty(RMSMidOverIn_NS_JetPtBin,fSystematicUncertInFile);

		ApplySystematicUncertainty(OutMinusIn_AS_JetPtBin,fSystematicUncertInFile);
		ApplySystematicUncertainty(MidMinusIn_AS_JetPtBin,fSystematicUncertInFile);
		ApplySystematicUncertainty(OutMinusIn_NS_JetPtBin,fSystematicUncertInFile);
		ApplySystematicUncertainty(MidMinusIn_NS_JetPtBin,fSystematicUncertInFile);

		ApplySystematicUncertainty(WidthsOutOverIn_AS_JetPtBin,fSystematicUncertInFile);
		ApplySystematicUncertainty(WidthsMidOverIn_AS_JetPtBin,fSystematicUncertInFile);
	}

	// Make Ratio Plots

	// Save time: build vector of all ratios, vector of all differences

	vector<vector<TGraphErrors *>> AllRatioGraphs = {OutOverIn_AS_JetPtBin,MidOverIn_AS_JetPtBin,OutOverInPar_AS_JetPtBin,MidOverInPar_AS_JetPtBin,OutOverIn_NS_JetPtBin,MidOverIn_NS_JetPtBin,RMSOutOverIn_AS_JetPtBin,RMSMidOverIn_AS_JetPtBin,RMSOutOverIn_NS_JetPtBin,RMSMidOverIn_NS_JetPtBin,WidthsOutOverIn_AS_JetPtBin,WidthsMidOverIn_AS_JetPtBin};
//,WidthsOutOverIn_NS_JetPtBin,WidthsMidOverIn_NS_JetPtBin};

	vector<vector<TGraphErrors *>> AllDifferenceGraphs = {OutMinusIn_AS_JetPtBin,MidMinusIn_AS_JetPtBin,OutMinusIn_NS_JetPtBin,MidMinusIn_NS_JetPtBin};

//	for(int i = 0; i < nTriggerPtBins; i++) {
	// Draw Ratios
	for (unsigned int i = 0; i < AllRatioGraphs.size(); i++) {
		for (unsigned int j = 0; j < AllRatioGraphs[i].size(); j++) {
			TGraphErrors * lGraph = AllRatioGraphs[i][j];
			canvas->Clear();
			lGraph->Draw("AP 3");
			constOne->Draw("SAME");
			lGraph->Draw("P 3 SAME");
			lGraph->GetYaxis()->SetRangeUser(RatioMin,RatioMax);
			canvas->Print(Form("EP_%s.pdf",lGraph->GetName()));
			canvas->Print(Form("EP_%s.C",lGraph->GetName()));
		}
	}
	// Draw Differences
	for (unsigned int i = 0; i < AllDifferenceGraphs.size(); i++) {
		for (unsigned int j = 0; j < AllDifferenceGraphs[i].size(); j++) {
			TGraphErrors * lGraph = AllDifferenceGraphs[i][j];
			canvas->Clear();
			lGraph->Draw("AP 3");
			constZero->Draw("SAME");
			lGraph->Draw("P 3 SAME");
			lGraph->GetYaxis()->SetRangeUser(DiffMin,DiffMax);
			canvas->Print(Form("EP_%s.pdf",lGraph->GetName()));
			canvas->Print(Form("EP_%s.C",lGraph->GetName()));
		}
	}

  	/*
		// Away Side
		lGraph = OutOverIn_AS_JetPtBin[i];
		canvas->Clear();
		lGraph->Draw("AP 3");
		constOne->Draw("SAME");
		lGraph->Draw("P 3 SAME");
		lGraph->GetYaxis()->SetRangeUser(RatioMin,RatioMax);a
		canvas->Print(Form("EP_%s.pdf",lGraph->GetName()));
		canvas->Print(Form("EP_%s.C",lGraph->GetName()));

		lGraph = MidOverIn_AS_JetPtBin[i];
		canvas->Clear();
		lGraph->Draw("AP 3");
		constOne->Draw("SAME");
		lGraph->Draw("P 3 SAME");
		lGraph->GetYaxis()->SetRangeUser(RatioMin,RatioMax);
		canvas->Print(Form("EP_%s.pdf",lGraph->GetName()));
		canvas->Print(Form("EP_%s.C",lGraph->GetName()));

		lGraph = OutMinusIn_AS_JetPtBin[i];
		canvas->Clear();
		lGraph->Draw("AP 3");
		constOne->Draw("SAME");
		lGraph->Draw("P 3 SAME");
		lGraph->GetYaxis()->SetRangeUser(RatioMin,RatioMax);
		canvas->Print(Form("EP_%s.pdf",lGraph->GetName()));
		canvas->Print(Form("EP_%s.C",lGraph->GetName()));

		// Near Side
	}
*/ 


	// Draw Plots of the AS Peak with each EP overlayed
	int n_y_SubDiv = TMath::FloorNint(TMath::Sqrt(nAssocParticleBins));
//	int n_x_SubDiv = TMath::CeilNint(TMath::Sqrt(nAssocParticleBins));
	int n_x_SubDiv = TMath::CeilNint(nAssocParticleBins*1.0/n_y_SubDiv);
	TLegend * legASPlotCmp = new TLegend(0.67,0.65,0.9,0.85);
	for (unsigned int i = 0; i < nTriggerPtBins; i++)
	{
	//	canvas->Clear();
//		canvas->Divide(n_x_SubDiv,n_y_SubDiv);
		index = 0;
		for (unsigned int j = 0; j < nAssocParticleBins; j++)
		{
			// Set appropraite part of canvas
	//		index++;
	//		canvas->cd(index);
			canvas->Clear();

			for (int k = 0; k < nEPBins; k++) {
				jetHProjection[k][i][j]->GetXaxis()->SetRangeUser(PI/2.,3.*PI/2);
				jetHProjection[k][i][j]->GetYaxis()->SetRangeUser(0,1.1*jetHProjection[k][i][j]->GetBinContent(jetHProjection[k][i][j]->GetMaximumBin()));
				jetHProjection[k][i][j]->SetLineColor(colorList[k]);
				jetHProjection[k][i][j]->SetMarkerColor(colorList[k]);
				jetHProjection[k][i][j]->SetMarkerStyle(markerList[k]);
				if (k==0) {
					jetHProjection[k][i][j]->Draw();
				} else {
					jetHProjection[k][i][j]->Draw("SAME");
				}
				if (i==0 && j==0) {
					legASPlotCmp->AddEntry(jetHProjection[k][i][j],titleList[k],"lp");
				}
			}
			legASPlotCmp->Draw("SAME");
			canvas->Print(Form("EP_AS_PlotCmp_" + jetClass + "_" + particleClass + ".pdf",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1),assocParticlePtBins.at(j),assocParticlePtBins.at(j+1)));
			canvas->Print(Form("EP_AS_PlotCmp_" + jetClass + "_" + particleClass + ".C",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1),assocParticlePtBins.at(j),assocParticlePtBins.at(j+1)));
		}
	}



	if (acceptanceCorrection) {
		fOut->Add(fAccCorr);
	}
	for (int k = 0; k < nEPBins; k++) {
		fOut->Add(leadingJetPtEP[k]);
	}


/*
	vector<TGraphErrors *> OutOverIn_AS_JetPtBin;
	vector<TGraphErrors *> MidOverIn_AS_JetPtBin;
	vector<TGraphErrors *> OutOverIn_NS_JetPtBin;
	vector<TGraphErrors *> MidOverIn_NS_JetPtBin;

	// Parameter Method
	vector<TGraphErrors *> OutOverInPar_AS_JetPtBin;
	vector<TGraphErrors *> MidOverInPar_AS_JetPtBin;*/

		for (int i = 0; i < nTriggerPtBins; i++) {
			fOut->Add(OutOverIn_AS_JetPtBin.at(i));
			fOut->Add(MidOverIn_AS_JetPtBin.at(i));
			fOut->Add(OutOverIn_NS_JetPtBin.at(i));
			fOut->Add(MidOverIn_NS_JetPtBin.at(i));

			fOut->Add(OutOverInPar_AS_JetPtBin.at(i));
			fOut->Add(MidOverInPar_AS_JetPtBin.at(i));

			fOut->Add(OutMinusIn_AS_JetPtBin.at(i));
			fOut->Add(MidMinusIn_AS_JetPtBin.at(i));
			fOut->Add(OutMinusIn_NS_JetPtBin.at(i));
			fOut->Add(MidMinusIn_NS_JetPtBin.at(i));

			fOut->Add(RMSOutOverIn_AS_JetPtBin.at(i));
			fOut->Add(RMSMidOverIn_AS_JetPtBin.at(i));
			fOut->Add(RMSOutOverIn_NS_JetPtBin.at(i));
			fOut->Add(RMSMidOverIn_NS_JetPtBin.at(i));

			fOut->Add(WidthsOutOverIn_AS_JetPtBin.at(i));
			fOut->Add(WidthsMidOverIn_AS_JetPtBin.at(i));

		}


	if (PartonPtByJetBinEP[0][0]) {
		for (int k = 0; k < nEPBins; k++) {
			for (int i = 0; i < nJetPtBins; i++) {
				normalizeTH1F(PartonPtByJetBinEP[k][i]);
				normalizeTH1F(qqbar_PartonPtByJetBinEP[k][i]);
				normalizeTH1F(qq_PartonPtByJetBinEP[k][i]);
				normalizeTH1F(gq_PartonPtByJetBinEP[k][i]);
				normalizeTH1F(qg_PartonPtByJetBinEP[k][i]);
				normalizeTH1F(gg_PartonPtByJetBinEP[k][i]);
		
				PartonPtByJetBinEP[k][i]->SetLineColor(colorList[k]);
				PartonPtByJetBinEP[k][i]->SetMarkerColor(colorList[k]);
				PartonPtByJetBinEP[k][i]->SetMarkerStyle(markerList[k]);

			}
		}
		TLegend * leg = new TLegend(0.5,0.65,0.90,0.85);
		for (int i = 0; i < nJetPtBins; i++) {
			canvas->Clear();
			leg->Clear();
			PartonPtByJetBinEP[1][i]->Draw();
			PartonPtByJetBinEP[2][i]->Draw("SAME");
			PartonPtByJetBinEP[3][i]->Draw("SAME");
			for (int k = 1; k < nEPBins; k++) leg->AddEntry(PartonPtByJetBinEP[k][i],titleList[k],"lp");
			leg->Draw("SAME");
			canvas->Print(Form("PartonByJetBin_EPCmp_%.0f_%.0f.pdf",jetPtBins.at(i),jetPtBins.at(i+1)));
			canvas->Print(Form("PartonByJetBin_EPCmp_%.0f_%.0f.C",jetPtBins.at(i),jetPtBins.at(i+1)));
		}
		for (int k = 0; k < nEPBins; k++) {
			for (int i = 0; i < nJetPtBins; i++) {
				fOut->Add(PartonPtByJetBinEP[k][i]);
				fOut->Add(qqbar_PartonPtByJetBinEP[k][i]);
				fOut->Add(qq_PartonPtByJetBinEP[k][i]);
				fOut->Add(gq_PartonPtByJetBinEP[k][i]);
				fOut->Add(qg_PartonPtByJetBinEP[k][i]);
				fOut->Add(gg_PartonPtByJetBinEP[k][i]);
			}
		}
	}
	if (ptLossLeadingJetPtBinEP[0][0]) {

		int loss_rebin = 4;

		for (int k = 0; k < nEPBins; k++) {
			for (int i = 0; i < nJetPtBins; i++) {
				normalizeTH1F(ptLossLeadingJetPtBinEP[k][i]);
				normalizeTH1F(energyLossLeadingJetPtBinEP[k][i]);
				normalizeTH1F(ptLossLeadingPartonPtBinEP[k][i]);
				normalizeTH1F(energyLossLeadingPartonPtBinEP[k][i]);
				
				ptLossLeadingJetPtBinEP[k][i]->Rebin(loss_rebin); 
				ptLossLeadingJetPtBinEP[k][i]->Scale(1./loss_rebin); 
				energyLossLeadingJetPtBinEP[k][i]->Rebin(loss_rebin);
				energyLossLeadingJetPtBinEP[k][i]->Scale(1./loss_rebin); 
				ptLossLeadingPartonPtBinEP[k][i]->Rebin(loss_rebin); 
				ptLossLeadingPartonPtBinEP[k][i]->Scale(1./loss_rebin); 
				energyLossLeadingPartonPtBinEP[k][i]->Rebin(loss_rebin);
				energyLossLeadingPartonPtBinEP[k][i]->Scale(1./loss_rebin); 

				ptLossLeadingJetPtBinEP[k][i]->SetLineColor(colorList[k]);
				ptLossLeadingJetPtBinEP[k][i]->SetMarkerColor(colorList[k]);
				ptLossLeadingJetPtBinEP[k][i]->SetMarkerStyle(markerList[k]);

				energyLossLeadingJetPtBinEP[k][i]->SetLineColor(colorList[k]);
				energyLossLeadingJetPtBinEP[k][i]->SetMarkerColor(colorList[k]);
				energyLossLeadingJetPtBinEP[k][i]->SetMarkerStyle(markerList[k]);

				ptLossLeadingPartonPtBinEP[k][i]->SetLineColor(colorList[k]);
				ptLossLeadingPartonPtBinEP[k][i]->SetMarkerColor(colorList[k]);
				ptLossLeadingPartonPtBinEP[k][i]->SetMarkerStyle(markerList[k]);

				energyLossLeadingPartonPtBinEP[k][i]->SetLineColor(colorList[k]);
				energyLossLeadingPartonPtBinEP[k][i]->SetMarkerColor(colorList[k]);
				energyLossLeadingPartonPtBinEP[k][i]->SetMarkerStyle(markerList[k]);
			}
		}
		TLegend * leg = new TLegend(0.5,0.65,0.90,0.85);


		for (int i = 0; i < nJetPtBins; i++) {
			canvas->Clear();
			leg->Clear();
			ptLossLeadingJetPtBinEP[1][i]->Draw();
			ptLossLeadingJetPtBinEP[2][i]->Draw("SAME");
			ptLossLeadingJetPtBinEP[3][i]->Draw("SAME");
			for (int k = 1; k < nEPBins; k++) leg->AddEntry(ptLossLeadingJetPtBinEP[k][i],titleList[k],"lp");
			leg->Draw("SAME");
			canvas->Print(Form("loss_pTLossLeadingJet_EPCmp_jetPt_%.0f_%.0f.pdf",jetPtBins.at(i),jetPtBins.at(i+1)));
			canvas->Print(Form("loss_pTLossLeadingJet_EPCmp_jetPt_%.0f_%.0f.C",jetPtBins.at(i),jetPtBins.at(i+1)));
		}
		for (int k = 0; k < nEPBins; k++) {
			for (int i = 0; i < nJetPtBins; i++) {
				fOut->Add(ptLossLeadingJetPtBinEP[k][i]);
			}
		}
		for (int i = 0; i < nJetPtBins; i++) {
			canvas->Clear();
			leg->Clear();
			energyLossLeadingJetPtBinEP[1][i]->Draw();
			energyLossLeadingJetPtBinEP[2][i]->Draw("SAME");
			energyLossLeadingJetPtBinEP[3][i]->Draw("SAME");
			for (int k = 1; k < nEPBins; k++) leg->AddEntry(energyLossLeadingJetPtBinEP[k][i],titleList[k],"lp");
			leg->Draw("SAME");
			canvas->Print(Form("loss_energyLossLeadingJet_EPCmp_jetPt_%.0f_%.0f.pdf",jetPtBins.at(i),jetPtBins.at(i+1)));
			canvas->Print(Form("loss_energyLossLeadingJet_EPCmp_jetPt_%.0f_%.0f.C",jetPtBins.at(i),jetPtBins.at(i+1)));
		}
		for (int k = 0; k < nEPBins; k++) {
			for (int i = 0; i < nJetPtBins; i++) {
				fOut->Add(energyLossLeadingJetPtBinEP[k][i]);
			}
		}
		for (int i = 0; i < nJetPtBins; i++) {
			canvas->Clear();
			leg->Clear();
			ptLossLeadingPartonPtBinEP[1][i]->Draw();
			ptLossLeadingPartonPtBinEP[2][i]->Draw("SAME");
			ptLossLeadingPartonPtBinEP[3][i]->Draw("SAME");
			for (int k = 1; k < nEPBins; k++) leg->AddEntry(ptLossLeadingPartonPtBinEP[k][i],titleList[k],"lp");
			leg->Draw("SAME");
			canvas->Print(Form("loss_pTLossLeadingJet_EPCmp_partonPt_%.0f_%.0f.pdf",jetPtBins.at(i),jetPtBins.at(i+1)));
			canvas->Print(Form("loss_pTLossLeadingJet_EPCmp_partonPt_%.0f_%.0f.C",jetPtBins.at(i),jetPtBins.at(i+1)));
		}
		for (int k = 0; k < nEPBins; k++) {
			for (int i = 0; i < nJetPtBins; i++) {
				fOut->Add(ptLossLeadingPartonPtBinEP[k][i]);
			}
		}
		for (int i = 0; i < nJetPtBins; i++) {
			canvas->Clear();
			leg->Clear();
			energyLossLeadingPartonPtBinEP[1][i]->Draw();
			energyLossLeadingPartonPtBinEP[2][i]->Draw("SAME");
			energyLossLeadingPartonPtBinEP[3][i]->Draw("SAME");
			for (int k = 1; k < nEPBins; k++) leg->AddEntry(energyLossLeadingPartonPtBinEP[k][i],titleList[k],"lp");
			leg->Draw("SAME");
			canvas->Print(Form("loss_energyLossLeadingJet_EPCmp_partonPt_%.0f_%.0f.pdf",jetPtBins.at(i),jetPtBins.at(i+1)));
			canvas->Print(Form("loss_energyLossLeadingJet_EPCmp_partonPt_%.0f_%.0f.C",jetPtBins.at(i),jetPtBins.at(i+1)));
		}
		for (int k = 0; k < nEPBins; k++) {
			for (int i = 0; i < nJetPtBins; i++) {
				fOut->Add(energyLossLeadingPartonPtBinEP[k][i]);
			}
		}





/*
	TH1F * ptLossLeadingJetPtBinEP[nEPBins][nJetPtBins];
	TH1F * energyLossLeadingJetPtBinEP[nEPBins][nJetPtBins];
	TH1F * ptLossLeadingPartonPtBinEP[nEPBins][nJetPtBins];
	TH1F * energyLossLeadingPartonPtBinEP[nEPBins][nJetPtBins];

	TH1F * ptLossSubleadingJetPtBinEP[nEPBins][nJetPtBins];
	TH1F * energyLossSubleadingJetPtBinEP[nEPBins][nJetPtBins];
	TH1F * ptLossSubleadingPartonPtBinEP[nEPBins][nJetPtBins];
	TH1F * energyLossSubleadingPartonPtBinEP[nEPBins][nJetPtBins];
*/


	}	
	}

  if (ptLossVertexRLeadingJet) {
    fOut->Add(ptLossVertexRLeadingJet);
    TH2F * ptLossVertexRLeadingJetNorm = normalizeHistogramColumns(ptLossVertexRLeadingJet);
    canvas->Clear();
    ptLossVertexRLeadingJetNorm->Draw("COLZ");
    canvas->Print("ptLossVertexRLeadingJetNorm.pdf");
  }
  if (ptLossVertexRSubleadingJet) {
    fOut->Add(ptLossVertexRSubleadingJet);
    TH2F * ptLossVertexRSubleadingJetNorm = normalizeHistogramColumns(ptLossVertexRSubleadingJet);
    canvas->Clear();
    ptLossVertexRSubleadingJetNorm->Draw("COLZ");
    canvas->Print("ptLossVertexRSubleadingJetNorm.pdf");
  }



  if (hsPtLeadingJetPt) {
    hsPtLeadingJetPt->Draw("COLZ");
    canvas->SetLogz();
    canvas->Print("hsPtLeadinJetPt.pdf");
    canvas->Clear();
    hsPtSubLeadingJetPt->Draw("COLZ");
    canvas->SetLogz();
    canvas->Print("hsPtSubLeadinJetPt.pdf");
    canvas->Clear();

    //renormalized version
    TH2F * hsPtLeadingJetPtNorm = normalizeHistogramColumns(hsPtLeadingJetPt);
    hsPtLeadingJetPtNorm->Draw("COLZ");
    canvas->SetLogz();
    canvas->Print("hsPtLeadingJetPtNorm.pdf");
    canvas->Clear();
    TH2F * hsPtSubLeadingJetPtNorm = normalizeHistogramColumns(hsPtSubLeadingJetPt);
    hsPtSubLeadingJetPtNorm->Draw("COLZ");
    canvas->SetLogz();
    canvas->Print("hsPtSubLeadingJetPtNorm.pdf");
    canvas->Clear();



    canvas->SetLogz(0);
    hsDEtaDPhiLeadingJet->Draw("COLZ");
    canvas->Print("hsDEtaDPhiLeadingJet.pdf");
    canvas->Clear();
    if(hsPtSubLeadingJetPt)  hsDEtaDPhiSubLeadingJet->Draw("COLZ");
    canvas->Print("hsDEtaDPhiSubLeadingJet.pdf");
    canvas->Clear();
  
	

    fOut->Add(hsPtLeadingJetPt);
    if(hsPtSubLeadingJetPt) fOut->Add(hsPtSubLeadingJetPt);
    fOut->Add(hsDEtaDPhiLeadingJet);
    fOut->Add(hsDEtaDPhiSubLeadingJet);
  }


	if (hs1hs2Pt) {
		TH1D * hsPt = hs1hs2Pt->ProjectionX();
		TH1D * hs2Pt = hs1hs2Pt->ProjectionY();
		hsPt->Add(hs2Pt);
		delete hs2Pt;
		hsPt->SetName("hsPt");
		hsPt->SetTitle("Hard Scatter p_{T}");
		hsPt->GetXaxis()->SetTitle("p_{T}^{hs} (GeV/c)");
//		fOut->Add(hsPt);
		fOut->Add(hs1hs2Pt);
		fOut->Add(hs1hs2dPhi);
	} 

  if(hJetPtZ && hJetPtZNorm) {
    fOut->Add(hJetPtZ);
    fOut->Add(hJetPtZNorm);
  }
	if(hJetPtDR && hJetPtDRNorm) { // new addition
		fOut->Add(hJetPtDR);
		fOut->Add(hJetPtDRNorm);
	}

	if (hJetPtZ_DRCut) fOut->Add(hJetPtZ_DRCut);
	if (hJetPtZ_DRCutNorm) fOut->Add(hJetPtZ_DRCutNorm);

	if (hJetPtSD2MassNorm) fOut->Add(hJetPtSD2MassNorm);

	printf("meep moop\n");

  if (hJetPtMass && hJetPtMassNorm) {
    fOut->Add(hJetPtMass);
    fOut->Add(hJetPtMassNorm);
 	 }  
	if (hSDJetPt) fOut->Add(hSDJetPt);
  if (hSDJetPtBkgSub) fOut->Add(hSDJetPtBkgSub);

	
	printf("Initiating Segmentation Fault ... \n");
  fOut->Write();
  fOut->Close();
}

int main(int argc, char * argv[])
{
  // Defined to allow rehlersi to use a different style for printing
  // It can be enabled by passing USER_DEFINED=-DRE_STYLE=1 to the makefile
  #ifdef RE_STYLE
  // Setup TStyle
  TStyle * readableStyle = initializeReadableStyle();
  #endif
  set_plot_style();
  // Parse input
  TString inputFilename = "";
  TString outputFilename = "";

  parseInput(argc, argv, inputFilename, outputFilename);
 
	if (outputFilename=="") {
		outputFilename="Final.root";
	}

  printf("Input File: %s \t\t Output File: %s \n",inputFilename.Data(),outputFilename.Data());
  printf("dPhidEta Background Method: %s\n",bkg2DMethodArray[bkg2DMethod]);
  printf("           dPhi Fit Method: %s\n",dPhiFitMethodArray[dPhiFitMethod]);
  printf("dPhi Fit Background Option: %s\n",dPhiBkgMethodArray[dPhiBkgMethod]);

	printf("Systematic Uncertainty Input File: %s\n",SystematicUncertInFilePath.Data());
	printf("Systematic Uncertainty Output File: %s\n",SystematicUncertOutFilePath.Data());

	if (SystematicUncertInFilePath.Length() && SystematicUncertOutFilePath.Length()) {
		printf("You are both using and creating systematic uncertainty files.  This is not advised.\n");
	}

//  exit(0);
 
  phase2(inputFilename,outputFilename);

  #ifdef RE_STYLE
  delete readableStyle;
  #endif

  return 0;
}

void parseInput(int argc, char * argv[], TString & inputFilename, TString & outputFilename)
{
  if (argc == 1)
  {
    std::cout << "Error: Must specify an option!" << std::endl;
    printHelp();
  }

  for (int i=1; i < argc; i++)
  {
    if ((argv[i] == std::string("--help")) || (argv[i] == std::string("-h")))
    {
      // Displays help
      printHelp();
    }
    if ((argv[i] == std::string("--inputFilename")) || (argv[i] == std::string("-i")) )
    {
      // Sets input filename
      if (argc > i+1)
      {
        inputFilename = argv[i+1];
        i++;
        continue;
      }
      else
      {
        std::cout << "An input filename must be passed!" << std::endl;
        printHelp();
      }
    }
    if ((argv[i] == std::string("--outputFilename")) || (argv[i] == std::string("-o")) )
    {
      // Sets output filename
      if (argc > i+1)
      {
        outputFilename = argv[i+1];
        i++;
        continue;
      }
      else
      {
        std::cout << "An output filename must be passed!" << std::endl;
        printHelp();
      }
    }
    if ((argv[i] == std::string("--outputDirPath")) || (argv[i] == std::string("-od")) )
    {
      // Sets output filename
      if (argc > i+1)
      {
        outputDirPath = argv[i+1];
        i++;
        continue;
      }
      else
      {
        std::cout << "An output directory must be passed!" << std::endl;
        printHelp();
      }
    }
    if ((argv[i] == std::string("--background2D")) || (argv[i] == std::string("-b")) )
    {
      // Sets the combinatorial background method for dPhiDEta
      if (argc > i+1)
      {
        char *t;
        long n = strtol(argv[i+1],&t,10);
        if (*t) { // entry is not a number
          int cmp = 1;
          for (int j = 0; j < Nbkg2DMethods; j++) {
            cmp = strcmp(argv[i+1],bkg2DMethodArray[j]);
            if (!cmp) {
              printf("Recognized option %s = %d\n",bkg2DMethodArray[j],j);
              bkg2DMethod = j;
              break;   
            } 
          }
          if (cmp) {
            fprintf(stderr,"Error: Unrecognized bkd2DMethod option: %s\n",argv[i+1]);
            exit(1);
          }
        } else { //entry is a number
          if (0 <= n && n < Nbkg2DMethods ) bkg2DMethod = n;
          else {
            fprintf(stderr,"Error: invalid choice for background2D Method\n");
            exit(1);
          }
        }
        i++;
        continue;
      }
      else
      {
        std::cout << "A choice of method must be passed!" << std::endl;
        printHelp();
      }
    }
    if ((argv[i] == std::string("--fitDPhi")) || (argv[i] == std::string("-f")) )
    {
      // Sets the fit function for the dPhi Projection
      if (argc > i+1)
      {
        char *t;
        long n = strtol(argv[i+1],&t,10);
        if (*t) { // entry is not a number
          int cmp = 1;
          for (int j = 0; j < NdPhiFitMethods; j++) {
            cmp = strcmp(argv[i+1],dPhiFitMethodArray[j]);
            if (!cmp) {
              printf("recognized option %s = %d\n",dPhiFitMethodArray[j],j);
              dPhiFitMethod = j;
              break;   
            } 
          }
          if (cmp) {
            fprintf(stderr,"Error: Unrecognized dPhi fit option: %s\n",argv[i+1]);
            exit(1);
          }
        } else { //entry is a number
          if (0 <= n && n < NdPhiFitMethods ) dPhiFitMethod = n;
          else {
            fprintf(stderr,"Error: invalid choice for dPhi fit option\n");
            exit(1);
          }
        }
        i++;
        continue;
      }
      else
      {
        std::cout << "A choice of dphi fit function must be passed!" << std::endl;
        printHelp();
      }
    }
    if ((argv[i] == std::string("--backgroundDPhi")) || (argv[i] == std::string("-B")) )
    {
      // Sets the choice for the background offset of the dPhi fit
      if (argc > i+1)
      {
        char *t;
        long n = strtol(argv[i+1],&t,10);
        if (*t) { // entry is not a number
          int cmp = 1;
          for (int j = 0; j < NdPhiBkgMethods; j++) {
            cmp = strcmp(argv[i+1],dPhiBkgMethodArray[j]);
            if (!cmp) {
              printf("Recognized option %s = %d\n",dPhiBkgMethodArray[j],j);
              dPhiBkgMethod = j;
              break;   
            } 
          }
          if (cmp) {
            fprintf(stderr,"Error: Unrecognized dPhi Bkg option: %s\n",argv[i+1]);
            exit(1);
          }
        } else { //entry is a number
          if (0 <= n && n < NdPhiBkgMethods ) dPhiBkgMethod = n;
          else {
            fprintf(stderr,"Error: invalid choice for dPhi Bkg Method\n");
            exit(1);
          }
        }
        i++;
        continue;
      }
      else
      {
        std::cout << "A choice of method must be passed!" << std::endl;
        printHelp();
      }
    }
    if ((argv[i] == std::string("--SysUncertOut")) || (argv[i] == std::string("-so")) ) {
      // Sets SysUncert input filename
      if (argc > i+1)
      {
        SystematicUncertOutFilePath = argv[i+1];
        i++;
        continue;
      }
      else
      {
        std::cout << "An systematic uncertainty input filename must be passed!" << std::endl;
        printHelp();
      }
		}
    if ((argv[i] == std::string("--SysUncertIn")) || (argv[i] == std::string("-si")) ) {
      // Sets SysUncert input filename
      if (argc > i+1)
      {
        SystematicUncertInFilePath = argv[i+1];
        i++;
        continue;
      }
      else
      {
        std::cout << "An systematic uncertainty output filename must be passed!" << std::endl;
        printHelp();
      }
		}
    if ((argv[i] == std::string("--special")) || (argv[i] == std::string("-s")) )
    {
      // Sets special options
      if (argc > i+1)
      {
        char *t;
        long n = strtol(argv[i+1],&t,10);
        if (*t) { // entry is not a number
          fprintf(stderr,"Error: Unrecognized special case option: %s\n",argv[i+1]);
          exit(1);
        } else { //entry is a number
          switch (n) {
            case 0:
              noNearSide = true;      
              break;
            case 1:
              noAwaySide = true;
              break;
						case 2:
							noFitting = true;
							break;
						case 3:
							acceptanceCorrection = true;
							fastMixed = false;
							break;
						case 4:
							acceptanceCorrection = false;
							fastMixed = false;
							break;
						case 5:
							fastMixingMode = 2;
							break;
						case 6:
							fastMixingMode = 3;
							break;
            default:
              fprintf(stderr,"Error: invalid choice for special case\n");
              exit(1);
          }
        }
        i++;
        continue;
      }
      else
      {
        std::cout << "A choice of method must be passed!" << std::endl;
        printHelp();
      }
    }



    else
    {
      printHelp(argv[i]);
    }
  }
}

void printHelp(std::string input)
{
  if (input != "")
  {
    std::cout << "Invalid option: " << input << std::endl;
  }

  // Prints help message
  // The formatting below is carefully written, but not nearly as fraigle as it used to be.
  std::cout << "Help: This output describes the possible parameters." << std::endl;
  std::cout << "Note: Any options that are not explicitly set are assumed to be default." << std::endl << std::endl;
  std::cout << "Valid options:" << std::endl
    << std::setw(5) << std::left << "\t-h" << "\t--help"
    << "\t\t\t\t-> Displays this help message" << std::endl
    << std::setw(5) << std::left << "\t-i" << "\t--intputFilename <filename>"
    << "\t-> Sets the input filename. Default: \"root/pp.root\"" << std::endl
    << std::setw(5) << std::left << "\t-o" << "\t--outputFilename <filename>"
    << "\t-> Sets the output filename. Default: \"root/final.root\"" << std::endl
    << std::setw(5) << std::left << "\t-od" << "\t--outputDirPath <path>"
    << "\t-> Sets the output directory for plots. Default: \"output\"" << std::endl
    << std::setw(5) << std::left << "\t-b" << "\t--background2D <algo>"
    << "\t-> Sets the method for combinatorial bkg subtraction." << std::endl 
    << "\t\t Choices: \"NoBkgSub = 0\" --> Do not subtract from dPhidEta" << std::endl
    << "\t\t          \"2DBkgSub = 1\" --> Fit a constant background + NS Peak, subtract constant background" << std::endl
    << "\t\t          \"dEtaBkgSub = 2\" --> Project dEta, fit constant + NS Peak, subtract constant" << std::endl
    << "\t\t          \"dEtaFarSub = 3\" --> Project dEta, use far dEta to get constant" << std::endl
    << "\t\t          \"dEtaGenGausFitSub = 4\" --> Project dEta, fit constant + NS Gen. Gaus. Peak, subtract constant" << std::endl
    << "\t\t          \"dEtaMixedGausFitSub = 5\" --> Project dEta, fit constant + NS Gen. Gaus. Peak for pT < cut, constant + 2 Gaus for pT > cut, subtract constant" << std::endl
    << "\t\t          \"dEtaFarZYAMSub = 6\" --> Project Far dEta region, rebin, use ZYAM" << std::endl
//    << "\t\t          \"fastMixed  = 4\" --> Use rapidity yields of jets and hadrons to estimate mixed event." << std::endl
    << "\t\t Default: \"NoBkgSub\"" << std::endl
    << std::setw(5) << std::left << "\t-f" << "\t--fitDPhi <algo>"
    << "\t-> Sets the fit function for the dPhi Projection." << std::endl 
    << "\t\t Choices: \"1G1G = 0\"  --> 1 Gaussian for each peak" << std::endl
    << "\t\t          \"1M1M = 1\"  --> 1 Modified Gaussian for each peak" << std::endl
    << "\t\t          \"2G1M = 2\"  --> 2 Gaus. for NS, 1 Mod. Gaus. for AS" << std::endl
    << "\t\t          \"2G1GG = 3\" --> 2 Gaus. for NS, 1 Gen. Gaus. for AS" << std::endl
    << "\t\t          \"1GG1GG = 4\" --> 1 Gen. Gaus. for NS, 1 Gen. Gaus. for AS" << std::endl
    << "\t\t          \"2G1SH = 5\" --> 2 Gaus. for NS, 1 SecH for AS" << std::endl
    << "\t\t          \"2G2G = 6\" --> 2 Gaus. for NS, 2 Gaus for AS" << std::endl
    << "\t\t Default: \"2G1GG\"" << std::endl
    << std::setw(5) << std::left << "\t-B" << "\t--backgroundDPhi <algo>"
    << "\t-> Sets the background method for the dPhi Projection." << std::endl 
    << "\t\t          \"NoDPhiBkg = 0\"  --> dPhi Background Fixed to 0." << std::endl
    << "\t\t          \"FreeDPhiBkg = 1\"  --> dPhi Background a free parameter." << std::endl
    << "\t\t          \"ZYAM = 2\" --> dPhi Background fixed via ZYAM" << std::endl
    << "\t\t Default: \"FreeDPhiBkg\"" << std::endl
    << std::setw(5) << std::left << "\t-si" << "\t--SysUncertIn <filename>"
    << "\t-> Sets the systematic uncertainty nput filename." << std::endl
    << "\t\t Default: none" << std::endl
    << std::setw(5) << std::left << "\t-so" << "\t--SysUncertOut <filename>"
    << "\t-> Sets the systematic uncertainty output filename.  Should only" << std::endl
    << std::setw(5) << "\t\t\t\t\tbe used if calculating the systematic uncertainies from event-plane independent data." << std::endl
    << "\t\t Default: none" << std::endl
    << std::setw(5) << std::left << "\t-s" << "\t--special <n>"
    << "\t-> Sets a special mode." << std::endl
    << "\t\t          \"0\"  --> NearSide Peaks fixed to 0" << std::endl
    << "\t\t          \"1\"  --> Awayside Peaks fixed to 0" << std::endl
    << "\t\t          \"2\"  --> No fitting of dPhi correlations; just show pre-fits" << std::endl
    << "\t\t          \"3\"  --> Do finite acceptance correction, not mixed" << std::endl
    << "\t\t          \"4\"  --> Do neither finite acceptance nor mixed correction" << std::endl
    << "\t\t          \"5\"  --> Set fastmixing to use measured associated yield" << std::endl
    << "\t\t          \"6\"  --> Set fastmixing to use scattering centers yield" << std::endl;
  std::exit(-1);
}




