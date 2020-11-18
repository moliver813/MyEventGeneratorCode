// Compares two outfiles from phase2, adds appropriate labels
//#include <boost/algorithm/string/replace.hpp>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TString.h>
#include <TRegexp.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TColor.h>
#include <TStyle.h>
#include <TKey.h>

#include "analysis_params.h"

int nTriggerPtBins = nJetPtBins;
std::vector<double> fTriggerPtBins = jetPtBins;

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



using namespace std;

struct fileLabel {
  fileLabel(TString _filepath, TString _label): filepath(_filepath),label(_label) {}
  TString filepath;
  TString label;
};

/** 
  * Calculates the systematic fractional uncertainty in the tgraphs in the input matrices by using k =1,2, nVersions
  * trials.  Saves them to the parameter output file
  */
void  CalculateSystematicUncertaintyN(vector <vector<TGraphErrors * >> tGraphs, TFile * f){
  // tGraphs expected to be indexed by k = 0,1, .., nVersions, then j = 0,1, .. nJetPtBins
	int nModes = tGraphs.size();

  for (int i = 0; i < nTriggerPtBins; i++) {
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
      for (int k = 0; k < nModes; k++) {
        y=0;
        tGraphs[k][i]->GetPoint(j,x,y);
        sumOfSquares += y*y;
        sum += y;
      }
			//y = sum / nModes;// New Fix
			tGraphs[0][i]->GetPoint(j,x,y); // Newer Fix, use the first value as the primary value to be reported.
      if (y != 0) SysErr->SetPoint(j,x,(sqrt((sumOfSquares - sum*sum/(nModes))/(nModes-1))/abs(y))); 
    }
    f->Add(SysErr);
  }

}



void cleanName(TString &str) {
		str.ReplaceAll(" ","_");
		str.ReplaceAll("%","_");
		str.ReplaceAll("#","_");
		str.ReplaceAll(",","_");
		str.ReplaceAll("-","_");
		str.ReplaceAll("=","_");
}

vector <vector<TGraphErrors * >> extractTGraphArray(TString name, vector<TFile *> * files) {

	int nModes = files->size();
	vector <vector<TGraphErrors * >> tGraphArray;// = new vector <vector<TGraphErrors * >>;

	for (int k = 0; k < nModes; k++) {
		TFile * f1 = files->at(k);

		vector<TGraphErrors *> tGraphSubArray;
		for (int i = 0; i < nTriggerPtBins; i++) {
			TGraphErrors * tGraph = (TGraphErrors *) f1->Get(Form(name,fTriggerPtBins.at(i), fTriggerPtBins.at(i+1)));
			if (!tGraph) {
				fprintf(stderr,"Error: Could not find %s\n",Form(name,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)));
				exit(1);
			}
//			printf("Found %s in file %s\n",tGraph->GetName(),f1->GetTitle());
			tGraphSubArray.push_back(tGraph);
		}
		tGraphArray.push_back(tGraphSubArray);
	}
	
	return tGraphArray;
}



const Int_t nTypeList = 9;

void calcSystUncert(vector<fileLabel> *fileLabels, char * fSystematicUncertOutFilePath)  {


  Int_t colorList[nTypeList] = {kBlack,kRed,kOrange-3,kGreen-3,kBlue,kViolet,kMagenta-5,kRed-10,kGray};
  Int_t markerList[nTypeList] = {kFullSquare,kFullCircle,kFullDiamond,kOpenSquare,kOpenCircle,kOpenDiamond,kFullStar,23,34};

  double c_width = 600;
  double c_height = 600;


	TFile * fSystematicUncertOutFile = TFile::Open(fSystematicUncertOutFilePath,"RECREATE");
	if (!fSystematicUncertOutFile) {
		fprintf(stderr,"Error: Unable to open systematic error output file!\n");
		exit(1);
	}


  vector<TFile *> *files = new vector<TFile *>;

  for (int i = 0 ; i < fileLabels->size(); i++ ) {
    fileLabel fl = fileLabels->at(i);
    printf("Opening File %s [%s] ... \n",fl.filepath.Data(),fl.label.Data());  
    TFile *f = TFile::Open(fl.filepath.Data(),"READ");
    if (!f) {
      fprintf(stderr,"Error: file %s not found!\nExiting.\n",fl.filepath.Data());
      exit(1);
    }
		TString nameNoSpace  = fl.label;
		cleanName(nameNoSpace);
//		printf("MHO nameNoSpace = %s\n",nameNoSpace.Data());
    f->SetTitle(fl.label.Data());
		f->SetName(nameNoSpace);
    files->push_back(f);
    
    printf("file %s, title %s\n",f->GetName(),f->GetTitle());

  }

  if ( !files->size()) exit(0);
  TFile *f0 = files->at(0);
   
  TCanvas *canvas = new TCanvas("canvas","canvas",c_width,c_height);
  gStyle->SetOptStat(0);

	// Need to construct arrays of objects, with an arbitrary first index,
  // second index is jetpt

	int iTrigger = 0; // 0 for jet, 1 for pi0Pt, 2 for pi0Zt

  TString jetClass = "jetPt_%.0f_%.0f";

	TString triggerClass = jetClass;

	if (iTrigger > 0) {
		triggerClass = "pi0Pt_%.0f_%.0f";
		nTriggerPtBins = nPi0PtBins;
		fTriggerPtBins = pi0PtBins;
	}

	vector <vector<TGraphErrors * >> ptBinASIntegralsGraphsEP_0 = extractTGraphArray(triggerClass + "_EP_0_AS_Integrals", files);
	vector <vector<TGraphErrors * >> ptBinASIntegralsGraphsEP_1 = extractTGraphArray(triggerClass + "_EP_1_AS_Integrals", files);
	vector <vector<TGraphErrors * >> ptBinASIntegralsGraphsEP_2 = extractTGraphArray(triggerClass + "_EP_2_AS_Integrals", files);
	vector <vector<TGraphErrors * >> ptBinASIntegralsGraphsEP_3 = extractTGraphArray(triggerClass + "_EP_3_AS_Integrals", files);

	CalculateSystematicUncertaintyN(ptBinASIntegralsGraphsEP_0,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinASIntegralsGraphsEP_1,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinASIntegralsGraphsEP_2,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinASIntegralsGraphsEP_3,fSystematicUncertOutFile);


	vector <vector<TGraphErrors * >> ptBinNSIntegralsGraphsEP_0 = extractTGraphArray(triggerClass + "_EP_0_NS_Integrals", files);
	vector <vector<TGraphErrors * >> ptBinNSIntegralsGraphsEP_1 = extractTGraphArray(triggerClass + "_EP_1_NS_Integrals", files);
	vector <vector<TGraphErrors * >> ptBinNSIntegralsGraphsEP_2 = extractTGraphArray(triggerClass + "_EP_2_NS_Integrals", files);
	vector <vector<TGraphErrors * >> ptBinNSIntegralsGraphsEP_3 = extractTGraphArray(triggerClass + "_EP_3_NS_Integrals", files);
	
	CalculateSystematicUncertaintyN(ptBinNSIntegralsGraphsEP_0,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinNSIntegralsGraphsEP_1,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinNSIntegralsGraphsEP_2,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinNSIntegralsGraphsEP_3,fSystematicUncertOutFile);


	vector <vector<TGraphErrors * >> ptBinASYieldsGraphsEP_0 = extractTGraphArray(triggerClass + "_EP_0_AS_Yields", files);
	vector <vector<TGraphErrors * >> ptBinASYieldsGraphsEP_1 = extractTGraphArray(triggerClass + "_EP_1_AS_Yields", files);
	vector <vector<TGraphErrors * >> ptBinASYieldsGraphsEP_2 = extractTGraphArray(triggerClass + "_EP_2_AS_Yields", files);
	vector <vector<TGraphErrors * >> ptBinASYieldsGraphsEP_3 = extractTGraphArray(triggerClass + "_EP_3_AS_Yields", files);

	CalculateSystematicUncertaintyN(ptBinASYieldsGraphsEP_0,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinASYieldsGraphsEP_1,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinASYieldsGraphsEP_2,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinASYieldsGraphsEP_3,fSystematicUncertOutFile);


	vector <vector<TGraphErrors * >> ptBinASRmsGraphsEP_0 = extractTGraphArray(triggerClass + "_EP_0_AS_Rms", files);
	vector <vector<TGraphErrors * >> ptBinASRmsGraphsEP_1 = extractTGraphArray(triggerClass + "_EP_1_AS_Rms", files);
	vector <vector<TGraphErrors * >> ptBinASRmsGraphsEP_2 = extractTGraphArray(triggerClass + "_EP_2_AS_Rms", files);
	vector <vector<TGraphErrors * >> ptBinASRmsGraphsEP_3 = extractTGraphArray(triggerClass + "_EP_3_AS_Rms", files);

	CalculateSystematicUncertaintyN(ptBinASRmsGraphsEP_0,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinASRmsGraphsEP_1,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinASRmsGraphsEP_2,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinASRmsGraphsEP_3,fSystematicUncertOutFile);


	vector <vector<TGraphErrors * >> ptBinNSRmsGraphsEP_0 = extractTGraphArray(triggerClass + "_EP_0_NS_Rms", files);
	vector <vector<TGraphErrors * >> ptBinNSRmsGraphsEP_1 = extractTGraphArray(triggerClass + "_EP_1_NS_Rms", files);
	vector <vector<TGraphErrors * >> ptBinNSRmsGraphsEP_2 = extractTGraphArray(triggerClass + "_EP_2_NS_Rms", files);
	vector <vector<TGraphErrors * >> ptBinNSRmsGraphsEP_3 = extractTGraphArray(triggerClass + "_EP_3_NS_Rms", files);
	
	CalculateSystematicUncertaintyN(ptBinNSRmsGraphsEP_0,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinNSRmsGraphsEP_1,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinNSRmsGraphsEP_2,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinNSRmsGraphsEP_3,fSystematicUncertOutFile);


	vector <vector<TGraphErrors * >> ptBinASWidthsGraphsEP_0 = extractTGraphArray(triggerClass + "_EP_0_AS_Widths", files);
	vector <vector<TGraphErrors * >> ptBinASWidthsGraphsEP_1 = extractTGraphArray(triggerClass + "_EP_1_AS_Widths", files);
	vector <vector<TGraphErrors * >> ptBinASWidthsGraphsEP_2 = extractTGraphArray(triggerClass + "_EP_2_AS_Widths", files);
	vector <vector<TGraphErrors * >> ptBinASWidthsGraphsEP_3 = extractTGraphArray(triggerClass + "_EP_3_AS_Widths", files);

	CalculateSystematicUncertaintyN(ptBinASWidthsGraphsEP_0,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinASWidthsGraphsEP_1,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinASWidthsGraphsEP_2,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(ptBinASWidthsGraphsEP_3,fSystematicUncertOutFile);



	// Now doing the Ratios 
	vector <vector<TGraphErrors * >> OutOverIn_AS = extractTGraphArray("OutOverIn_AS_" + triggerClass, files);
	vector <vector<TGraphErrors * >> MidOverIn_AS = extractTGraphArray("MidOverIn_AS_" + triggerClass, files);
	vector <vector<TGraphErrors * >> OutOverIn_NS = extractTGraphArray("OutOverIn_NS_" + triggerClass, files);
	vector <vector<TGraphErrors * >> MidOverIn_NS = extractTGraphArray("MidOverIn_NS_" + triggerClass, files);
	vector <vector<TGraphErrors * >> OutOverInPar_AS = extractTGraphArray("OutOverInPar_AS_" + triggerClass, files);
	vector <vector<TGraphErrors * >> MidOverInPar_AS = extractTGraphArray("MidOverInPar_AS_" + triggerClass, files);

	CalculateSystematicUncertaintyN(OutOverIn_AS,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(MidOverIn_AS,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(OutOverIn_NS,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(MidOverIn_NS,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(OutOverInPar_AS,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(MidOverInPar_AS,fSystematicUncertOutFile);


	vector <vector<TGraphErrors * >> OutMinusIn_AS = extractTGraphArray("OutMinusIn_AS_" + triggerClass, files);
	vector <vector<TGraphErrors * >> MidMinusIn_AS = extractTGraphArray("MidMinusIn_AS_" + triggerClass, files);
	vector <vector<TGraphErrors * >> OutMinusIn_NS = extractTGraphArray("OutMinusIn_NS_" + triggerClass, files);
	vector <vector<TGraphErrors * >> MidMinusIn_NS = extractTGraphArray("MidMinusIn_NS_" + triggerClass, files);

	CalculateSystematicUncertaintyN(OutMinusIn_AS,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(MidMinusIn_AS,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(OutMinusIn_NS,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(MidMinusIn_NS,fSystematicUncertOutFile);

	vector <vector<TGraphErrors * >> RMSOutOverIn_AS = extractTGraphArray("RMSOutOverIn_AS_" + triggerClass, files);
	vector <vector<TGraphErrors * >> RMSMidOverIn_AS = extractTGraphArray("RMSMidOverIn_AS_" + triggerClass, files);
	vector <vector<TGraphErrors * >> RMSOutOverIn_NS = extractTGraphArray("RMSOutOverIn_NS_" + triggerClass, files);
	vector <vector<TGraphErrors * >> RMSMidOverIn_NS = extractTGraphArray("RMSMidOverIn_NS_" + triggerClass, files);

	CalculateSystematicUncertaintyN(RMSOutOverIn_AS,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(RMSMidOverIn_AS,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(RMSOutOverIn_NS,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(RMSMidOverIn_NS,fSystematicUncertOutFile);


	// FIXME add width ratios?
	vector <vector<TGraphErrors * >> WidthsOutOverIn_AS = extractTGraphArray("WidthsOutOverIn_AS_" + triggerClass, files);
	vector <vector<TGraphErrors * >> WidthsMidOverIn_AS = extractTGraphArray("WidthsMidOverIn_AS_" + triggerClass, files);

	CalculateSystematicUncertaintyN(WidthsOutOverIn_AS,fSystematicUncertOutFile);
	CalculateSystematicUncertaintyN(WidthsMidOverIn_AS,fSystematicUncertOutFile);

// MidOverIn_AS_pi0Pt

  /*
  vector <vector<TGraphErrors *> > ptBinASWidthsGraphsEP;
  vector <vector<TGraphErrors *> > ptBinASRmsGraphsEP;
  vector <vector<TGraphErrors *> > ptBinASIntegralsGraphsEP;
  vector <vector<TGraphErrors *> > ptBinASFwhmGraphsEP;
  vector <vector<TGraphErrors *> > ptBinASFwhmRescaledGraphsEP;
  vector <vector<TGraphErrors *> > ptBinASYieldsGraphsEP;
*/
 
	fSystematicUncertOutFile->Write();

	printf("Finished calculating systematic uncertainties for integrals,yields,rms,width\n");

}

int main(int argc, char * argv[]) {

  if (argc < 3) {
    printf("Usage: %s <SystematicUncertaintyRootFile> [Filepath_1] [Label_1] [Filepath_2] [Label_2] ... \n",argv[0]);
//    printf("    Or %s -p [PATTERN] [Filepath_1] [Label_1] [Filepath_2] [Label_2] ... \n",argv[0]);
    exit(0);
  }



  int num = ((argc) / 2) ;
  vector<fileLabel> *fileLabels = new vector<fileLabel>;
  for (int i = 1; i < num; i++) {
    fileLabels->push_back(fileLabel(TString(argv[2*i]),TString(argv[2*i+1])));
  }
	printf("Starting calculation of systematic uncertainties using files:");
	for (int i = 0; i < num - 1 ; i++) {
		printf("%s ",fileLabels->at(i).filepath.Data());
	}
	printf("\n");
  calcSystUncert(fileLabels,argv[1]);
}




