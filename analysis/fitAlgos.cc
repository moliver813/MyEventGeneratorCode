
#ifdef RE_STYLE
#include <readableStyle.h>
#endif

#define PI TMath::Pi()

#include <vector>
#include <math.h>
#include <cmath>
#include <iostream>
#include <cstdlib>

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

#include "phase2.h"
//#include "analysis_params.h"


extern bool noAwaySide;
extern bool noNearSide;
extern bool noFitting;

using namespace std;

/** Method to print all of the parameter limits, with or without errors
  * 
  */


void printParameters(TF1 * func, bool includeErrors, std::string label = "Parameters") {

  int nPar = func->GetNpar();
  double parMin;
  double parMax;
  printf("======================================================================\n");
  printf("||  %s\n",label.c_str());
  printf("======================================================================\n");
  const int nBuff = 15;

  for (int i = 0; i < nPar; i++) {
    func->GetParLimits(i,parMin,parMax);
    const char *parName = func->GetParName(i);
    char parNameBuffer[nBuff] = "              ";
    strncpy(parNameBuffer,parName,min((int) strlen(parName),nBuff));
    int parNameLen = strlen(parName);
    if (includeErrors) {
      printf("%s:\t[%f| \t%f \\pm %f \t|%f]  \n",parNameBuffer,parMin,func->GetParameter(i),func->GetParError(i),parMax);
    //  printf("%s:\t\t\t[%.2f, %.2f \\pm %.2f  ,%.2f]  \n",func->GetParName(i),parMin,func->GetParameter(i),func->GetParError(i),parMax);
    } else {
      printf("%s:\t[%f| \t%f \t|%f]  \n",parNameBuffer,parMin,func->GetParameter(i),parMax);
    }
  }
	if (includeErrors) {
			// If errors are included, also show ChiSq, NDF
			printf("Chi^{2} =  %f\t\tNDF = %d\n",func->GetChisquare(),func->GetNDF());	
	}
  printf("======================================================================\n");
  return;
}


double getDPhiBackground(TH1D * hist){
  double B_g = 0;
// another method to add
/*
  //Method 2 fit pi/4,3pi/4 with v2, take minimum
  hist_clone->GetXaxis()->SetRangeUser(0.6,2.2);
  //TF1 *v2_func = new TF1("v2_func","[0]*TMath::Cos(x - TMath::Pi()/2.0) + [1]",PI/4,3.0*PI/4);
  TF1 *v2_func = new TF1("v2_func","[0]*TMath::Cos(x - [2]) + [1]",PI/4,3.0*PI/4);
  v2_func->SetParLimits(2,0,PI);
  v2_func->SetParameter(0,-1);
  v2_func->SetParameter(1,1);
  v2_func->SetParameter(2,hist_clone->GetXaxis()->GetBinCenter(hist_clone->GetMinimumBin()));
  // v2_func->SetParameter(2,PI/2);
  hist_clone->Fit(v2_func,"R0Q WW");
  //  B_g = max((Double_t) v2_func->GetMinimum(),0.);// minimum value, but not below zero
  B_g = 0 ;// FIXME
  printf("Subtracting background ...\n");
  TF1 *constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
  // Subtracting constant background
  constant_function->SetParameter(0,-B_g);
  hist_clone->Add(constant_function);
*/


  switch (dPhiBkgMethod) {
    case aNoDPhiBkg :
      B_g = 0 ;
      break;
    case aFreeDPhiBkg :
      // use ZYAM method for estimating other parameters
      B_g = hist->GetBinContent(hist->GetMinimumBin());
      break;
    case aZYAM:
      B_g = hist->GetBinContent(hist->GetMinimumBin());
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }
  return B_g;
}

// Fitting functions which returns the width
TF1 * fitPeak(TH1D * hist, minMax peak, std::string name)
{
  // Set name
  // Fit to the gain curves
  TF1 * peakFit = new TF1(name.c_str(), "gaus", peak.min, peak.max);

  // Fit to histogram
  // R uses the range defined in the fit function
  // 0 ensures that the fit isn't drawn
  // Q ensures minimum printing
  // + adds the fit to function list to ensure that it is not deleted on the creation of a new fit
  printParameters(peakFit,false,"Parameter Guesses");
  //FIT HERE
  hist->Fit(peakFit, "RQ0+");
  printParameters(peakFit,true,"Final Parameters");
  //hist->Fit(peakFit, "R0+");

  //double peakMean = peakFit->GetParameter("Mean");
  //double width = peakFit->GetParameter("Sigma");

  return peakFit;
}


// Fit correlation data to function used in STAR paper:
//  C(dPhi) = 
//    (Y_{NS}/sqrt(2*pi*sigma_{NS}^2)e^(-dPhi^2/(2*pi*sigma_{NS}^2)) + 
//    (Y_{AS}/sqrt(2*pi*sigma_{AS}^2)e^(-dPhi^2/(2*pi*sigma_{AS}^2)) + 
//    

// can use TMath::Gaus(x,mean,sigma,kTrue)  for 1/sqrt(2*pi*sigma^2) * e ^( ...
TF1 * phiCorrFit(TH1D * hist, std::string name)
{
  TF1 * fit = new TF1(name.c_str(),"[0]*TMath::Gaus(x,0,[1],1) + [2]*TMath::Gaus(x,TMath::Pi(),[3],1) + [4]*(1 + 2*[5]*cos(2*x) + 2*[6]*cos(3*x))",-PI/2,3.0*PI/2);
  //name parameters
  fit->SetParName(0,"Y_NS");
  fit->SetParName(1,"Sigma_NS");
  fit->SetParName(2,"Y_AS");
  fit->SetParName(3,"Sigma_AS");
  fit->SetParName(4,"B");
  fit->SetParName(5,"v2");
  fit->SetParName(6,"v3");

  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  //
  // 3.94 \aprox (pi**(3/2))/sqrt(2) 
  fit->SetParLimits(0,0.,3.94*max_val);
  fit->SetParLimits(1,0.,PI/2);
  fit->SetParLimits(2,0.,3.94*max_val);
  fit->SetParLimits(3,0.,PI/2);
  fit->SetParLimits(4,0.,max_val);
  fit->SetParLimits(5,-1,1);
  fit->SetParLimits(6,-1,1);

  double sigma_ns_g = 0.8;
  double sigma_as_g = 0.8;
  double v2_g = 0;
  double v3_g = 0;
  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();
  //Method 1 for B_g 
  //to get good RMS guess: remove constant background via ZYAM
  //  B_g = hist->GetBinContent(hist->GetMinimumBin());  

  //Method 2 fit pi/4,3pi/4 with v2, take minimum
  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/8,5.0*PI/8);
  TF1 *v2_func = new TF1("v2_func","[0]*TMath::Cos(x - 3.1415927/2.0) + [1]",PI/4,3.0*PI/4);
  hist_clone->Fit(v2_func,"R0Q");
  B_g = v2_func->Eval(PI/2);

  TF1 *constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
  // Subtracting constant background
  constant_function->SetParameter(0,-B_g);
  hist_clone->Add(constant_function);

  hist_clone->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  sigma_ns_g = hist_clone->GetRMS() /  1.3;
  hist_clone->GetXaxis()->SetRangeUser(-PI/8,PI/8);
  //  double Y_NS_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_ns_g;
  // new method: use the A*cos(x) + B to avoid fluctuation effects.
  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);
  smooth_func->SetParameter(0,1);
  smooth_func->SetParameter(1,sigma_ns_g);
  smooth_func->SetParLimits(0,0,1e5);
  smooth_func->SetParLimits(1,0,PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  double Y_NS_g = smooth_func->Eval(0)*TMath::Sqrt(2*PI)*sigma_ns_g;


  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  sigma_as_g = hist_clone->GetRMS() / 1.5;
  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  delete smooth_func;
  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",-PI/2,3.0*PI/2);
  smooth_func->SetParameter(0,1);
  smooth_func->SetParameter(1,sigma_as_g);
  smooth_func->SetParLimits(0,0,1e5);
  smooth_func->SetParLimits(1,0,PI/2);
  // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
  hist_clone->Fit(smooth_func,"R0Q");
  double Y_as_g = smooth_func->Eval(PI) *TMath::Sqrt(2*PI)*sigma_as_g;

  //	printf("Guess NS %f \t\t AS %f \n",Y_NS_g,Y_as_g);

  // double Y_NS_g = hist_clone->GetBinContent(hist->FindBin(0))*TMath::Sqrt(2*PI)*sigma_ns_g;
  //double Y_as_g = hist_clone->GetBinContent(hist->FindBin(PI))*TMath::Sqrt(2*PI)*sigma_as_g;

  delete hist_clone;

  //  printf("--------------------------------\n");
  //  printf("sigma_ns_g = %f\t\t sigma_as_g = %f\n",sigma_ns_g,sigma_as_g);
  //  printf("--------------------------------\n");
  // I could also integrate the histogram times 1, cos(2x) cos(3x) 
  // to get guesses for the background

  // will need a function that takes a histogram, a minmax, a mean, and then calculates 
  // the standard deviation, assuming the given mean is the correct one.

  //set parameters (guessing)
  fit->SetParameter(0,Y_NS_g);
  //  printf("Guess Y_NS = %f\n",fit->GetParameter(0));
  fit->SetParameter(1,sigma_ns_g);
  fit->SetParameter(2,Y_as_g);
  fit->SetParameter(3,sigma_as_g);
  fit->SetParameter(4,B_g);
  fit->SetParameter(5,v2_g);
  fit->SetParameter(6,v3_g);

  //FIXME trying turning off the v2,v3
  fit->FixParameter(5,0);
  fit->FixParameter(6,0);


	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(2,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
    fit->FixParameter(2,0);
  }


  // Fit to histogram
  // R uses the range defined in the fit function
  // 0 ensures that the fit isn't drawn
  // Q ensures minimum printing
  // + adds the fit to function list to ensure that it is not deleted on the creation of a new fit

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");

  //	hist->Fit(fit, "R0+");
  return fit;
}


// Simple, Two Gaussians + Constant background
TF1 * phiCorrFit_1Gaus_1Gaus(TH1D * hist, std::string name)
{
  TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2] * (TMath::Gaus(x,-TMath::Pi(),[3],1) + TMath::Gaus(x,TMath::Pi(),[3],1) + TMath::Gaus(x,3*TMath::Pi(),[3],1)) + [4]",-PI/2,3.0*PI/2);
  //name parameters
  fit->SetParName(0,"Y_NS");
  fit->SetParName(1,"Sigma_NS");
  fit->SetParName(2,"Y_AS");
  fit->SetParName(3,"Sigma_AS");
  fit->SetParName(4,"B");

  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  //
  // 3.94 \aprox (pi**(3/2))/sqrt(2) 
  double sigma_ns_g = 0.8;
  double sigma_as_g = 0.8;
  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();
  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }

  //Near side prefit analysis
  // hist_clone->GetXaxis()->SetRangeUser(-PI/2,PI/2);
//  hist_clone->GetXaxis()->SetRangeUser(-R,R);
//  sigma_ns_g = TMath::Max(hist_clone->GetRMS(),0.2);

  double Y_ns_g;

  //hist_clone->GetXaxis()->SetRangeUser(-PI/8,PI/8);
  //  double Y_NS_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_ns_g;
  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);
  smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  smooth_func->SetParameter(1,sigma_ns_g);

  smooth_func->SetParLimits(0,0,2*max_val);
  smooth_func->SetParLimits(1,0.05,2.0*PI/3);

  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_g = smooth_func->GetParameter(0);
  sigma_ns_g = smooth_func->GetParameter(1);


  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  sigma_as_g = TMath::Max(hist_clone->GetRMS() / 1.5,0.5);
  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  delete smooth_func;
  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",-PI/2,3.0*PI/2);
  smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  smooth_func->SetParameter(1,sigma_as_g);
  smooth_func->SetParLimits(0,0,2*max_val);
  smooth_func->SetParLimits(1,0,PI);
  // // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
  hist_clone->Fit(smooth_func,"R0Q");
  sigma_as_g = smooth_func->GetParameter(1);
  double Y_as_g = smooth_func->GetParameter(0);


  delete hist_clone;

  //  printf("--------------------------------\n");
  //  printf("sigma_ns_g = %f\t\t sigma_as_g = %f\n",sigma_ns_g,sigma_as_g);
  //  printf("--------------------------------\n");
  // I could also integrate the histogram times 1, cos(2x) cos(3x) 
  // to get guesses for the background

  // will need a function that takes a histogram, a minmax, a mean, and then calculates 
  // the standard deviation, assuming the given mean is the correct one.


  fit->SetParLimits(0,0.,3.94*(max_val - B_g));
  fit->SetParLimits(1,0.,PI);
  fit->SetParLimits(2,0.,3.94*(max_val - B_g));
  fit->SetParLimits(3,0.,PI);
  fit->SetParLimits(4,0,max_val);

  //set parameters (guessing)
  fit->SetParameter(0,Y_ns_g);
  fit->SetParameter(1,sigma_ns_g);
  fit->SetParameter(2,Y_as_g);
  fit->SetParameter(3,sigma_as_g);
  //FIXME??



	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(2,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
  }



//  fit->SetParameter(6,B_g);
  // fit->FixParameter(6,B_g);
  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(4,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(4,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(4,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }

  // Fit to histogram
  // R uses the range defined in the fit function
  // 0 ensures that the fit isn't drawn
  // Q ensures minimum printing
  // + adds the fit to function list to ensure that it is not deleted on the creation of a new fit

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");
  return fit;
}





// Simple, Two Gaussians + Constant background, using FWHM for AS width
TF1 * phiCorrFit_1Gaus_1Gaus_Fwhm(TH1D * hist, std::string name)
{
  TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2] * (2./[3]) * TMath::Sqrt(TMath::Log(2) / (2 * TMath::Pi() )) * (TMath::Exp(-TMath::Power(x+TMath::Pi(),2)*4*TMath::Log(2)/TMath::Power([3],2)) + TMath::Exp(-TMath::Power(x-TMath::Pi(),2)*4*TMath::Log(2)/TMath::Power([3],2)) + TMath::Exp(-TMath::Power(x-3.*TMath::Pi(),2)*4*TMath::Log(2)/TMath::Power([3],2)) )+ [4]",-PI/2,3.0*PI/2);
  //name parameters 
  fit->SetParName(0,"Y_NS");
  fit->SetParName(1,"Sigma_NS");
  fit->SetParName(2,"Y_AS");
  fit->SetParName(3,"Omega_AS");
  fit->SetParName(4,"B");

  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  //
  // 3.94 \aprox (pi**(3/2))/sqrt(2) 
  double sigma_ns_g = 0.8;
  double omega_as_g = 0.8;
  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();
  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }

  //Near side prefit analysis
  // hist_clone->GetXaxis()->SetRangeUser(-PI/2,PI/2);
//  hist_clone->GetXaxis()->SetRangeUser(-R,R);
//  sigma_ns_g = TMath::Max(hist_clone->GetRMS(),0.2);

  double Y_ns_g;

  //hist_clone->GetXaxis()->SetRangeUser(-PI/8,PI/8);
  //  double Y_NS_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_ns_g;
  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);
	hist_clone->GetXaxis()->SetRangeUser(-PI/4,PI/4);
	Y_ns_g = hist_clone->Integral("width");
  smooth_func->SetParameter(0,Y_ns_g);
 // smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  smooth_func->SetParameter(1,sigma_ns_g);

  smooth_func->SetParLimits(0,0,2*max_val);
  smooth_func->SetParLimits(1,0.05,2.0*PI/3);

  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_g = smooth_func->GetParameter(0);
  sigma_ns_g = smooth_func->GetParameter(1);

	double Y_as_g;

  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
	Y_as_g = hist_clone->Integral("width");
  omega_as_g = TMath::Max(hist_clone->GetRMS() / 1.5,0.5);
  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  delete smooth_func;
  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",-PI/2,3.0*PI/2);
  smooth_func->SetParameter(0,Y_as_g);
//  smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*omega_as_g);
  smooth_func->SetParameter(1,omega_as_g);
  smooth_func->SetParLimits(0,0,2*max_val);
  smooth_func->SetParLimits(1,0,PI);
  // // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
  hist_clone->Fit(smooth_func,"R0Q");
  omega_as_g = smooth_func->GetParameter(1);
  Y_as_g = smooth_func->GetParameter(0);
//  double Y_as_g = smooth_func->Eval(PI) *TMath::Sqrt(2*PI)*omega_as_g;
// correct omega to FWHM
	omega_as_g = 2 * omega_as_g * TMath::Sqrt(2.*TMath::Log(2));


  delete hist_clone;

  //  printf("--------------------------------\n");
  //  printf("sigma_ns_g = %f\t\t sigma_as_g = %f\n",sigma_ns_g,sigma_as_g);
  //  printf("--------------------------------\n");
  // I could also integrate the histogram times 1, cos(2x) cos(3x) 
  // to get guesses for the background

  // will need a function that takes a histogram, a minmax, a mean, and then calculates 
  // the standard deviation, assuming the given mean is the correct one.


 // fit->SetParLimits(0,0.,3.94*(max_val - B_g));
  fit->SetParLimits(0,0.,10*(max_val - B_g));
  fit->SetParLimits(1,0.,PI);
  fit->SetParLimits(2,0.,3.94*(max_val - B_g));
  fit->SetParLimits(3,0.,PI);
  fit->SetParLimits(4,0,max_val);

  //set parameters (guessing)
  fit->SetParameter(0,Y_ns_g);
  fit->SetParameter(1,sigma_ns_g);
  fit->SetParameter(2,Y_as_g);
  fit->SetParameter(3,omega_as_g);

	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(2,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
  }

  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(4,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(4,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(4,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }

  // Fit to histogram
  // R uses the range defined in the fit function
  // 0 ensures that the fit isn't drawn
  // Q ensures minimum printing
  // + adds the fit to function list to ensure that it is not deleted on the creation of a new fit

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");
  return fit;
}











// second method: use two gaussians for near side, one gaussian for away side
// + constant backgorund
// DEV note: try generalized gaussian for AS?  --> will make widths, yields, 
// slightly more complicated

// DEV note: running into wraparond, where a gauusian is wide enough to have a small 
// contribution to the nearside.  This currently only is reflected in the PI/2 area,
// not in the cut -PI/2, 3PI/2 area
// possible solutions: symmetrize and fit just in [0,PI]
// 										 modify function such that the gaussian 
// 										   is properly wrapped around -> tricky, since this is 
//										   an infinite series.
//  									   self-response: i'm an idiot
// 										   just add another AS at -PI and 2 PI
// DEV note: could change parameters in widths for nearside such that one is
// always wider than thother: [1] = first width, [1] + [3] = second widths


TF1 * phiCorrFit_2Gaus(TH1D * hist, std::string name)
{
  TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * (TMath::Gaus(x,-TMath::Pi(),[5],1) + TMath::Gaus(x,TMath::Pi(),[5],1) + TMath::Gaus(x,3*TMath::Pi(),[5],1)) + [6]",-PI/2,3.0*PI/2);
  //name parameters
  fit->SetParName(0,"Y_NS_1");
  fit->SetParName(1,"Sigma_1_NS");
  fit->SetParName(2,"Y_NS_2");
  // NOTE: sigma_2_NS is the difference between sigma_2 and sigma_1
  fit->SetParName(3,"Sigma_2_NS");
  fit->SetParName(4,"Y_AS");
  fit->SetParName(5,"Sigma_AS");
  fit->SetParName(6,"B");

  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  //
  // 3.94 \aprox (pi**(3/2))/sqrt(2) 

  double sigma_ns_g = 0.8;
  double sigma_ns_1_g = 0.3;
  double sigma_ns_2_g = 0.3;
  double sigma_as_g = 0.8;
  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();
  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }


  //Near side prefit analysis
  // hist_clone->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  hist_clone->GetXaxis()->SetRangeUser(-R,R);
  sigma_ns_g = TMath::Max(hist_clone->GetRMS(),0.2);

  //	sigma_ns_1_g = sigma_ns_g/1.3;
  //	sigma_ns_2_g = sigma_ns_g/3;



  double Y_ns_1_g,Y_ns_2_g;

  //hist_clone->GetXaxis()->SetRangeUser(-PI/8,PI/8);
  //  double Y_NS_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_ns_g;
  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);
  smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  smooth_func->SetParameter(1,sigma_ns_g);

  smooth_func->SetParLimits(0,0,2*max_val);
  smooth_func->SetParLimits(1,0.05,2.0*PI/3);

  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_1_g = smooth_func->GetParameter(0);
  sigma_ns_1_g = smooth_func->GetParameter(1);

  // new idea: use separate fits: one for peak < R, on for peak > R
  // since root doesn't let me fit using two ranges, I will use each
  // range separately, and average

  smooth_func->SetParameter(1,R);
  smooth_func->SetParLimits(1,R/1.2,2.0*PI/3);
  hist_clone->GetXaxis()->SetRangeUser(-PI/2,-R);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = smooth_func->GetParameter(0);
  sigma_ns_2_g = smooth_func->GetParameter(1);

  hist_clone->GetXaxis()->SetRangeUser(R,PI/2);
  hist_clone->Fit(smooth_func,"R0Q");\
    Y_ns_2_g = (Y_ns_2_g + smooth_func->GetParameter(0)) / 2.;
  sigma_ns_2_g = (sigma_ns_2_g + smooth_func->GetParameter(1)) / 2.;


  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  sigma_as_g = TMath::Max(hist_clone->GetRMS() / 1.5,0.5);
  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  delete smooth_func;
  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",-PI/2,3.0*PI/2);
  smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  smooth_func->SetParameter(1,sigma_as_g);
  smooth_func->SetParLimits(0,0,2*max_val);
  smooth_func->SetParLimits(1,0,PI);
  // // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
  hist_clone->Fit(smooth_func,"R0Q");
  sigma_as_g = smooth_func->GetParameter(1);
  double Y_as_g = smooth_func->Eval(PI) *TMath::Sqrt(2*PI)*sigma_as_g;


  delete hist_clone;

  //  printf("--------------------------------\n");
  //  printf("sigma_ns_g = %f\t\t sigma_as_g = %f\n",sigma_ns_g,sigma_as_g);
  //  printf("--------------------------------\n");
  // I could also integrate the histogram times 1, cos(2x) cos(3x) 
  // to get guesses for the background

  // will need a function that takes a histogram, a minmax, a mean, and then calculates 
  // the standard deviation, assuming the given mean is the correct one.


  fit->SetParLimits(0,0.,3.94*(max_val - B_g));
  fit->SetParLimits(1,0.,PI);
  fit->SetParLimits(2,0.,3.94*(max_val - B_g));
  fit->SetParLimits(3,0.,PI);
  fit->SetParLimits(4,0.,3.94*(max_val - B_g));
  fit->SetParLimits(5,0.,PI);
  fit->SetParLimits(6,0,max_val);

  //set parameters (guessing)
  fit->SetParameter(0,Y_ns_1_g);
  fit->SetParameter(1,sigma_ns_1_g);
  fit->SetParameter(2,Y_ns_2_g);
  fit->SetParameter(3,sigma_ns_2_g);
  fit->SetParameter(4,Y_as_g);
  fit->SetParameter(5,sigma_as_g);

	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(4,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
    fit->FixParameter(2,0);
  }

  //FIXME??
//  fit->SetParameter(6,B_g);
  // fit->FixParameter(6,B_g);
  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(6,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(6,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(6,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }

  // Fit to histogram
  // R uses the range defined in the fit function
  // 0 ensures that the fit isn't drawn
  // Q ensures minimum printing
  // + adds the fit to function list to ensure that it is not deleted on the creation of a new fit

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");
  return fit;
}

/** New method: 2 gaussians for nearside, and a generalized gaussian for the 
 * awayside
 *
 */

TF1 * phiCorrFit_2Gaus_1GenGaus(TH1D * hist, std::string name)
{
  TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * (2.*[6] * pow(TMath::Gamma(3./[6]),0.5))/ (2.*[5]*pow(TMath::Gamma(1./[6]),1.5)) * (TMath::Exp(-pow(abs(x-TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])) + TMath::Exp(-pow(abs(x+TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])) +TMath::Exp(-pow(abs(x-3*TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6]))  ) + [7]",-PI/2,3.0*PI/2);

  //name parameters
  fit->SetParName(0,"Y_NS_1");
  fit->SetParName(1,"Sigma_1_NS");
  fit->SetParName(2,"Y_NS_2");
  // NOTE: sigma_2_NS is the difference between sigma_2 and sigma_1
  fit->SetParName(3,"Sigma_2_NS");
  fit->SetParName(4,"Y_AS");
  fit->SetParName(5,"Sigma_AS");
  fit->SetParName(6,"Beta_AS");
  fit->SetParName(7,"B");

  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  double sigma_ns_g = 0.8;
  double sigma_ns_1_g = 0.3;
  double sigma_ns_2_g = 0.3;
  double sigma_as_g = 0.8;
  double beta_as_g = 1.5;
  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();

  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }


  double near_yield_max,away_yield_max;
  hist->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  //  near_yield_max = hist->Integral("width") - 0.9 * B_g * PI; //0.8 in case bad B_g
  near_yield_max = hist->Integral("width") - 1.0 * B_g * PI; //0.8 in case bad B_g
  hist->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  //  away_yield_max = hist->Integral("width") - 0.9 * B_g * PI;
  away_yield_max = hist->Integral("width") - 1.0 * B_g * PI;
  hist->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);

  printf("near_yield_max = %f \t away_yield_max = %f\n",near_yield_max,away_yield_max);

  //Near side prefit analysis
  sigma_ns_g = TMath::Max(hist_clone->GetRMS(),0.2);

  double Y_ns_1_g,Y_ns_2_g;

  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);


  printf("Starting near side ...\n");

  // new idea: use separate fits: one for peak < R, on for peak > R
  // since root doesn't let me fit using two ranges, I will use each
  // range separately, and average

  // better idea: do the outer fit first
  // then subtract it, do the second fit

  // Trying to make this independent of jet resolution parameter
  //	double r_ring = 0.2;  // estimated area where ring effect 
  double r_ring = R; //FIXME This might break everything
  // from jet reconstruction affects near side



  printf("Wide peak pre-fit ...\n");
  smooth_func->SetParLimits(0,0,2*near_yield_max);
  //peak2 (leftside)
  smooth_func->SetParameter(1,r_ring);
  smooth_func->SetParLimits(1,0,PI);
  hist_clone->GetXaxis()->SetRangeUser(-PI/2,-r_ring/1.25);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = smooth_func->GetParameter(0);
  sigma_ns_2_g = smooth_func->GetParameter(1);

  //peak2 (rightside)
  smooth_func->SetParameter(1,r_ring);
  hist_clone->GetXaxis()->SetRangeUser(r_ring/1.25,PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = (Y_ns_2_g + smooth_func->GetParameter(0)) / 2.;
  sigma_ns_2_g = (sigma_ns_2_g + smooth_func->GetParameter(1)) / 2.;

  Y_ns_2_g = min(near_yield_max,Y_ns_2_g/5);

  //subtract guess for wider peak from hist_clone
  /*  smooth_func->SetParLimits(0,0,near_yield_max);
    smooth_func->SetParLimits(1,0,PI);
    smooth_func->SetParameter(0,0.3*Y_ns_2_g); //playing it safe
    smooth_func->SetParameter(1,sigma_ns_2_g);
    */
  TH1F *hist_clone_clone = (TH1F *) hist_clone->Clone();
  // hist_clone_clone->Add(smooth_func,-1.);

  printf("Sharp peak pre-fit ...\n");
  //peak 1

  //fixing the maximum point:  use unnormalized gaussian

  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],0)",-PI/2,PI/2);

  hist_clone_clone->GetXaxis()->SetRangeUser(-r_ring,r_ring);

  smooth_func->SetParameter(1,hist_clone_clone->GetRMS());


  // smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  // smooth_func->SetParameter(0,near_yield_max/3);
  double clone_clone_max_value = hist_clone_clone->GetBinContent(hist_clone_clone->GetMaximumBin());

  smooth_func->FixParameter(0,clone_clone_max_value);
  //	smooth_func->SetParameter(1,0.8*sigma_ns_2_g);
  smooth_func->SetParLimits(0,0,1.2*clone_clone_max_value);
  //	smooth_func->SetParLimits(1,0,1.2*sigma_ns_2_g);
  smooth_func->SetParLimits(1,0,2.2*sigma_ns_2_g);

  hist_clone_clone->Fit(smooth_func,"R0Q");
  //	Y_ns_1_g = smooth_func->GetParameter(0);
  sigma_ns_1_g = smooth_func->GetParameter(1);
  Y_ns_1_g = smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g;

  delete hist_clone_clone; // unless we want to use this for awayside?

  // Redefining sigma_ns_2_g to be the difference between the two 
  // sigmas, thus maintaining the ordering in the final function

  printf("Found Y_ns_1_g     = %f\t Y_ns_2_g     = %f\n",Y_ns_1_g,Y_ns_2_g);
  printf("Found sigma_ns_1_g = %f\t sigma_ns_2_g = %f\n",sigma_ns_1_g,sigma_ns_2_g);

  if (sigma_ns_2_g > sigma_ns_1_g) sigma_ns_2_g = sigma_ns_2_g - sigma_ns_1_g;
  else sigma_ns_2_g = r_ring;


  printf("Starting away side ...\n");

  // Away Side
  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  sigma_as_g = TMath::Max(hist_clone->GetRMS() ,0.1);
  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  delete smooth_func;
  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",+PI/2,3.0*PI/2);
  double yield_as = hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2));
  // double yield_as = abs(hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2)));
  smooth_func->SetParLimits(0,0,away_yield_max);
  smooth_func->SetParLimits(1,0,4.0*PI/3);

  smooth_func->SetParameter(0,0.9*yield_as);
  //	smooth_func->SetParameter(0,(max_val_as - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  smooth_func->SetParameter(1,sigma_as_g);
  // // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
  hist_clone->Fit(smooth_func,"R0Q");
  //	sigma_as_g = smooth_func->GetParameter(1);
  //double Y_as_g = smooth_func->Eval(PI) *TMath::Sqrt(2*PI)*sigma_as_g;
  double Y_as_g = smooth_func->GetParameter(0);
  //Y_as_g = yield_as;
  //FIXME

  delete hist_clone;


  fit->SetParLimits(0,0.,2.5*near_yield_max);
  fit->SetParLimits(1,0.,PI);
  fit->SetParLimits(2,0.,1.0*near_yield_max);
  fit->SetParLimits(3,0.,PI/4);
  fit->SetParLimits(4,0.,away_yield_max);
  fit->SetParLimits(5,0.,PI);
  fit->SetParLimits(6,0.1,4.);
  fit->SetParLimits(7,0,max_val);

  //set parameters 
  fit->SetParameter(0,Y_ns_1_g);
  fit->SetParameter(1,sigma_ns_1_g);
  fit->SetParameter(2,Y_ns_2_g);
  fit->SetParameter(3,sigma_ns_2_g);
  fit->SetParameter(4,Y_as_g);
  fit->SetParameter(5,sigma_as_g);
  fit->SetParameter(6,beta_as_g);

	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(4,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
    fit->FixParameter(2,0);
  }


  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(7,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(7,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(7,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }

  // Fit to histogram
  // R uses the range defined in the fit function
  // 0 ensures that the fit isn't drawn
  // Q ensures minimum printing
  // + adds the fit to function list to ensure that it is not deleted on the creation of a new fit

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");

  return fit;
}


/** 
 * Same method as above, but with the FWHM as a parameter instead of sigma
 */
TF1 * phiCorrFit_2Gaus_1GenGaus_Fwhm(TH1D * hist, std::string name)
{
   TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * (2.* [6] * pow(TMath::Log(2),1./[6]))/ ([5]*TMath::Gamma(1./[6])) * (TMath::Exp(-TMath::Log(2) * pow(2*abs(x-TMath::Pi())/([5]),[6])) + TMath::Exp(-TMath::Log(2) * pow(2*abs(x+TMath::Pi())/([5]),[6])) +  TMath::Exp(-TMath::Log(2) * pow(2*abs(x-3.*TMath::Pi())/([5]),[6]))) + [7]",-PI/2,3.0*PI/2);
  // FIXME this might break something
  // this uses an attempt to get the right shape from the ring caused by the jet 
  // algorithm
//  TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1)*((TMath::Abs(x) >= [8]) + (TMath::Abs(x) < [8])*(TMath::Erfc(TMath::Sqrt(TMath::Abs([8]*[8] - x*x))))) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * ([6] * pow(TMath::Log(2),1./[6]))/ ([5]*TMath::Gamma(1./[6])) * (TMath::Exp(-TMath::Log(2) * pow(2*abs(x-TMath::Pi())/([5]),[6])) + TMath::Exp(-TMath::Log(2) * pow(2*abs(x+TMath::Pi())/([5]),[6]))) + [7]",-PI/2,3.0*PI/2);
  //name parameters
  fit->SetParName(0,"Y_NS_1");
  fit->SetParName(1,"Sigma_1_NS");
  fit->SetParName(2,"Y_NS_2");
  // NOTE: sigma_2_NS is the difference between sigma_2 and sigma_1
  fit->SetParName(3,"Sigma_2_NS");
  fit->SetParName(4,"Y_AS");
  fit->SetParName(5,"Omega_AS");
  fit->SetParName(6,"Beta_AS");
  fit->SetParName(7,"B");

  fit->SetParName(8,"R");
  fit->FixParameter(8,R);


  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  double sigma_ns_g = 0.8;
  double sigma_ns_1_g = 0.3;
  double sigma_ns_2_g = 0.3;
  double sigma_as_g = 0.8;
  double omega_as_g = 1.9;
  double beta_as_g = 1.5;
  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();
  //Method 1 for B_g 
  //to get good RMS guess: remove constant background via ZYAM
  //  B_g = hist->GetBinContent(hist->GetMinimumBin());  

  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }

  //Method 2 fit pi/4,3pi/4 with v2, take minimum
//  hist_clone->GetXaxis()->SetRangeUser(0.6,2.2);
  //TF1 *v2_func = new TF1("v2_func","[0]*TMath::Cos(x - TMath::Pi()/2.0) + [1]",PI/4,3.0*PI/4);
 /*
  TF1 *v2_func = new TF1("v2_func","[0]*TMath::Cos(x - [2]) + [1]",PI/4,3.0*PI/4);
  v2_func->SetParLimits(2,0,PI);
  v2_func->SetParameter(0,-1);
  v2_func->SetParameter(1,1);
  v2_func->SetParameter(2,hist_clone->GetXaxis()->GetBinCenter(hist_clone->GetMinimumBin()));
  // v2_func->SetParameter(2,PI/2);
  hist_clone->Fit(v2_func,"R0Q WW");
  //  B_g = max((Double_t) v2_func->GetMinimum(),0.);// minimum value, but not below zero
*/





  double near_yield_max,away_yield_max;
  hist->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  near_yield_max = abs(2*hist->Integral("width") - 1.0 * B_g);
  hist->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  away_yield_max = abs(5*hist->Integral("width") - 1.0 * B_g);
  hist->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);

  printf("near_yield_max = %f \t away_yield_max = %f\n",near_yield_max,away_yield_max);

  //Near side prefit analysis
  sigma_ns_g = TMath::Max(hist_clone->GetRMS(),0.2);

  double Y_ns_1_g,Y_ns_2_g;

  //TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);
  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)+[2]",-PI/2,3.0*PI/2);


  printf("Starting near side ...\n");

  // new idea: use separate fits: one for peak < R, on for peak > R
  // since root doesn't let me fit using two ranges, I will use each
  // range separately, and average

  // better idea: do the outer fit first
  // then subtract it, do the second fit

  // Trying to make this independent of jet resolution parameter
  //	double r_ring = 0.2;  // estimated area where ring effect 
  double r_ring = R*1.2; //FIXME This might break everything
  // from jet reconstruction affects near side



  printf("Wide peak pre-fit ...\n");
  smooth_func->SetParLimits(0,0,2*near_yield_max/PI);
  //peak2 (leftside)
  smooth_func->SetParameter(1,r_ring);
  smooth_func->SetParLimits(1,0,PI);
  smooth_func->SetParameter(2,0.1*near_yield_max/PI);
  smooth_func->SetParLimits(2,0,2.1*near_yield_max/PI);

  hist_clone->GetXaxis()->SetRangeUser(-PI/2,-r_ring/1.25);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = smooth_func->GetParameter(0);
  sigma_ns_2_g = smooth_func->GetParameter(1);

  //peak2 (rightside)
  smooth_func->SetParameter(1,r_ring);
  hist_clone->GetXaxis()->SetRangeUser(r_ring/1.25,PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = (Y_ns_2_g + smooth_func->GetParameter(0)) / 2.;
  sigma_ns_2_g = (sigma_ns_2_g + smooth_func->GetParameter(1)) / 2.;

  //Y_ns_2_g = min(near_yield_max,Y_ns_2_g/5);
  Y_ns_2_g = min(near_yield_max,Y_ns_2_g/2.);
  if (sigma_ns_2_g > PI/2) sigma_ns_2_g = PI/4; // checking for bullshit

  //subtract guess for wider peak from hist_clone
  /*  smooth_func->SetParLimits(0,0,near_yield_max);
    smooth_func->SetParLimits(1,0,PI);
    smooth_func->SetParameter(0,0.3*Y_ns_2_g); //playing it safe
    smooth_func->SetParameter(1,sigma_ns_2_g);
    */
  TH1F *hist_clone_clone = (TH1F *) hist_clone->Clone();
  // hist_clone_clone->Add(smooth_func,-1.);

  printf("Sharp peak pre-fit ...\n");
  //peak 1

  //fixing the maximum point:  use unnormalized gaussian

  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],0)",-PI/2,PI/2);

  hist_clone_clone->GetXaxis()->SetRangeUser(-r_ring,r_ring);

  smooth_func->SetParameter(1,hist_clone_clone->GetRMS());


  // smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  // smooth_func->SetParameter(0,near_yield_max/3);
  double clone_clone_max_value = hist_clone_clone->GetBinContent(hist_clone_clone->GetMaximumBin());

  smooth_func->FixParameter(0,clone_clone_max_value);
  //	smooth_func->SetParameter(1,0.8*sigma_ns_2_g);
  smooth_func->SetParLimits(0,0,1.2*clone_clone_max_value);
  //	smooth_func->SetParLimits(1,0,1.2*sigma_ns_2_g);
  //	smooth_func->SetParLimits(1,0,2.2*sigma_ns_2_g);
  smooth_func->SetParLimits(1,0,1.5);

  hist_clone_clone->Fit(smooth_func,"R0Q");
  //	Y_ns_1_g = smooth_func->GetParameter(0);
  sigma_ns_1_g = max(smooth_func->GetParameter(1),0.001);
  Y_ns_1_g = min(smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g,near_yield_max);
  //Y_ns_1_g = 0.5* (smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g + near_yield_max);
  //	Y_ns_1_g = near_yield_max;



  delete hist_clone_clone; // unless we want to use this for awayside?

  // Redefining sigma_ns_2_g to be the difference between the two 
  // sigmas, thus maintaining the ordering in the final function

  printf("Found Y_ns_1_g     = %f\t Y_ns_2_g     = %f\n",Y_ns_1_g,Y_ns_2_g);
  printf("Found sigma_ns_1_g = %f\t sigma_ns_2_g = %f\n",sigma_ns_1_g,sigma_ns_2_g);

  if (sigma_ns_2_g > sigma_ns_1_g) sigma_ns_2_g = sigma_ns_2_g - sigma_ns_1_g;
  else sigma_ns_2_g = r_ring;


  printf("Starting away side ...\n");

  // Away Side
  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  sigma_as_g = TMath::Max(hist_clone->GetRMS() ,0.1);
  omega_as_g = sigma_as_g * 2.35;
//  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  delete smooth_func;
  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",+PI/2,3.0*PI/2);
  double yield_as = hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2),"width");
  // double yield_as = abs(hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2)));
  smooth_func->SetParLimits(0,0,away_yield_max);
  smooth_func->SetParLimits(1,0,4.0*PI/3);

  smooth_func->SetParameter(0,0.9*yield_as);
  //	smooth_func->SetParameter(0,(max_val_as - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  smooth_func->SetParameter(1,sigma_as_g);
  // // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
 // hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  //	sigma_as_g = smooth_func->GetParameter(1);
  //double Y_as_g = smooth_func->Eval(PI) *TMath::Sqrt(2*PI)*sigma_as_g;
  //  double Y_as_g = smooth_func->GetParameter(0);
 // double Y_as_g = 0.5 * (smooth_func->GetParameter(0) + yield_as);
  double Y_as_g = TMath::Min(smooth_func->GetParameter(0),yield_as);
  

  // Experiment: trying first step of wikipedia estimation
  // getting m1, the first statistical moment of the absolute values
  double m1,m1_1,m1_2,m2;
  m2 = hist_clone->GetStdDev();
  hist_clone->GetXaxis()->SetRangeUser(PI/2,PI);
  m1_1 = PI - hist_clone->GetMean();
  hist_clone->GetXaxis()->SetRangeUser(PI,3.0*PI/2);
  m1_2 = PI + hist_clone->GetMean();
  m1 = m1_1 + m1_2;
  printf("m1_1 = %f \tm1_2 = %f \tm1 = %f \t m2 = %f\n",m1_1,m1_2,m1,m2);

  beta_as_g = m1 / sqrt(m2) ;
  if (beta_as_g > 4.) beta_as_g = 2.0;

  delete hist_clone;

  fit->SetParLimits(0,0.,2.5*near_yield_max);
  fit->SetParLimits(1,0.,PI);
  fit->SetParLimits(2,0.,1.0*near_yield_max);
  fit->SetParLimits(3,0.,PI/4);
  fit->SetParLimits(4,0.,away_yield_max);
  fit->SetParLimits(5,0.,PI);
  fit->SetParLimits(6,0.1,4.);
  fit->SetParLimits(7,0,max_val);

  //set parameters 

  fit->SetParameter(0,Y_ns_1_g);
  fit->SetParameter(1,sigma_ns_1_g);
  fit->SetParameter(2,Y_ns_2_g);
  fit->SetParameter(3,sigma_ns_2_g);
  fit->SetParameter(4,Y_as_g);
  fit->SetParameter(5,omega_as_g);
  fit->SetParameter(6,beta_as_g);

	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(4,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
    fit->FixParameter(2,0);
  }


  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(7,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(7,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(7,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }
  // Fit to histogram
  // R uses the range defined in the fit function
  // 0 ensures that the fit isn't drawn
  // Q ensures minimum printing
  // + adds the fit to function list to ensure that it is not deleted on the creation of a new fit
  // M Improves fit results

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");
 




  return fit;
}




/**
 * 2 Gaus for near side, 2 Gaus for Awayside
 */
TF1 * phiCorrFit_2Gaus_2Gaus(TH1D * hist, std::string name)
{
 //TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * (2.*[6] * pow(TMath::Gamma(3./[6]),0.5))/ (2.*[5]*pow(TMath::Gamma(1./[6]),1.5)) * (TMath::Exp(-pow(abs(x-TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])) + TMath::Exp(-pow(abs(x+TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])) +TMath::Exp(-pow(abs(x-3*TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6]))  ) + [7]",-PI/2,3.0*PI/2);
  TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * (TMath::Gaus(x,-TMath::Pi(),[5],1) + TMath::Gaus(x,TMath::Pi(),[5],1) + TMath::Gaus(x,3.*TMath::Pi(),[5],1)) + [4]*[6]*(TMath::Gaus(x,-TMath::Pi(),[5]+[7],1) + TMath::Gaus(x,TMath::Pi(),[5]+[7],1) + TMath::Gaus(x,3*TMath::Pi(),[5]+[7],1))+ [8]",-PI/2,3.0*PI/2);

  //name parameters
  fit->SetParName(0,"Y_NS_1");
  fit->SetParName(1,"Sigma_1_NS");
  fit->SetParName(2,"Y_NS_2");
  // NOTE: sigma_2_NS is the difference between sigma_2 and sigma_1
  fit->SetParName(3,"Sigma_2_NS");

  fit->SetParName(4,"Y_AS_1");
	fit->SetParName(5,"Sigma_1_AS");

	// Possible next improvement to make: Y_AS_2 is the ratio between Y_2 and Y_1 in the awayside
  fit->SetParName(6,"Y_AS_2");
  // NOTE: sigma_2_NS is the difference between sigma_2 and sigma_1
  fit->SetParName(7,"Sigma_2_AS");
  fit->SetParName(8,"B");

  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  double sigma_ns_g = 0.8;
  double sigma_ns_1_g = 0.3;
  double sigma_ns_2_g = 0.3;
  double sigma_as_g = 0.8;
  double sigma_as_g_1, sigma_as_g_2;
  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();

  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }


  double near_yield_max,away_yield_max;
//  hist->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  //  near_yield_max = hist->Integral("width") - 0.9 * B_g * PI; //0.8 in case bad B_g
//  near_yield_max = hist->Integral("width") - 1.0 * B_g * PI; //0.8 in case bad B_g
  hist_clone->GetXaxis()->SetRangeUser(-PI/2,PI/2);
	near_yield_max = hist_clone->Integral("width");
//  hist->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  //  away_yield_max = hist->Integral("width") - 0.9 * B_g * PI;
//  away_yield_max = hist->Integral("width") - 1.0 * B_g * PI;
  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  away_yield_max = hist_clone->Integral("width");


  printf("near_yield_max = %f \t away_yield_max = %f\n",near_yield_max,away_yield_max);

  hist_clone->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);
///  hist->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);

  printf("Starting near side ...\n");

	// NearSide prefit analysis
  sigma_ns_g = TMath::Max(hist_clone->GetRMS(),0.2);

  double Y_ns_1_g,Y_ns_2_g;

  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);

  // Trying to make this independent of jet resolution parameter
  //	double r_ring = 0.2;  // estimated area where ring effect 
  double r_ring = R; //FIXME This might break everything
	r_ring = 0.4; 
  // from jet reconstruction affects near side

  printf("Wide peak pre-fit ...\n");
  smooth_func->SetParLimits(0,0,2*near_yield_max);
  //peak2 (leftside)
  smooth_func->SetParameter(1,r_ring);
  smooth_func->SetParLimits(1,0,PI);
  hist_clone->GetXaxis()->SetRangeUser(-PI/2,-r_ring);
 // hist_clone->GetXaxis()->SetRangeUser(-PI/2,-r_ring/1.25);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = smooth_func->GetParameter(0);
  sigma_ns_2_g = smooth_func->GetParameter(1);

  //peak2 (rightside)
  smooth_func->SetParameter(1,r_ring);
  hist_clone->GetXaxis()->SetRangeUser(r_ring,PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = (Y_ns_2_g + smooth_func->GetParameter(0)) / 2.;
  sigma_ns_2_g = (sigma_ns_2_g + smooth_func->GetParameter(1)) / 2.;

  Y_ns_2_g = min(near_yield_max,Y_ns_2_g/5);

  //subtract guess for wider peak from hist_clone
  /*  smooth_func->SetParLimits(0,0,near_yield_max);
    smooth_func->SetParLimits(1,0,PI);
    smooth_func->SetParameter(0,0.3*Y_ns_2_g); //playing it safe
    smooth_func->SetParameter(1,sigma_ns_2_g);
    */
  TH1F *hist_clone_clone = (TH1F *) hist_clone->Clone();
  // hist_clone_clone->Add(smooth_func,-1.);

  printf("Sharp peak pre-fit ...\n");
  //peak 1

  //fixing the maximum point:  use unnormalized gaussian

  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],0)",-PI/2,PI/2);

  hist_clone_clone->GetXaxis()->SetRangeUser(-r_ring,r_ring);

  smooth_func->SetParameter(1,hist_clone_clone->GetRMS());


  // smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  // smooth_func->SetParameter(0,near_yield_max/3);
  double clone_clone_max_value = hist_clone_clone->GetBinContent(hist_clone_clone->GetMaximumBin());

  smooth_func->FixParameter(0,clone_clone_max_value);
  //	smooth_func->SetParameter(1,0.8*sigma_ns_2_g);
  smooth_func->SetParLimits(0,0,1.2*clone_clone_max_value);
  //	smooth_func->SetParLimits(1,0,1.2*sigma_ns_2_g);
  smooth_func->SetParLimits(1,0,2.2*sigma_ns_2_g);

  hist_clone_clone->Fit(smooth_func,"R0Q");
  //	Y_ns_1_g = smooth_func->GetParameter(0);
  sigma_ns_1_g = smooth_func->GetParameter(1);
  Y_ns_1_g = smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g;

  delete hist_clone_clone; // unless we want to use this for awayside?

  // Redefining sigma_ns_2_g to be the difference between the two 
  // sigmas, thus maintaining the ordering in the final function

  printf("Found Y_ns_1_g     = %f\t Y_ns_2_g     = %f\n",Y_ns_1_g,Y_ns_2_g);
  printf("Found sigma_ns_1_g = %f\t sigma_ns_2_g = %f\n",sigma_ns_1_g,sigma_ns_2_g);

  if (sigma_ns_2_g > sigma_ns_1_g) sigma_ns_2_g = sigma_ns_2_g - sigma_ns_1_g;
  else sigma_ns_2_g = r_ring;


  printf("Starting away side ...\n");

  // Away Side
//  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
 // sigma_as_g = TMath::Max(hist_clone->GetRMS() ,0.1);
  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  sigma_as_g_1 = TMath::Max(hist_clone->GetRMS() ,0.1);
  delete smooth_func;

  sigma_as_g_1 = 0.8;
  sigma_as_g_2 = 0.4;


  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",+PI/2,3.0*PI/2);
  double yield_as = hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2));
  // double yield_as = abs(hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2)));
  smooth_func->SetParLimits(0,0,away_yield_max);
  smooth_func->SetParLimits(1,0,4.0*PI/3);

  smooth_func->SetParameter(0,0.9*yield_as);
  //	smooth_func->SetParameter(0,(max_val_as - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
 // smooth_func->SetParameter(1,sigma_as_g);
  smooth_func->SetParameter(1,sigma_as_g_1);
  // // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
  hist_clone->Fit(smooth_func,"R0Q");
  //	sigma_as_g = smooth_func->GetParameter(1);
  //double Y_as_g = smooth_func->Eval(PI) *TMath::Sqrt(2*PI)*sigma_as_g;
  double Y_as_g = smooth_func->GetParameter(0);
	
	double Y_as_g_1 = 0.8 * Y_as_g;
	double Y_as_g_2 = 0.2;
	
  //Y_as_g = yield_as;
  //FIXME

  delete hist_clone;


  fit->SetParLimits(0,0.,2.5*near_yield_max);
  fit->SetParLimits(1,0.,PI);
  fit->SetParLimits(2,0.,1.0*near_yield_max);
  fit->SetParLimits(3,0.,PI/4);
  fit->SetParLimits(4,0.,away_yield_max);
  fit->SetParLimits(5,0.,PI);
  fit->SetParLimits(6,0.,away_yield_max);
  fit->SetParLimits(7,0.,PI/3);
  fit->SetParLimits(8,0,max_val);

  //set parameters 
  fit->SetParameter(0,Y_ns_1_g);
  fit->SetParameter(1,sigma_ns_1_g);
  fit->SetParameter(2,Y_ns_2_g);
  fit->SetParameter(3,sigma_ns_2_g);
  fit->SetParameter(4,Y_as_g_1);
  fit->SetParameter(5,sigma_as_g_1);
  fit->SetParameter(6,Y_as_g_2);
  fit->SetParameter(7,sigma_as_g_2);
  fit->SetParameter(8,B_g);

	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(4,0);
		fit->FixParameter(6,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
    fit->FixParameter(2,0);
  }


  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(8,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(8,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(8,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }

  // Fit to histogram
  // R uses the range defined in the fit function
  // 0 ensures that the fit isn't drawn
  // Q ensures minimum printing
  // + adds the fit to function list to ensure that it is not deleted on the creation of a new fit

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");

  return fit;
}


/** 
 * Same method as above, but with the FWHM as a parameter instead of sigma
 */
TF1 * phiCorrFit_2Gaus_2Gaus_Fwhm(TH1D * hist, std::string name)
{
   TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * (2.* [6] * pow(TMath::Log(2),1./[6]))/ ([5]*TMath::Gamma(1./[6])) * (TMath::Exp(-TMath::Log(2) * pow(2*abs(x-TMath::Pi())/([5]),[6])) + TMath::Exp(-TMath::Log(2) * pow(2*abs(x+TMath::Pi())/([5]),[6])) +  TMath::Exp(-TMath::Log(2) * pow(2*abs(x-3.*TMath::Pi())/([5]),[6]))) + [7]",-PI/2,3.0*PI/2);
  // FIXME this might break something
  // this uses an attempt to get the right shape from the ring caused by the jet 
  // algorithm
//  TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1)*((TMath::Abs(x) >= [8]) + (TMath::Abs(x) < [8])*(TMath::Erfc(TMath::Sqrt(TMath::Abs([8]*[8] - x*x))))) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * ([6] * pow(TMath::Log(2),1./[6]))/ ([5]*TMath::Gamma(1./[6])) * (TMath::Exp(-TMath::Log(2) * pow(2*abs(x-TMath::Pi())/([5]),[6])) + TMath::Exp(-TMath::Log(2) * pow(2*abs(x+TMath::Pi())/([5]),[6]))) + [7]",-PI/2,3.0*PI/2);
  //name parameters
  fit->SetParName(0,"Y_NS_1");
  fit->SetParName(1,"Sigma_1_NS");
  fit->SetParName(2,"Y_NS_2");
  // NOTE: sigma_2_NS is the difference between sigma_2 and sigma_1
  fit->SetParName(3,"Sigma_2_NS");
  fit->SetParName(4,"Y_AS");
  fit->SetParName(5,"Omega_AS");
  fit->SetParName(6,"Beta_AS");
  fit->SetParName(7,"B");

  fit->SetParName(8,"R");
  fit->FixParameter(8,R);


  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  double sigma_ns_g = 0.8;
  double sigma_ns_1_g = 0.3;
  double sigma_ns_2_g = 0.3;
  double sigma_as_g = 0.8;
  double omega_as_g = 1.9;
  double beta_as_g = 1.5;
  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();
  //Method 1 for B_g 
  //to get good RMS guess: remove constant background via ZYAM
  //  B_g = hist->GetBinContent(hist->GetMinimumBin());  

  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }

  //Method 2 fit pi/4,3pi/4 with v2, take minimum
//  hist_clone->GetXaxis()->SetRangeUser(0.6,2.2);
  //TF1 *v2_func = new TF1("v2_func","[0]*TMath::Cos(x - TMath::Pi()/2.0) + [1]",PI/4,3.0*PI/4);
 /*
  TF1 *v2_func = new TF1("v2_func","[0]*TMath::Cos(x - [2]) + [1]",PI/4,3.0*PI/4);
  v2_func->SetParLimits(2,0,PI);
  v2_func->SetParameter(0,-1);
  v2_func->SetParameter(1,1);
  v2_func->SetParameter(2,hist_clone->GetXaxis()->GetBinCenter(hist_clone->GetMinimumBin()));
  // v2_func->SetParameter(2,PI/2);
  hist_clone->Fit(v2_func,"R0Q WW");
  //  B_g = max((Double_t) v2_func->GetMinimum(),0.);// minimum value, but not below zero
*/





  double near_yield_max,away_yield_max;
  hist->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  near_yield_max = abs(2*hist->Integral("width") - 1.0 * B_g);
  hist->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  away_yield_max = abs(5*hist->Integral("width") - 1.0 * B_g);
  hist->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);

  printf("near_yield_max = %f \t away_yield_max = %f\n",near_yield_max,away_yield_max);

  //Near side prefit analysis
  sigma_ns_g = TMath::Max(hist_clone->GetRMS(),0.2);

  double Y_ns_1_g,Y_ns_2_g;

  //TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);
  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)+[2]",-PI/2,3.0*PI/2);


  printf("Starting near side ...\n");

  // new idea: use separate fits: one for peak < R, on for peak > R
  // since root doesn't let me fit using two ranges, I will use each
  // range separately, and average

  // better idea: do the outer fit first
  // then subtract it, do the second fit

  // Trying to make this independent of jet resolution parameter
  //	double r_ring = 0.2;  // estimated area where ring effect 
  double r_ring = R*1.2; //FIXME This might break everything
  // from jet reconstruction affects near side



  printf("Wide peak pre-fit ...\n");
  smooth_func->SetParLimits(0,0,2*near_yield_max/PI);
  //peak2 (leftside)
  smooth_func->SetParameter(1,r_ring);
  smooth_func->SetParLimits(1,0,PI);
  smooth_func->SetParameter(2,0.1*near_yield_max/PI);
  smooth_func->SetParLimits(2,0,2.1*near_yield_max/PI);

  hist_clone->GetXaxis()->SetRangeUser(-PI/2,-r_ring/1.25);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = smooth_func->GetParameter(0);
  sigma_ns_2_g = smooth_func->GetParameter(1);

  //peak2 (rightside)
  smooth_func->SetParameter(1,r_ring);
  hist_clone->GetXaxis()->SetRangeUser(r_ring/1.25,PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = (Y_ns_2_g + smooth_func->GetParameter(0)) / 2.;
  sigma_ns_2_g = (sigma_ns_2_g + smooth_func->GetParameter(1)) / 2.;

  //Y_ns_2_g = min(near_yield_max,Y_ns_2_g/5);
  Y_ns_2_g = min(near_yield_max,Y_ns_2_g/2.);
  if (sigma_ns_2_g > PI/2) sigma_ns_2_g = PI/4; // checking for bullshit

  //subtract guess for wider peak from hist_clone
  /*  smooth_func->SetParLimits(0,0,near_yield_max);
    smooth_func->SetParLimits(1,0,PI);
    smooth_func->SetParameter(0,0.3*Y_ns_2_g); //playing it safe
    smooth_func->SetParameter(1,sigma_ns_2_g);
    */
  TH1F *hist_clone_clone = (TH1F *) hist_clone->Clone();
  // hist_clone_clone->Add(smooth_func,-1.);

  printf("Sharp peak pre-fit ...\n");
  //peak 1

  //fixing the maximum point:  use unnormalized gaussian

  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],0)",-PI/2,PI/2);

  hist_clone_clone->GetXaxis()->SetRangeUser(-r_ring,r_ring);

  smooth_func->SetParameter(1,hist_clone_clone->GetRMS());


  // smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  // smooth_func->SetParameter(0,near_yield_max/3);
  double clone_clone_max_value = hist_clone_clone->GetBinContent(hist_clone_clone->GetMaximumBin());

  smooth_func->FixParameter(0,clone_clone_max_value);
  //	smooth_func->SetParameter(1,0.8*sigma_ns_2_g);
  smooth_func->SetParLimits(0,0,1.2*clone_clone_max_value);
  //	smooth_func->SetParLimits(1,0,1.2*sigma_ns_2_g);
  //	smooth_func->SetParLimits(1,0,2.2*sigma_ns_2_g);
  smooth_func->SetParLimits(1,0,1.5);

  hist_clone_clone->Fit(smooth_func,"R0Q");
  //	Y_ns_1_g = smooth_func->GetParameter(0);
  sigma_ns_1_g = max(smooth_func->GetParameter(1),0.001);
  Y_ns_1_g = min(smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g,near_yield_max);
  //Y_ns_1_g = 0.5* (smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g + near_yield_max);
  //	Y_ns_1_g = near_yield_max;



  delete hist_clone_clone; // unless we want to use this for awayside?

  // Redefining sigma_ns_2_g to be the difference between the two 
  // sigmas, thus maintaining the ordering in the final function

  printf("Found Y_ns_1_g     = %f\t Y_ns_2_g     = %f\n",Y_ns_1_g,Y_ns_2_g);
  printf("Found sigma_ns_1_g = %f\t sigma_ns_2_g = %f\n",sigma_ns_1_g,sigma_ns_2_g);

  if (sigma_ns_2_g > sigma_ns_1_g) sigma_ns_2_g = sigma_ns_2_g - sigma_ns_1_g;
  else sigma_ns_2_g = r_ring;


  printf("Starting away side ...\n");

  // Away Side
  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  sigma_as_g = TMath::Max(hist_clone->GetRMS() ,0.1);
  omega_as_g = sigma_as_g * 2.35;
//  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  delete smooth_func;
  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",+PI/2,3.0*PI/2);
  double yield_as = hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2),"width");
  // double yield_as = abs(hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2)));
  smooth_func->SetParLimits(0,0,away_yield_max);
  smooth_func->SetParLimits(1,0,4.0*PI/3);

  smooth_func->SetParameter(0,0.9*yield_as);
  //	smooth_func->SetParameter(0,(max_val_as - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  smooth_func->SetParameter(1,sigma_as_g);
  // // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
 // hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  //	sigma_as_g = smooth_func->GetParameter(1);
  //double Y_as_g = smooth_func->Eval(PI) *TMath::Sqrt(2*PI)*sigma_as_g;
  //  double Y_as_g = smooth_func->GetParameter(0);
 // double Y_as_g = 0.5 * (smooth_func->GetParameter(0) + yield_as);
  double Y_as_g = TMath::Min(smooth_func->GetParameter(0),yield_as);
  

  // Experiment: trying first step of wikipedia estimation
  // getting m1, the first statistical moment of the absolute values
  double m1,m1_1,m1_2,m2;
  m2 = hist_clone->GetStdDev();
  hist_clone->GetXaxis()->SetRangeUser(PI/2,PI);
  m1_1 = PI - hist_clone->GetMean();
  hist_clone->GetXaxis()->SetRangeUser(PI,3.0*PI/2);
  m1_2 = PI + hist_clone->GetMean();
  m1 = m1_1 + m1_2;
  printf("m1_1 = %f \tm1_2 = %f \tm1 = %f \t m2 = %f\n",m1_1,m1_2,m1,m2);

  beta_as_g = m1 / sqrt(m2) ;
  if (beta_as_g > 4.) beta_as_g = 2.0;

  delete hist_clone;

  fit->SetParLimits(0,0.,2.5*near_yield_max);
  fit->SetParLimits(1,0.,PI);
  fit->SetParLimits(2,0.,1.0*near_yield_max);
  fit->SetParLimits(3,0.,PI/4);
  fit->SetParLimits(4,0.,away_yield_max);
  fit->SetParLimits(5,0.,PI);
  fit->SetParLimits(6,0.1,4.);
  fit->SetParLimits(7,0,max_val);

  //set parameters 

  fit->SetParameter(0,Y_ns_1_g);
  fit->SetParameter(1,sigma_ns_1_g);
  fit->SetParameter(2,Y_ns_2_g);
  fit->SetParameter(3,sigma_ns_2_g);
  fit->SetParameter(4,Y_as_g);
  fit->SetParameter(5,omega_as_g);
  fit->SetParameter(6,beta_as_g);

	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(4,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
    fit->FixParameter(2,0);
  }


  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(7,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(7,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(7,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }
  // Fit to histogram
  // R uses the range defined in the fit function
  // 0 ensures that the fit isn't drawn
  // Q ensures minimum printing
  // + adds the fit to function list to ensure that it is not deleted on the creation of a new fit
  // M Improves fit results

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");
 




  return fit;
}



/*  1 Generalized Gaussian for each peak
 *
 */

TF1 * phiCorrFit_1GenGaus_1GenGaus(TH1D * hist, std::string name)
{
  TF1 * fit = new TF1(name.c_str(),"[0] * (2.*[2] * pow(TMath::Gamma(3./[2]),0.5))/ (2.*[1]*pow(TMath::Gamma(1./[2]),1.5)) * (TMath::Exp(-pow(abs(x+2.*TMath::Pi())/([1]*pow(TMath::Gamma(1./[2])/TMath::Gamma(3./[2]),0.5)),[2])) + TMath::Exp(-pow(abs(x)/([1]*pow(TMath::Gamma(1./[2])/TMath::Gamma(3./[2]),0.5)),[2])) +TMath::Exp(-pow(abs(x-2.*TMath::Pi())/([1]*pow(TMath::Gamma(1./[2])/TMath::Gamma(3./[2]),0.5)),[2]))) + [3] * (2.*[5] * pow(TMath::Gamma(3./[5]),0.5))/ (2.*[4]*pow(TMath::Gamma(1./[5]),1.5)) * (TMath::Exp(-pow(abs(x-TMath::Pi())/([4]*pow(TMath::Gamma(1./[5])/TMath::Gamma(3./[5]),0.5)),[5])) + TMath::Exp(-pow(abs(x+TMath::Pi())/([4]*pow(TMath::Gamma(1./[5])/TMath::Gamma(3./[5]),0.5)),[5])) +TMath::Exp(-pow(abs(x-3*TMath::Pi())/([4]*pow(TMath::Gamma(1./[5])/TMath::Gamma(3./[5]),0.5)),[5]))  ) + [6]",-PI/2,3.0*PI/2);

  //name parameters
  fit->SetParName(0,"Y_NS");
  fit->SetParName(1,"Sigma_NS");
  fit->SetParName(2,"Beta_NS");

  fit->SetParName(3,"Y_AS");
  fit->SetParName(4,"Sigma_AS");
  fit->SetParName(5,"Beta_AS");
  fit->SetParName(6,"B");

  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  double sigma_ns_g = 0.8;
	double beta_ns_g = 1.5;
  double sigma_as_g = 0.8;
  double beta_as_g = 1.5;
  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();

  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }


  double near_yield_max,away_yield_max;
  hist->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  //  near_yield_max = hist->Integral("width") - 0.9 * B_g * PI; //0.8 in case bad B_g
  near_yield_max = hist->Integral("width") - 1.0 * B_g * PI; //0.8 in case bad B_g
  hist->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  //  away_yield_max = hist->Integral("width") - 0.9 * B_g * PI;
  away_yield_max = hist->Integral("width") - 1.0 * B_g * PI;
  hist->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);

  printf("near_yield_max = %f \t away_yield_max = %f\n",near_yield_max,away_yield_max);

  //Near side prefit analysis
  sigma_ns_g = TMath::Max(hist_clone->GetRMS(),0.2);


  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);
	
	double Y_ns_g = 0.8 * near_yield_max;


  printf("Starting near side ...\n");

  // new idea: use separate fits: one for peak < R, on for peak > R
  // since root doesn't let me fit using two ranges, I will use each
  // range separately, and average

  // better idea: do the outer fit first
  // then subtract it, do the second fit

  // Trying to make this independent of jet resolution parameter
  //	double r_ring = 0.2;  // estimated area where ring effect 
//  double r_ring = R; //FIXME This might break everything
  // from jet reconstruction affects near side



/*  printf("Wide peak pre-fit ...\n");
  smooth_func->SetParLimits(0,0,2*near_yield_max);
  //peak2 (leftside)
  smooth_func->SetParameter(1,r_ring);
  smooth_func->SetParLimits(1,0,PI);
  hist_clone->GetXaxis()->SetRangeUser(-PI/2,-r_ring/1.25);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = smooth_func->GetParameter(0);
  sigma_ns_2_g = smooth_func->GetParameter(1);

  //peak2 (rightside)
  smooth_func->SetParameter(1,r_ring);
  hist_clone->GetXaxis()->SetRangeUser(r_ring/1.25,PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = (Y_ns_2_g + smooth_func->GetParameter(0)) / 2.;
  sigma_ns_2_g = (sigma_ns_2_g + smooth_func->GetParameter(1)) / 2.;

  Y_ns_2_g = min(near_yield_max,Y_ns_2_g/5);
*/
  //subtract guess for wider peak from hist_clone
  /*  smooth_func->SetParLimits(0,0,near_yield_max);
    smooth_func->SetParLimits(1,0,PI);
    smooth_func->SetParameter(0,0.3*Y_ns_2_g); //playing it safe
    smooth_func->SetParameter(1,sigma_ns_2_g);
    */
  TH1F *hist_clone_clone = (TH1F *) hist_clone->Clone();
  // hist_clone_clone->Add(smooth_func,-1.);

  printf("Sharp peak pre-fit ...\n");
  //peak 1

  //fixing the maximum point:  use unnormalized gaussian

/*  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],0)",-PI/2,PI/2);

  hist_clone_clone->GetXaxis()->SetRangeUser(-r_ring,r_ring);

  smooth_func->SetParameter(1,hist_clone_clone->GetRMS());


  // smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  // smooth_func->SetParameter(0,near_yield_max/3);
  double clone_clone_max_value = hist_clone_clone->GetBinContent(hist_clone_clone->GetMaximumBin());

  smooth_func->FixParameter(0,clone_clone_max_value);
  //	smooth_func->SetParameter(1,0.8*sigma_ns_2_g);
  smooth_func->SetParLimits(0,0,1.2*clone_clone_max_value);
  //	smooth_func->SetParLimits(1,0,1.2*sigma_ns_2_g);
  smooth_func->SetParLimits(1,0,2.2*sigma_ns_2_g);

  hist_clone_clone->Fit(smooth_func,"R0Q");
  //	Y_ns_1_g = smooth_func->GetParameter(0);
  sigma_ns_1_g = smooth_func->GetParameter(1);
  Y_ns_1_g = smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g;
*/
  delete hist_clone_clone; // unless we want to use this for awayside?



  printf("Starting away side ...\n");

  // Away Side
  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  sigma_as_g = TMath::Max(hist_clone->GetRMS() ,0.1);
  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  delete smooth_func;
  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",+PI/2,3.0*PI/2);
  double yield_as = hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2));
  // double yield_as = abs(hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2)));
  smooth_func->SetParLimits(0,0,away_yield_max);
  smooth_func->SetParLimits(1,0,4.0*PI/3);

  smooth_func->SetParameter(0,0.9*yield_as);
  //	smooth_func->SetParameter(0,(max_val_as - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  smooth_func->SetParameter(1,sigma_as_g);
  // // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
  hist_clone->Fit(smooth_func,"R0Q");
  //	sigma_as_g = smooth_func->GetParameter(1);
  //double Y_as_g = smooth_func->Eval(PI) *TMath::Sqrt(2*PI)*sigma_as_g;
  double Y_as_g = smooth_func->GetParameter(0);
  //Y_as_g = yield_as;
  //FIXME

  delete hist_clone;


  fit->SetParLimits(0,0.,2.5*near_yield_max);
  fit->SetParLimits(1,0.,PI);
  fit->SetParLimits(2,0.,4.);
  fit->SetParLimits(3,0.,away_yield_max);
  fit->SetParLimits(4,0.,PI);
  fit->SetParLimits(5,0.1,4.);
  fit->SetParLimits(6,0,max_val);

  //set parameters 
  fit->SetParameter(0,Y_ns_g);
  fit->SetParameter(1,sigma_ns_g);
  fit->SetParameter(2,beta_ns_g);
  fit->SetParameter(3,Y_as_g);
  fit->SetParameter(4,sigma_as_g);
  fit->SetParameter(5,beta_as_g);

	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(3,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
  }


  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(6,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(6,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(6,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }

  // Fit to histogram
  // R uses the range defined in the fit function
  // 0 ensures that the fit isn't drawn
  // Q ensures minimum printing
  // + adds the fit to function list to ensure that it is not deleted on the creation of a new fit

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");

  return fit;
}


/** 
 * Same method as above, but with the FWHM as a parameter instead of sigma
 *
 */
TF1 * phiCorrFit_1GenGaus_1GenGaus_Fwhm(TH1D * hist, std::string name)
{
  TF1 * fit = new TF1(name.c_str(),"[0] * (2.*[2] * pow(TMath::Gamma(3./[2]),0.5))/ (2.*[1]*pow(TMath::Gamma(1./[2]),1.5)) * (TMath::Exp(-pow(abs(x+2.*TMath::Pi())/([1]*pow(TMath::Gamma(1./[2])/TMath::Gamma(3./[2]),0.5)),[2])) + TMath::Exp(-pow(abs(x)/([1]*pow(TMath::Gamma(1./[2])/TMath::Gamma(3./[2]),0.5)),[2])) +TMath::Exp(-pow(abs(x-2.*TMath::Pi())/([1]*pow(TMath::Gamma(1./[2])/TMath::Gamma(3./[2]),0.5)),[2]))) + [3] * (2.* [5] * pow(TMath::Log(2),1./[5]))/ ([4]*TMath::Gamma(1./[5])) * (TMath::Exp(-TMath::Log(2) * pow(2*abs(x-TMath::Pi())/([4]),[5])) + TMath::Exp(-TMath::Log(2) * pow(2*abs(x+TMath::Pi())/([4]),[5])) +  TMath::Exp(-TMath::Log(2) * pow(2*abs(x-3.*TMath::Pi())/([4]),[5]))) + [6]",-PI/2,3.0*PI/2);
  // TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * (2.* [6] * pow(TMath::Log(2),1./[6]))/ ([5]*TMath::Gamma(1./[6])) * (TMath::Exp(-TMath::Log(2) * pow(2*abs(x-TMath::Pi())/([5]),[6])) + TMath::Exp(-TMath::Log(2) * pow(2*abs(x+TMath::Pi())/([5]),[6])) +  TMath::Exp(-TMath::Log(2) * pow(2*abs(x-3.*TMath::Pi())/([5]),[6]))) + [7]",-PI/2,3.0*PI/2);
  // FIXME this might break something
  // this uses an attempt to get the right shape from the ring caused by the jet 
  // algorithm
//  TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1)*((TMath::Abs(x) >= [8]) + (TMath::Abs(x) < [8])*(TMath::Erfc(TMath::Sqrt(TMath::Abs([8]*[8] - x*x))))) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * ([6] * pow(TMath::Log(2),1./[6]))/ ([5]*TMath::Gamma(1./[6])) * (TMath::Exp(-TMath::Log(2) * pow(2*abs(x-TMath::Pi())/([5]),[6])) + TMath::Exp(-TMath::Log(2) * pow(2*abs(x+TMath::Pi())/([5]),[6]))) + [7]",-PI/2,3.0*PI/2);
  //name parameters
  fit->SetParName(0,"Y_NS");
  fit->SetParName(1,"Sigma_NS");
  fit->SetParName(2,"Beta_NS");
  fit->SetParName(3,"Y_AS");
  fit->SetParName(4,"Omega_AS");
  fit->SetParName(5,"Beta_AS");
  fit->SetParName(6,"B");

//  fit->SetParName(8,"R");
//  fit->FixParameter(8,R);


  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  double min_val = hist->GetBinContent(hist->GetMinimumBin());
  double sigma_ns_g = 0.8;
	double beta_ns_g = 1.5;
  double sigma_as_g = 0.8;
  double omega_as_g = 1.9;
  double beta_as_g = 1.5;
  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();
  //Method 1 for B_g 
  //to get good RMS guess: remove constant background via ZYAM
  //  B_g = hist->GetBinContent(hist->GetMinimumBin());  

  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }

  //Method 2 fit pi/4,3pi/4 with v2, take minimum
//  hist_clone->GetXaxis()->SetRangeUser(0.6,2.2);
  //TF1 *v2_func = new TF1("v2_func","[0]*TMath::Cos(x - TMath::Pi()/2.0) + [1]",PI/4,3.0*PI/4);
 /*
  TF1 *v2_func = new TF1("v2_func","[0]*TMath::Cos(x - [2]) + [1]",PI/4,3.0*PI/4);
  v2_func->SetParLimits(2,0,PI);
  v2_func->SetParameter(0,-1);
  v2_func->SetParameter(1,1);
  v2_func->SetParameter(2,hist_clone->GetXaxis()->GetBinCenter(hist_clone->GetMinimumBin()));
  // v2_func->SetParameter(2,PI/2);
  hist_clone->Fit(v2_func,"R0Q WW");
  //  B_g = max((Double_t) v2_func->GetMinimum(),0.);// minimum value, but not below zero
*/





  double near_yield_max,away_yield_max,Y_ns_g;
  hist->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  near_yield_max = abs(2*hist->Integral("width") - 1.0 * B_g);
	//Y_ns_g = abs(hist->Integral("width") - 1.0 * B_g);
		
  hist->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  away_yield_max = abs(5*hist->Integral("width") - 1.0 * B_g);
  hist->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);

  printf("near_yield_max = %f \t away_yield_max = %f\n",near_yield_max,away_yield_max);

  printf("Starting near side ...\n");
  //Near side prefit analysis
	
	hist_clone->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  sigma_ns_g = TMath::Max(hist_clone->GetRMS(),0.01);
	Y_ns_g = hist_clone->Integral("width");

  double m1,m1_1,m1_2,m2;
  m2 = hist_clone->GetStdDev();
  hist_clone->GetXaxis()->SetRangeUser(-PI/2,0);
  m1_1 = PI - hist_clone->GetMean();
  hist_clone->GetXaxis()->SetRangeUser(0,PI/2);
  m1_2 = PI + hist_clone->GetMean();
  m1 = m1_1 + m1_2;
  printf("NearSide: m1_1 = %f \tm1_2 = %f \tm1 = %f \t m2 = %f\n",m1_1,m1_2,m1,m2);

  beta_as_g = m1 / sqrt(m2) ;
  if (beta_as_g > 4.) beta_as_g = 2.0;



	hist_clone->GetXaxis()->SetRangeUser(-PI/2,3.0 * PI/2);

	//estimating beta

//  double Y_ns_1_g,Y_ns_2_g;

  //TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);
  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)+[2]",-PI/2,3.0*PI/2);



	


  // new idea: use separate fits: one for peak < R, on for peak > R
  // since root doesn't let me fit using two ranges, I will use each
  // range separately, and average

  // better idea: do the outer fit first
  // then subtract it, do the second fit

  // Trying to make this independent of jet resolution parameter
  //	double r_ring = 0.2;  // estimated area where ring effect 
//  double r_ring = R*1.2; //FIXME This might break everything
  // from jet reconstruction affects near side



  //subtract guess for wider peak from hist_clone
  /*  smooth_func->SetParLimits(0,0,near_yield_max);
    smooth_func->SetParLimits(1,0,PI);
    smooth_func->SetParameter(0,0.3*Y_ns_2_g); //playing it safe
    smooth_func->SetParameter(1,sigma_ns_2_g);
    */
//  TH1F *hist_clone_clone = (TH1F *) hist_clone->Clone();
  // hist_clone_clone->Add(smooth_func,-1.);

  //peak 1

  //fixing the maximum point:  use unnormalized gaussian

 // smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],0)",-PI/2,PI/2);

//  hist_clone_clone->GetXaxis()->SetRangeUser(-r_ring,r_ring);

//	sigma_ns_g = hist_clone_clone->GetRMS();

//  smooth_func->SetParameter(1,hist_clone_clone->GetRMS());


  // smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
//  // smooth_func->SetParameter(0,near_yield_max/3);
//  double clone_clone_max_value = hist_clone_clone->GetBinContent(hist_clone_clone->GetMaximumBin());

//	smooth_func->SetParameter(clone_clone_max_value);
//	smooth_func->FixParameter(0,clone_clone_max_value);
//  smooth_func->SetParameter(0,clone_clone_max_value);
  //	smooth_func->SetParameter(1,0.8*sigma_ns_2_g);
//  smooth_func->SetParLimits(0,0,1.2*clone_clone_max_value);
  //	smooth_func->SetParLimits(1,0,1.2*sigma_ns_2_g);
  //	smooth_func->SetParLimits(1,0,2.2*sigma_ns_2_g);
  //smooth_func->SetParLimits(1,0,1.5);

//  hist_clone_clone->Fit(smooth_func,"R0Q");
  //	Y_ns_1_g = smooth_func->GetParameter(0);
//  sigma_ns_g = max(smooth_func->GetParameter(1),0.001);
//  Y_ns_1_g = min(smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g,near_yield_max);
  //Y_ns_1_g = 0.5* (smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g + near_yield_max);
  //	Y_ns_1_g = near_yield_max;

//	sigma_ns_g = smooth_func->GetParameter(1);


//  delete hist_clone_clone; // unless we want to use this for awayside?

  // Redefining sigma_ns_2_g to be the difference between the two 
  // sigmas, thus maintaining the ordering in the final function



  printf("Starting away side ...\n");

  // Away Side
  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  sigma_as_g = TMath::Max(hist_clone->GetRMS() ,0.1);
  omega_as_g = sigma_as_g * 2.35;
//  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  delete smooth_func;
  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",+PI/2,3.0*PI/2);
  double yield_as = hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2),"width");
  // double yield_as = abs(hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2)));
  smooth_func->SetParLimits(0,0,away_yield_max);
  smooth_func->SetParLimits(1,0,4.0*PI/3);

  smooth_func->SetParameter(0,0.9*yield_as);
  //	smooth_func->SetParameter(0,(max_val_as - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  smooth_func->SetParameter(1,sigma_as_g);
  // // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
 // hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  //	sigma_as_g = smooth_func->GetParameter(1);
  //double Y_as_g = smooth_func->Eval(PI) *TMath::Sqrt(2*PI)*sigma_as_g;
  //  double Y_as_g = smooth_func->GetParameter(0);
 // double Y_as_g = 0.5 * (smooth_func->GetParameter(0) + yield_as);
  double Y_as_g = TMath::Min(smooth_func->GetParameter(0),yield_as);
  

  // Experiment: trying first step of wikipedia estimation
  // getting m1, the first statistical moment of the absolute values
  m1,m1_1,m1_2,m2;
  m2 = hist_clone->GetStdDev();
  hist_clone->GetXaxis()->SetRangeUser(PI/2,PI);
  m1_1 = PI - hist_clone->GetMean();
  hist_clone->GetXaxis()->SetRangeUser(PI,3.0*PI/2);
  m1_2 = PI + hist_clone->GetMean();
  m1 = m1_1 + m1_2;
  printf("AwaySide m1_1 = %f \tm1_2 = %f \tm1 = %f \t m2 = %f\n",m1_1,m1_2,m1,m2);

  beta_as_g = m1 / sqrt(m2) ;
  if (beta_as_g > 4.) beta_as_g = 2.0;

  delete hist_clone;

  fit->SetParLimits(0,0.,2.5*near_yield_max);
  fit->SetParLimits(1,0.,PI);
  fit->SetParLimits(2,0.,4.);
  fit->SetParLimits(3,0.,away_yield_max);
  fit->SetParLimits(4,0.,PI);
  fit->SetParLimits(5,0.1,4.);
  fit->SetParLimits(6,0,1.1*min_val);

  //set parameters 

  fit->SetParameter(0,Y_ns_g);
  fit->SetParameter(1,sigma_ns_g);
  fit->SetParameter(2,beta_ns_g);
  fit->SetParameter(3,Y_as_g);
  fit->SetParameter(4,omega_as_g);
  fit->SetParameter(5,beta_as_g);

	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(3,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
  }


  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(6,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(6,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(6,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }
  // Fit to histogram
  // R uses the range defined in the fit function
  // 0 ensures that the fit isn't drawn
  // Q ensures minimum printing
  // + adds the fit to function list to ensure that it is not deleted on the creation of a new fit
  // M Improves fit results

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "RQ0+M");
  printParameters(fit,true,"Final Parameters");
 




  return fit;
}







/**
  * 2 Gaussians for Nearside, 1 Modified Gaussian for awayside
  * G(X) = A(1+(1/(2n))((x/sigma)/2)^2)^(-n)
  * see https://arxiv.org/pdf/1509.04732v2.pdf
	* is that link right? https://arxiv.org/pdf/0904.1733.pdf
  * Standard deviation sigma as width parameter
  */
TF1 * phiCorrFit_2Gaus_1ModGaus(TH1D * hist, std::string name)
{
/*
  TString functionString = "[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1))"; 
  functionString += "+ [4] * (TMath::Power(1 + TMath::Power((x+TMath::Pi())/[5],2)/(2*[6]),-[6]) + TMath::Power(1 + TMath::Power((x-TMath::Pi())/[5],2)/(2*[6]),-[6])";
  functionString += " + TMath::Power(1 + TMath::Power((x-3*TMath::Pi())/[5],2)/(2*[6]),-[6])) + [7]";
*/
  TString functionString = "[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1))"; 
  functionString += "+ [4] * (1. / ([5] *TMath::Sqrt(2 * TMath::Pi() * [6]))) * (TMath::Gamma([6])/TMath::Gamma([6]-0.5)) * (TMath::Power(1 + TMath::Power((x+TMath::Pi())/[5],2)/(2*[6]),-[6]) + TMath::Power(1 + TMath::Power((x-TMath::Pi())/[5],2)/(2*[6]),-[6])";
  functionString += " + TMath::Power(1 + TMath::Power((x-3*TMath::Pi())/[5],2)/(2*[6]),-[6])) + [7]";

/*
  TString functionString = "[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1))";
  functionString += "+ [4] * (2*TMath::Sqrt(TMath::Power(2,1/[6])-1)/(TMath::Sqrt(TMath::Pi())*[5])) * (TMath::Gamma([6])/TMath::Gamma([6]-0.5)) * (TMath::Power(1 + (TMath::Power(2,1/[6])-1)*TMath::Power(2*(x+TMath::Pi())/[5],2)/(2*[6]),-[6]) + TMath::Power(1 + (TMath::Power(2,1/[6])-1)*TMath::Power(2*(x-TMath::Pi())/[5],2)/(2*[6]),-[6])";
  functionString += " + TMath::Power(1 + (TMath::Power(2,1/[6])-1)*TMath::Power(2*(x-3*TMath::Pi())/[5],2)/(2*[6]),-[6])) + [7]";
*/

/*

  TString functionString = "[0] * (TMath::Gamma([2])/(TMath::Sqrt(2*TMath::Pi()*([2]-1.5))*[1]*TMath::Gamma([2]-0.5))) *(TMath::Power(1 + TMath::Power((x+2.*TMath::Pi())/[1],2)/(2*([2]-1.5)),-[2]) + TMath::Power(1 + TMath::Power((x)/[1],2)/(2*([2]-1.5)),-[2])";
  functionString += " + TMath::Power(1 + TMath::Power((x-2*TMath::Pi())/[1],2)/(2*([2]-1.5)),-[2]))";
  functionString += "+ [3] *  (TMath::Gamma([2])/(TMath::Sqrt(2*TMath::Pi()*([5]-1.5))*[4]*TMath::Gamma([5]-0.5))) * (TMath::Power(1 + TMath::Power((x+TMath::Pi())/[4],2)/(2*([5]-1.5)),-[5]) + TMath::Power(1 + TMath::Power((x-TMath::Pi())/[4],2)/(2*([5]-1.5)),-[5])";
  functionString += " + TMath::Power(1 + TMath::Power((x-3*TMath::Pi())/[4],2)/(2*([5]-1.5)),-[5])) + [6]";




*/



  TF1 * fit = new TF1(name.c_str(),functionString,-PI/2,3.0*PI/2);
  
  // FIXME treating sigma as omega for now
  // treating n as beta for now

  //name parameters
  fit->SetParName(0,"Y_NS_1");
  fit->SetParName(1,"Sigma_1_NS");
  fit->SetParName(2,"Y_NS_2");
  // NOTE: sigma_2_NS is the difference between sigma_2 and sigma_1
  fit->SetParName(3,"Sigma_2_NS");
  fit->SetParName(4,"Y_AS");
  fit->SetParName(5,"Sigma_AS");
  fit->SetParName(6,"Beta_AS");
  fit->SetParName(7,"B");

//  fit->SetParName(8,"R");
//  fit->FixParameter(8,R);


  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  double sigma_ns_g = 0.8;
  double sigma_ns_1_g = 0.3;
  double sigma_ns_2_g = 0.3;
  double sigma_as_g = 0.8;
  double omega_as_g = 1.9;
  double beta_as_g = 1.5;
  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();
  //Method 1 for B_g 
  //to get good RMS guess: remove constant background via ZYAM
  //  B_g = hist->GetBinContent(hist->GetMinimumBin());  

  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }

  //Method 2 fit pi/4,3pi/4 with v2, take minimum
  hist_clone->GetXaxis()->SetRangeUser(0.6,2.2);
  //TF1 *v2_func = new TF1("v2_func","[0]*TMath::Cos(x - TMath::Pi()/2.0) + [1]",PI/4,3.0*PI/4);
 /*
  TF1 *v2_func = new TF1("v2_func","[0]*TMath::Cos(x - [2]) + [1]",PI/4,3.0*PI/4);
  v2_func->SetParLimits(2,0,PI);
  v2_func->SetParameter(0,-1);
  v2_func->SetParameter(1,1);
  v2_func->SetParameter(2,hist_clone->GetXaxis()->GetBinCenter(hist_clone->GetMinimumBin()));
  // v2_func->SetParameter(2,PI/2);
  hist_clone->Fit(v2_func,"R0Q WW");
  //  B_g = max((Double_t) v2_func->GetMinimum(),0.);// minimum value, but not below zero
*/


  double near_yield_max,away_yield_max;
  hist->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  //  near_yield_max = hist->Integral("width") - 0.9 * B_g * PI; //0.8 in case bad B_g
  near_yield_max = abs(1.2*hist->Integral("width") - 1.0 * B_g); //0.8 in case bad B_g
  hist->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  //  away_yield_max = hist->Integral("width") - 0.9 * B_g * PI;
  away_yield_max = abs(1.3*hist->Integral("width") - 0.9 * B_g);
  hist->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);



  printf("near_yield_max = %f \t away_yield_max = %f\n",near_yield_max,away_yield_max);

  //Near side prefit analysis
  sigma_ns_g = TMath::Max(hist_clone->GetRMS(),0.2);

  double Y_ns_1_g,Y_ns_2_g;

  //TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);
  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)+[2]",-PI/2,3.0*PI/2);


  printf("Starting near side ...\n");

  // new idea: use separate fits: one for peak < R, on for peak > R
  // since root doesn't let me fit using two ranges, I will use each
  // range separately, and average

  // better idea: do the outer fit first
  // then subtract it, do the second fit

  // Trying to make this independent of jet resolution parameter
  //	double r_ring = 0.2;  // estimated area where ring effect 
  double r_ring = R*1.2; //FIXME This might break everything
  // from jet reconstruction affects near side



  printf("Wide peak pre-fit ...\n");
  smooth_func->SetParLimits(0,0,2*near_yield_max/PI);
  //peak2 (leftside)
  smooth_func->SetParameter(1,r_ring);
  smooth_func->SetParLimits(1,0,PI);
  smooth_func->SetParameter(2,0.1*near_yield_max/PI);
  smooth_func->SetParLimits(2,0,2.1*near_yield_max/PI);

  hist_clone->GetXaxis()->SetRangeUser(-PI/2,-r_ring/1.25);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = smooth_func->GetParameter(0);
  sigma_ns_2_g = smooth_func->GetParameter(1);

  //peak2 (rightside)
  smooth_func->SetParameter(1,r_ring);
  hist_clone->GetXaxis()->SetRangeUser(r_ring/1.25,PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = (Y_ns_2_g + smooth_func->GetParameter(0)) / 2.;
  sigma_ns_2_g = (sigma_ns_2_g + smooth_func->GetParameter(1)) / 2.;

  //Y_ns_2_g = min(near_yield_max,Y_ns_2_g/5);
  Y_ns_2_g = min(near_yield_max,Y_ns_2_g/2.);
  if (sigma_ns_2_g > PI/2) sigma_ns_2_g = PI/4; // checking for bullshit

  //subtract guess for wider peak from hist_clone
  /*  smooth_func->SetParLimits(0,0,near_yield_max);
    smooth_func->SetParLimits(1,0,PI);
    smooth_func->SetParameter(0,0.3*Y_ns_2_g); //playing it safe
    smooth_func->SetParameter(1,sigma_ns_2_g);
    */
  TH1F *hist_clone_clone = (TH1F *) hist_clone->Clone();
  // hist_clone_clone->Add(smooth_func,-1.);

  printf("Sharp peak pre-fit ...\n");
  //peak 1

  //fixing the maximum point:  use unnormalized gaussian

  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],0)",-PI/2,PI/2);

  hist_clone_clone->GetXaxis()->SetRangeUser(-r_ring,r_ring);

  smooth_func->SetParameter(1,hist_clone_clone->GetRMS());


  // smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  // smooth_func->SetParameter(0,near_yield_max/3);
  double clone_clone_max_value = hist_clone_clone->GetBinContent(hist_clone_clone->GetMaximumBin());

  smooth_func->FixParameter(0,clone_clone_max_value);
  //	smooth_func->SetParameter(1,0.8*sigma_ns_2_g);
  smooth_func->SetParLimits(0,0,1.2*clone_clone_max_value);
  //	smooth_func->SetParLimits(1,0,1.2*sigma_ns_2_g);
  //	smooth_func->SetParLimits(1,0,2.2*sigma_ns_2_g);
  smooth_func->SetParLimits(1,0,1.5);

  hist_clone_clone->Fit(smooth_func,"R0Q");
  //	Y_ns_1_g = smooth_func->GetParameter(0);
  sigma_ns_1_g = max(smooth_func->GetParameter(1),0.001);
  Y_ns_1_g = min(smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g,near_yield_max);
  //Y_ns_1_g = 0.5* (smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g + near_yield_max);
  //	Y_ns_1_g = near_yield_max;



  delete hist_clone_clone; // unless we want to use this for awayside?

  // Redefining sigma_ns_2_g to be the difference between the two 
  // sigmas, thus maintaining the ordering in the final function

  printf("Found Y_ns_1_g     = %f\t Y_ns_2_g     = %f\n",Y_ns_1_g,Y_ns_2_g);
  printf("Found sigma_ns_1_g = %f\t sigma_ns_2_g = %f\n",sigma_ns_1_g,sigma_ns_2_g);

  if (sigma_ns_2_g > sigma_ns_1_g) sigma_ns_2_g = sigma_ns_2_g - sigma_ns_1_g;
  else sigma_ns_2_g = r_ring;


  printf("Starting away side ...\n");

  // Away Side
  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  sigma_as_g = TMath::Max(hist_clone->GetRMS() ,0.1);
  omega_as_g = sigma_as_g * 2.35;
//  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  delete smooth_func;
  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",+PI/2,3.0*PI/2);
  double yield_as = hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2),"width");
  // double yield_as = abs(hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2)));
  smooth_func->SetParLimits(0,0,away_yield_max);
  smooth_func->SetParLimits(1,0,4.0*PI/3);

  smooth_func->SetParameter(0,0.9*yield_as);
  //	smooth_func->SetParameter(0,(max_val_as - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  smooth_func->SetParameter(1,sigma_as_g);
  // // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
 // hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  //	sigma_as_g = smooth_func->GetParameter(1);
  //double Y_as_g = smooth_func->Eval(PI) *TMath::Sqrt(2*PI)*sigma_as_g;
//    double Y_as_g = smooth_func->GetParameter(0);
    double Y_as_g = yield_as;
 // double Y_as_g = 0.5 * (smooth_func->GetParameter(0) + yield_as);
 // double Y_as_g = TMath::Min(smooth_func->GetParameter(0),yield_as);
  	

  // Experiment: trying first step of wikipedia estimation
  // getting m1, the first statistical moment of the absolute values
/*  double m1,m1_1,m1_2,m2;
  m2 = hist_clone->GetStdDev();
  hist_clone->GetXaxis()->SetRangeUser(PI/2,PI);
  m1_1 = PI - hist_clone->GetMean();
  hist_clone->GetXaxis()->SetRangeUser(PI,3.0*PI/2);
  m1_2 = PI + hist_clone->GetMean();
  m1 = m1_1 + m1_2;
  printf("m1_1 = %f \tm1_2 = %f \tm1 = %f \t m2 = %f\n",m1_1,m1_2,m1,m2);

  beta_as_g = m1 / sqrt(m2) ;
  if (beta_as_g > 4.) beta_as_g = 2.0;
*/
  double beta_as_min  = 0.1;
  beta_as_g = 3;
  double beta_as_max = 20;

  delete hist_clone;

  fit->SetParLimits(0,0.,2.5*near_yield_max);
  fit->SetParLimits(1,0.,PI);
//  fit->SetParLimits(2,0.,2.5*near_yield_max);
  fit->SetParLimits(2,0.,1.0*near_yield_max);
  fit->SetParLimits(3,0.,PI/4);
  fit->SetParLimits(4,0.,away_yield_max);
  fit->SetParLimits(5,0.,PI);
  fit->SetParLimits(6,beta_as_min,beta_as_max);
  fit->SetParLimits(7,0,max_val);

  //set parameters 

  fit->SetParameter(0,Y_ns_1_g);
  fit->SetParameter(1,sigma_ns_1_g);
  fit->SetParameter(2,Y_ns_2_g);
  fit->SetParameter(3,sigma_ns_2_g);
  fit->SetParameter(4,Y_as_g);
  fit->SetParameter(5,sigma_as_g);
  fit->SetParameter(6,beta_as_g);

	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(4,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
    fit->FixParameter(2,0);
  }

  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(7,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(7,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(7,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }
  printf("B_g:          [%f \t %f \t %f]\n",0.,B_g,max_val);

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
//  if (!noFitting) hist->Fit(fit, "RQ0+ MM");
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");




  return fit;
}



/**
  * 2 Gaussians for Nearside, 1 Modified Gaussian for awayside
  * G(X) = A(1+(1/(2n))((x/sigma)/2)^2)^(-n)
  * see https://arxiv.org/pdf/1509.04732v2.pdf
  * FWHM omega as width parameter
  */
TF1 * phiCorrFit_2Gaus_1ModGaus_Fwhm(TH1D * hist, std::string name)
{
  TString functionString = "[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1))";
  functionString += "+ [4] * (2*TMath::Sqrt(TMath::Power(2,1/[6])-1)/(TMath::Sqrt(TMath::Pi())*[5])) *(TMath::Gamma([6])/TMath::Gamma([6]-0.5)) * (TMath::Power(1 + (TMath::Power(2,1/[6])-1)*TMath::Power(2*(x+TMath::Pi())/[5],2)/(2*[6]),-[6]) + TMath::Power(1 + (TMath::Power(2,1/[6])-1)*TMath::Power(2*(x-TMath::Pi())/[5],2)/(2*[6]),-[6])";
  functionString += " + TMath::Power(1 + (TMath::Power(2,1/[6])-1)*TMath::Power(2*(x-3*TMath::Pi())/[5],2)/(2*[6]),-[6])) + [7]";

  TF1 * fit = new TF1(name.c_str(),functionString,-PI/2,3.0*PI/2);
  
  // FIXME treating sigma as omega for now
  // treating n as beta for now

  //name parameters
  fit->SetParName(0,"Y_NS_1");
  fit->SetParName(1,"Sigma_1_NS");
  fit->SetParName(2,"Y_NS_2");
  // NOTE: sigma_2_NS is the difference between sigma_2 and sigma_1
  fit->SetParName(3,"Sigma_2_NS");
  fit->SetParName(4,"Y_AS");
  fit->SetParName(5,"Omega_AS");
 // fit->SetParName(6,"N_AS");
  fit->SetParName(6,"Beta_AS");
  fit->SetParName(7,"B");

//  fit->SetParName(8,"R");
//  fit->FixParameter(8,R);

  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  double sigma_ns_g = 0.8;
  double sigma_ns_1_g = 0.3;
  double sigma_ns_2_g = 0.3;
  double sigma_as_g = 0.8;
  double omega_as_g = 1.9;
//  double beta_as_g = 3.0;

  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();
  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }


  double near_yield_max,away_yield_max;
  hist->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  //  near_yield_max = hist->Integral("width") - 0.9 * B_g * PI; //0.8 in case bad B_g
  near_yield_max = abs(2*hist->Integral("width") - 1.0 * B_g); //0.8 in case bad B_g
  hist->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  //  away_yield_max = hist->Integral("width") - 0.9 * B_g * PI;
  away_yield_max = abs(5*hist->Integral("width") - 1.0 * B_g);
  hist->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);

  printf("near_yield_max = %f \t away_yield_max = %f\n",near_yield_max,away_yield_max);

  //Near side prefit analysis
  sigma_ns_g = TMath::Max(hist_clone->GetRMS(),0.2);
  double Y_ns_1_g,Y_ns_2_g;


  //TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);
  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)+[2]",-PI/2,3.0*PI/2);

  printf("Starting near side ...\n");

  double r_ring = R; //FIXME This might break everything
  // from jet reconstruction affects near side

  printf("Wide peak pre-fit ...\n");
  smooth_func->SetParLimits(0,0,2*near_yield_max);
 // peak2 (leftside)
  smooth_func->SetParameter(1,r_ring);
  smooth_func->SetParLimits(1,0,PI);
  hist_clone->GetXaxis()->SetRangeUser(-PI/2,-r_ring/1.25);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = smooth_func->GetParameter(0);
  sigma_ns_2_g = smooth_func->GetParameter(1);

  //peak2 (rightside)
  smooth_func->SetParameter(1,r_ring);
  hist_clone->GetXaxis()->SetRangeUser(r_ring/1.25,PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = (Y_ns_2_g + smooth_func->GetParameter(0)) / 2.;
  sigma_ns_2_g = (sigma_ns_2_g + smooth_func->GetParameter(1)) / 2.;

  Y_ns_2_g = min(near_yield_max,Y_ns_2_g/5);

  Y_ns_2_g = min(near_yield_max,Y_ns_2_g/2.);
  if (sigma_ns_2_g > PI/2) sigma_ns_2_g = PI/4; // checking for bullshit

  //subtract guess for wider peak from hist_clone
  /*  smooth_func->SetParLimits(0,0,near_yield_max);
    smooth_func->SetParLimits(1,0,PI);
    smooth_func->SetParameter(0,0.3*Y_ns_2_g); //playing it safe
    smooth_func->SetParameter(1,sigma_ns_2_g);
    */
  TH1F *hist_clone_clone = (TH1F *) hist_clone->Clone();
  // hist_clone_clone->Add(smooth_func,-1.);

  printf("Sharp peak pre-fit ...\n");
  //peak 1

  //fixing the maximum point:  use unnormalized gaussian

  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],0)",-PI/2,PI/2);

  hist_clone_clone->GetXaxis()->SetRangeUser(-r_ring,r_ring);

  smooth_func->SetParameter(1,hist_clone_clone->GetRMS());

  // smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  // smooth_func->SetParameter(0,near_yield_max/3);
  double clone_clone_max_value = hist_clone_clone->GetBinContent(hist_clone_clone->GetMaximumBin());

  smooth_func->FixParameter(0,clone_clone_max_value);
  //  smooth_func->SetParameter(1,0.8*sigma_ns_2_g);
  smooth_func->SetParLimits(0,0,1.2*clone_clone_max_value);
  //  smooth_func->SetParLimits(1,0,1.2*sigma_ns_2_g);
  //  smooth_func->SetParLimits(1,0,2.2*sigma_ns_2_g);
  smooth_func->SetParLimits(1,0,1.5);

  hist_clone_clone->Fit(smooth_func,"R0Q");
  //  Y_ns_1_g = smooth_func->GetParameter(0);
  sigma_ns_1_g = max(smooth_func->GetParameter(1),0.001);
  Y_ns_1_g = min(smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g,near_yield_max);
  //Y_ns_1_g = 0.5* (smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g + near_yield_max);
  //  Y_ns_1_g = near_yield_max;



  delete hist_clone_clone; // unless we want to use this for awayside?

  // Redefining sigma_ns_2_g to be the difference between the two 
  // sigmas, thus maintaining the ordering in the final function

  printf("Found Y_ns_1_g     = %f\t Y_ns_2_g     = %f\n",Y_ns_1_g,Y_ns_2_g);
  printf("Found sigma_ns_1_g = %f\t sigma_ns_2_g = %f\n",sigma_ns_1_g,sigma_ns_2_g);

  if (sigma_ns_2_g > sigma_ns_1_g) sigma_ns_2_g = sigma_ns_2_g - sigma_ns_1_g;
  else sigma_ns_2_g = r_ring;

  printf("Starting away side ...\n");

  // Away Side
  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  sigma_as_g = TMath::Max(hist_clone->GetRMS() ,0.1);
  omega_as_g = sigma_as_g * 2.35;
//  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  delete smooth_func;
  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",+PI/2,3.0*PI/2);
  double yield_as = hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2),"width");
  // double yield_as = abs(hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2)));
  smooth_func->SetParLimits(0,0,away_yield_max);
  smooth_func->SetParLimits(1,0,4.0*PI/3);

  smooth_func->SetParameter(0,0.9*yield_as);
  //  smooth_func->SetParameter(0,(max_val_as - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  smooth_func->SetParameter(1,sigma_as_g);
  // // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
 // hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  //  sigma_as_g = smooth_func->GetParameter(1);
  //double Y_as_g = smooth_func->Eval(PI) *TMath::Sqrt(2*PI)*sigma_as_g;
  //  double Y_as_g = smooth_func->GetParameter(0);
 // double Y_as_g = 0.5 * (smooth_func->GetParameter(0) + yield_as);
  double Y_as_g = TMath::Min(smooth_func->GetParameter(0),yield_as);

  double beta_as_min  = 0.1;
  double beta_as_g = 3;
  double beta_as_max = 20;
  

  delete hist_clone;

  fit->SetParLimits(0,0.,2.5*near_yield_max);
  fit->SetParLimits(1,0.,PI);
  fit->SetParLimits(2,0.,1.0*near_yield_max);
  fit->SetParLimits(3,0.,PI/4);
  fit->SetParLimits(4,0.,away_yield_max);
  fit->SetParLimits(5,0.,PI);
  fit->SetParLimits(6,beta_as_min,beta_as_max);
  fit->SetParLimits(7,0,max_val);
  
  //set parameters 
  fit->SetParameter(0,Y_ns_1_g);
  fit->SetParameter(1,sigma_ns_1_g);
  fit->SetParameter(2,Y_ns_2_g);
  fit->SetParameter(3,sigma_ns_2_g);
  fit->SetParameter(4,Y_as_g);
  fit->SetParameter(5,sigma_as_g);
  fit->SetParameter(6,beta_as_g);

	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(4,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
    fit->FixParameter(2,0);
  }

   switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(7,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(7,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(7,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "R0+");
  printParameters(fit,true,"Final Parameters");

  return fit;
}



/**
  * 1 Modified Gaussian for nearside, 1 Modified Gaussian for awayside
  * G(X) = A(1+(1/(2n))((x/sigma)/2)^2)^(-n)
  * see https://arxiv.org/pdf/1509.04732v2.pdf
  */
TF1 * phiCorrFit_1ModGaus_1ModGaus(TH1D * hist, std::string name)
{
  //TString functionString = "[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1))"; 
  TString functionString = "[0] * (TMath::Gamma([2])/(TMath::Sqrt(2*TMath::Pi()*([2]-1.5))*[1]*TMath::Gamma([2]-0.5))) *(TMath::Power(1 + TMath::Power((x+2.*TMath::Pi())/[1],2)/(2*([2]-1.5)),-[2]) + TMath::Power(1 + TMath::Power((x)/[1],2)/(2*([2]-1.5)),-[2])";
  functionString += " + TMath::Power(1 + TMath::Power((x-2*TMath::Pi())/[1],2)/(2*([2]-1.5)),-[2]))";
  functionString += "+ [3] *  (TMath::Gamma([2])/(TMath::Sqrt(2*TMath::Pi()*([5]-1.5))*[4]*TMath::Gamma([5]-0.5))) * (TMath::Power(1 + TMath::Power((x+TMath::Pi())/[4],2)/(2*([5]-1.5)),-[5]) + TMath::Power(1 + TMath::Power((x-TMath::Pi())/[4],2)/(2*([5]-1.5)),-[5])";
  functionString += " + TMath::Power(1 + TMath::Power((x-3*TMath::Pi())/[4],2)/(2*([5]-1.5)),-[5])) + [6]";

  TF1 * fit = new TF1(name.c_str(),functionString,-PI/2,3.0*PI/2);
  
  // FIXME treating sigma as omega for now
  // treating n as beta for now

  //name parameters
  fit->SetParName(0,"Y_NS");
  fit->SetParName(1,"Sigma_NS");
  fit->SetParName(2,"Beta_NS");
//  fit->SetParName(3,"Sigma_2_NS");
  fit->SetParName(3,"Y_AS");
  fit->SetParName(4,"Sigma_AS");
  fit->SetParName(5,"Beta_AS");
  fit->SetParName(6,"B");

//  fit->SetParName(8,"R");
//  fit->FixParameter(8,R);


  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  double B_g = 0;
  double Y_ns_g = max_val / 2.;
  double omega_ns_g = 0.4;
  double beta_ns_g = 3.;
  double Y_as_g = max_val / 2.;
  double omega_as_g = 0.4;
  double beta_as_g = 3.;
  
  double beta_min = 0.5;
//  double beta_min = 1.5;
  double beta_max = 100;  // damn you, VHS


  TH1D *hist_clone = (TH1D *) hist->Clone();
  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }

  double near_yield_max,away_yield_max;
  hist_clone->GetXaxis()->SetRangeUser(-PI/2,PI/2);

  near_yield_max = abs(1.2*hist_clone->Integral("width"));

  int zero_bin = hist->FindFixBin(0.);
	hist_clone->GetXaxis()->SetRangeUser(-PI/2,PI);
//	int minBin = hist->GetMiniumBin();

  int bin_diff = abs(hist_clone->GetMinimumBin() - zero_bin);
  hist_clone->GetXaxis()->SetRange(zero_bin - bin_diff, zero_bin + bin_diff);
  Y_ns_g = hist_clone->Integral("width");
  omega_ns_g = TMath::Max(hist_clone->GetRMS(),0.04);

  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  //  away_yield_max = hist->Integral("width") - 0.9 * B_g * PI;
  // approximate method: ignoring the wraparound 
  hist_clone->GetXaxis()->SetRange(zero_bin + bin_diff,hist_clone->FindFixBin(3.*PI/2));

  omega_as_g = hist_clone->GetRMS();

  away_yield_max = abs(1.2*hist_clone->Integral("width"));
//  away_yield_max = abs(1.2*hist->Integral("width") - 1.0 * B_g);

  hist_clone->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);

 // Y_ns_g = near_yield_max / 2;
//  Y_as_g = away_yield_max / 2;
  Y_as_g = hist_clone->Integral("width") - Y_ns_g;



/*  //  near_yield_max = hist->Integral("width") - 0.9 * B_g * PI; //0.8 in case bad B_g
//  near_yield_max = abs(2*hist->Integral("width") - 1.0 * B_g); //0.8 in case bad B_g
  Y_ns_g = abs(hist->Integral("width") - 1.0 * B_g);
  near_yield_max = 2* Y_ns_g;
  hist->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  //  away_yield_max = hist->Integral("width") - 0.9 * B_g * PI;
  Y_as_g = abs(2*hist->Integral("width") - 1.0 * B_g);
  away_yield_max = 2 * Y_as_g;
  hist->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);
//  Y_ns_g = near_yield_max;
//  Y_as_g = away_yield_max;
  printf("near_yield_max = %f \t away_yield_max = %f\n",near_yield_max,away_yield_max);
  //Near side prefit analysis
  omega_ns_g = TMath::Max(hist_clone->GetRMS(),0.2);
*/


  fit->SetParLimits(0,0.,2.5*near_yield_max);
  fit->SetParLimits(1,0.,5*PI); //FIXME
  fit->SetParLimits(2,beta_min,beta_max);
  fit->SetParLimits(3,0.,away_yield_max);
  fit->SetParLimits(4,0.,5*PI);
  fit->SetParLimits(5,beta_min,beta_max);
  fit->SetParLimits(6,0,max_val);

  //set parameters 

  fit->SetParameter(0,Y_ns_g);
  fit->SetParameter(1,omega_ns_g);
  fit->SetParameter(2,beta_ns_g);
  fit->SetParameter(3,Y_as_g);
  fit->SetParameter(4,omega_as_g);
  fit->SetParameter(5,beta_as_g);

	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(3,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
  }

  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(6,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(6,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(6,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }
  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");

  delete hist_clone;
  return fit;
}




/**
  * 1 Modified Gaussian for nearside, 1 Modified Gaussian for awayside
  * G(X) = A(1+(1/(2n))((x/sigma)/2)^2)^(-n)
  * see https://arxiv.org/pdf/1509.04732v2.pdf
  * FWHM as width parameter
  */

TF1 * phiCorrFit_1ModGaus_1ModGaus_Fwhm(TH1D * hist, std::string name)
{
  //TString functionString = "[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1))"; 

//  TString functionString = "[0] * (TMath::Gamma([2])/(TMath::Sqrt(2*TMath::Pi()*[2])*[1]*TMath::Gamma([2]-0.5))) *(TMath::Power(1 + TMath::Power((x+2.*TMath::Pi())/[1],2)/(2*[2]),-[2]) + TMath::Power(1 + TMath::Power((x)/[1],2)/(2*[2]),-[2])";
//  functionString += " + TMath::Power(1 + TMath::Power((x-2*TMath::Pi())/[1],2)/(2*[2]),-[2]))";
//  functionString += "+ [3] *  (TMath::Gamma([2])/(TMath::Sqrt(2*TMath::Pi()*[5])*[1]*TMath::Gamma([5]-0.5))) * (TMath::Power(1 + TMath::Power((x+TMath::Pi())/[4],2)/(2*[5]),-[5]) + TMath::Power(1 + TMath::Power((x-TMath::Pi())/[4],2)/(2*[5]),-[5])";
//  functionString += " + TMath::Power(1 + TMath::Power((x-3*TMath::Pi())/[4],2)/(2*[5]),-[5])) + [6]";




 // TString functionString = "[0] * (2*TMath::Sqrt(TMath::Power(2,1/[2])-1)/(TMath::Sqrt(TMath::Pi())*[1])) * (TMath::Gamma([2])/TMath::Gamma([2]-0.5)) * (TMath::Power(1 + (TMath::Power(2,1/[2])-1)*TMath::Power(2*(x+2*TMath::Pi())/[1],2)/(2*[2]),-[2]) + TMath::Power(1 + (TMath::Power(2,1/[2])-1)*TMath::Power(2*(x)/[1],2)/(2*[2]),-[2])";
  //functionString += " + TMath::Power(1 + (TMath::Power(2,1/[2])-1)*TMath::Power(2*(x-2*TMath::Pi())/[1],2)/(2*[2]),-[2])) ";

  //functionString += "+ [3] * (2*TMath::Sqrt(TMath::Power(2,1/[5])-1)/(TMath::Sqrt(TMath::Pi())*[4])) * (TMath::Gamma([5])/TMath::Gamma([5]-0.5)) * (TMath::Power(1 + (TMath::Power(2,1/[5])-1)*TMath::Power(2*(x+TMath::Pi())/[4],2)/(2*[5]),-[5]) + TMath::Power(1 + (TMath::Power(2,1/[5])-1)*TMath::Power(2*(x-TMath::Pi())/[4],2)/(2*[5]),-[5])";
 // functionString += " + TMath::Power(1 + (TMath::Power(2,1/[5])-1)*TMath::Power(2*(x-3*TMath::Pi())/[4],2)/(2*[5]),-[5])) + [6]";


// checking if factor of 2 is in error
  
  TString functionString = "[0] * (2*TMath::Sqrt(TMath::Power(2,1/[2])-1)/(TMath::Sqrt(TMath::Pi())*[1])) * (TMath::Gamma([2])/TMath::Gamma([2]-0.5)) * (TMath::Power(1 + (TMath::Power(2,1/[2])-1)*TMath::Power(2*(x+2*TMath::Pi())/[1],2)/(2*[2]),-[2]) + TMath::Power(1 + (TMath::Power(2,1/[2])-1)*TMath::Power(2*(x)/[1],2)/(2*[2]),-[2])";
  functionString += " + TMath::Power(1 + (TMath::Power(2,1/[2])-1)*TMath::Power(2*(x-2*TMath::Pi())/[1],2)/(2*[2]),-[2])) ";

  functionString += "+ [3] * (2*TMath::Sqrt(TMath::Power(2,1/[5])-1)/(TMath::Sqrt(TMath::Pi())*[4])) * (TMath::Gamma([5])/TMath::Gamma([5]-0.5)) * (TMath::Power(1 + (TMath::Power(2,1/[5])-1)*TMath::Power(2*(x+TMath::Pi())/[4],2)/(2*[5]),-[5]) + TMath::Power(1 + (TMath::Power(2,1/[5])-1)*TMath::Power(2*(x-TMath::Pi())/[4],2)/(2*[5]),-[5])";
  functionString += " + TMath::Power(1 + (TMath::Power(2,1/[5])-1)*TMath::Power(2*(x-3*TMath::Pi())/[4],2)/(2*[5]),-[5])) + [6]";

//  functionString += "+ [4] * (2*TMath::Sqrt(TMath::Power(2,1/[6])-1)/(TMath::Sqrt(TMath::Pi())*[5])) * (TMath::Gamma([6])/TMath::Gamma([6]-0.5)) * (TMath::Power(1 + (TMath::Power(2,1/[6])-1)*TMath::Power(2*(x+TMath::Pi())/[5],2)/(2*[6]),-[6]) + TMath::Power(1 + (TMath::Power(2,1/[6])-1)*TMath::Power(2*(x-TMath::Pi())/[5],2)/(2*[6]),-[6])";
//  functionString += " + TMath::Power(1 + (TMath::Power(2,1/[6])-1)*TMath::Power(2*(x-3*TMath::Pi())/[5],2)/(2*[6]),-[6])) + [7]";

  TF1 * fit = new TF1(name.c_str(),functionString,-PI/2,3.0*PI/2);
  
  // FIXME treating sigma as omega for now
  // treating n as beta for now

  //name parameters
  fit->SetParName(0,"Y_NS");
  fit->SetParName(1,"Omega_NS");
  fit->SetParName(2,"Beta_NS");
  fit->SetParName(3,"Y_AS");
  fit->SetParName(4,"Omega_AS");
  fit->SetParName(5,"Beta_AS");
  fit->SetParName(6,"B");


  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  double B_g = 0;
  double Y_ns_g = max_val / 2.;
  double omega_ns_g = 0.4;
  double beta_ns_g = 3.;
  double Y_as_g = max_val / 2.;
  double omega_as_g = 0.4;
  double beta_as_g = 3.;
  
  double beta_min = 0.5;
  double beta_max = 100;  // damn you, VHS




  TH1D *hist_clone = (TH1D *) hist->Clone();
  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }

  double near_yield_max,away_yield_max;
  // near_side
 // hist->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  hist_clone->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  //  near_yield_max = hist->Integral("width") - 0.9 * B_g * PI; //0.8 in case bad B_g
  //near_yield_max = abs(1.2*hist->Integral("width") - 1.0 * B_g); //0.8 in case bad B_g
  near_yield_max = abs(1.2*hist_clone->Integral("width"));
  
  int zero_bin = hist->FindFixBin(0.);
  int bin_diff = abs(hist->GetMinimumBin() - zero_bin);
  hist_clone->GetXaxis()->SetRange(zero_bin - bin_diff, zero_bin + bin_diff);
  Y_ns_g = hist_clone->Integral("width");
  omega_ns_g = TMath::Max(hist_clone->GetRMS(),0.04);

  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  //  away_yield_max = hist->Integral("width") - 0.9 * B_g * PI;
  // approximate method: ignoring the wraparound 
  hist_clone->GetXaxis()->SetRange(zero_bin + bin_diff,hist_clone->FindFixBin(3.*PI/2));
  
  omega_as_g = hist_clone->GetRMS();

  away_yield_max = abs(1.2*hist_clone->Integral("width"));
//  away_yield_max = abs(1.2*hist->Integral("width") - 1.0 * B_g);

  hist_clone->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);

 // Y_ns_g = near_yield_max / 2;
//  Y_as_g = away_yield_max / 2;
  Y_as_g = hist_clone->Integral("width") - Y_ns_g;


  printf("near_yield_max = %f \t away_yield_max = %f\n",near_yield_max,away_yield_max);




  fit->SetParLimits(0,0.,2.5*near_yield_max);
  fit->SetParLimits(1,0.,PI);
  fit->SetParLimits(2,beta_min,beta_max);
  fit->SetParLimits(3,0.,away_yield_max);
  fit->SetParLimits(4,0.,PI);
  fit->SetParLimits(5,beta_min,beta_max);
  fit->SetParLimits(6,0,max_val);

  //set parameters 

  fit->SetParameter(0,Y_ns_g);
  fit->SetParameter(1,omega_ns_g);
  fit->SetParameter(2,beta_ns_g);
  fit->SetParameter(3,Y_as_g);
  fit->SetParameter(4,omega_as_g);
  fit->SetParameter(5,beta_as_g);

	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(3,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
  }

  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(6,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(6,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(6,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }
  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
//  hist->Fit(fit, "RQ0");
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");

  delete hist_clone;

  return fit;
}






/** 
 * 2 Gaussians for nearside, sech(|x|^beta) for awayside
 * This is slightly worse than a normalized gaussian
 * Standard deviation sigma as width for awayside
 */

TF1 * phiCorrFit_2Gaus_SecH(TH1D * hist, std::string name)
{
  //  TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * ([6] * pow(TMath::Log(2),1./[6]))/ ([5]*TMath::Gamma(1./[6])) * (TMath::Exp(-TMath::Log(2) * pow(2*abs(x-TMath::Pi())/([5]),[6])) + TMath::Exp(-TMath::Log(2) * pow(2*abs(x+TMath::Pi())/([5]),[6]))) + [7]",-PI/2,3.0*PI/2);
  TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * 2*TMath::Power(TMath::ACosH(2),1/[6])/(TMath::Pi()*[5]) * (1./TMath::CosH(-TMath::ACosH(2)*TMath::Power(TMath::Abs(2*(x+TMath::Pi())/[5]),[6])) + 1./TMath::CosH(-TMath::ACosH(2)*TMath::Power(TMath::Abs(2*(x-TMath::Pi())/[5]),[6] ))) + [7]",-PI/2,3.0*PI/2);

  fprintf(stderr,"Warning, this doesn't work yet");
  // FIXME: this doesn't actually use sigma yet


  //name parameters
  fit->SetParName(0,"Y_NS_1");
  fit->SetParName(1,"Sigma_1_NS");
  fit->SetParName(2,"Y_NS_2");
  // NOTE: sigma_2_NS is the difference between sigma_2 and sigma_1
  fit->SetParName(3,"Sigma_2_NS");
  fit->SetParName(4,"Y_AS");
  fit->SetParName(5,"Sigma_AS");
//  fit->SetParName(5,"Omega_AS");
  fit->SetParName(6,"Beta_AS");
  fit->SetParName(7,"B");

  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  double min_val = hist->GetBinContent(hist->GetMinimumBin());
  double sigma_ns_g = 0.8;
  double sigma_ns_1_g = 0.3;
  double sigma_ns_2_g = 0.3;
  double sigma_as_g = 0.8;
  double omega_as_g = 1.9;
  double beta_as_g = 1.5;
  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();

  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }



  double near_yield_max,away_yield_max;
  hist->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  //  near_yield_max = hist->Integral("width") - 0.9 * B_g * PI; //0.8 in case bad B_g
  near_yield_max = hist->Integral("width") - 1.0 * B_g * PI; //0.8 in case bad B_g
  hist->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  //  away_yield_max = hist->Integral("width") - 0.9 * B_g * PI;
  away_yield_max = hist->Integral("width") - 1.0 * B_g * PI;
  hist->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);

  printf("near_yield_max = %f \t away_yield_max = %f\n",near_yield_max,away_yield_max);

  //Near side prefit analysis
  sigma_ns_g = TMath::Max(hist_clone->GetRMS(),0.2);

  double Y_ns_1_g,Y_ns_2_g;

  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);


  printf("Starting near side ...\n");

  // new idea: use separate fits: one for peak < R, on for peak > R
  // since root doesn't let me fit using two ranges, I will use each
  // range separately, and average

  // better idea: do the outer fit first
  // then subtract it, do the second fit

  // Trying to make this independent of jet resolution parameter
  // double r_ring = 0.2;  // estimated area where ring effect 
  double r_ring = R; //FIXME This might break everything
  // from jet reconstruction affects near side



  printf("Wide peak pre-fit ...\n");
  smooth_func->SetParLimits(0,0,2*near_yield_max);
  //peak2 (leftside)
  smooth_func->SetParameter(0,near_yield_max / 10.);
  smooth_func->SetParameter(1,r_ring);
  smooth_func->SetParLimits(1,0,PI);
  hist_clone->GetXaxis()->SetRangeUser(-PI/2,-r_ring/1.25);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = smooth_func->GetParameter(0);
  sigma_ns_2_g = smooth_func->GetParameter(1);

  //peak2 (rightside)
  smooth_func->SetParameter(0,near_yield_max / 10.);
  smooth_func->SetParameter(1,r_ring);
  hist_clone->GetXaxis()->SetRangeUser(r_ring/1.25,PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = (Y_ns_2_g + smooth_func->GetParameter(0)) / 2.;
  sigma_ns_2_g = (sigma_ns_2_g + smooth_func->GetParameter(1)) / 2.;

  Y_ns_2_g = min(near_yield_max,Y_ns_2_g/5);

  //subtract guess for wider peak from hist_clone
  TH1F *hist_clone_clone = (TH1F *) hist_clone->Clone();
  // hist_clone_clone->Add(smooth_func,-1.);

  printf("Sharp peak pre-fit ...\n");
  //peak 1

  //fixing the maximum point:  use unnormalized gaussian
  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],0)",-PI/2,PI/2);

  hist_clone_clone->GetXaxis()->SetRangeUser(-r_ring,r_ring);

  smooth_func->SetParameter(1,hist_clone_clone->GetRMS());
  // smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  // smooth_func->SetParameter(0,near_yield_max/3);
  double clone_clone_max_value = hist_clone_clone->GetBinContent(hist_clone_clone->GetMaximumBin());

  smooth_func->FixParameter(0,clone_clone_max_value);
  //  smooth_func->SetParameter(1,0.8*sigma_ns_2_g);
  smooth_func->SetParLimits(0,0,1.2*clone_clone_max_value);
  //  smooth_func->SetParLimits(1,0,1.2*sigma_ns_2_g);
  smooth_func->SetParLimits(1,0,2.2*sigma_ns_2_g);

  hist_clone_clone->Fit(smooth_func,"R0Q");
  //  Y_ns_1_g = smooth_func->GetParameter(0);
  sigma_ns_1_g = smooth_func->GetParameter(1);
  Y_ns_1_g = smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g;

  delete hist_clone_clone; // unless we want to use this for awayside?

  // Redefining sigma_ns_2_g to be the difference between the two 
  // sigmas, thus maintaining the ordering in the final function

  printf("Found Y_ns_1_g     = %f\t Y_ns_2_g     = %f\n",Y_ns_1_g,Y_ns_2_g);
  printf("Found sigma_ns_1_g = %f\t sigma_ns_2_g = %f\n",sigma_ns_1_g,sigma_ns_2_g);

  if (sigma_ns_2_g > sigma_ns_1_g) sigma_ns_2_g = sigma_ns_2_g - sigma_ns_1_g;
  else sigma_ns_2_g = r_ring;


  printf("Starting away side ...\n");
  // Away Side
  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  sigma_as_g = TMath::Max(hist_clone->GetRMS() ,0.1);
  omega_as_g = sigma_as_g * 2.35;
  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  delete smooth_func;
  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",+PI/2,3.0*PI/2);
  double yield_as = hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2));
  // double yield_as = abs(hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2)));
  smooth_func->SetParLimits(0,0,away_yield_max);
  smooth_func->SetParLimits(1,0,4.0*PI/3);

  smooth_func->SetParameter(0,0.9*yield_as);
  //  smooth_func->SetParameter(0,(max_val_as - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  smooth_func->SetParameter(1,sigma_as_g);
  // // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
  hist_clone->Fit(smooth_func,"R0Q");
  //  sigma_as_g = smooth_func->GetParameter(1);
  //double Y_as_g = smooth_func->Eval(PI) *TMath::Sqrt(2*PI)*sigma_as_g;
  double Y_as_g = smooth_func->GetParameter(0);
  //Y_as_g = yield_as;
  //FIXME

  delete hist_clone;

  fit->SetParLimits(0,0.,2.5*near_yield_max);
  fit->SetParLimits(1,0.,PI);
  fit->SetParLimits(2,0.,1.0*near_yield_max);
  fit->SetParLimits(3,0.,PI/4);
  fit->SetParLimits(4,0.,away_yield_max);
  fit->SetParLimits(5,0.,PI);
  fit->SetParLimits(6,0.5,4.);
  fit->SetParLimits(7,0,1.1*min_val);

  //set parameters 

  fit->SetParameter(0,Y_ns_1_g);
  fit->SetParameter(1,sigma_ns_1_g);
  fit->SetParameter(2,Y_ns_2_g);
  fit->SetParameter(3,sigma_ns_2_g);
  fit->SetParameter(4,Y_as_g);
  fit->SetParameter(5,omega_as_g);
  fit->SetParameter(6,beta_as_g);

	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(4,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
    fit->FixParameter(2,0);
  }

  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(7,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(7,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(7,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }


  // Fit to histogram
  // R uses the range defined in the fit function
  // 0 ensures that the fit isn't drawn
  // Q ensures minimum printing
  // + adds the fit to function list to ensure that it is not deleted on the creation of a new fit

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");


  return fit;
}



/** 
 * 2 Gaussians for nearside, sech(|x|^beta) for awayside
 * This is slightly worse than a normalized gaussian
 * FWHM omega as width for awayside
 */

TF1 * phiCorrFit_2Gaus_SecH_Fwhm(TH1D * hist, std::string name)
{
  //  TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * ([6] * pow(TMath::Log(2),1./[6]))/ ([5]*TMath::Gamma(1./[6])) * (TMath::Exp(-TMath::Log(2) * pow(2*abs(x-TMath::Pi())/([5]),[6])) + TMath::Exp(-TMath::Log(2) * pow(2*abs(x+TMath::Pi())/([5]),[6]))) + [7]",-PI/2,3.0*PI/2);
  TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[1]+[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * 2*TMath::Power(TMath::ACosH(2),1/[6])/(TMath::Pi()*[5]) * (1./TMath::CosH(-TMath::ACosH(2)*TMath::Power(TMath::Abs(2*(x+TMath::Pi())/[5]),[6])) + 1./TMath::CosH(-TMath::ACosH(2)*TMath::Power(TMath::Abs(2*(x-TMath::Pi())/[5]),[6] ))) + [7]",-PI/2,3.0*PI/2);

  //name parameters
  fit->SetParName(0,"Y_NS_1");
  fit->SetParName(1,"Sigma_1_NS");
  fit->SetParName(2,"Y_NS_2");
  // NOTE: sigma_2_NS is the difference between sigma_2 and sigma_1
  fit->SetParName(3,"Sigma_2_NS");
  fit->SetParName(4,"Y_AS");
  fit->SetParName(5,"Omega_AS");
  fit->SetParName(6,"Beta_AS");
  fit->SetParName(7,"B");

  //setting limits
  double max_val = hist->GetBinContent(hist->GetMaximumBin());
  double min_val = hist->GetBinContent(hist->GetMinimumBin());
  double sigma_ns_g = 0.8;
  double sigma_ns_1_g = 0.3;
  double sigma_ns_2_g = 0.3;
  double sigma_as_g = 0.8;
  double omega_as_g = 1.9;
  double beta_as_g = 1.5;
  double B_g = 0;

  TH1D *hist_clone = (TH1D *) hist->Clone();
  printf("Estimating background ...\n");
  B_g = getDPhiBackground(hist);
  if (B_g != 0) {
    TF1 * constant_function = new TF1("constant","[0]",-PI/2,3.0*PI/2);
   //  Subtracting constant background
    constant_function->SetParameter(0,-B_g);
    hist_clone->Add(constant_function);
  }


  double near_yield_max,away_yield_max;
  hist->GetXaxis()->SetRangeUser(-PI/2,PI/2);
  //  near_yield_max = hist->Integral("width") - 0.9 * B_g * PI; //0.8 in case bad B_g
  near_yield_max = hist->Integral("width") - 1.0 * B_g * PI; //0.8 in case bad B_g
  hist->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  //  away_yield_max = hist->Integral("width") - 0.9 * B_g * PI;
  away_yield_max = hist->Integral("width") - 1.0 * B_g * PI;
  hist->GetXaxis()->SetRangeUser(-PI/2,3.*PI/2);

  printf("near_yield_max = %f \t away_yield_max = %f\n",near_yield_max,away_yield_max);

  //Near side prefit analysis
  sigma_ns_g = TMath::Max(hist_clone->GetRMS(),0.2);

  double Y_ns_1_g,Y_ns_2_g;

  TF1 *smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,0,[1],1)",-PI/2,3.0*PI/2);


  printf("Starting near side ...\n");

  // new idea: use separate fits: one for peak < R, on for peak > R
  // since root doesn't let me fit using two ranges, I will use each
  // range separately, and average

  // better idea: do the outer fit first
  // then subtract it, do the second fit

  // Trying to make this independent of jet resolution parameter
  // double r_ring = 0.2;  // estimated area where ring effect 
  double r_ring = R; //FIXME This might break everything
  // from jet reconstruction affects near side



  printf("Wide peak pre-fit ...\n");
  smooth_func->SetParLimits(0,0,2*near_yield_max);
  //peak2 (leftside)
  smooth_func->SetParameter(0,near_yield_max / 10.);
  smooth_func->SetParameter(1,r_ring);
  smooth_func->SetParLimits(1,0,PI);
  hist_clone->GetXaxis()->SetRangeUser(-PI/2,-r_ring/1.25);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = smooth_func->GetParameter(0);
  sigma_ns_2_g = smooth_func->GetParameter(1);

  //peak2 (rightside)
  smooth_func->SetParameter(0,near_yield_max / 10.);
  smooth_func->SetParameter(1,r_ring);
  hist_clone->GetXaxis()->SetRangeUser(r_ring/1.25,PI/2);
  hist_clone->Fit(smooth_func,"R0Q");
  Y_ns_2_g = (Y_ns_2_g + smooth_func->GetParameter(0)) / 2.;
  sigma_ns_2_g = (sigma_ns_2_g + smooth_func->GetParameter(1)) / 2.;

  Y_ns_2_g = min(near_yield_max,Y_ns_2_g/5);

  //subtract guess for wider peak from hist_clone
  /*  smooth_func->SetParLimits(0,0,near_yield_max);
    smooth_func->SetParLimits(1,0,PI);
    smooth_func->SetParameter(0,0.3*Y_ns_2_g); //playing it safe
    smooth_func->SetParameter(1,sigma_ns_2_g);
    */
  TH1F *hist_clone_clone = (TH1F *) hist_clone->Clone();
  // hist_clone_clone->Add(smooth_func,-1.);

  printf("Sharp peak pre-fit ...\n");
  //peak 1

  //fixing the maximum point:  use unnormalized gaussian

  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],0)",-PI/2,PI/2);

  hist_clone_clone->GetXaxis()->SetRangeUser(-r_ring,r_ring);

  smooth_func->SetParameter(1,hist_clone_clone->GetRMS());
  // smooth_func->SetParameter(0,(max_val - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
  // smooth_func->SetParameter(0,near_yield_max/3);
  double clone_clone_max_value = hist_clone_clone->GetBinContent(hist_clone_clone->GetMaximumBin());

  smooth_func->FixParameter(0,clone_clone_max_value);
  //  smooth_func->SetParameter(1,0.8*sigma_ns_2_g);
  smooth_func->SetParLimits(0,0,1.2*clone_clone_max_value);
  //  smooth_func->SetParLimits(1,0,1.2*sigma_ns_2_g);
  smooth_func->SetParLimits(1,0,2.2*sigma_ns_2_g);

  hist_clone_clone->Fit(smooth_func,"R0Q");
  //  Y_ns_1_g = smooth_func->GetParameter(0);
  sigma_ns_1_g = smooth_func->GetParameter(1);
  Y_ns_1_g = smooth_func->GetParameter(0) * TMath::Sqrt(2*PI) * sigma_ns_1_g;

  delete hist_clone_clone; // unless we want to use this for awayside?

  // Redefining sigma_ns_2_g to be the difference between the two 
  // sigmas, thus maintaining the ordering in the final function

  printf("Found Y_ns_1_g     = %f\t Y_ns_2_g     = %f\n",Y_ns_1_g,Y_ns_2_g);
  printf("Found sigma_ns_1_g = %f\t sigma_ns_2_g = %f\n",sigma_ns_1_g,sigma_ns_2_g);

  if (sigma_ns_2_g > sigma_ns_1_g) sigma_ns_2_g = sigma_ns_2_g - sigma_ns_1_g;
  else sigma_ns_2_g = r_ring;


  printf("Starting away side ...\n");
  // Away Side
  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
  sigma_as_g = TMath::Max(hist_clone->GetRMS() ,0.1);
  omega_as_g = sigma_as_g * 2.35;
  hist_clone->GetXaxis()->SetRangeUser(3.0*PI/4,5.0*PI/4);
  delete smooth_func;
//  smooth_func = new TF1("2nd_smooth_func","[0]*TMath::Gaus(x,TMath::Pi(),[1],1)",+PI/2,3.0*PI/2);
//  double yield_as = hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2),"width");
  // double yield_as = abs(hist_clone->Integral(hist_clone->FindBin(PI/2),hist_clone->FindBin(3*PI/2)));
//  smooth_func->SetParLimits(0,0,away_yield_max);
//  smooth_func->SetParLimits(1,0,4.0*PI/3);

 // smooth_func->SetParameter(0,0.9*yield_as);
  //  smooth_func->SetParameter(0,(max_val_as - B_g)* TMath::Sqrt(2*PI)*sigma_ns_g);
//  smooth_func->SetParameter(1,sigma_as_g);
  // // double Y_as_g = hist_clone->GetBinContent(hist_clone->GetMaximumBin())*TMath::Sqrt(2*PI)*sigma_as_g;
//  hist_clone->Fit(smooth_func,"R0Q");
  //  sigma_as_g = smooth_func->GetParameter(1);
  //double Y_as_g = smooth_func->Eval(PI) *TMath::Sqrt(2*PI)*sigma_as_g;
//  double Y_as_g = smooth_func->GetParameter(0);
  //Y_as_g = yield_as;
  //FIXME
	
	hist_clone->GetXaxis()->SetRangeUser(PI/2,3.0*PI/2);
	double Y_as_g = hist_clone->Integral("width");
	hist_clone->GetXaxis()->SetRangeUser(-PI/2,3.0*PI/2);


  double m1,m1_1,m1_2,m2;
  hist_clone->GetXaxis()->SetRangeUser(PI/2,3.*PI/2);
  m2 = hist_clone->GetStdDev();
  hist_clone->GetXaxis()->SetRangeUser(PI/2,PI);
//  m1_1 = PI - hist_clone->GetMean();
  m1_1 = hist_clone->GetMean();
  hist_clone->GetXaxis()->SetRangeUser(PI,3.*PI/2);
 // m1_2 = PI + hist_clone->GetMean();
  m1_2 =hist_clone->GetMean();
  m1 = m1_1 + m1_2;
  printf("AwaySide: m1_1 = %f \tm1_2 = %f \tm1 = %f \t m2 = %f\n",m1_1,m1_2,m1,m2);
  beta_as_g = max(1 + 6*m2/m1,0.5) ;
	printf("AwaySide: beta = m1/sqrt(m2) = %f\n",beta_as_g);
//  if (beta_as_g > 4.) beta_as_g = 2.0;



  delete hist_clone;

  fit->SetParLimits(0,0.,2.5*near_yield_max);
  fit->SetParLimits(1,0.,PI);
  fit->SetParLimits(2,0.,1.0*near_yield_max);
  fit->SetParLimits(3,0.,PI/4);
  fit->SetParLimits(4,0.,away_yield_max);
  fit->SetParLimits(5,0.,PI);
  fit->SetParLimits(6,0.5,4);
  fit->SetParLimits(7,0,1.1*min_val);


  fit->SetParameter(0,Y_ns_1_g);
  fit->SetParameter(1,sigma_ns_1_g);
  fit->SetParameter(2,Y_ns_2_g);
  fit->SetParameter(3,sigma_ns_2_g);
  fit->SetParameter(4,Y_as_g);
  fit->SetParameter(5,omega_as_g);
  fit->SetParameter(6,beta_as_g);



	// Switch for single-jet tests
	if (noAwaySide) {
		fit->FixParameter(4,0);
	}
  if (noNearSide) {
    fit->FixParameter(0,0);
    fit->FixParameter(2,0);
  }


  //set parameters 
  switch (dPhiBkgMethod) {
    case aNoDPhiBkg: 
      fit->FixParameter(7,B_g);
      break;
    case aFreeDPhiBkg:
      fit->SetParameter(7,B_g); 
      break;
    case aZYAM:
      fit->FixParameter(7,B_g);
      break;
    default:
      fprintf(stderr,"Error: Invalid dPhiBkgMethod: %d\n",dPhiBkgMethod);
      exit(1);
      break;
  }

  // Fit to histogram
  // R uses the range defined in the fit function
  // 0 ensures that the fit isn't drawn
  // Q ensures minimum printing
  // + adds the fit to function list to ensure that it is not deleted on the creation of a new fit

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hist->Fit(fit, "RQ0+");
  printParameters(fit,true,"Final Parameters");

  return fit;
}








/** Newer method: 2 gaussians for nearside, and a Voigt profile for the 
 * awayside
 *  To be coded
 */






TF2 * etaPhiCorr_fit(TH2F * hist, std::string name)
{
  // TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * ([6] * pow(TMath::Gamma(3./[6]),0.5))/ (2.*[5]*pow(TMath::Gamma(1./[6]),1.5)) * (TMath::Exp(-pow(abs(x-TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])) + TMath::Exp(-pow(abs(x+TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6]))) + [7]",-PI/2,3.0*PI/2);


  printf("======================================================================\n");
  printf("||  Attempting 2D Fit                                               ||\n");
  printf("======================================================================\n");


  TH2F * hClone = (TH2F *) hist->Clone();

  double local_delta_eta_range = hist->GetXaxis()->GetXmax();

  //rebin?
  int rebinX = 2; 
  int rebinY = 2; 
  hClone->Rebin2D(rebinX,rebinY);
  hClone->Scale(1./(rebinX * rebinY));


  TString function_string = "[0]*(TMath::Gaus(x,0,[1],1)*TMath::Gaus(y,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1)*TMath::Gaus(y,0,[1],1) + TMath::Gaus(x,0,[1],1)*TMath::Gaus(y,2*TMath::Pi(),[1],1))" ;
  function_string += " + [2]*(TMath::Gaus(x,0,[3],1)*TMath::Gaus(y,-2*TMath::Pi(),[3],1) + TMath::Gaus(x,0,[1]+[3],1)*TMath::Gaus(y,0,[1]+[3],1)  + TMath::Gaus(x,0,[1]+[3],1)*TMath::Gaus(y,2*TMath::Pi(),[1]+[3],1))";
  function_string += " + [4] * pow( ([6] * pow(TMath::Gamma(3./[7]),0.5))/ (2.*[5]*pow(TMath::Gamma(1./[7]),1.5)),2) * (TMath::Exp(-pow(abs(x)/([5]*pow(TMath::Gamma(1./[7])/TMath::Gamma(3./[7]),0.5)),[7]) - pow(abs(y-TMath::Pi())/([5]*pow(TMath::Gamma(1./[7])/TMath::Gamma(3./[7]),0.5)),[7]) ) + TMath::Exp(-pow(abs(x)/([5]*pow(TMath::Gamma(1./[7])/TMath::Gamma(3./[7]),0.5)),[7]) -pow(abs(y+TMath::Pi())/([5]*pow(TMath::Gamma(1./[7])/TMath::Gamma(3./[7]),0.5)),[7]) ))";
  // function_string += " + [4] * ([6] * pow(TMath::Gamma(3./[6]),0.5))/ (2.*[5]*pow(TMath::Gamma(1./[6]),1.5)) * (TMath::Exp(-pow(abs(y-TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])) + TMath::Exp(-pow(abs(y+TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])))";
 // function_string += "+ [8]*(1-TMath::Abs(x/2))";
 // function_string += "+ [8] * (1 - ((TMath::Abs(x))+[10])/(2*[9])) * (0 < 1 - ((TMath::Abs(x))+[10])/(2*[9]))";
  function_string += "+ [8] * (1 > ((TMath::Abs(x))+[10])/(2*[9]))";


// The old trapezoid code:
//  function_string += "*((TMath::Abs(x)<=[10])";
//  function_string += "+(TMath::Abs(x)>[10])*(2.*[9]-[10] - TMath::Abs(x))/(2.*[9]-2*[10])  )";

// [9] == eta_cut
// [10] == R

// reference 

//  TString function_string = "[0] * (1 > ((TMath::Abs(x))+[6])/(2*[5]))";
//  function_string += "*((TMath::Abs(x)<=[6])";
//  function_string += "+(TMath::Abs(x)>[6])*(2.*[5]-[6] - TMath::Abs(x))/(2.*[5]-2*[6])  )";


  TF2 * fit = new TF2(name.c_str(),function_string,-local_delta_eta_range,local_delta_eta_range,-PI/2,3.0*PI/2);

  double near_yield_max,away_yield_max,Y_ns_1_g,sigma_ns_1_g,Y_ns_2_g,sigma_ns_2_g;
  double B_g = 0.3;

  double max_val = hClone->GetBinContent(hClone->GetMaximumBin());

  fit->SetParName(0,"Y_NS_1");
  fit->SetParName(1,"Sigma_1_NS");
  fit->SetParName(2,"Y_NS_2");
  // NOTE: sigma_2_NS is the difference between sigma_2 and sigma_1
  fit->SetParName(3,"Sigma_2_NS");
  fit->SetParName(4,"Y_AS");
  fit->SetParName(5,"Sigma_Eta_AS");
  fit->SetParName(6,"Sigma_Phi_AS");
  fit->SetParName(7,"Beta_AS");
  fit->SetParName(8,"B");
  fit->SetParName(9,"EtaCut");
  fit->SetParName(10,"R");

  fit->FixParameter(9,eta_cut);
  fit->FixParameter(10,R);


  hClone->GetYaxis()->SetRangeUser(-PI/2,PI/2);
  hClone->GetXaxis()->SetRangeUser(1,2);
  B_g = hClone->Integral("width");  
  hClone->GetXaxis()->SetRangeUser(-2,-1);
  B_g += hClone->Integral("width");  
  B_g = B_g / (2 * PI); //  is now mean value in outer region.
  B_g = 4 * B_g;  // is now the appropriate value.
  B_g = TMath::Abs(B_g); 


  // Returning to normal boundaries
  hClone->GetXaxis()->SetRangeUser(-2.,2.);
  hClone->GetYaxis()->SetRangeUser(-PI/2,3*PI/2);

  // Finding other guesses and ranges
  // subtract current tent, take some values:
  TH2F * hClone2 = (TH2F *) hClone->Clone();
  fit->SetParameter("Y_NS_1",0.);
  fit->SetParameter("Sigma_1_NS",1.);
  fit->SetParameter("Y_NS_2",0.);
  fit->SetParameter("Sigma_2_NS",0.1);
  fit->SetParameter("Y_AS",0.);
  fit->SetParameter("Sigma_Eta_AS",1.);
  fit->SetParameter("Sigma_Phi_AS",0.5);
  fit->SetParameter("Beta_AS",2);
  fit->SetParameter("B",B_g * 0.8); //assuming overestimation
  hClone2->Add(fit,-1);

  // Nearside
  hClone2->GetYaxis()->SetRangeUser(-PI/2,PI/2);
  Y_ns_1_g = hClone2->Integral("width");
  // checking for oversubtraction
  if (Y_ns_1_g < 0) {
    hClone2->Add(fit,1);
    Y_ns_1_g = hClone2->Integral("width");
  }
  near_yield_max = 2 * Y_ns_1_g;
  Y_ns_2_g = Y_ns_1_g * 0.111;
  Y_ns_1_g = Y_ns_1_g * 0.889;
  sigma_ns_1_g = hClone2->GetRMS()/1.4;
  sigma_ns_2_g = 0.1; 


  // Awayside
  hClone2->GetYaxis()->SetRangeUser(PI/2,3.*PI/2);

  double Y_as_g = hClone2->Integral("width");
  away_yield_max = 2 * Y_as_g;
  double sigma_as_eta_g = 1;
  double sigma_as_phi_g = 0.3;

  double beta_as_g = 2;

  delete hClone2;

  fit->SetParameter(0,Y_ns_1_g);      fit->SetParLimits(0,0.,near_yield_max);
  fit->SetParameter(1,sigma_ns_1_g);  fit->SetParLimits(1,0.,PI);
  fit->SetParameter(2,Y_ns_2_g);      fit->SetParLimits(2,0.,near_yield_max);
  fit->SetParameter(3,sigma_ns_2_g);  fit->SetParLimits(3,0.,PI/4);
  fit->SetParameter(4,Y_as_g);        fit->SetParLimits(4,0.,away_yield_max);
  fit->SetParameter(5,sigma_as_eta_g);  fit->SetParLimits(5,0.,PI);
  fit->SetParameter(6,sigma_as_phi_g);  fit->SetParLimits(6,0.,PI);
  fit->SetParameter(7,beta_as_g);       fit->SetParLimits(7,0.1,5);
  //  fit->FixParameter(8,B_g);             fit->SetParLimits(8,0.,2*B_g);
  fit->SetParameter(8,B_g);             fit->SetParLimits(8,0.,max_val);
  // fit->SetParameter(7,B_g);
  //  fit->SetParameter(8,B_slope_g);


  // FIXME: check if B_g is very small.  If so, fix it to 0


  // Tests
  hClone->GetYaxis()->SetRangeUser(-PI/2,PI/2);
  fit->FixParameter(4,0);
  fit->FixParameter(5,1);
  fit->FixParameter(6,1);
  fit->FixParameter(7,2);
    

  printParameters(fit,false,"Parameter Guesses");
  //FIT HERE
  if (!noFitting) hClone->Fit(fit, "R0+ WW");
  printParameters(fit,true,"Final Parameters");
  //	hist->Fit(fit, "R0+");
  delete hClone;

  return fit;
}



/**
 * Method that takes a dPhi projection and returns the TF2 that matches the 
 * combinatorial background of the unprojected fit.
 * Assumes finite acceptance has been corrected for
 */
TF1 * dEta_fit(TH1D * hist, std::string name, double deltaEtaRange, bool doFit)
{
  // TF1 * fit = new TF1(name.c_str(),"[0]*(TMath::Gaus(x,-2*TMath::Pi(),[1],1) + TMath::Gaus(x,0,[1],1) + TMath::Gaus(x,2*TMath::Pi(),[1],1)) + [2]*(TMath::Gaus(x,-2*TMath::Pi(),[3],1) + TMath::Gaus(x,0,[1]+[3],1) + TMath::Gaus(x,2*TMath::Pi(),[1]+[3],1)) + [4] * ([6] * pow(TMath::Gamma(3./[6]),0.5))/ (2.*[5]*pow(TMath::Gamma(1./[6]),1.5)) * (TMath::Exp(-pow(abs(x-TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])) + TMath::Exp(-pow(abs(x+TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6]))) + [7]",-PI/2,3.0*PI/2);

  //double local_delta_eta_range = hist->GetXaxis()->GetXmax();
	double local_delta_eta_range = deltaEtaRange;
	printf("dEta_fit found hist with DeltaEtaRange = %f\n",local_delta_eta_range);
  //  TH1F * hClone = (TH1F *) hist->Clone();
  // [5] == eta_cut
  // [6] == R
 // TString function_string = "[0] * (1 - ((TMath::Abs(x))+[6])/(2*[5])) ";
 // TString function_string = "[0] * (1 - ((TMath::Abs(x))+[6])/(2*[5])) * (0 < 1 - ((TMath::Abs(x))+[6])/(2*[5]))";

// FIXME this needs to be corrected to not use a trapezoid, but a constant, since 
// acceptance effect corrections have been implemented.

 // TString function_string = "[0] * (1 > ((TMath::Abs(x))+[6])/(2*[5]))";
  TString function_string = "[0]";

// old trapezoid code
//  function_string += "*((TMath::Abs(x)<=[6])";
//  function_string += "+(TMath::Abs(x)>[6])*(2.*[5]-[6] - TMath::Abs(x))/(2.*[5]-2*[6])  )";


//  TString function_string = "[0] * (1 - TMath::Abs(x/2))";
  function_string += "+ [1]*TMath::Gaus(x,0,[2],1)";
  function_string += "+ [3]*TMath::Gaus(x,0,[2] + [4],1)";

  // function_string += " + [4] * ([6] * pow(TMath::Gamma(3./[6]),0.5))/ (2.*[5]*pow(TMath::Gamma(1./[6]),1.5)) * (TMath::Exp(-pow(abs(y-TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])) + TMath::Exp(-pow(abs(y+TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])))";

  double A_g, sigma_g, A_2_g, sigma_delta_g, B_g;
 // TF1 * fit1 = new TF1("1d_fit",function_string,-2*eta_cut,2*eta_cut);
  TF1 * fit1 = new TF1("1d_fit",function_string,-local_delta_eta_range,local_delta_eta_range);

	
	// Defining the range to use for estimating the background
	double cut_from_edge = 0.1;
  double min_dEta = -(local_delta_eta_range - cut_from_edge);
  double max_dEta = -local_delta_eta_range/2;

	printf("DEBUG: min = %f, max = %f\n",min_dEta,max_dEta);

  double max_val = 1.2*hist->GetBinContent(hist->GetMaximumBin());

  B_g = hist->Integral(hist->FindFixBin(min_dEta),hist->FindFixBin(max_dEta),"width");
//  min_dEta = -max_dEta;
//  max_dEta = eta_range/1.1;
//	min_dEta = 2;
//	max_dEta = 3.7;

  min_dEta = local_delta_eta_range/2.;
  max_dEta = local_delta_eta_range - cut_from_edge;

	printf("DEBUG: min = %f, max = %f\n",min_dEta,max_dEta);

  B_g += hist->Integral(hist->FindFixBin(min_dEta),hist->FindFixBin(max_dEta),"width");
  B_g = B_g / 2 ; // mean of integrals
//	printf("B_g_2 = %f\n",B_g);
  // next, set B_g to solution of integral of fit function equal to integral 
  // of data
//  B_g = B_g * (2*eta_cut - 2*R) / (3*eta_cut - R);

	// old formula (for trapezoid)
//  B_g = 2*B_g * (eta_cut - R) / (2*eta_cut - R-(min_dEta + max_dEta)/2.);

	B_g = B_g / (max_dEta - min_dEta);


  A_g = hist->GetBinContent(hist->FindFixBin(0)) - B_g;
 // if (A_g < 0) A_g = B_g * 0.2;
  if (A_g < 0) A_g = max_val / 2.;
  sigma_g = 0.2;
  A_2_g= A_g/4;
  sigma_delta_g = 0.4;




 // fit1->SetParameter(0,B_g);     fit1->SetParLimits(2,0,max_val);
  fit1->SetParameter(1,A_g);     fit1->SetParLimits(1,0,max_val);
  fit1->SetParameter(2,sigma_g); fit1->SetParLimits(2,0,0.8);
  fit1->SetParameter(3,A_2_g);     fit1->SetParLimits(3,0,max_val);
  fit1->SetParameter(4,sigma_delta_g); fit1->SetParLimits(4,0,0.6);
//  fit1->FixParameter(5,eta_cut);
//    fit1->FixParameter(6,R);


	// FIXME experiment
//	fit1->FixParameter(3,0);
	

  if (doFit) {
    fit1->SetParameter(0,B_g);     
    fit1->SetParLimits(0,0,(B_g+A_g));     
  } else {
    fit1->FixParameter(0,B_g);     
  }
  if (!noFitting) hist->Fit(fit1,"NRQ0 WW");
  //  function_string += "[0] * (1 - TMath::Abs(x/2))";
  //  TF2 *fit2 = new TF2(name.c_str(),function_string,-2,2,-PI/2,3*PI/2);

  //rescale tent parameter to undo the scaling from projecting on to dEta
  //  fit2->SetParameter(0,fit1->GetParameter(2)/dPhi_range);
  return fit1;
  // return fit2;

}


TF1 * dEta_GenGaus_fit(TH1D * hist, std::string name, bool doFit) {

  double local_delta_eta_range = hist->GetXaxis()->GetXmax();
  TString function_string = "[0]";

  function_string += "+ [1]*TMath::Exp(-TMath::Power(TMath::Abs(x/[2]),[3]))";

  // function_string += " + [4] * ([6] * pow(TMath::Gamma(3./[6]),0.5))/ (2.*[5]*pow(TMath::Gamma(1./[6]),1.5)) * (TMath::Exp(-pow(abs(y-TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])) + TMath::Exp(-pow(abs(y+TMath::Pi())/([5]*pow(TMath::Gamma(1./[6])/TMath::Gamma(3./[6]),0.5)),[6])))";

  double A_g, sigma_g, beta_g, B_g;
 // TF1 * fit1 = new TF1("1d_fit",function_string,-2*eta_cut,2*eta_cut);
  TF1 * fit1 = new TF1("1d_fit",function_string,-local_delta_eta_range,local_delta_eta_range);

	

  double min_dEta = -(local_delta_eta_range - 0.1);
  double max_dEta = -local_delta_eta_range/2;

	printf("DEBUG: min = %f, max = %f\n",min_dEta,max_dEta);

  double max_val = 1.2*hist->GetBinContent(hist->GetMaximumBin());

  B_g = hist->Integral(hist->FindFixBin(min_dEta),hist->FindFixBin(max_dEta),"width");
//  min_dEta = -max_dEta;
//  max_dEta = eta_range/1.1;
//	min_dEta = 2;
//	max_dEta = 3.7;

  min_dEta = local_delta_eta_range/2.;
  max_dEta = local_delta_eta_range - 0.1;

	printf("DEBUG: min = %f, max = %f\n",min_dEta,max_dEta);

  B_g += hist->Integral(hist->FindFixBin(min_dEta),hist->FindFixBin(max_dEta),"width");
  B_g = B_g / 2 ; // mean of integrals
//	printf("B_g_2 = %f\n",B_g);
  // next, set B_g to solution of integral of fit function equal to integral 
  // of data
//  B_g = B_g * (2*eta_cut - 2*R) / (3*eta_cut - R);

	// old formula (for trapezoid)
//  B_g = 2*B_g * (eta_cut - R) / (2*eta_cut - R-(min_dEta + max_dEta)/2.);

	B_g = B_g / (max_dEta - min_dEta);


  A_g = hist->GetBinContent(hist->FindFixBin(0)) - B_g;
 // if (A_g < 0) A_g = B_g * 0.2;
  if (A_g < 0) A_g = max_val / 2.;
 
	// Estimating Sigma from a truncated rms
	TH1F * hClone = (TH1F *) hist->Clone();
	TF1 * fConst = new TF1("fConst","[0]",-local_delta_eta_range,local_delta_eta_range);
	fConst->SetParameter(0,B_g);
	hClone->Add(fConst);
	double dEta_RMS_Range = local_delta_eta_range / 3.;
	hClone->GetXaxis()->SetRangeUser(-dEta_RMS_Range,dEta_RMS_Range);
	sigma_g = hClone->GetRMS();


//	sigma_g = 0.2;
	beta_g = 2;

 // fit1->SetParameter(0,B_g);     fit1->SetParLimits(2,0,max_val);
  fit1->SetParameter(1,A_g);     fit1->SetParLimits(1,0,max_val);
  fit1->SetParameter(2,sigma_g); fit1->SetParLimits(2,0,1.0);
  fit1->SetParameter(3,beta_g);     fit1->SetParLimits(3,1.,10.);
//  fit1->FixParameter(5,eta_cut);
//    fit1->FixParameter(6,R);

	// FIXME experiment
//	fit1->FixParameter(3,0);
	

  if (doFit) {
    fit1->SetParameter(0,B_g);     
    fit1->SetParLimits(0,0,(B_g+A_g));     
  } else {
    fit1->FixParameter(0,B_g);     
  }
  if (!noFitting) hist->Fit(fit1,"NRQ0 WW");
  //  function_string += "[0] * (1 - TMath::Abs(x/2))";
  //  TF2 *fit2 = new TF2(name.c_str(),function_string,-2,2,-PI/2,3*PI/2);

  //rescale tent parameter to undo the scaling from projecting on to dEta
  //  fit2->SetParameter(0,fit1->GetParameter(2)/dPhi_range);
  return fit1;
  // return fit2;






}





TF1 * swiftJetFit(TH1F * hist, std::string name) {
 
  TF1 * fit = new TF1(name.c_str(),"[0]*TMath::Exp(-TMath::Power(TMath::Abs(x/[1]),[2])) + [3]",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  fit->SetParName(0,"Y");
  fit->SetParName(1,"alpha");
  fit->SetParName(2,"beta");
  fit->SetParName(3,"C");

  double max_val = hist->GetBinContent(hist->GetMaximumBin());

  fit->SetParLimits(0,0,1.3 * max_val);
  fit->SetParameter(0,max_val);
  fit->SetParLimits(1,0.01,2);
  fit->SetParameter(1,0.6);	
  fit->SetParLimits(2,2,12);
  fit->SetParameter(2,5.5);
  fit->SetParLimits(3,0,max_val/5);
  fit->SetParameter(3,0);

  hist->Fit(fit, "R0Q+");

  return fit;



    //old version

/*
 ///	TF1 * fit = new TF1(name.c_str(),"[0]*TMath::Gaus(x,0,[1],1) ",-(eta_cut-R),eta_cut-R);
  //	TF1 * fit = new TF1(name.c_str(),"[0]*TMath::Gaus(x,0,[1],1) ",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  TF1 * fit = new TF1(name.c_str(),"[0] ",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  fit->SetParName(0,"A");
  fit->SetParLimits(0,0,hist->GetBinContent(hist->GetMaximumBin()));		


  hist->Fit(fit, "R0Q+");

  return fit;
*/
}

TF1 * swiftParticleFit(TH1F * hist, std::string name) {
  ///	TF1 * fit = new TF1(name.c_str(),"[0]*TMath::Gaus(x,0,[1],1) ",-(eta_cut-R),eta_cut-R);
 // TF1 * fit = new TF1(name.c_str(),"[0]*TMath::Exp(-TMath::Power(TMath::Abs(x/[1]),[2])) + [3]",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  TF1 * fit = new TF1(name.c_str(),"[0]*TMath::Exp(-TMath::Power(TMath::Abs((x-[4])/[1]),[2])) + [3]",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  fit->SetParName(0,"Y");
  fit->SetParName(1,"alpha");
  fit->SetParName(2,"beta");
  fit->SetParName(3,"C");
  fit->SetParName(3,"mu");

  double max_val = hist->GetBinContent(hist->GetMaximumBin());

  fit->SetParLimits(0,0,1.3 * max_val);
  fit->SetParameter(0,max_val);
  fit->SetParLimits(1,0.01,2);
  fit->SetParameter(1,0.6);	
  fit->SetParLimits(2,2,12);
  fit->SetParameter(2,5.5);
  fit->SetParLimits(3,0,max_val/5);
  fit->SetParameter(3,0);
  fit->SetParLimits(4,-eta_cut,eta_cut);
  fit->SetParameter(4,0);

  hist->Fit(fit, "R0Q+");

  return fit;
}


TF2 * createSwiftBackground(TF1 *jetFit, TF1 * particleFit, std::string name) {

  TF2 * fit = new TF2(name.c_str(),"[0]*TMath::Exp(-TMath::Power(TMath::Abs(x/[1]),[2])) + [3]",jetFit->GetXmin() + particleFit->GetXmin(),jetFit->GetXmax());

  
  return fit;
}


