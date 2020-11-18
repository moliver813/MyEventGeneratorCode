// Compares two outfiles from phase2, adds appropriate labels
//#include <boost/algorithm/string/replace.hpp>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLine.h>
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
#include <TError.h>

bool ratioMode = false;
// FIXME find a way to turn ratio mode off for histograms that won't divide
bool enableCFiles = true;
bool enableRootFiles = false;
bool proj2D = false;

Int_t LineStyle = 2;
Int_t LineWidth = 4;
Int_t LineColor = kGray;	




const bool useCustomColor = true; // whether to use the rotating color palette, or just fixed list of colors.

// setting for space between pads
const float small = 1e-5;

void set_plot_style() {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);

}

TGraphErrors * DivideTGraphErrors(TGraphErrors * num, TGraphErrors * denom, TString name) {
  int nPoints = TMath::Min(num->GetN(),denom->GetN());
  int nSkipPoints = 0;
  TGraphErrors * newTGraph = new TGraphErrors(nPoints);
  newTGraph->SetName(name);
  for (int i = 0; i < nPoints; i++) {
    double x = num->GetX()[i];
    double x_err = num->GetEX()[i];
    double y1 = num->GetY()[i];
    double y2 = denom->GetY()[i];
    double y,y_err;
    if (y1 != 0 && y2 != 0) { // or maybe skip the point completely
      y = y1/y2;
      double y1_err = num->GetEY()[i];
      double y2_err = denom->GetEY()[i];
      y_err = y * TMath::Sqrt(TMath::Power(y1_err/y1,2)+TMath::Power(y2_err/y2,2));    
      newTGraph->SetPoint(i-nSkipPoints,x,y);
      newTGraph->SetPointError(i-nSkipPoints,x_err,y_err);
    } else {
      newTGraph->RemovePoint(i-nSkipPoints);
      nSkipPoints++;
      y = 0;
      y_err = 0;
    }

  }
  return newTGraph;

}



using namespace std;

struct fileLabel {
  fileLabel(TString _filepath, TString _label): filepath(_filepath),label(_label) {}
  TString filepath;
  TString label;
};

void cleanName(TString &str) {
		str.ReplaceAll(" ","_");
		str.ReplaceAll(".","_");
		str.ReplaceAll("/","_");
		str.ReplaceAll("%","_");
		str.ReplaceAll("+","_");
		str.ReplaceAll("#","_");
		str.ReplaceAll(",","_");
		str.ReplaceAll("-","_");
		str.ReplaceAll("(","");
		str.ReplaceAll(")","");
		str.ReplaceAll("<=","lte");
		str.ReplaceAll(">=","gte");
		str.ReplaceAll("<","lt");
		str.ReplaceAll(">","gt");
		str.ReplaceAll("=","_");
}

int GetCustomColor(int x) {
  // r,g,b in [0,255]
  //int r = (13 * x) % 256 ;
  //would be simpler to use GetColor(Float_t r, Float_t g, Float_t b)
	//that probably rounds any way.
	const float omega_r = 2*0.12;  const float phi_r = -3.1415;
	const float omega_g = 2*0.08;  const float phi_g = 1.57;
	const float omega_b = 1.5*0.16;  const float phi_b = -3.1415/2.;
// Original values:
//	omega_r = 0.06;  phi_r = -3.1415;
//	omega_g = 0.04;  phi_g = 1.57;
//	omega_b = 0.02;  phi_b = -3.1415/2.;

  float theta_r = 3.1415 * omega_r * x + phi_r;
  float theta_g = 3.1415 * omega_g * x + phi_g;
  float theta_b = 3.1415 * omega_b * x + phi_b;
  int r = 128 + round(127.*TMath::Cos(theta_r)); 
  int g = 128 + round(127.*TMath::Cos(theta_g));
  int b = 128 + round(127.*TMath::Cos(theta_b));
//  printf("r,g,b = (%d,%d,%d)\n",r,g,b);
  return TColor::GetColor(r,g,b);
}

const Int_t nTypeList = 11;// 9;

void compare(vector<fileLabel> *fileLabels, char * outputDirPath)  {

	gErrorIgnoreLevel = kWarning;

  Int_t colorList[nTypeList] = {kBlack,kRed,kOrange-3,kGreen-3,kBlue,kViolet,kMagenta-5,kRed-10,kGray,kCyan,kAzure};
  Int_t markerList[nTypeList] = {kFullSquare,kFullCircle,kFullDiamond,kOpenSquare,kOpenCircle,kOpenDiamond,kFullStar,23,34,46,42};

  double c_width = 500;
  double c_height = 500;

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
    f->SetTitle(fl.label.Data());
		f->SetName(nameNoSpace);
    files->push_back(f);
    
    printf("file %s, title %s\n",f->GetName(),f->GetTitle());

  }

  if ( !files->size()) exit(0);
  TFile *f0 = files->at(0);
	   
	// Moving to output directory
	gSystem->ChangeDirectory(outputDirPath);

  TCanvas *canvas = new TCanvas("canvas","canvas",c_width,c_height);
  gStyle->SetOptStat(0);
	gStyle->SetTitleSize(0.04,"xy");
 
  TList * keys = f0->GetListOfKeys();

  TObject *obj,*obj2;
  TIter next(keys);
  while ((obj = next())) {
    if (!obj->InheritsFrom("TNamed")) continue;
    if (!obj->InheritsFrom("TKey")) continue;

    TKey *kobj = dynamic_cast<TKey*>(obj);
    obj2 = kobj->ReadObj();

    vector<TObject *> * otherObj = new vector<TObject *>;

		TString objNameNoSpace;
	
    // Look for object in other files
    for (int j = 1 ; j < files->size(); j++) {
      TObject *o = files->at(j)->Get(obj2->GetName());
      if (o) {
        // FIXME testing
        if (!o->InheritsFrom("TF1") ) ((TNamed *) o)->SetTitle(files->at(j)->GetTitle()); 
        else ((TF1 *) o)->GetYaxis()->SetTitle(files->at(j)->GetTitle());
				objNameNoSpace = ((TNamed *) o)->GetName();
				objNameNoSpace = TString::Format("%s_%s",objNameNoSpace.Data(),files->at(j)->GetTitle());
				cleanName(objNameNoSpace);
				((TNamed *) o)->SetName(objNameNoSpace);
        otherObj->push_back(o);
      } else {
        printf("Object %s not found in file %s.\n",obj2->GetName(),files->at(j)->GetName());
      }

    }
    
    // May need to do something about legend location
    // Hoping for the best for now
    TLegend * legend = new TLegend(0.65,0.75,0.90,0.95);


		objNameNoSpace = obj->GetName();
		cleanName(objNameNoSpace);

		TFile *fOut = 0;
	  if (enableRootFiles) fOut = TFile::Open(Form("%s.root",objNameNoSpace.Data()),"RECREATE");

		TString objname = obj2->GetName();

		objNameNoSpace = obj2->GetName();
		objNameNoSpace = TString::Format("%s_%s",objNameNoSpace.Data(),files->at(0)->GetTitle());
		cleanName(objNameNoSpace);

		if (obj2->InheritsFrom("TH3")) continue; // don't do anything for 

    if (obj2->InheritsFrom("TH2")) {
      // Nothing obvious to do.  Maybe calculate differences/ratios ??
      TH2 * h2obj = (TH2 *) obj2;
 //     printf("TH2 %s\n",h2obj->GetName());

      double max_y, min_y;    
			int max_bin, min_bin;

			int nBinsX = h2obj->GetNbinsX();
			int nBinsY = h2obj->GetNbinsY();
			// If there is a small number of bins, project the histogram finitely many times
			if (nBinsX <= nTypeList) {
				printf("Projecting and Printing TH2F %s\n",h2obj->GetName());
				vector<vector<TH1D *>> arrProj;
				vector<TH1D *> arrProjFile;
				for (int i = 0; i < nBinsX; i++) {
					TH1D * projY = h2obj->ProjectionY(Form("%s_%d",h2obj->GetName(),i),i+1,i+1);
					arrProjFile.push_back(projY);
//					fOut->Add(projY);
				}
				arrProj.push_back(arrProjFile);

				//Do the same for the other files
        for (int j = 0; j < otherObj->size(); j++) {  
					TH2F *otherH = (TH2F *) otherObj->at(j);
					objNameNoSpace = otherH->GetName();
					objNameNoSpace = TString::Format("%s_%s",objNameNoSpace.Data(),files->at(0)->GetTitle());
					cleanName(objNameNoSpace);
					otherH->SetName(objNameNoSpace);
					if (enableRootFiles) fOut->Add(otherH);

					vector <TH1D *> arrProjOtherFile;
					for(int i = 0; i < nBinsX; i++) {
						TH1D * projY = otherH->ProjectionY(Form("%s_%d",otherH->GetName(),i),i+1,i+1);
						arrProjOtherFile.push_back(projY);
					}
					arrProj.push_back(arrProjOtherFile);
				}
				// Now, create the canvases
				//idea: dynamically scale tlegend with nobj
				TLegend * legend = new TLegend(0.55,0.45,0.90,0.85);
				//TLegend * legend = new TLegend(0.55,0.65,0.90,0.85);
			//	TLegend * legend = new TLegend(0.65,0.65,0.90,0.85);
				for (int i = 0; i < nBinsX; i++) {
					legend->Clear();
				  canvas->Clear();
					TH1D * h1 = arrProj[0][i];
					h1->SetLineColor(1);
					h1->SetMarkerColor(1);
					h1->SetMarkerStyle(markerList[0]);
					TAxis * h2objXaxis = h2obj->GetXaxis();
					h1->SetTitle(Form("%s (%.2f < %s < %.2f)",h2obj->GetYaxis()->GetTitle(),h2objXaxis->GetBinLowEdge(i+1),h2objXaxis->GetTitle(),h2objXaxis->GetBinUpEdge(i+1)));
					legend->AddEntry(h1,f0->GetTitle(),"lp");

					max_bin = h1->GetMaximumBin();
					min_bin = h1->GetMinimumBin();
					max_y = h1->GetBinContent(max_bin) + h1->GetBinErrorUp(max_bin);
					min_y = h1->GetBinContent(min_bin) - h1->GetBinErrorLow(min_bin);


					h1->Draw();

					for (int j = 1; j <= otherObj->size(); j++) {							
						TH1D * hN = arrProj[j][i];
//						int color = colorList[j];
						int color = 1;
						if (useCustomColor) color = GetCustomColor(j);
						else color = colorList[j];
						hN->SetLineColor(color);        
						hN->SetMarkerColor(color);        
						hN->SetMarkerStyle(markerList[j]);

						max_bin = hN->GetMaximumBin();
						min_bin = hN->GetMinimumBin();
						max_y = TMath::Max(max_y,hN->GetBinContent(max_bin) + hN->GetBinErrorUp(max_bin));
						min_y = TMath::Min(min_y,hN->GetBinContent(min_bin) - hN->GetBinErrorLow(min_bin));
						hN->Draw("SAME");							
						legend->AddEntry(hN,hN->GetTitle(),"lp");

					}	
					h1->GetYaxis()->SetRangeUser(min_y,max_y);
					legend->Draw("SAME");
					canvas->Print(TString::Format("%s_%d.pdf",objname.Data(),i));
					if (enableCFiles) canvas->Print(TString::Format("%s_%d.C",objname.Data(),i));
				}
//				fOut->Write();
//				fOut->Close();
			} else {
			// FIXME add in a mode where we project and compare the x-axis and y-axis 
			// separately
				vector<TH1D*> xArrProj = {};
				vector<TH1D*> yArrProj = {};
				if (proj2D) {
					printf("Projecting separately the X and Y axes for %s\n",h2obj->GetName());

					TH1D * projX = h2obj->ProjectionX(Form("%s_ProjX",h2obj->GetName()));
					TH1D * projY = h2obj->ProjectionY(Form("%s_ProjY",h2obj->GetName()));

					xArrProj.push_back(projX);
					yArrProj.push_back(projY);

					// reference: nBinsX, nBinsY.  TH2 h2obj
					for (int j = 0; j < otherObj->size(); j++) {  
						TH2F *otherH = (TH2F *) otherObj->at(j);
						objNameNoSpace = otherH->GetName();
						objNameNoSpace = TString::Format("%s_%s",objNameNoSpace.Data(),files->at(0)->GetTitle());
						cleanName(objNameNoSpace);
						otherH->SetName(objNameNoSpace);
						//	if (enableRootFiles) fOut->Add(otherH);

						TH1D * projX_o = otherH->ProjectionX(Form("%s_ProjX",otherH->GetName()));
						TH1D * projY_o = otherH->ProjectionY(Form("%s_ProjY",otherH->GetName()));

						xArrProj.push_back(projX_o);
						yArrProj.push_back(projY_o);
					}

					TLegend * legend = new TLegend(0.55,0.65,0.90,0.85);
					canvas->Clear();
					projX->SetLineColor(1);
					projX->SetMarkerColor(1);
					projX->SetMarkerStyle(markerList[0]);
					projX->DrawNormalized("",1);
					legend->AddEntry(projX,f0->GetTitle(),"lp");
					for (int i = 1; i <= otherObj->size(); i++) {
						TH1D * hX = xArrProj[i];
						//int color = colorList[i];
						int color = 1;
						if (useCustomColor) color = GetCustomColor(i);
						else color = colorList[i];
						hX->SetLineColor(color);
						hX->SetMarkerColor(color);
						hX->SetMarkerStyle(markerList[i]);
						hX->DrawNormalized("SAME",1);
						legend->AddEntry(hX,hX->GetTitle(),"lp");
					}
		
					legend->Draw("SAME");
					canvas->Print(TString::Format("%s_ProjX.pdf",objname.Data()));
					if (enableCFiles) canvas->Print(TString::Format("%s_ProjX.C",objname.Data()));
					// Now for Y
					canvas->Clear();
					legend->Clear();
			// temp fix
					projY->GetXaxis()->SetTitleSize(0.05);
					projY->GetYaxis()->SetTitleSize(0.05);

					projY->SetLineColor(1);
					projY->SetMarkerColor(1);
					projY->SetMarkerStyle(markerList[0]);
					projY->DrawNormalized("",1);
					legend->AddEntry(projY,f0->GetTitle(),"lp");
					for (int i = 1; i <= otherObj->size(); i++) {
						TH1D * hY = yArrProj[i];
					//	int color = colorList[i];
						int color = 1;
						if (useCustomColor) color = GetCustomColor(i);
						else color = colorList[i];
						hY->SetLineColor(color);
						hY->SetMarkerColor(color);
						hY->SetMarkerStyle(markerList[i]);
						hY->DrawNormalized("SAME",1);
						legend->AddEntry(hY,hY->GetTitle(),"lp");
					}
	
					legend->Draw("SAME");
					canvas->Print(TString::Format("%s_ProjY.pdf",objname.Data()));
					if (enableCFiles) canvas->Print(TString::Format("%s_ProjY.C",objname.Data()));
				} // end of projection section
		
				




			}

//      h2obj->Draw("COLZ"); 
//      canvas->Print(TString::Format("%s.pdf",obj2->GetName()));

    } else if(obj2->InheritsFrom("TProfile")) {
			TProfile * pobj = (TProfile *) obj2;
			printf("TRYING TO PRINT PROFILE\n");
			pobj->SetName(objNameNoSpace); 


      double max_y, min_y;    
			int max_bin, min_bin;

      canvas->SetLogy(0);

      pobj->SetLineColor(1);
      pobj->SetMarkerColor(1);
      pobj->SetMarkerStyle(markerList[0]);
  //    hobj->SetMinimum(0.0000000001);
      pobj->Draw("PE");
    //  max_y = pobj->GetMaximum();
   //   min_y = pobj->GetMinimum();
			max_bin = pobj->GetMaximumBin();
			min_bin = pobj->GetMinimumBin();

      max_y = pobj->GetBinContent(max_bin) + pobj->GetBinErrorUp(max_bin);
      min_y = pobj->GetBinContent(min_bin) - pobj->GetBinErrorLow(min_bin);
      legend->AddEntry(pobj,f0->GetTitle(),"lp");

	//		printf("outfile name = %s, title = %s\n",fOut->GetName(),fOut->GetTitle());
			if (enableRootFiles) fOut->Add(pobj);

      for (int j = 0; j < otherObj->size(); j++) {  
        TH1 *otherP = (TProfile *) otherObj->at(j);
				objNameNoSpace = otherP->GetName();
				objNameNoSpace = TString::Format("%s_%s",objNameNoSpace.Data(),files->at(0)->GetTitle());
				cleanName(objNameNoSpace);
				otherP->SetName(objNameNoSpace);// FIXME
				if (enableRootFiles) fOut->Add(otherP);
//        int color = colorList[j+1];
				int color = 1;
				if (useCustomColor) color = GetCustomColor(j);
				else color = colorList[j+1];
//        int color = j+2;
        otherP->SetLineColor(color);        
        otherP->SetMarkerColor(color);        
        otherP->SetMarkerStyle(markerList[j+1]);
//        otherH->SetMinimum(0.0000000001);
        otherP->Draw("SAME PE");
				max_bin = otherP->GetMaximumBin();
				min_bin = otherP->GetMinimumBin();
        max_y = TMath::Max(max_y,otherP->GetBinContent(max_bin) + otherP->GetBinErrorUp(max_bin));
        min_y = TMath::Min(min_y,otherP->GetBinContent(min_bin) - otherP->GetBinErrorLow(min_bin));
				printf("MHO: min %f, max %f\n", otherP->GetBinErrorUp(otherP->GetMaximumBin()),otherP->GetBinErrorLow(otherP->GetMinimumBin()));
		

        legend->AddEntry(otherP,otherP->GetTitle(),"lp");
      }
      // slight adjust
     // max_y += 0.46714 * (max_y - min_y);
//      max_y += 0.2 * (max_y - min_y);

      pobj->GetYaxis()->SetRangeUser(min_y,max_y);
      legend->Draw("SAME");

		//	fOut->Write();
      canvas->Print(TString::Format("%s.pdf",objname.Data()));
 //     canvas->Print(TString::Format("%s.eps",objname.Data()));
      if (enableCFiles) canvas->Print(TString::Format("%s.C",objname.Data()));
			canvas->SetLogy(0);








		} else if (obj2->InheritsFrom("TH1" )) {
      TH1 * hobj = (TH1 *) obj2;
      printf("TH1 %s\n",hobj->GetName());   

			//printf("MHO objNameNoSpace = %s\n",objNameNoSpace.Data());
			hobj->SetName(objNameNoSpace);

      //checking for specific histograms of interest
      if ( !strcmp(hobj->GetName(),"jetR_AA")) {
        hobj->GetXaxis()->SetRangeUser(10,40);
        hobj->SetTitle("Jet R_{AA}");
      }

			bool isDeltaEta = (strstr(hobj->GetName(), "DEta") != NULL || strstr(hobj->GetName(), "Deta") != NULL); 

      double max_y, min_y;    
			int max_bin, min_bin;

			double max_ratio_y=1.,min_ratio_y=1;

      canvas->SetLogy(0);

      hobj->SetLineColor(1);
      hobj->SetMarkerColor(1);
      hobj->SetMarkerStyle(markerList[0]);
  //    hobj->SetMinimum(0.0000000001);


			max_bin = hobj->GetMaximumBin();
			min_bin = hobj->GetMinimumBin();
      max_y = hobj->GetBinContent(max_bin) + hobj->GetBinErrorUp(max_bin);
      min_y = hobj->GetBinContent(min_bin) - hobj->GetBinErrorLow(min_bin);


      double maxEntry = hobj->GetBinContent(hobj->GetMaximumBin());
      double max_rms = hobj->GetStdDev();
      double xRange = hobj->GetXaxis()->GetXmax() - hobj->GetXaxis()->GetXmin();

			bool bDivisionPossible = false;

      if (ratioMode) {
				// FIXME check that division is actually possible
				int nBinsFirst = hobj->GetNbinsX();
				for (int j = 0; j < otherObj->size(); j++) {
					TH1 * otherH = (TH1 *) otherObj->at(j);
					if (nBinsFirst == otherH->GetNbinsX()) bDivisionPossible = true;
				}
				if (bDivisionPossible) {
					canvas->Divide(1,2,small,small);
					canvas->cd(1);
					gPad->SetBottomMargin(small);
				}
      }

      hobj->Draw("PE");
      legend->AddEntry(hobj,f0->GetTitle(),"lp");
			if (enableRootFiles) fOut->Add(hobj);

			std::vector<TH1 *> RatioArray = {} ;

      for (int j = 0; j < otherObj->size(); j++) {  
        TH1 *otherH = (TH1 *) otherObj->at(j);
				objNameNoSpace = otherH->GetName();
				objNameNoSpace = TString::Format("%s_%s",objNameNoSpace.Data(),files->at(0)->GetTitle());
				cleanName(objNameNoSpace);
				otherH->SetName(objNameNoSpace);
				if (enableRootFiles) fOut->Add(otherH);
//        int color = colorList[j+1];
				int color = 1;
				if (useCustomColor) color = GetCustomColor(j);
				else color = colorList[j+1];
        otherH->SetLineColor(color);        
        otherH->SetMarkerColor(color);        
        otherH->SetMarkerStyle(markerList[j+1]);
		
				// Ratio
				if (ratioMode && bDivisionPossible ) {
					TH1 * LocalRatio = (TH1 *) otherH->Clone(Form("%s_Ratio",otherH->GetName()));
					LocalRatio->Divide(hobj);
					double localRatioMax = LocalRatio->GetBinContent(LocalRatio->GetMaximumBin());
					double localRatioMin = LocalRatio->GetBinContent(LocalRatio->GetMinimumBin());
					if (localRatioMax > max_ratio_y) max_ratio_y = localRatioMax;
					if (localRatioMin < min_ratio_y) min_ratio_y = localRatioMin;

					RatioArray.push_back(LocalRatio);
				}

        otherH->Draw("SAME PE");
   //     max_y = TMath::Max(max_y,otherH->GetMaximum());
   //     min_y = TMath::Min(min_y,otherH->GetMinimum());
				max_bin = otherH->GetMaximumBin();
				min_bin = otherH->GetMinimumBin();
        max_y = TMath::Max(max_y,otherH->GetBinContent(max_bin) + otherH->GetBinErrorUp(max_bin));
        min_y = TMath::Min(min_y,otherH->GetBinContent(min_bin) - otherH->GetBinErrorLow(min_bin));
				max_rms = TMath::Max(max_rms,otherH->GetStdDev());

        legend->AddEntry(otherH,otherH->GetTitle(),"lp");
      }
      // slight adjust
      max_y += 0.2 * (max_y - min_y);

      hobj->GetYaxis()->SetRangeUser(min_y,max_y);
	    if (max_rms < 0.005*xRange) { 
		//		hobj->GetYaxis()->SetRangeUser()
      	hobj->GetYaxis()->SetRangeUser(max_y*1e-9,max_y);
				canvas->SetLogy(1);
			}
			if (isDeltaEta) {
				TLine * line1 = new TLine(-1.3,min_y,-1.3,max_y);
				TLine * line2 = new TLine(1.3,min_y,1.3,max_y);

				line1->SetLineStyle(LineStyle);
				line1->SetLineWidth(LineWidth);
				line1->SetLineColor(LineColor);

				line2->SetLineStyle(LineStyle);
				line2->SetLineWidth(LineWidth);
				line2->SetLineColor(LineColor);
	
				line1->Draw("SAME");
				line2->Draw("SAME");
			}
      legend->Draw("SAME");

			if (ratioMode && bDivisionPossible) {
				canvas->cd(2);
				gPad->SetTopMargin(small);
				gPad->SetTickx();
				gPad->SetBottomMargin(0.2);

				// Increasing range via magic
				if (max_ratio_y > 1./min_ratio_y) min_ratio_y = 1./max_ratio_y;
				else max_ratio_y = 1./min_ratio_y;
			
				// If ratios exploded, try a default range
				if (max_ratio_y > 30) {
					max_ratio_y = 1.3;
					min_ratio_y = 1.0/1.3;
				}

				min_ratio_y = TMath::Power(min_ratio_y,1.6);
				max_ratio_y = TMath::Power(max_ratio_y,1.6);

				RatioArray[0]->SetTitle("");
				RatioArray[0]->GetYaxis()->SetTitle(Form("Ratio Over %s",f0->GetTitle()));
				RatioArray[0]->GetYaxis()->SetTitleOffset(0.7);
				RatioArray[0]->GetYaxis()->SetRangeUser(min_ratio_y,max_ratio_y);
				RatioArray[0]->Draw("P");

				TF1 * unityFunction = new TF1("unity","1",hobj->GetXaxis()->GetXmin(),hobj->GetXaxis()->GetXmax());
				unityFunction->SetLineStyle(2);
				unityFunction->SetLineWidth(4);
				unityFunction->SetLineColor(kGray);	
				unityFunction->Draw("SAME");

				for (Int_t i = 1; i < otherObj->size(); i++) {
					RatioArray[i]->Draw("SAME P");
				}

				if (isDeltaEta) {
					TLine * line3 = new TLine(-1.3,min_ratio_y,-1.3,max_ratio_y);
					TLine * line4 = new TLine(1.3,min_ratio_y,1.3,max_ratio_y);

					line3->SetLineStyle(LineStyle);
					line3->SetLineWidth(LineWidth);
					line3->SetLineColor(LineColor);

					line4->SetLineStyle(LineStyle);
					line4->SetLineWidth(LineWidth);
					line4->SetLineColor(LineColor);
		
					line3->Draw("SAME");
					line4->Draw("SAME");
				}

			}

//      canvas->SetLogy(1);
    // experiment: check if logY is appropriate

//			fOut->Write();
//			fOut->Close();
//			delete fOut;
 //     canvas->Print(TString::Format("%s.eps",objname.Data()));
      canvas->Print(TString::Format("%s.pdf",objname.Data()));
      if (enableCFiles) canvas->Print(TString::Format("%s.C",objname.Data()));
     // canvas->Print(TString::Format("%s.C",objNameNoSpace.Data()));
//      canvas->Print(TString::Format("%s.pdf",hobj->GetName()));
//      canvas->Print(TString::Format("%s.C",hobj->GetName()));
 //     canvas->Print(TString::Format("%s.pdf",obj2->GetName()));
//      canvas->Print(TString::Format("%s.C",obj2->GetName()));
			canvas->SetLogy(0);
 //     canvas->Print(TString::Format("%s.root",obj2->GetName()));
    //  canvas->Print(TString::Format("%s.C",obj2->GetName()));
    } else if (obj2->InheritsFrom("TGraphErrors")) {  // FIXME nothing done for ordinary TGraphs
      TGraph * gobj = (TGraph *) obj2;
      printf("TGraphErrors %s\n",gobj->GetName());
			TString gobjStringClean = gobj->GetName();
			cleanName(gobjStringClean);
			gobj->SetName(gobjStringClean);
      double max_y, min_y;
      gobj->SetLineColor(1);  
			gobj->SetLineStyle(1);
      gobj->SetMarkerColor(1);
      gobj->SetMarkerStyle(markerList[0]);

			if (enableRootFiles) fOut->Add(gobj);

      TMultiGraph *mg = new TMultiGraph();
      mg->Add(gobj);
      TMultiGraph *ratioMG = new TMultiGraph();
			
      std::vector<TGraphErrors> ratioArray;
    
      if (ratioMode) {
//        canvas->Divide(1,2,0.02,0.00);
        canvas->Divide(1,2,small,small);
        canvas->cd(1);
				gPad->SetBottomMargin(small);
      }

     // gobj->Draw("ALP");
      max_y = gobj->GetMaximum();
      min_y = gobj->GetMinimum();
      legend->AddEntry(gobj,f0->GetTitle(),"lp");
      for (int j = 0; j < otherObj->size(); j++) {  
        TGraph *otherG = (TGraph *) otherObj->at(j);
        
        // Creating ratios
        TGraphErrors * localRatio = DivideTGraphErrors((TGraphErrors *) otherG, (TGraphErrors *) gobj, Form("%s_ratio",otherG->GetName()));

        //int color = colorList[j+1];
				int color = 1;
				if (useCustomColor) color = GetCustomColor(j);
				else color = colorList[j+1];
        otherG->SetLineColor(color);        
				otherG->SetLineStyle(1);
        otherG->SetMarkerColor(color);        
        otherG->SetMarkerStyle(markerList[j+1]);
  //      otherG->Draw("SAME LP");
				if (enableRootFiles) fOut->Add(otherG);
        mg->Add(otherG);

        localRatio->SetLineColor(color);
				localRatio->SetLineStyle(1);
        localRatio->SetMarkerColor(color);        
        localRatio->SetMarkerStyle(markerList[j+1]);
  
        ratioMG->Add(localRatio);

		
        max_y = TMath::Max(max_y,otherG->GetMaximum());
        min_y = TMath::Min(min_y,otherG->GetMinimum());

        legend->AddEntry(otherG,otherG->GetTitle(),"lp");
      }
      // slight adjust
      max_y += 0.46714 * (max_y - min_y);
    //  max_y += 0.3 * (max_y - min_y);
      //not sure if this is working
	
			if (objname.Contains("OverIn")) {
				min_y = 0.5;
				max_y = 1.5;
			} else if (objname.Contains("Minus")) {
				min_y = -0.8;
				max_y = 1.45;
			}

//      gobj->SetMinimum(min_y);
//      gobj->SetMaximum(max_y);

 //     mg->GetYaxis()->SetRangeUser(min_y,max_y);
 //     mg->SetMaximum(0.6);
  
      mg->Draw("ALP");
      mg->GetXaxis()->SetTitle(gobj->GetXaxis()->GetTitle());
      mg->GetYaxis()->SetTitle(gobj->GetYaxis()->GetTitle());
			mg->GetYaxis()->SetRangeUser(min_y,max_y);

			if (objname.Contains("OverIn")) {
				TF1 * unityFunction = new TF1("unity","1",0,15);
				unityFunction->SetLineStyle(LineStyle);
				unityFunction->SetLineWidth(LineWidth);
				unityFunction->SetLineColor(LineColor);	
				unityFunction->Draw("SAME");
			} else if (objname.Contains("Minus")) {
				TF1 * zeroFunction = new TF1("zero","0",0,15);
				zeroFunction->SetLineStyle(LineStyle);
				zeroFunction->SetLineWidth(LineWidth);
				zeroFunction->SetLineColor(LineColor);	
				zeroFunction->Draw("SAME");

			}
	
			if (objname.Contains("Spectrum")) {
        if (!ratioMode) canvas->cd(1);
				canvas->SetLogy();
			}

      legend->Draw("SAME");

      if (ratioMode) {
        canvas->cd(2);
				gPad->SetTopMargin(small);
				gPad->SetTickx();
				gPad->SetBottomMargin(0.2);
        ratioMG->Draw("ALP");
      //  ratioMG->GetYaxis()->SetRangeUser(0,2);
				ratioMG->GetYaxis()->SetTitle(Form("Ratio Over %s",f0->GetTitle()));
        ratioMG->GetXaxis()->SetTitle(gobj->GetXaxis()->GetTitle());
      }


      canvas->Print(TString::Format("%s.pdf",gobj->GetName()));
 //     canvas->Print(TString::Format("%s.eps",gobj->GetName()));
     // canvas->Print(TString::Format("%s.C",gobj->GetName()));
     // canvas->Print(TString::Format("%s.C",objNameNoSpace.Data()));
      if (enableCFiles) canvas->Print(TString::Format("%s.C",gobjStringClean.Data()));
//      canvas->Print(TString::Format("%s.pdf",obj2->GetName()));
//      canvas->Print(TString::Format("%s.C",obj2->GetName()));
			canvas->SetLogy(0);
    } else if (obj2->InheritsFrom("TF1")) {

      TF1 * fobj = (TF1 *) obj2;
      printf("TF1 %s\n",fobj->GetName());
			if (useCustomColor) fobj->SetLineColor(GetCustomColor(0));
			else fobj->SetLineColor(colorList[0]);
			fobj->Draw();

      if (enableRootFiles) fOut->Add(fobj);

      legend->AddEntry(fobj,f0->GetTitle(),"lp");
      for (int j = 0; j < otherObj->size(); j++) {  
        TF1 *otherG = (TF1 *) otherObj->at(j);
				int color = 1;
				if (useCustomColor) color = GetCustomColor(j);
				else color = colorList[j];
        otherG->SetLineColor(color);        
//        otherG->SetMarkerColor(color);        
//        otherG->SetMarkerStyle(22);
  //      otherG->Draw("SAME LP");
				otherG->Draw("SAME");
  //      mg->Add(otherG);
      //  max_y = TMath::Max(max_y,otherG->GetMaximum());
     //   min_y = TMath::Min(min_y,otherG->GetMinimum());

        if (enableRootFiles) fOut->Add(otherG);

      //  legend->AddEntry(otherG,otherG->GetTitle(),"l"); // FIXME
        legend->AddEntry(otherG,otherG->GetYaxis()->GetTitle(),"l");
      }


			legend->Draw("SAME");
	//		canvas->SetLogy(1);
      canvas->Print(TString::Format("%s.pdf",fobj->GetName()));
      if (enableCFiles) canvas->Print(TString::Format("%s.C",fobj->GetName()));
 //     canvas->Print(TString::Format("%s.eps",fobj->GetName()));
  //    canvas->Print(TString::Format("%s.pdf",obj2->GetName()));
		}
		if (enableRootFiles) {
			fOut->Write();
			fOut->Close();
		}
	//	delete fOut;
    legend->Clear();
    canvas->Clear();
  }
}

int main(int argc, char * argv[]) {

  if (argc == 1) {
    printf("Usage: %s <OutputDirectory> [Filepath_1] [Label_1] [Filepath_2] [Label_2] ... \n",argv[0]);
//    printf("    Or %s -p [PATTERN] [Filepath_1] [Label_1] [Filepath_2] [Label_2] ... \n",argv[0]);
    exit(0);
  }


	
  int num = argc/ 2;
  vector<fileLabel> *fileLabels = new vector<fileLabel>;
  for (int i = 1; i < num; i++) {
    fileLabels->push_back(fileLabel(TString(argv[2*i]),TString(argv[2*i+1])));
  }
/*a  int num = (argc - 1) / 2;
  vector<fileLabel> *fileLabels = new vector<fileLabel>;
  for (int i = 0; i < num; i++) {
    fileLabels->push_back(fileLabel(TString(argv[2*i + 1]),TString(argv[2*(i+1)])));
  }*/
  compare(fileLabels,argv[1]);
}




