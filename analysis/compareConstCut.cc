
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
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TStyle.h>
#include <TColor.h>


#include "analysis_params.h"

using namespace std;


struct file_cut {
  double c_value;
  TFile *file;
};

bool file_cut_cmp (file_cut f1, file_cut f2) {
  return f1.c_value < f2.c_value;
}

vector<vector<double> > transpose(vector<vector<double> > inMat) {
  vector<vector<double> > outMat(inMat[0].size(),vector<double>(inMat.size()));
  for (size_t i = 0; i < inMat.size(); ++i) 
    for (size_t j = 0; j < inMat[0].size(); ++j)
      outMat[j][i] = inMat[i][j];
  return outMat;
}

TH2F * normalizeHistogramColumns(TH2F *inHist) {

  TH2F * normHist = (TH2F *) inHist->Clone();
  normHist->SetName(Form("%sNorm",inHist->GetName()));
  normHist->SetTitle(Form("%s (Normalized)",inHist->GetTitle()));
  TH2F * scaleHist = (TH2F *) inHist->Clone();
  scaleHist->SetName(Form("%sScale",inHist->GetName()));
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




vector<vector<TH1F *> > transpose(vector<vector<TH1F *> > inMat) {
  vector<vector<TH1F *> > outMat(inMat[0].size(),vector<TH1F *>(inMat.size()));
//  vector<vector<TH1 *> > outMat(inMat[0].size(),vector<TH1 *>(inMat.size()));
  for (size_t i = 0; i < inMat.size(); ++i) 
    for (size_t j = 0; j < inMat[0].size(); ++j)
      outMat[j][i] = inMat[i][j];
  return outMat;
}



void compareConstCut(TString AA_final_pattern = "output/PbPb_2760GeV/PbPb_W_2760GeV_CENT_0_5_1_CCUT_") {

// outline:
/**
 input a pattern for pp files "output/pp_W_200GeV_CCUT_%d.final.root"
 input a pattern for AA files "output/AuAu_W_200GeV_CENT_0_5_CCUT_%d.final.root"
  For each histogram of interest, create vectors of histograms.
    Add %d to the title?  or create array of %d?
    vertexXRotJetPt
    vertexXRotJetPtHardCore
*/

  gStyle->SetOptStat(0);
  set_plot_style();


  double c_width = 600;
  double c_height = 600;

  vector<double> ccut;
  vector<file_cut> file_cuts;

  // Get the directory from the path
  int index_slash = AA_final_pattern.Last('/');
  TString dirname = AA_final_pattern(0,index_slash+1);
  TString search_pattern = AA_final_pattern(index_slash+1,AA_final_pattern.Length()-index_slash-1);

  printf("%s\t\t%s\n",dirname.Data(),search_pattern.Data());

  TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (!file->IsDirectory() && fname.BeginsWith(search_pattern)) {
            cout << fname.Data() << endl;
            //getting the list of ccut values
            TString sub1 = fname(TRegexp("CCUT_[^.]*"),0);   
            sub1 = sub1(5,sub1.Length() - 5);
            double ccut_value = sub1.Atof();
            ccut.push_back(ccut_value);
            TString filename = dirname + fname;

            TFile *f = new TFile(filename.Data(),"READ","");
            file_cut fc = {.c_value = ccut_value, .file = f};
            file_cuts.push_back(fc);
         }
      }
   }

  std::sort(ccut.begin(),ccut.end());
  std::sort(file_cuts.begin(),file_cuts.end(),file_cut_cmp);

  for (vector<double>::iterator it=ccut.begin(); it!=ccut.end(); ++it) {
    cout << Form("%f",*it) << endl;
  }  

  vector<TH2F> vertexXRotJetPtArray;
  vector<TH2F> vertexXRotJetPtHardCoreArray;


// In the beginning  
//   first index: constCutPtBin
//   second index: jetPtBin
// Then transpose the matrix before creating TGraphs
  vector<vector<double> > meanX;
  vector<vector<double> > meanX_err;
  vector<vector<double> > meanXHardCore;
  vector<vector<double> > meanXHardCore_err;
  vector<vector<double> > s_values;
  vector<vector<double> > s_values_err;
  vector<vector<double> > sHardCore_values;
  vector<vector<double> > sHardCore_values_err;


  TCanvas *c = new TCanvas("c","c",c_width,c_height);
  // Go through the file_cuts
  for (vector<file_cut>::iterator it=file_cuts.begin(); it!=file_cuts.end(); ++it) {
    cout << Form("%s",it->file->GetName()) << endl;

    TH2F * vertexXRotJetPt = (TH2F *) it->file->Get("vertexXRotJetPt");
    TH2F * vertexXRotJetPtHardCore = (TH2F *) it->file->Get("vertexXRotJetPtHardCore");
  
    if (!vertexXRotJetPt || !vertexXRotJetPtHardCore) {
      printf("Error: vertexXRotJetPt or vertexXRotJetPtHardCore missing.\n");
      exit(1);
    }
    vertexXRotJetPt->Draw("COLZ");
    c->Print(Form("output/ccut_vertexXRotJetPt_ccut_%.0f.pdf",it->c_value));
    c->Clear();
    vertexXRotJetPtHardCore->Draw("COLZ");
    c->Print(Form("output/ccut_vertexXRotJetPt_HardCore_ccut_%.0f.pdf",it->c_value));

    vector<double> meanX_jetPtBin;
    vector<double> meanX_jetPtBin_err;
    vector<double> meanXHardCore_jetPtBin;
    vector<double> meanXHardCore_jetPtBin_err;
    vector<double> s_jetPtBin;
    vector<double> s_jetPtBin_err;
    vector<double> sHardCore_jetPtBin;
    vector<double> sHardCore_jetPtBin_err;

    c->Clear();
    TLegend *leg = new TLegend(0.5,0.65,0.9,0.85);

    for (int i = 0; i < jetPtBins.size()-1; i++) {
      TH1D * proj = vertexXRotJetPt->ProjectionY(Form("%s_px_%d",vertexXRotJetPt->GetName(),i),i+1,i+1);
      TH1D * projHardCore = vertexXRotJetPtHardCore->ProjectionY(Form("%s_px_%d",vertexXRotJetPtHardCore->GetName(),i),i+1,i+1);

// normalize for TH1 plots
      TH1D *proj_clone = (TH1D *) proj->Clone();
      proj_clone->Scale(1./proj_clone->Integral());
      proj_clone->SetLineColor(i+1);
      if(i) proj_clone->Draw("SAME"); else proj_clone->Draw();
      leg->AddEntry(proj_clone,Form("%.1f < p^{jet}_{t} < %.1f",jetPtBins.at(i),jetPtBins.at(i+1)),"lp");

      //calculating means, s_values 
      double x_m = proj->GetMean(1);
      double x_m_err = proj->GetMeanError(1); 
    //  double x_m_err = proj->GetStdDev(1); 
      Double_t nearside_err, awayside_err;
      double nearside = proj->IntegralAndError(0,proj->GetXaxis()->FindBin(0.), nearside_err);
     // double nearside = proj->IntegralAndError(0,proj->GetXaxis()->FindBin(0.),(Double_t&) &nearside_err);
      double awayside = proj->IntegralAndError(proj->GetXaxis()->FindBin(0.),proj->GetNbinsX()+1,awayside_err);
      double s = awayside ? nearside / awayside : 0;
      double s_err = awayside ? s * TMath::Sqrt(TMath::Power(nearside_err/nearside,2) + TMath::Power(awayside_err/awayside,2)) : 0;
      double x_m_HardCore = projHardCore->GetMean(1);      
   //   double x_m_HardCore_err = projHardCore->GetStdDev(1);    
      double x_m_HardCore_err = projHardCore->GetMeanError(1);    
      nearside = projHardCore->IntegralAndError(0,projHardCore->GetXaxis()->FindBin(0.),nearside_err);
      awayside = projHardCore->IntegralAndError(projHardCore->GetXaxis()->FindBin(0.),projHardCore->GetNbinsX()+1,awayside_err);
      double s_HardCore  = awayside ? nearside / awayside : 0;
      double s_HardCore_err = awayside ? s_HardCore * TMath::Sqrt(TMath::Power(nearside_err/nearside,2) + TMath::Power(awayside_err/awayside,2)) : 0;


      printf("<x> = %.3f \\pm %.3f \t s = %.3f\t HC: <x> =%.3f \\pm %.3f \t s = %.3f\n",x_m,x_m_err,s,x_m_HardCore,x_m_HardCore_err,s_HardCore);
      
      meanX_jetPtBin.push_back(x_m);
      meanX_jetPtBin_err.push_back(x_m_err);
      s_jetPtBin.push_back(s);
      s_jetPtBin_err.push_back(s_err);
      meanXHardCore_jetPtBin.push_back(x_m_HardCore);
      meanXHardCore_jetPtBin_err.push_back( x_m_HardCore_err);
      sHardCore_jetPtBin.push_back(s_HardCore);
      sHardCore_jetPtBin_err.push_back(s_HardCore_err);

    } 
    meanX.push_back(meanX_jetPtBin);
    meanX_err.push_back(meanX_jetPtBin_err);
    s_values.push_back(s_jetPtBin);
    s_values_err.push_back(s_jetPtBin_err);
    meanXHardCore.push_back(meanXHardCore_jetPtBin);
    meanXHardCore_err.push_back(meanXHardCore_jetPtBin_err);
    sHardCore_values.push_back(sHardCore_jetPtBin);
    sHardCore_values_err.push_back(sHardCore_jetPtBin_err);
    
    leg->Draw("SAME");
    c->Print(Form("output/ccut_vertexXRot_JetPt_c_%.0f.pdf",it->c_value));
  }
 
 
  meanX = transpose(meanX);
  meanX_err = transpose(meanX_err);
  meanXHardCore = transpose(meanXHardCore);
  meanXHardCore_err = transpose(meanXHardCore_err);
  s_values = transpose(s_values);
  s_values_err = transpose(s_values_err);
  sHardCore_values = transpose(sHardCore_values);
  sHardCore_values_err = transpose(sHardCore_values_err);
//   first index: jetPtBin
//   second index: constCutPtBin

  
//  for (vector<double>::iterator jpit=jetPtBins.begin(); jpit!=jetPtBins.end(); ++jpit) {
  TLegend *legMeanX = new TLegend(0.65,0.35,0.95,0.45);
  TLegend *legSValue = new TLegend(0.65,0.45,0.95,0.65);
  TCanvas *c_MeanX_all = new TCanvas("c_MeanX_all","c_MeanX_all",c_width,c_height);
  TCanvas *c_SValue_all = new TCanvas("c_SValue_all","c_SValue_all",c_width,c_height);
  c_MeanX_all->Divide(2,2);
  c_SValue_all->Divide(2,2);


  TGraphErrors *gMeanX, *gMeanXHardCore, *gSValue, *gSValueHardCore;

  TFile * fOut = new TFile("output/ccut.root","RECREATE");
	// FIXME 

	float meanX_max_y = 0;
	float meanX_min_y = -2;

	float sValue_max_y = 6;
	float sValue_min_y = 0;

  for (int i = 0; i < jetPtBins.size() - 1; i++) {
    c->Clear(); 
    legMeanX->Clear();
    c->cd();

   // TGraph *gMeanX = new TGraph((int) ccut.size(),&ccut[0],&meanX[i][0]);
    gMeanX = new TGraphErrors((int) ccut.size(),&ccut[0],&meanX[i][0],0,&meanX_err[i][0]);
    gMeanX->SetName(Form("x_value_AA_nohc_jetPt_%.0f_%.0f",jetPtBins.at(i),jetPtBins.at(i+1)));
    gMeanX->SetTitle("Mean X_{vertex,rot} vs Constituent Cut");
    gMeanX->GetXaxis()->SetTitle("Const Cut (GeV/c)");      
    gMeanX->GetYaxis()->SetTitle("<X>");      
    gMeanX->SetMarkerStyle(25);

   // TGraph *gMeanXHardCore = new TGraph((int) ccut.size(),&ccut[0],&meanXHardCore[i][0]);
    gMeanXHardCore = new TGraphErrors((int) ccut.size(),&ccut[0],&meanXHardCore[i][0],0,&meanXHardCore_err[i][0]);
 //   TGraphErrors *gMeanXHardCore = new TGraphErrors((int) ccut.size(),&ccut[0],&meanXHardCore[i][0],0,0);
    gMeanXHardCore->SetTitle(Form("Mean X_{vertex,rot} vs Constituent Cut (%.0f #leq p_{T}^{jet} < %.0f)",jetPtBins.at(i),jetPtBins.at(i+1)));
    gMeanXHardCore->SetName(Form("x_value_AA_hc_jetPt_%.0f_%.0f",jetPtBins.at(i),jetPtBins.at(i+1)));
    gMeanXHardCore->GetXaxis()->SetTitle("Const Cut (GeV/c)");      
    gMeanXHardCore->GetYaxis()->SetTitle("<X>");      
    gMeanXHardCore->SetMarkerStyle(21);
    gMeanXHardCore->SetMarkerColor(kRed);
    gMeanXHardCore->SetLineColor(kRed);


    gMeanXHardCore->Draw("ALP");
    gMeanX->Draw("LP SAME");
    legMeanX->AddEntry(gMeanX,"No Core Cut","lp");  
    legMeanX->AddEntry(gMeanXHardCore,"6 GeV/c Core Cut","lp");  

    float max_y = max(gMeanX->GetYaxis()->GetXmax(),gMeanXHardCore->GetYaxis()->GetXmax());     
    float min_y = min(gMeanX->GetYaxis()->GetXmin(),gMeanXHardCore->GetYaxis()->GetXmin());     


    //gMeanX->GetYaxis()->SetRangeUser(min_y-0.1,max_y+0.1);
    gMeanXHardCore->GetYaxis()->SetRangeUser(meanX_min_y,meanX_max_y);
    legMeanX->Draw("SAME");

 //   c->Print(Form("output/ccut_JetPt_%.0f_%.0f.pdf",*jpit,*(jpit+1)));
    c->Print(Form("output/ccut_MeanX_JetPt_%.0f_%.0f.pdf",jetPtBins.at(i),jetPtBins.at(i+1)));
    c->Clear();
    // Add to plot with all of them
    c_MeanX_all->cd(i+1);
    gMeanXHardCore->Draw("ALP");
    gMeanX->Draw("LP SAME");
    legMeanX->Draw("SAME");
  //will this work???
    // S value
    legSValue->Clear();
    c->cd();    

    gSValue = new TGraphErrors((int) ccut.size(),&ccut[0],&s_values[i][0],0,&s_values_err[i][0]);
    gSValue->SetName(Form("s_value_AA_nohc_jetPt_%.0f_%.0f",jetPtBins.at(i),jetPtBins.at(i+1)));
    gSValue->SetTitle("S Value vs Constituent Cut");
    gSValue->GetXaxis()->SetTitle("Const Cut (GeV/c)");
    gSValue->GetYaxis()->SetTitle("s");
    gSValue->SetMarkerStyle(25);
    
    gSValueHardCore = new TGraphErrors((int) ccut.size(),&ccut[0],&sHardCore_values[i][0],0,&sHardCore_values_err[i][0]);
    gSValueHardCore->SetTitle(Form("S Value vs Constituent Cut (%.0f #leq p_{T}^{jet} < %.0f)",jetPtBins.at(i),jetPtBins.at(i+1)));
    gSValueHardCore->SetName(Form("s_value_AA_hc_jetPt_%.0f_%.0f",jetPtBins.at(i),jetPtBins.at(i+1)));
    gSValueHardCore->GetXaxis()->SetTitle("Const Cut (GeV/c)");
    gSValueHardCore->GetYaxis()->SetTitle("s");
    gSValueHardCore->SetLineColor(kRed);
    gSValueHardCore->SetMarkerColor(kRed);
    gSValueHardCore->SetMarkerStyle(21);
    


    gSValueHardCore->Draw("ALP");
    gSValue->Draw("LP SAME");
    legSValue->AddEntry(gSValue,"No Core Cut","lp");
    legSValue->AddEntry(gSValueHardCore,"6 GeV/c Core Cut","lp");
    
    max_y = max(gSValue->GetYaxis()->GetXmax(),gSValueHardCore->GetYaxis()->GetXmax());     
    min_y = min(gSValue->GetYaxis()->GetXmin(),gSValueHardCore->GetYaxis()->GetXmin());     

  //  gSValue->GetYaxis()->SetRangeUser(min_y-0.1,max_y+0.1);
    gSValueHardCore->GetYaxis()->SetRangeUser(sValue_min_y,sValue_max_y);
    legSValue->Draw("SAME");
    c->Print(Form("output/ccut_S_JetPt_%.0f_%.0f.pdf",jetPtBins.at(i),jetPtBins.at(i+1)));
    c->Clear();
    c_SValue_all->cd(i+1);
    gSValueHardCore->Draw("ALP");
    gSValue->Draw("LP SAME");
    legSValue->Draw("SAME");
    
    fOut->Add(gMeanX);  
    fOut->Add(gMeanXHardCore);  
    fOut->Add(gSValue);
    fOut->Add(gSValueHardCore);
    
  }

  c_MeanX_all->Print("output/ccut_MeanX_all.pdf");
  c_SValue_all->Print("output/ccut_SValue_all.pdf");

  // Parton Pt by JetBin, type
  cout << "Starting parton by pt part" << endl;
  
  const int numType = 5;
  string typeArr[numType] = {"qqbar","qq","gq","qg","gg"};
  string titleArr[numType] = {"q#bar{q}","qq","gq","qg","gg"};
  

  // Making appropriate consituent cut bins 
  vector<double > ccutBins;
  for (int i = 0; i < ccut.size(); i++) {
    ccutBins.push_back(ccut.at(i));
  }
  ccutBins.push_back(1.+ccut.at(ccut.size()-1));


  TCanvas *cParton = new TCanvas("cParton","cParton",c_width,c_height);
 
  bool partonPtByJetBin_absent = false;

  vector<vector<TH1F *> > allTypes_parton_PartonPtByJetBin_BinArr; 
  for (vector<file_cut>::iterator it=file_cuts.begin(); it!=file_cuts.end(); ++it) {
    cout << Form("%s",it->file->GetName()) << endl;

  // Differentially by hard scatter type
  // Create combined PartonPt along the way.
  for (int z = 0; z < numType; z++){

    vector<vector<TH1F *> > parton_PartonPtByJetBin_BinArr; 
    for (vector<file_cut>::iterator it=file_cuts.begin(); it!=file_cuts.end(); ++it) {
      cout << Form("%s",it->file->GetName()) << endl;

      vector<TH1F *> parton_PartonPtByJetBin;
      for (int i = 0; i < jetPtBins.size()-1; i++) {
        TH1F * temp_parton = (TH1F *) it->file->Get(Form("%s_PartonPtByJetBin_EP_0_jetPt_%.0f_%.0f",typeArr[z].c_str(),jetPtBins.at(i),jetPtBins.at(i+1)));   
        if (!temp_parton) {
          fprintf(stderr,"Error: %s not found\n",Form("%s_PartonPtByJetBin_EP_0_jetPt_%.0f_%.0f",typeArr[z].c_str(),jetPtBins.at(i),jetPtBins.at(i+1)));
          partonPtByJetBin_absent = true;
          break;
        }
        temp_parton->SetName(Form("%s_ccut_%.1f",temp_parton->GetName(),it->c_value));
        parton_PartonPtByJetBin.push_back(temp_parton);
      }
      parton_PartonPtByJetBin_BinArr.push_back(parton_PartonPtByJetBin);
    }
    
    // Transpose, such that the first index is jet pt bin, and the second is 
    // constituent cut
  //  parton_PartonPtByJetBin_BinArr = (vector<vector<TH1F *> >) transpose((vector<vector<TH1 *> >) parton_PartonPtByJetBin_BinArr);





    if(!partonPtByJetBin_absent) {
      
      parton_PartonPtByJetBin_BinArr = transpose(parton_PartonPtByJetBin_BinArr);
    
      if (parton_PartonPtByJetBin_BinArr.size() != jetPtBins.size() - 1) {
        fprintf(stderr,"There was an error in the transposition!\n");
        fprintf(stderr,"%d != %d\n",(int) parton_PartonPtByJetBin_BinArr.size(),(int) jetPtBins.size());
        exit(1);
      }
      for (int i = 0; i < parton_PartonPtByJetBin_BinArr.size(); i++) {
        vector<TH1F *> thisRow = parton_PartonPtByJetBin_BinArr.at(i);
        int nBinsX = thisRow.size();
        // FIXME need ccut 'bins' 0,1,2,3,4,6   + 1?
        // 0-1 is 0, 1-2 is 1, ... 4-6 is 4, 6-7 is 6???
        int nBinsY = thisRow.at(i)->GetNbinsX();  
        double minY = thisRow.at(i)->GetXaxis()->GetXmin();
        double maxY = thisRow.at(i)->GetXaxis()->GetXmax();

        TString title = Form("%s_jetPt_%.0f_%.0f",typeArr[z].c_str(),jetPtBins.at(i),jetPtBins.at(i+1));
        TH2F * partonPtVsCCut = new TH2F(title,Form("%s, %.0f < p_{T}^{jet} < %.0f;Const Cut (GeV/c);p_{T}^{parton} (GeV/c)",titleArr[z].c_str(),jetPtBins.at(i),jetPtBins.at(i+1)),nBinsX,&ccutBins[0],nBinsY,minY,maxY);

        // remember: for x-axis, want to avoid  underflow bin on ccut axis
        // bin = 0  -> underflow
      
        // bin = nbins + 1 -> overflow

        for (int j = 0; j < nBinsX; j++) {  //cycle over ccuts 
          for (int k = 0; k < thisRow.at(i)->GetNbinsX(); k++) { // cycle over parton pt
            partonPtVsCCut->SetBinContent(j+1,k,thisRow.at(j)->GetBinContent(k));
          }
        }
        partonPtVsCCut = normalizeHistogramColumns(partonPtVsCCut);
        partonPtVsCCut->Draw("COLZ");
        cParton->Print(Form("output/%s.pdf",partonPtVsCCut->GetName()));
        cParton->Print(Form("output/%s.C",partonPtVsCCut->GetName()));

      }
    }
  }
  }



  fOut->Write();

}

int main(int argc, char * argv[]) {

  if (argc == 1) {
    printf("Usage: %s root/path/title_CCUT_\n",argv[0]);
    exit(0);
  }

  compareConstCut(TString(argv[1]));
  
}




