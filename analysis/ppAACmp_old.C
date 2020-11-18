//ppAACmp(TString pp_final_name = "output/pp_final.root", TString AA_final_name = "output/AA_final.root") {


ppAACmp(TString pp_final_name = "output/pp.final.root", TString AA_final_name = "output/AuAu.final.root", TString outFile_path = "output/cmp/cmp.root") {

	bool do_JetPt_things = true;

  double c_width = 600;
  double c_height = 600;


	// I don't know how to include the following straight out of the analysis_params.h file yet
	// ---------------------------------------------------------------------------------
  static enum { nJetPtBins = 4, nParticlePtBins = 11 };
	double jetPtBins[nJetPtBins+1] = {10,15,20,40,60};
	double particlePtBins[nParticlePtBins+1] = {0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 17};

	double particle_pt_min = 0.2;
	double jet_constituent_cut = 6.0;
	double eta_cut = 1.0;

// ---------------------------------------------------------------------------------


  printf("Creating comparision plots using %s for pp and %s for AA.\n",pp_final_name.Data(),AA_final_name.Data());
  TFile *ppf = new TFile(pp_final_name,"READ");
  TFile *AAf = new TFile(AA_final_name,"READ");
  
  if (!ppf) {
    fprintf(stderr,"Error: %s not loaded!\n",pp_final_name.Data());
    return;
  }
  if (!AAf) {
    fprintf(stderr,"Error: %s not loaded!\n",AA_final_name.Data());
    return;
  }

  TString jetClass = "jetPt_%.0f_%.0f";
  TString particleClass = "particlePt_%.2f_%.2f";
	TString combinedClass = "";


  // Extracting Projections from each file, ploting comparisons.
  // May also include code to look at differences, more useful 
  // for looking at keeping recoils vs not keeping recoils.
//	  TH2F * jetHadron[nJetPtBins][nParticlePtBins];
  
  TH1D * jetHProjection_pp[nJetPtBins][nParticlePtBins];
  TH1D * jetHProjection_AA[nJetPtBins][nParticlePtBins];

  TCanvas *cJetHProjection = new TCanvas("cJetHProjection","cJetHProjection",c_width,c_height);
  TCanvas *cJetHProjectionAS = new TCanvas("cJetHProjectionAS","cJetHProjectionAS",c_width,c_height);
  TCanvas *cJetHProjectionDiff = new TCanvas("cJetHProjectionDiff","cJetHProjectionDiff",c_width,c_height);

	for (int i = 0; i < nJetPtBins; i++) {
    int index = 1;


    cJetHProjection->Divide(TMath::CeilNint(TMath::Sqrt(nParticlePtBins)),
            TMath::FloorNint(TMath::Sqrt(nParticlePtBins)));
    cJetHProjectionAS->Divide(TMath::CeilNint(TMath::Sqrt(nParticlePtBins)),
            TMath::FloorNint(TMath::Sqrt(nParticlePtBins)));
    cJetHProjectionDiff->Divide(TMath::CeilNint(TMath::Sqrt(nParticlePtBins)),
            TMath::FloorNint(TMath::Sqrt(nParticlePtBins)));

    for (int j = 0; j < nParticlePtBins; j++) {
      
      cJetHProjection->cd(index);


      combinedClass = Form("jetHadron_jetPt_%.0f_%.0f_particlePt_%.2f_%.2f_Projection", jetPtBins[i], jetPtBins[i+1], particlePtBins[j], particlePtBins[j+1]);
      jetHProjection_pp[i][j] = (TH1D *) ppf->Get(combinedClass);
      jetHProjection_pp[i][j]->SetName(Form("%s_pp",jetHProjection_pp[i][j]->GetName()));
  
      jetHProjection_AA[i][j] = (TH1D *) AAf->Get(combinedClass);
      jetHProjection_AA[i][j]->SetName(Form("%s_AA",jetHProjection_pp[i][j]->GetName()));

      // Appearance
      jetHProjection_pp[i][j]->SetLineColor(kBlack);
      jetHProjection_pp[i][j]->SetMarkerColor(kBlack);
      jetHProjection_AA[i][j]->SetLineColor(kRed);
      jetHProjection_AA[i][j]->SetMarkerColor(kRed);

      jetHProjection_AA[i][j]->GetXaxis()->SetRangeUser(-TMath::Pi()/2,3*TMath::Pi()/2);

      jetHProjection_AA[i][j]->Draw();
      jetHProjection_pp[i][j]->Draw("SAME");



      cJetHProjectionDiff->cd(index);
      TH1D * clone_diff = (TH1D * ) jetHProjection_AA[i][j]->Clone();
      clone_diff->SetName(Form("%s_diff",clone_diff->GetName()));
      clone_diff->Add(jetHProjection_pp[i][j],-1);
      clone_diff->Rebin(2);
      clone_diff->Draw();

      //Awayside only
      // not the most memory efficient way of doing things, but easy to work with
      TH1D * pp_clone = (TH1D *) jetHProjection_pp[i][j]->Clone();
      TH1D * AA_clone = (TH1D *) jetHProjection_AA[i][j]->Clone();
      cJetHProjectionAS->cd(index);
      pp_clone->GetXaxis()->SetRangeUser(TMath::Pi()/2,3.*TMath::Pi()/2);
      AA_clone->GetXaxis()->SetRangeUser(TMath::Pi()/2,3.*TMath::Pi()/2);
      double y_min = 0;
      int y_max_bin = pp_clone->GetMaximumBin();
      double y_max = pp_clone->GetBinContent(y_max_bin) + pp_clone->GetBinError(y_max_bin);
      y_max_bin = AA_clone->GetMaximumBin();
      y_max = TMath::Max(y_max,AA_clone->GetBinContent(y_max_bin) + AA_clone->GetBinError(y_max_bin));
      pp_clone->GetYaxis()->SetRangeUser(y_min,y_max);    



      pp_clone->Draw();
      AA_clone->Draw("SAME");


//      if (clone_diff->GetYaxis()->GetXmin() > 0){ clone_diff->GetYaxis()->SetRangeUser(0,clone_diff->GetY->GetXmax());
//}

//      if (index ==
//            TMath::CeilNint(TMath::Sqrt(nParticlePtBins))*
//            TMath::FloorNint(TMath::Sqrt(nParticlePtBins)))
//            index++;

      index++;
    }
    combinedClass = Form("jetHadron_jetPt_%.0f_%.0f_Projection", jetPtBins[i], jetPtBins[i+1]);

    cJetHProjection->Print(Form("output/cmp/%s_cmp.pdf",combinedClass.Data()));
    cJetHProjectionAS->Print(Form("output/cmp/%sAS_cmp.pdf",combinedClass.Data()));
    cJetHProjectionDiff->Print(Form("output/cmp/%s_diff.pdf",combinedClass.Data()));

    cJetHProjection->Clear();
    cJetHProjectionAS->Clear();
    cJetHProjectionDiff->Clear();
  }


// the original code for this doesn't seem to work with CINT or whatever root 
// uses now.
  TGraphErrors *ptBinYields_pp[nJetPtBins];
  TGraphErrors *ptBinYields_AA[nJetPtBins];
  TGraphErrors *ptBinMeans_pp[nJetPtBins];
  TGraphErrors *ptBinMeans_AA[nJetPtBins];
  TGraphErrors *ptBinWidths_pp[nJetPtBins];
  TGraphErrors *ptBinWidths_AA[nJetPtBins];
  TGraphErrors *ptBinRms_pp[nJetPtBins];
  TGraphErrors *ptBinRms_AA[nJetPtBins];
  TGraphErrors *ptBinFwhm_pp[nJetPtBins];
  TGraphErrors *ptBinFwhm_AA[nJetPtBins];
  TGraphErrors *ptBinBetas_pp[nJetPtBins];
  TGraphErrors *ptBinBetas_AA[nJetPtBins];



  double ptBin_Y_pp[nJetPtBins][nParticlePtBins],ptBin_Y_pp_err[nJetPtBins][nParticlePtBins];
  double ptBin_Mean_pp[nJetPtBins][nParticlePtBins],ptBin_Mean_pp_err[nJetPtBins][nParticlePtBins];
  double ptBin_Y_AA[nJetPtBins][nParticlePtBins],ptBin_Y_AA_err[nJetPtBins][nParticlePtBins];
  double ptBin_Mean_AA[nJetPtBins][nParticlePtBins],ptBin_Mean_AA_err[nJetPtBins][nParticlePtBins];

// ---------------------------------------------------------------------------------
	for (int i = 0; i < nJetPtBins; i++) {


		//Get Yields
		combinedClass.Form("jetPt_%.0f_%.0f_Yields",jetPtBins[i],jetPtBins[i+1]);
		ptBinYields_pp[i] = (TGraphErrors *) ppf->Get(combinedClass);
		if(!ptBinYields_pp[i]) {
			printf("ERROR.  Can't find %s\n",combinedClass.Data());
			exit(1);
		}
  	ptBinYields_pp[i]->SetName(TString::Format("%s_pp",ptBinYields_pp[i]->GetName()));
		ptBinYields_AA[i] = (TGraphErrors *) AAf->Get(combinedClass);
		if(!ptBinYields_pp[i]) {
			printf("ERROR.  Can't find %s\n",combinedClass.Data());
			exit(1);
		}
  	ptBinYields_AA[i]->SetName(TString::Format("%s_AA",ptBinYields_AA[i]->GetName()));

    //Get Means
		combinedClass.Form("jetPt_%.0f_%.0f_Means",jetPtBins[i],jetPtBins[i+1]);
		ptBinMeans_pp[i] = (TGraphErrors *) ppf->Get(combinedClass);
		if(!ptBinMeans_pp[i]) {
			printf("ERROR.  Can't find %s\n",combinedClass.Data());
			exit(1);
		}
  	ptBinMeans_pp[i]->SetName(TString::Format("%s_pp",ptBinMeans_pp[i]->GetName()));
		ptBinMeans_AA[i] = (TGraphErrors *) AAf->Get(combinedClass);
		if(!ptBinMeans_pp[i]) {
			printf("ERROR.  Can't find %s\n",combinedClass.Data());
			exit(1);
		}
  	ptBinMeans_AA[i]->SetName(TString::Format("%s_AA",ptBinMeans_AA[i]->GetName()));



		//Get Widths
		combinedClass.Form("jetPt_%.0f_%.0f_Widths",jetPtBins[i],jetPtBins[i+1]);
		ptBinWidths_pp[i] = (TGraphErrors *) ppf->Get(combinedClass);
		if(!ptBinWidths_pp[i]) {
			printf("ERROR.  Can't find %s\n",combinedClass.Data());
			exit(1);
		}
  	ptBinWidths_pp[i]->SetName(TString::Format("%s_pp",ptBinWidths_pp[i]->GetName()));
		ptBinWidths_AA[i] = (TGraphErrors *) AAf->Get(combinedClass);
		if(!ptBinWidths_pp[i]) {
			printf("ERROR.  Can't find %s\n",combinedClass.Data());
			exit(1);
		}
  	ptBinWidths_AA[i]->SetName(TString::Format("%s_AA",ptBinWidths_AA[i]->GetName()));


		//Get RMS
		combinedClass.Form("jetPt_%.0f_%.0f_Rms",jetPtBins[i],jetPtBins[i+1]);
		ptBinRms_pp[i] = (TGraphErrors *) ppf->Get(combinedClass);
		if(!ptBinRms_pp[i]) {
			printf("ERROR.  Can't find %s\n",combinedClass.Data());
			exit(1);
		}
  	ptBinRms_pp[i]->SetName(TString::Format("%s_pp",ptBinRms_pp[i]->GetName()));
		ptBinRms_AA[i] = (TGraphErrors *) AAf->Get(combinedClass);
		if(!ptBinRms_pp[i]) {
			printf("ERROR.  Can't find %s\n",combinedClass.Data());
			exit(1);
		}
  	ptBinRms_AA[i]->SetName(TString::Format("%s_AA",ptBinRms_AA[i]->GetName()));


		//Get FWHM
		combinedClass.Form("jetPt_%.0f_%.0f_Fwhm",jetPtBins[i],jetPtBins[i+1]);
		ptBinFwhm_pp[i] = (TGraphErrors *) ppf->Get(combinedClass);
		if(!ptBinFwhm_pp[i]) {
			printf("ERROR.  Can't find %s\n",combinedClass.Data());
			exit(1);
		}
  	ptBinFwhm_pp[i]->SetName(TString::Format("%s_pp",ptBinFwhm_pp[i]->GetName()));
		ptBinFwhm_AA[i] = (TGraphErrors *) AAf->Get(combinedClass);
		if(!ptBinFwhm_pp[i]) {
			printf("ERROR.  Can't find %s\n",combinedClass.Data());
			exit(1);
		}
  	ptBinFwhm_AA[i]->SetName(TString::Format("%s_AA",ptBinFwhm_AA[i]->GetName()));


    //Get Betas
		combinedClass.Form("jetPt_%.0f_%.0f_Betas",jetPtBins[i],jetPtBins[i+1]);
		ptBinBetas_pp[i] = (TGraphErrors *) ppf->Get(combinedClass);
		if(!ptBinBetas_pp[i]) {
			printf("ERROR.  Can't find %s\n",combinedClass.Data());
			exit(1);
		}
  	ptBinBetas_pp[i]->SetName(TString::Format("%s_pp",ptBinBetas_pp[i]->GetName()));
		ptBinBetas_AA[i] = (TGraphErrors *) AAf->Get(combinedClass);
		if(!ptBinBetas_pp[i]) {
			printf("ERROR.  Can't find %s\n",combinedClass.Data());
			exit(1);
		}
  	ptBinBetas_AA[i]->SetName(TString::Format("%s_AA",ptBinBetas_AA[i]->GetName()));

    //pp
		double *tempPt_Y = ptBinYields_pp[i]->GetY();
		double *tempPt_EY = ptBinYields_pp[i]->GetEY();
		double *tempPt_M = ptBinMeans_pp[i]->GetY();
		double *tempPt_EM = ptBinMeans_pp[i]->GetEY();
		for (int j = 0; j < nParticlePtBins; j++ ) {
			ptBin_Y_pp[i][j] = tempPt_Y[j];
			ptBin_Y_pp_err[i][j] = tempPt_EY[j];
			ptBin_Mean_pp[i][j] = tempPt_M[j];
			ptBin_Mean_pp_err[i][j] = tempPt_EM[j];
		}

		//AA
		tempPt_Y = ptBinYields_AA[i]->GetY();
		tempPt_EY = ptBinYields_AA[i]->GetEY();
		tempPt_M = ptBinMeans_AA[i]->GetY();
		tempPt_EM = ptBinMeans_AA[i]->GetEY();
		for (int j = 0; j < nParticlePtBins; j++ ) {
			ptBin_Y_AA[i][j] = tempPt_Y[j];
			ptBin_Y_AA_err[i][j] = tempPt_EY[j];
			ptBin_Mean_AA[i][j] = tempPt_M[j];
			ptBin_Mean_AA_err[i][j] = tempPt_EM[j];
		}




		//Basic Settings
		ptBinYields_pp[i]->SetMarkerStyle(kOpenSquare);
		ptBinYields_pp[i]->SetFillColor(0);
		ptBinYields_pp[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		ptBinYields_pp[i]->GetYaxis()->SetTitle("#sigma_{AS}");
  	ptBinYields_pp[i]->SetLineStyle(2);
		
    ptBinYields_AA[i]->SetMarkerStyle(kFullSquare);
		ptBinYields_AA[i]->SetFillColor(0);
		ptBinYields_AA[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		ptBinYields_AA[i]->GetYaxis()->SetTitle("#sigma_{AS}");

		ptBinMeans_pp[i]->SetMarkerStyle(kOpenSquare);
		ptBinMeans_pp[i]->SetFillColor(0);
		ptBinMeans_pp[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		ptBinMeans_pp[i]->GetYaxis()->SetTitle("#beta_{AS}");
		ptBinMeans_pp[i]->SetLineStyle(2);
		
    ptBinMeans_AA[i]->SetMarkerStyle(kFullSquare);
		ptBinMeans_AA[i]->SetFillColor(0);
		ptBinMeans_AA[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		ptBinMeans_AA[i]->GetYaxis()->SetTitle("#beta_{AS}");


		ptBinWidths_pp[i]->SetMarkerStyle(kOpenSquare);
		ptBinWidths_pp[i]->SetFillColor(0);
		ptBinWidths_pp[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		ptBinWidths_pp[i]->GetYaxis()->SetTitle("#sigma_{AS}");
  	ptBinWidths_pp[i]->SetLineStyle(2);
		
    ptBinWidths_AA[i]->SetMarkerStyle(kFullSquare);
		ptBinWidths_AA[i]->SetFillColor(0);
		ptBinWidths_AA[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		ptBinWidths_AA[i]->GetYaxis()->SetTitle("#sigma_{AS}");

		ptBinRms_pp[i]->SetMarkerStyle(kOpenSquare);
		ptBinRms_pp[i]->SetFillColor(0);
		ptBinRms_pp[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		ptBinRms_pp[i]->GetYaxis()->SetTitle("RMS");
  	ptBinRms_pp[i]->SetLineStyle(2);
		
    ptBinRms_AA[i]->SetMarkerStyle(kFullSquare);
		ptBinRms_AA[i]->SetFillColor(0);
		ptBinRms_AA[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		ptBinRms_AA[i]->GetYaxis()->SetTitle("RMS");

		ptBinFwhm_pp[i]->SetMarkerStyle(kOpenSquare);
		ptBinFwhm_pp[i]->SetFillColor(0);
		ptBinFwhm_pp[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		ptBinFwhm_pp[i]->GetYaxis()->SetTitle("FWHM");
  	ptBinFwhm_pp[i]->SetLineStyle(2);
		
    ptBinFwhm_AA[i]->SetMarkerStyle(kFullSquare);
		ptBinFwhm_AA[i]->SetFillColor(0);
		ptBinFwhm_AA[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		ptBinFwhm_AA[i]->GetYaxis()->SetTitle("FWHM");

		ptBinBetas_pp[i]->SetMarkerStyle(kOpenSquare);
		ptBinBetas_pp[i]->SetFillColor(0);
		ptBinBetas_pp[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		ptBinBetas_pp[i]->GetYaxis()->SetTitle("#beta_{AS}");
		ptBinBetas_pp[i]->SetLineStyle(2);
		
    ptBinBetas_AA[i]->SetMarkerStyle(kFullSquare);
		ptBinBetas_AA[i]->SetFillColor(0);
		ptBinBetas_AA[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		ptBinBetas_AA[i]->GetYaxis()->SetTitle("#beta_{AS}");

	}



// ------------------------------------------------------------
//Settings for STAR PLOT

  ptBinWidths_pp[0]->SetTitle("10<p_{T}^{jet} < 15 GeV/c");
  ptBinWidths_pp[1]->SetTitle("15<p_{T}^{jet} < 20 GeV/c");
  ptBinWidths_pp[2]->SetTitle("20<p_{T}^{jet} < 40 GeV/c");
  ptBinWidths_pp[3]->SetTitle("40<p_{T}^{jet} < 60 GeV/c");

	ptBinWidths_pp[0]->SetLineColor(kRed);
	ptBinWidths_pp[0]->SetMarkerColor(kRed);
	ptBinWidths_pp[0]->SetMarkerStyle(kOpenCircle);
	ptBinWidths_pp[0]->SetLineStyle(2);
	ptBinWidths_pp[2]->SetLineColor(kBlack);
	ptBinWidths_pp[2]->SetMarkerColor(kBlack);
	ptBinWidths_pp[2]->SetMarkerStyle(kOpenSquare);
	ptBinWidths_pp[2]->SetLineStyle(2);
	


  ptBinWidths_AA[0]->SetTitle("10<p_{T}^{jet} < 15 GeV/c");
  ptBinWidths_AA[1]->SetTitle("15<p_{T}^{jet} < 20 GeV/c");
  ptBinWidths_AA[2]->SetTitle("20<p_{T}^{jet} < 40 GeV/c");
  ptBinWidths_AA[3]->SetTitle("40<p_{T}^{jet} < 60 GeV/c");
	

	ptBinWidths_AA[0]->SetLineColor(kRed);
	ptBinWidths_AA[0]->SetMarkerColor(kRed);
	ptBinWidths_AA[0]->SetMarkerStyle(kFullCircle);
	ptBinWidths_AA[0]->SetLineStyle(1);
	ptBinWidths_AA[2]->SetLineColor(kBlack);
	ptBinWidths_AA[2]->SetMarkerColor(kBlack);
	ptBinWidths_AA[2]->SetMarkerStyle(kFullSquare);
	ptBinWidths_AA[2]->SetLineStyle(1);


// cloning to make 
//	TGraphErrors *pt10pp = (TGraphErrors ptBinWidths_pp[0]->Clone();


	TLegend * leg = new TLegend(0.60, 0.65, 0.80, 0.85);
	leg->AddEntry(ptBinWidths_pp[0], "p+p 10 #leq p_{T}^{jet} #leq 15", "p");
	leg->AddEntry(ptBinWidths_pp[2], "p+p 20 #leq p_{T}^{jet} #leq 40", "p");
	leg->AddEntry(ptBinWidths_AA[0], "A+A 10 #leq p_{T}^{jet} #leq 15", "p");
	leg->AddEntry(ptBinWidths_AA[2], "A+A 20 #leq p_{T}^{jet} #leq 40", "p");
 
  TCanvas *c1 = new TCanvas("c1","c1",c_width,c_height);
  c1->cd();
  ptBinWidths_pp[0]->Draw("ALP");
  ptBinWidths_AA[0]->Draw("LP SAME");
  ptBinWidths_pp[2]->Draw("LP SAME");
  ptBinWidths_AA[2]->Draw("LP SAME");
	leg->Draw("SAME");
	c1->SetTitle("STAR PLOT");
	c1->SetLogy();
//	ptBinWidths_pp[0]->GetYaxis()->SetRangeUser(0.01,1.3);
	c1->Print("output/cmp/star_plot.pdf");
  c1->SetLogy(0);

// ------------------------------------------------------------



// ------------------------------------------------------------
//Color stuff

  int colorArray[4] = {kBlack,kRed,kGreen,kBlue};

	for (int i = 0; i < nJetPtBins; i++) {
    ptBinYields_pp[i]->SetMarkerColor(colorArray[i]);
    ptBinYields_pp[i]->SetLineColor(colorArray[i]);
    ptBinYields_AA[i]->SetMarkerColor(colorArray[i]);
    ptBinYields_AA[i]->SetLineColor(colorArray[i]);
    ptBinMeans_pp[i]->SetMarkerColor(colorArray[i]);
    ptBinMeans_pp[i]->SetLineColor(colorArray[i]);
    ptBinMeans_AA[i]->SetMarkerColor(colorArray[i]);
    ptBinMeans_AA[i]->SetLineColor(colorArray[i]);

  // FIXME doing stuff for poster
    ptBinWidths_pp[i]->SetMarkerColor(kBlack);
    ptBinWidths_pp[i]->SetLineColor(kBlack);
    ptBinWidths_pp[i]->SetMarkerStyle(kFullSquare);
    ptBinWidths_AA[i]->SetMarkerColor(kRed);
    ptBinWidths_AA[i]->SetLineColor(kRed);
    
    ptBinRms_pp[i]->SetMarkerColor(kBlack);
    ptBinRms_pp[i]->SetLineColor(kBlack);
    ptBinRms_pp[i]->SetMarkerStyle(kFullSquare);
    ptBinRms_AA[i]->SetMarkerColor(kRed);
    ptBinRms_AA[i]->SetLineColor(kRed);

    ptBinFwhm_pp[i]->SetMarkerColor(kBlack);
    ptBinFwhm_pp[i]->SetLineColor(kBlack);
    ptBinFwhm_pp[i]->SetMarkerStyle(kFullSquare);
    ptBinFwhm_AA[i]->SetMarkerColor(kRed);
    ptBinFwhm_AA[i]->SetLineColor(kRed);

   /* ptBinWidths_pp[i]->SetMarkerColor(colorArray[i]);
    ptBinWidths_pp[i]->SetLineColor(colorArray[i]);
    ptBinWidths_AA[i]->SetMarkerColor(colorArray[i]);
    ptBinWidths_AA[i]->SetLineColor(colorArray[i]);
     */ 
    ptBinBetas_pp[i]->SetMarkerColor(colorArray[i]);
    ptBinBetas_pp[i]->SetLineColor(colorArray[i]);
    ptBinBetas_AA[i]->SetMarkerColor(colorArray[i]);
    ptBinBetas_AA[i]->SetLineColor(colorArray[i]);
  
  }    




// Plots for each jet bin: AA vs pp


	TCanvas *cPtBinYields[nJetPtBins];
	TLegend *legPtBinYields[nJetPtBins];

	for (int i = 0; i < nJetPtBins; i++) {
		combinedClass.Form("cJetPt_%.0f_%.0f_Yields",jetPtBins[i],jetPtBins[i+1]);
		cPtBinYields[i] = new TCanvas(combinedClass,combinedClass,c_width,c_height);
		legPtBinYields[i] = new TLegend(0.60, 0.65, 0.80, 0.85);
		legPtBinYields[i]->AddEntry(ptBinYields_pp[i],"p+p");
		legPtBinYields[i]->AddEntry(ptBinYields_AA[i],"A+A");
		cPtBinYields[i]->cd();
	  double max_y = TMath::Max(ptBinYields_AA[i]->GetYaxis()->GetXmax(),ptBinYields_pp[i]->GetYaxis()->GetXmax());
	  double min_y = TMath::Min(ptBinYields_AA[i]->GetYaxis()->GetXmin(),ptBinYields_pp[i]->GetYaxis()->GetXmin());
    ptBinYields_AA[i]->GetYaxis()->SetRangeUser(min_y,max_y);
	  ptBinYields_AA[i]->Draw("ALP");
		ptBinYields_pp[i]->Draw("LP SAME");
		legPtBinYields[i]->Draw("SAME");
	 
//    cPtBinYields[i]->SetLogy();
		cPtBinYields[i]->Print(Form("output/cmp/%s.pdf",combinedClass.Data()));

	}

	TCanvas *cPtBinMeans[nJetPtBins];
	TLegend *legPtBinMeans[nJetPtBins];

	for (int i = 0; i < nJetPtBins; i++) {
		combinedClass.Form("cJetPt_%.0f_%.0f_Means",jetPtBins[i],jetPtBins[i+1]);
		cPtBinMeans[i] = new TCanvas(combinedClass,combinedClass,c_width,c_height);
		legPtBinMeans[i] = new TLegend(0.60, 0.65, 0.80, 0.85);
		legPtBinMeans[i]->AddEntry(ptBinMeans_pp[i],"p+p");
		legPtBinMeans[i]->AddEntry(ptBinMeans_AA[i],"A+A");
		cPtBinMeans[i]->cd();
	  double max_y = TMath::Max(ptBinMeans_AA[i]->GetYaxis()->GetXmax(),ptBinMeans_pp[i]->GetYaxis()->GetXmax());
	  double min_y = TMath::Min(ptBinMeans_AA[i]->GetYaxis()->GetXmin(),ptBinMeans_pp[i]->GetYaxis()->GetXmin());
    ptBinMeans_AA[i]->GetYaxis()->SetRangeUser(min_y,max_y);

		ptBinMeans_AA[i]->Draw("ALP");
		ptBinMeans_pp[i]->Draw("LP SAME");
		legPtBinMeans[i]->Draw("SAME");
	 
		cPtBinMeans[i]->Print(Form("output/cmp/%s.pdf",combinedClass.Data()));
	}



	TCanvas *cPtBinWidths[nJetPtBins];
	TLegend *legPtBinWidths[nJetPtBins];

	for (int i = 0; i < nJetPtBins; i++) {
		combinedClass.Form("cJetPt_%.0f_%.0f_Widths",jetPtBins[i],jetPtBins[i+1]);
		cPtBinWidths[i] = new TCanvas(combinedClass,combinedClass,c_width,c_height);
		legPtBinWidths[i] = new TLegend(0.60, 0.65, 0.80, 0.85);
		legPtBinWidths[i]->AddEntry(ptBinWidths_pp[i],"p+p");
		legPtBinWidths[i]->AddEntry(ptBinWidths_AA[i],"A+A");
		cPtBinWidths[i]->cd();
	  double max_y = TMath::Max(ptBinWidths_AA[i]->GetYaxis()->GetXmax(),ptBinWidths_pp[i]->GetYaxis()->GetXmax());
	  double min_y = TMath::Min(ptBinWidths_AA[i]->GetYaxis()->GetXmin(),ptBinWidths_pp[i]->GetYaxis()->GetXmin());
    ptBinWidths_AA[i]->GetYaxis()->SetRangeUser(min_y,max_y);
	  ptBinWidths_AA[i]->Draw("AP");
		ptBinWidths_pp[i]->Draw("P SAME");
    //FIXME
	//  ptBinWidths_AA[i]->Draw("ALP");
//		ptBinWidths_pp[i]->Draw("LP SAME");
		legPtBinWidths[i]->Draw("SAME");
 //   cPtBinWidths[i]->SetLogy();
		cPtBinWidths[i]->Print(Form("output/cmp/%s.pdf",combinedClass.Data()));

	}
	
  // RMS
  TCanvas *cPtBinRms[nJetPtBins];
	TLegend *legPtBinRms[nJetPtBins];

	for (int i = 0; i < nJetPtBins; i++) {
		combinedClass.Form("cJetPt_%.0f_%.0f_Rms",jetPtBins[i],jetPtBins[i+1]);
		cPtBinRms[i] = new TCanvas(combinedClass,combinedClass,c_width,c_height);
		legPtBinRms[i] = new TLegend(0.60, 0.65, 0.80, 0.85);
		legPtBinRms[i]->AddEntry(ptBinRms_pp[i],"p+p");
		legPtBinRms[i]->AddEntry(ptBinRms_AA[i],"A+A");
		cPtBinRms[i]->cd();
	  double max_y = TMath::Max(ptBinRms_AA[i]->GetYaxis()->GetXmax(),ptBinRms_pp[i]->GetYaxis()->GetXmax());
	  double min_y = TMath::Min(ptBinRms_AA[i]->GetYaxis()->GetXmin(),ptBinRms_pp[i]->GetYaxis()->GetXmin());
    ptBinRms_AA[i]->GetYaxis()->SetRangeUser(min_y,max_y);
	  ptBinRms_AA[i]->Draw("AP");
		ptBinRms_pp[i]->Draw("P SAME");
    //FIXME
	//  ptBinWidths_AA[i]->Draw("ALP");
//		ptBinWidths_pp[i]->Draw("LP SAME");
		legPtBinRms[i]->Draw("SAME");
 //   cPtBinWidths[i]->SetLogy();
		cPtBinRms[i]->Print(Form("output/cmp/%s.pdf",combinedClass.Data()));

	}

  // FWHM
	TCanvas *cPtBinFwhm[nJetPtBins];
	TLegend *legPtBinFwhm[nJetPtBins];

	for (int i = 0; i < nJetPtBins; i++) {
		combinedClass.Form("cJetPt_%.0f_%.0f_Fwhm",jetPtBins[i],jetPtBins[i+1]);
		cPtBinFwhm[i] = new TCanvas(combinedClass,combinedClass,c_width,c_height);
		legPtBinFwhm[i] = new TLegend(0.60, 0.65, 0.80, 0.85);
		legPtBinFwhm[i]->AddEntry(ptBinFwhm_pp[i],"p+p");
		legPtBinFwhm[i]->AddEntry(ptBinFwhm_AA[i],"A+A");
		cPtBinFwhm[i]->cd();
	  double max_y = TMath::Max(ptBinFwhm_AA[i]->GetYaxis()->GetXmax(),ptBinFwhm_pp[i]->GetYaxis()->GetXmax());
	  double min_y = TMath::Min(ptBinFwhm_AA[i]->GetYaxis()->GetXmin(),ptBinFwhm_pp[i]->GetYaxis()->GetXmin());
    ptBinFwhm_AA[i]->GetYaxis()->SetRangeUser(min_y,max_y);
	  ptBinFwhm_AA[i]->Draw("AP");
		ptBinFwhm_pp[i]->Draw("P SAME");
    //FIXME
	//  ptBinFwhm_AA[i]->Draw("ALP");
//		ptBinFwhm_pp[i]->Draw("LP SAME");
		legPtBinFwhm[i]->Draw("SAME");
 //   cPtBinFwhm[i]->SetLogy();
		cPtBinFwhm[i]->Print(Form("output/cmp/%s.pdf",combinedClass.Data()));

	}





	TCanvas *cPtBinBetas[nJetPtBins];
	TLegend *legPtBinBetas[nJetPtBins];

	for (int i = 0; i < nJetPtBins; i++) {
		combinedClass.Form("cJetPt_%.0f_%.0f_Betas",jetPtBins[i],jetPtBins[i+1]);
		cPtBinBetas[i] = new TCanvas(combinedClass,combinedClass,c_width,c_height);
		legPtBinBetas[i] = new TLegend(0.60, 0.65, 0.80, 0.85);
		legPtBinBetas[i]->AddEntry(ptBinBetas_pp[i],"p+p");
		legPtBinBetas[i]->AddEntry(ptBinBetas_AA[i],"A+A");
		cPtBinBetas[i]->cd();
	  double max_y = TMath::Max(ptBinBetas_AA[i]->GetYaxis()->GetXmax(),ptBinBetas_pp[i]->GetYaxis()->GetXmax());
	  double min_y = TMath::Min(ptBinBetas_AA[i]->GetYaxis()->GetXmin(),ptBinBetas_pp[i]->GetYaxis()->GetXmin());
    ptBinBetas_AA[i]->GetYaxis()->SetRangeUser(min_y,max_y);

		ptBinBetas_AA[i]->Draw("ALP");
		ptBinBetas_pp[i]->Draw("LP SAME");
		legPtBinBetas[i]->Draw("SAME");
	 
		cPtBinBetas[i]->Print(Form("output/cmp/%s.pdf",combinedClass.Data()));
	}


 



// R_AA	
  // should change to using the pt bins, find yield, and make tgraphs.
  // or create separate bins for R_AA, but that sounds like more work.
    
  // i don't remember what i was thinking with the previous comment
  
	if(do_JetPt_things) {
 
  TH1D * jetPt_pp = (TH1D *) ppf->Get("recJetPt");
  jetPt_pp->SetName(TString::Format("%s_pp",jetPt_pp->GetName()));
  TH1D * jetPt_AA = (TH1D *) AAf->Get("recJetPt");
  if (!jetPt_AA) {
    fprintf(stderr,"Error: recJetPt for AA missing.\n");
    exit();
  }
  jetPt_AA->SetName(TString::Format("%s_AA",jetPt_AA->GetName()));

  TCanvas *cCarl = new TCanvas("cCarl","cCarl",c_width,c_height);
  cCarl->cd();
  jetPt_pp->SetTitle("Reconstructed Jet Pt");
  jetPt_pp->GetXaxis()->SetTitle("p_{T}^{jet} GeV/c");
  jetPt_pp->GetYaxis()->SetTitle("1/N_{ev} dN/dp_{T}^{jet} (GeV/c)^{-1}");
  jetPt_pp->SetMarkerStyle(kOpenSquare);
  jetPt_AA->SetMarkerStyle(kFullSquare);
  jetPt_pp->Draw("P");
  jetPt_AA->Draw("P SAME");
	TLegend * legCarl = new TLegend(0.60, 0.65, 0.80, 0.85);
	legCarl->AddEntry(jetPt_pp, "pp", "p");
	legCarl->AddEntry(jetPt_AA, "AA", "p");
	legCarl->Draw("same");
  cCarl->SetLogy();
  cCarl->Print("output/cmp/RecJetCmp.pdf");


	int jet_rebin_n = 1;
	jetPt_AA->Rebin(jet_rebin_n);
	jetPt_AA->Scale(1./jet_rebin_n);
	jetPt_pp->Rebin(jet_rebin_n);
	jetPt_pp->Scale(1./jet_rebin_n);
	


  TH1D * jetR_AA = jetPt_AA->Clone();  
  jetR_AA->SetName("jetR_AA");
  jetR_AA->SetTitle("Jet R_{AA} (Background Subtracted)");
  jetR_AA->GetXaxis()->SetTitle("p_{T}^{jet} GeV/c");
  jetR_AA->GetYaxis()->SetTitle("R_{AA}");
  jetR_AA->Divide(jetPt_pp);
  jetR_AA->SetMarkerStyle(kFullTriangleUp);

  TCanvas *c2 = new TCanvas("c2","c2",c_width,c_height);
  c2->cd();
  jetR_AA->GetXaxis()->SetRangeUser(10,100);
  jetR_AA->GetYaxis()->SetRangeUser(jetR_AA->GetYaxis()->GetXmin(),1.4);
	jetR_AA->Draw("P");

  c2->Print("output/cmp/jetR_AA.pdf");

//raw stuff
  TH1D * raw_jetPt_pp = (TH1D *) ppf->Get("jetPt");
  raw_jetPt_pp->SetName(TString::Format("%s_pp",raw_jetPt_pp->GetName()));
  TH1D * raw_jetPt_AA = (TH1D *) AAf->Get("jetPt");
  raw_jetPt_pp->SetName(TString::Format("%s_AA",raw_jetPt_AA->GetName()));

  TH1D * raw_jetR_AA = (TH1D *) raw_jetPt_AA->Clone();  
  raw_jetR_AA->SetName("raw_jetR_AA");
  raw_jetR_AA->SetTitle("Jet R_{AA} (Raw)");
  raw_jetR_AA->GetXaxis()->SetTitle("p_{T}^{jet} GeV/c");
  raw_jetR_AA->GetYaxis()->SetTitle("R_{AA}");
  raw_jetR_AA->Divide(raw_jetPt_pp);
  raw_jetR_AA->SetMarkerStyle(kFullTriangleUp);

  TCanvas *c35 = new TCanvas("c35","c35",c_width,c_height);
  c35->cd();
  raw_jetR_AA->GetYaxis()->SetRangeUser(0,1.0);
//  raw_jetR_AA->GetXaxis()->SetRangeUser();
  raw_jetR_AA->Draw("P");

  c35->Print("output/cmp/Jet_R_AA_raw.pdf");  


	TH1D * trackPt_pp = (TH1D *) ppf->Get("trackPt");
	trackPt_pp->SetName(TString::Format("%s_pp",trackPt_pp->GetName()));
	TH1D * trackPt_AA = (TH1D *) AAf->Get("trackPt");
	trackPt_AA->SetName(TString::Format("%s_AA",trackPt_AA->GetName()));

	trackPt_pp->SetMarkerStyle(kOpenCircle);
	trackPt_AA->SetMarkerStyle(kFullCircle);
	


	TCanvas *cTrackCmp = new TCanvas("TrackCmp","TrackCmp",c_width,c_height);
	cTrackCmp->cd();
	cTrackCmp->SetLogy();
	trackPt_AA->Draw("P");
	trackPt_pp->Draw("P SAME");
	cTrackCmp->Print("output/cmp/TrackCmp.pdf");

	TH1D *trackR_AA = (TH1D *) trackPt_AA->Clone();
	trackR_AA->Divide(trackPt_pp);
	trackR_AA->SetName("trackR_AA");
	trackR_AA->SetTitle("R_{AA}^{Track}");
	TCanvas *cTrackR_AA = new TCanvas("TrackR_AA","TrackR_AA",c_width,c_height);
	cTrackR_AA->cd();
	trackR_AA->GetYaxis()->SetRangeUser(0.,1.6);
	trackR_AA->Draw();
	cTrackR_AA->Print("output/cmp/TrackRAA.pdf");

	//Gamma R_AA
/*	TH1D *gammaEt_pp = (TH1D *) ppf->Get("gammaEt");
  gammaEt_pp->SetName(TString::Format("%s_pp",gammaEt_pp->GetName()));
	TH1D *gammaEt_AA = (TH1D *) AAf->Get("gammaEt");
  gammaEt_AA->SetName(TString::Format("%s_pp",gammaEt_AA->GetName()));
	gammaEt_AA->GetXaxis()->SetTitle("E_{T}^{#gamma} GeV");
	int rebin_n = 4;
	gammaEt_AA->Rebin(rebin_n);
	gammaEt_AA->Scale(1.0/rebin_n);
	gammaEt_pp->Rebin(rebin_n);
	gammaEt_pp->Scale(1.0/rebin_n);

	TH1D *gammaR_AA = (TH1D *) gammaEt_AA->Clone();
	gammaR_AA->Divide(gammaEt_pp);
	gammaR_AA->SetName("gammaR_AA");
	gammaR_AA->SetTitle("Gamma R_{AA}");
	gammaR_AA->GetYaxis()->SetTitle("Gamma R_{AA}");
	gammaR_AA->SetMarkerStyle(kFullTriangleUp);
	gammaR_AA->SetMarkerColor(kRed);
	gammaR_AA->SetLineColor(kRed);
	TCanvas *cGammaRAA = new TCanvas("cGammaRAA");
	cGammaRAA->cd();
	gammaR_AA->Draw("P");
	*/
	}

  //D_AA 
  double *pt10_Y_pp,*pt10_Y_pp_err;
  double *pt20_Y_pp,*pt20_Y_pp_err;
  
  double *pt10_mean_pp,*pt10_mean_pp_err;
  double *pt20_mean_pp,*pt20_mean_pp_err;

  double *pt10_Y_AA,*pt10_Y_AA_err;
  double *pt20_Y_AA,*pt20_Y_AA_err;
  
  double *pt10_mean_AA,*pt10_mean_AA_err;
  double *pt20_mean_AA,*pt20_mean_AA_err;

  std::vector <double> ptBin;
  std::vector <double> ptBin_err;   

/*  std::vector <double> pt10D;
  std::vector <double> pt10D_err;
  std::vector <double> pt20D;
  std::vector <double> pt20D_err;

*/
	double ptBinD_AA[nJetPtBins][nParticlePtBins];
	double ptBinD_AA_err[nJetPtBins][nParticlePtBins];

  double ptBinSigmaD_AA[nJetPtBins];
  double ptBinSigmaD_AA_err[nJetPtBins];

	double ptBinI_AA[nJetPtBins][nParticlePtBins];
	double ptBinI_AA_err[nJetPtBins][nParticlePtBins];


// Making x axis.  Assume shared.
  TGraphErrors *ptBin_pp_1 = (TGraphErrors *) ppf->Get(Form("jetPt_%.0f_%.0f_Widths",jetPtBins[0],jetPtBins[1]));

  for (int i = 0; i < ptBin_pp_1->GetN(); i++) {
		// this is extremely inelegant, but it's just so unimportant
    ptBin.push_back(ptBin_pp_1->GetX()[i]);
    ptBin_err.push_back(ptBin_pp_1->GetErrorX(i));    
	}

	for (int i = 0; i < nJetPtBins; i++) {
    ptBinSigmaD_AA[i] = 0;
    ptBinSigmaD_AA_err[i] = 0;

		for (int j = 0; j < nParticlePtBins; j++) {
			ptBinD_AA[i][j] = ptBin_Y_AA[i][j] * ptBin_Mean_AA[i][j] -
							ptBin_Y_pp[i][j] * ptBin_Mean_pp[i][j] ;

    	ptBinD_AA_err[i][j] = TMath::Sqrt(
        ptBin_Y_AA[i][j] * ptBin_Mean_AA[i][j] * TMath::Sqrt(
          TMath::Power(ptBin_Y_AA_err[i][j] / ptBin_Y_AA[i][j],2) 
          + TMath::Power(ptBin_Mean_AA_err[i][j]/ptBin_Mean_AA[i][j],2)
          )
        +
        ptBin_Y_pp[i][j] * ptBin_Mean_pp[i][j] * TMath::Sqrt(
          TMath::Power(ptBin_Y_pp_err[i][j] / ptBin_Y_pp[i][j],2) 
          + TMath::Power(ptBin_Mean_pp_err[i][j]/ptBin_Mean_pp[i][j],2)
          )
      );  
	
			ptBinI_AA[i][j] = ptBin_Y_AA[i][j] / ptBin_Y_pp[i][j] ;
      ptBinI_AA_err[i][j] = ptBinI_AA[i][j] * TMath::Sqrt(
        TMath::Power(ptBin_Y_AA_err[i][j]/ptBin_Y_AA[i][j],2)
        + TMath::Power(ptBin_Y_pp_err[i][j]/ptBin_Y_pp[i][j],2) 
      );
			if (ptBinI_AA[i][j] < 0) ptBinI_AA[i][j] = 0.;

      ptBinSigmaD_AA[i] += ptBinD_AA[i][j];
      ptBinSigmaD_AA_err[i] += TMath::Power(ptBinD_AA_err[i][j],2);
		}
    ptBinSigmaD_AA_err[i] = TMath::Sqrt(ptBinSigmaD_AA_err[i]);
	}




  

	double max_y= 1;
	double min_y = -1;

	TGraphErrors *gPtBin_D_AA[nJetPtBins];
  TGraphErrors *gPtBin_SigmaD_AA = new TGraphErrors(nJetPtBins);
	for (int i = 0; i < nJetPtBins; i++) {
    gPtBin_SigmaD_AA->SetPoint(i,(jetPtBins[i+1]+jetPtBins[i])/2.,ptBinSigmaD_AA[i]);
    gPtBin_SigmaD_AA->SetPointError(i,(jetPtBins[i+1]-jetPtBins[i])/2.,ptBinSigmaD_AA_err[i]);
		gPtBin_D_AA[i] = new TGraphErrors(nParticlePtBins);
		for (int j = 0; j < nParticlePtBins; j++) {
			gPtBin_D_AA[i]->SetPoint(j,ptBin[j],ptBinD_AA[i][j]);
			gPtBin_D_AA[i]->SetPointError(j,ptBin_err[j],ptBinD_AA_err[i][j]);
		}
	
		if (gPtBin_D_AA[i]->GetYaxis()->GetXmax() > max_y) max_y = gPtBin_D_AA[i]->GetYaxis()->GetXmax();
		if (gPtBin_D_AA[i]->GetYaxis()->GetXmin() < min_y) min_y = gPtBin_D_AA[i]->GetYaxis()->GetXmin();

		gPtBin_D_AA[i]->SetTitle("D_{AA}"); //FIXME

    
	  gPtBin_D_AA[i]->SetName(Form("D_AA_jetPt_%.0f_%.0f",jetPtBins[i],jetPtBins[i+1]));
	//  gPtBin_D_AA[i]->SetTitle(Form("%.0f < p_{T}^{jet}  < %.0f GeV/c",jetPtBins[i],jetPtBins[i+1]));
		gPtBin_D_AA[i]->GetXaxis()->SetTitle("p_{T}^{assoc}");
		gPtBin_D_AA[i]->GetYaxis()->SetTitle("D_{AA} (GeV/c)");
		
		gPtBin_D_AA[i]->SetMarkerStyle(kFullSquare);
    gPtBin_D_AA[i]->SetMarkerColor(colorArray[i]);
    gPtBin_D_AA[i]->SetLineColor(colorArray[i]);

	}
  gPtBin_SigmaD_AA->SetName("SigmaD_AA");
  gPtBin_SigmaD_AA->SetTitle("#SigmaD_{AA}");
  gPtBin_SigmaD_AA->GetXaxis()->SetTitle("p_T^{jet} (GeV/c)");
  gPtBin_SigmaD_AA->GetYaxis()->SetTitle("#Sigma D_{AA} (GeV/c)");
  gPtBin_SigmaD_AA->SetMarkerStyle(kFullSquare);
  gPtBin_SigmaD_AA->SetMarkerColor(kBlack);
  gPtBin_SigmaD_AA->SetLineColor(kBlack);
  TCanvas *cSigmaD_AA = new TCanvas("cSigmaD_AA","cSigmaD_AA",c_width,c_height);
  gPtBin_SigmaD_AA->Draw("");
  cSigmaD_AA->SetGridy();
  cSigmaD_AA->Print("output/cmp/SigmaD_AA.pdf");


  //Setting y axis to appropriate size
	gPtBin_D_AA[0]->GetYaxis()->SetRangeUser(min_y,max_y);
  TCanvas *cD_AA = new TCanvas("cD_AA","cD_AA",c_width,c_height);
  cD_AA->cd();
  cD_AA->SetGrid();
  gPtBin_D_AA[0]->GetYaxis()->SetRangeUser(-15,15);
	gPtBin_D_AA[0]->Draw("ALP");
	for (int i = 1; i < nJetPtBins; i++) {
		gPtBin_D_AA[i]->Draw("SAME LP");
	}
	TLegend * legD_AA = new TLegend(0.60, 0.65, 0.85, 0.85);
	for (int i = 0; i < nJetPtBins; i++) {
	  legD_AA->AddEntry(gPtBin_D_AA[i], Form("%.0f < p_{T}^{jet}  < %.0f GeV/c",jetPtBins[i],jetPtBins[i+1]), "p");
	}
	legD_AA->Draw("same");
	cD_AA->Print("output/cmp/D_AA.pdf");


	max_y= 1;
	min_y = 0;

	TGraphErrors *gPtBin_I_AA[nJetPtBins];
	for (int i = 0; i < nJetPtBins; i++) {
		gPtBin_I_AA[i] = new TGraphErrors(nParticlePtBins);
		for (int j = 0; j < nParticlePtBins; j++) {
			gPtBin_I_AA[i]->SetPoint(j,ptBin[j],ptBinI_AA[i][j]);
			gPtBin_I_AA[i]->SetPointError(j,ptBin_err[j],ptBinI_AA_err[i][j]);
		}
	
		if (gPtBin_I_AA[i]->GetYaxis()->GetXmax() > max_y) max_y = gPtBin_I_AA[i]->GetYaxis()->GetXmax();
		if (gPtBin_I_AA[i]->GetYaxis()->GetXmin() < min_y) min_y = gPtBin_I_AA[i]->GetYaxis()->GetXmin();

		gPtBin_I_AA[i]->SetTitle("I_{AA}"); //FIXME

	  gPtBin_I_AA[i]->SetName(Form("I_AA_jetPt_%.0f_%.0f",jetPtBins[i],jetPtBins[i+1]));
//	  gPtBin_I_AA[i]->SetTitle(Form("%.0f < p_{T}^{jet}  < %.0f GeV/c",jetPtBins[i],jetPtBins[i+1]));
		gPtBin_I_AA[i]->GetXaxis()->SetTitle("p_{T}^{assoc}");
		gPtBin_I_AA[i]->GetYaxis()->SetTitle("I_{AA} (GeV/c)");
		
		gPtBin_I_AA[i]->SetMarkerStyle(kFullSquare);
    gPtBin_I_AA[i]->SetMarkerColor(colorArray[i]);
    gPtBin_I_AA[i]->SetLineColor(colorArray[i]);

	}

  //Setting y axis to appropriate size
	gPtBin_I_AA[0]->GetYaxis()->SetRangeUser(min_y,max_y);
  TCanvas *cI_AA = new TCanvas("cI_AA","cI_AA",c_width,c_height);
  cI_AA->cd();
  cI_AA->SetGrid();
  gPtBin_I_AA[0]->GetYaxis()->SetRangeUser(-0.1,2);
//  gPtBin_I_AA[0]->GetYaxis()->SetRangeUser(-15,15);
	gPtBin_I_AA[0]->Draw("ALP");
	for (int i = 1; i < nJetPtBins; i++) {
		gPtBin_I_AA[i]->Draw("SAME LP");
	}
	TLegend * legI_AA = new TLegend(0.60, 0.65, 0.85, 0.85);
	for (int i = 0; i < nJetPtBins; i++) {
	  legI_AA->AddEntry(gPtBin_I_AA[i], Form("%.0f < p_{T}^{jet}  < %.0f GeV/c",jetPtBins[i],jetPtBins[i+1]), "p");
	}
	legI_AA->Draw("same");
	cI_AA->Print("output/cmp/I_AA.pdf");




  // Dijet Comparisons
  TH1F *dijetAj_pp = (TH1F *) ppf->Get("dijetAj");
  if (dijetAj_pp) dijetAj_pp->SetName(Form("%s_pp",dijetAj_pp->GetName()));
  TH1F *dijetAj_AA = (TH1F *) AAf->Get("dijetAj");
  if (dijetAj_AA) dijetAj_AA->SetName(Form("%s_AA",dijetAj_AA->GetName()));
  if (dijetAj_pp && dijetAj_AA) {
//    TH1F *dijetAj_cmp = (TH1F *) dijetAj_AA->Clone();
//    dijetAj_cmp->Divide(dijetAj_pp);
    TCanvas *cAj = new TCanvas("cA_j","cA_j",c_width,c_height);  
    cAj->cd();
	  TLegend * legAj = new TLegend(0.60, 0.65, 0.85, 0.85);
//    dijetAj_cmp->GetYaxis()->SetRangeUser(0.,3.);
//    dijetAj_cmp->Draw();
    dijetAj_AA->SetLineColor(kRed);
    dijetAj_AA->SetMarkerColor(kRed);
    dijetAj_AA->SetMarkerStyle(kFullCircle);
    dijetAj_pp->SetLineColor(kBlack);
    dijetAj_pp->SetMarkerColor(kBlack);
    dijetAj_pp->SetMarkerStyle(kOpenDiamond);
    dijetAj_AA->Draw("LP");
    dijetAj_pp->Draw("LP SAME");
    legAj->AddEntry(dijetAj_pp,"pp","p");
    legAj->AddEntry(dijetAj_AA,"AA","p");

    legAj->Draw("SAME");
    cAj->Print("output/cmp/dijetAj_cmp.pdf");
  }


  TH1F *dijetXj_pp = (TH1F *) ppf->Get("dijetXj");
  if (dijetXj_pp) dijetXj_pp->SetName(Form("%s_pp",dijetXj_pp->GetName()));
  TH1F *dijetXj_AA = (TH1F *) AAf->Get("dijetXj");
  if (dijetXj_AA) dijetXj_AA->SetName(Form("%s_AA",dijetXj_AA->GetName()));
  if (dijetXj_pp && dijetXj_AA) {
  //  TH1F *dijetXj_cmp = (TH1F *) dijetXj_AA->Clone();
 //   dijetXj_cmp->Divide(dijetXj_pp);
    TCanvas *cXj = new TCanvas("cX_j","cX_j",c_width,c_height);  
	  TLegend * legXj = new TLegend(0.60, 0.65, 0.85, 0.85);
    cXj->cd();

    dijetXj_AA->SetLineColor(kRed);
    dijetXj_AA->SetMarkerColor(kRed);
    dijetXj_AA->SetMarkerStyle(kFullCircle);
    dijetXj_pp->SetLineColor(kBlack);
    dijetXj_pp->SetMarkerColor(kBlack);
    dijetXj_pp->SetMarkerStyle(kOpenDiamond);
    dijetXj_AA->Draw("LP");
    dijetXj_pp->Draw("LP SAME");
    legXj->AddEntry(dijetXj_pp,"pp","p");
    legXj->AddEntry(dijetXj_AA,"AA","p");

    legXj->Draw("SAME");

    
//    dijetXj_cmp->GetYaxis()->SetRangeUser();
   // dijetXj_cmp->Draw();
    cXj->Print("output/cmp/dijetXj_cmp.pdf");
  }


  // Saving some things to output:
  TFile *outFile = new TFile(outFile_path,"RECREATE");
// FIXME

  outFile->Add(jetR_AA); 
  outFile->Add(raw_jetR_AA); 
  outFile->Add(trackR_AA);
//  if (dijetAj_pp && dijetAj_AA) outFile->Add(dijetAj_cmp);
 // if (dijetXj_pp && dijetXj_AA) outFile->Add(dijetXj_cmp);

	for (int i = 0; i < nJetPtBins; i++) {
    outFile->Add(gPtBin_D_AA[i]);
    outFile->Add(gPtBin_I_AA[i]);

		outFile->Add(ptBinWidths_pp[i]);
	  outFile->Add(ptBinWidths_AA[i]);
		outFile->Add(ptBinRms_pp[i]);
	  outFile->Add(ptBinRms_AA[i]);
		outFile->Add(ptBinFwhm_pp[i]);
	  outFile->Add(ptBinFwhm_AA[i]);
	}

  outFile->Write();

}



