// basicAnalysis file
// 
// Can compile with:
//  make mBasicAnalysis
//
//  or
//
//  g++ -Wl,--no-as-needed -std=c++11 basicAnalysis.cc -o basicAnalysis -I`root-config --incdir` `root-config --libs` `fastjet-config --cxxflags --libs --plugins`
//
//  or (if you have readableStyle.h)
//
//  make USER_DEFINED=-DRE_STYLE=1 basicAnalysis

#ifdef RE_STYLE
#include <readableStyle.h>
#endif

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>

#include <fastjet/ClusterSequence.hh>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>

// Put deltaPhi in the range of [-Pi/2, 3*Pi/2]
double calculateDeltaPhi(const double phiOne, const double phiTwo)
{
	double deltaPhi = phiOne - phiTwo;
	if (deltaPhi < -TMath::Pi()/2)
	{
		deltaPhi = deltaPhi + 2*TMath::Pi();
	}
	if (deltaPhi > 3*TMath::Pi()/2)
	{
		deltaPhi = deltaPhi - 2*TMath::Pi();
	}

	return deltaPhi;
}

// Needed for peak fitting
struct minMax
{
	minMax(double _min, double _max): min(_min), max(_max) {}
	double min;
	double max;
};

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
	hist->Fit(peakFit, "RQ0+");
	//hist->Fit(peakFit, "R0+");

	//double peakMean = peakFit->GetParameter("Mean");
	//double width = peakFit->GetParameter("Sigma");

	return peakFit;
}

void basicAnalysis(TString inputFilename = "yajemDefaultProfile", TString
        outputFilename = "defaultOut")
{
	// Open input file
	inputFilename = Form("output/%s.root", inputFilename.Data());
	TFile * fIn = TFile::Open(inputFilename.Data());

	// Check if file is open
	if (fIn->IsOpen() == kFALSE)
	{
		Printf("inputFilename \"%s\" is not found!", inputFilename.Data());
		std::exit(1);
	}

	// The name of the tree is "T"
	TTree * tree = (TTree *) fIn->Get("T");

	// Define our basic variables
	// Basic types don't need to be pointers, but classes do
	//Int_t eventNumber = 0;
	//Int_t nParticles = 0;
	std::vector <double> * particleID = 0;
	std::vector <double> * px = 0;
	std::vector <double> * py = 0;
	std::vector <double> * pz = 0;
	std::vector <double> * energy = 0;
	std::vector <double> * mass = 0;

	// Set Brnahces to our variables
	//tree->SetBranchAddress("eventNumber", &eventNumber);
	//tree->SetBranchAddress("nParticles", &nParticles);
	tree->SetBranchAddress("particleID", &particleID);
	tree->SetBranchAddress("px", &px);
	tree->SetBranchAddress("py", &py);
	tree->SetBranchAddress("pz", &pz);
	tree->SetBranchAddress("energy", &energy);
	tree->SetBranchAddress("mass", &mass);

	// Use either TLorentzVector or fastjet::PseudoJet for the particles
	// Need to use fastjet::PesudoJet for jets regardless
	//std::vector <TLorentzVector> particles;
	std::vector <fastjet::PseudoJet> particles;
	std::vector <fastjet::PseudoJet> jets;
	
	// Define the pt bins
	// Size of the bins is from STAR Jet-h paper
	// Desired bins are 0 - 2 based on the convention below.
	std::vector <double> jetPtBins = {10, 15, 20, 40};
	// Desired bins are 0 - 10 based on the convention below.
	std::vector <double> particlePtBins = {0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 17};
	
	// Define histograms
	TH2F * dEtaDPhi = new TH2F("dEtaDPhi", "dEtaDPhi", 100, -2, 2, 100, -TMath::Pi()/2, 3.0*TMath::Pi()/2);
	TH2F * etaPhi = new TH2F("etaPhi", "etaPhi", 100, -5, 5, 100, 0, 2.0*TMath::Pi());
	TH1D * jetPt = new TH1D("jetPt", "jetPt", 100, 0, 100);
	TH1D * particleEnergy = new TH1D("particleEnergy", "Total Energy of Outgoing Particles", 40, 180, 220);
//MHO
	TH2F * jetEtaPhi = new TH2F("jetEtaPhi", "jetEtaPhi", 100, -5, 5, 100, 0, 2.0*TMath::Pi());


	// Jet-hadron is done more differentially
	enum { nJetPtBins = 3, nParticlePtBins = 11 };
	TH2F * jetHadron[nJetPtBins][nParticlePtBins];
	TString jetClass = "jetPt_%.0f_%.0f";
	TString particleClass = "particlePt_%.2f_%.2f";
	TString combinedClass = "";
	for (unsigned int i = 0; i < nJetPtBins; i++)
	{
		for (unsigned int j = 0; j < nParticlePtBins; j++)
		{
			combinedClass = Form("jetHadron_" + jetClass + "_" + particleClass, jetPtBins.at(i), jetPtBins.at(i+1), particlePtBins.at(j), particlePtBins.at(j+1));
			//Printf("combinedClass: %s", combinedClass.Data());
			jetHadron[i][j] = new TH2F(combinedClass, combinedClass, 100, -2, 2, 100, -TMath::Pi()/2, 3.0*TMath::Pi()/2);
		}
	}

	// Event loop
	for (Int_t eventNumber = 0; eventNumber < tree->GetEntries(); eventNumber++)
	{
		tree->GetEntry(eventNumber);

		// Fill particles ito array
		for (unsigned int particleNumber = 0; particleNumber < px->size(); particleNumber++)
		{
			// For TLorentzVector
			/*TLorentzVector particle;
			particle.SetPxPyPzE(px->at(particleNumber), py->at(particleNumber),
					pz->at(particleNumber), energy->at(particleNumber) );
			particles.push_back(particle);*/
			// For Fastjet
			particles.push_back( fastjet::PseudoJet( px->at(particleNumber), py->at(particleNumber),
					pz->at(particleNumber), energy->at(particleNumber) ));
		}

		// Refactor from here so that it can be used by everyone
		// Define jet algorithm
		double R = 0.4;
		fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R);

		// Perform jet finding
		fastjet::ClusterSequence cs(particles, jetDef);
		jets = fastjet::sorted_by_pt(cs.inclusive_jets());

		// Begin the actual analysis 
		// Jet-h correlations
		for (unsigned int i = 0; i < jets.size(); i++)
		{
			// Fill jet spectra
			jetPt->Fill(jets.at(i).pt());
      jetEtaPhi->Fill(jets.at(i).eta(),jets.at(i).phi());


			std::vector<double>::iterator jetPtBinIter;
			jetPtBinIter = std::lower_bound(jetPtBins.begin(), jetPtBins.end(), jets.at(i).pt());
			// The -1 is because lower_bound iterates to the first value greater than the given (jetPt) value
			Int_t jetPtBin = jetPtBinIter - jetPtBins.begin() - 1;

			if ((jetPtBin >= 0) && (jetPtBin < nJetPtBins))
			{
				// Reassign jetPtBin from 2 -> 1 to allow the array to be continuous
				//if (jetPtBin == 2) { jetPtBin = 1; }

				//Printf("jetPt: %f, Bin: %d", jets.at(i).pt(), jetPtBin);
				std::vector <double>::iterator particlePtBinIter;
				Int_t particlePtBin = -1;
				for (unsigned int j = 0; j < particles.size(); j++)
				{
					particlePtBinIter = std::lower_bound(particlePtBins.begin(), particlePtBins.end(), particles.at(j).pt());
					// The -1 is because lower_bound counts below the lowest value as defining the 0 index.
					// We are interested in everything above that index.
					particlePtBin = particlePtBinIter - particlePtBins.begin() - 1;
					//Printf("particlePt: %f, Bin: %d", particles.at(j).pt(), particlePtBin);

					if ((particlePtBin >= 0) && (particlePtBin < nParticlePtBins))
					{
						jetHadron[jetPtBin][particlePtBin]->Fill(jets.at(i).eta()-particles.at(j).eta(), calculateDeltaPhi(jets.at(i).phi(), particles.at(j).phi()));
					}
					else
					{
						//Printf("\t\tWarning: Particle pt out of range!! particlePt: %f, Bin: %d", particles.at(j).pt(), particlePtBin);
					}
				}
			}
			else
			{
				// We don't care about anything with less than 10 GeV. This is condition will only be met with a jet over 40
				if (jets.at(i).pt() > 10)
				{
					//Printf("\t\tOutside range:jetPt: %f, Bin: %d", jets.at(i).pt(), jetPtBin);
				}
			}
		}


		// Do angular correlations
		double particleEnergyTemp = 0;
		for (unsigned int i = 0; i < particles.size(); i++)
		{
			// Eta-Phi of each particle
			etaPhi->Fill(particles.at(i).eta(), particles.at(i).phi());

//MHO fix: we need to double count here 
// to avoid biases in the ordering of particles
//      for (unsigned int j = i+1; j < particles.size(); j++)
			for (unsigned int j = 0; j < particles.size(); j++)
			{
// MHO
        if (i == j) continue;
				// deltaEta-deltaPhi or each particle pair
				// For TLorentzVector
				//dEtaDPhi->Fill(particles.at(i).Eta()-particles.at(j).Eta(), calculateDeltaPhi(particles.at(i).Phi(), particles.at(j).Phi()));
				// For Fastjet
				dEtaDPhi->Fill(particles.at(i).eta()-particles.at(j).eta(), calculateDeltaPhi(particles.at(i).phi(), particles.at(j).phi()) );
			}
			particleEnergyTemp += particles.at(i).E();
		}

		// Fill total particle energy
		particleEnergy->Fill(particleEnergyTemp);

		// Clear vector of PseudoJets after each event
		particles.clear();
		jets.clear();

		if ((eventNumber % 1000 == 0) && (eventNumber != 0)) { printf("eventNumber = %d\n", eventNumber); }

		// For testing purposes
		//if (eventNumber > 99) {break;}
	}

 // fIn->Close();

	outputFilename = Form("output/%s.root", outputFilename.Data());
	TFile * fOut = TFile::Open(outputFilename.Data(),"RECREATE");
	if (fOut->IsOpen() == kFALSE)
	{
		Printf("outputFilename \"%s\" failed to open!", outputFilename.Data());
		std::exit(1);
	}



	// Draw and save histograms
	TCanvas canvas;
	// Plot dEtaDPhi for multiple types of drawings
	dEtaDPhi->GetXaxis()->SetTitle("#Delta#eta");
	dEtaDPhi->GetYaxis()->SetTitle("#Delta#phi");
	dEtaDPhi->Draw("surf1");
	canvas.Print("output/dEtaDPhi.surf.pdf");
	dEtaDPhi->Draw("colz");
	canvas.Print("output/dEtaDPhi.pdf");
  
  fOut->Add(dEtaDPhi);

	/*jetHadron->GetXaxis()->SetTitle("#Delta#eta");
	jetHadron->GetYaxis()->SetTitle("#Delta#phi");
	jetHadron->Draw("surf1");
	canvas.Print("output/jetHadron.surf.pdf");
	jetHadron->Draw("colz");
	canvas.Print("output/jetHadron.pdf");*/
	// etaPhi plot
	etaPhi->GetXaxis()->SetTitle("#eta");
	etaPhi->GetYaxis()->SetTitle("#phi");
	etaPhi->Draw("colz");
	canvas.Print("output/etaPhi.pdf");
  fOut->Add(etaPhi);
	// Overall energy of particle
	particleEnergy->GetXaxis()->SetTitle("E (GeV)");
	particleEnergy->Draw();
	canvas.Print("output/particleEnergy.pdf");
  fOut->Add(particleEnergy);
	// Jet spectra
	jetPt->GetXaxis()->SetTitle("p_{T} (GeV)");
	jetPt->Draw();
	canvas.SetLogy(1);
	canvas.Print("output/jetPt.pdf");
	canvas.SetLogy(0);
  fOut->Add(jetPt);

//MHO JetEtaPhi
  jetEtaPhi->GetXaxis()->SetTitle("#eta");
  jetEtaPhi->GetYaxis()->SetTitle("#phi");
  jetEtaPhi->Draw("colz");
  canvas.Print("output/jetEtaPhi.pdf");
  fOut->Add(jetEtaPhi);


	// JetH
	canvas.Clear();
	// -1 since we are ignoring 15-20 to recreate the STAR plots
	canvas.Divide(nParticlePtBins, nJetPtBins-1);
	Int_t index = 1;
	// Refactor this into a function later
	// Surf seems to be far to computationally intensive to be worth it!
	/*for (unsigned int i = 0; i < nJetPtBins; i++)
	{
		for (unsigned int j = 0; j < nParticlePtBins; j++)
		{
			// The canvas is indexed starting at 1
			canvas.cd(index);
			index++;
			jetHadron[i][j]->Draw("surf1");
		}
	}
	canvas.Print("output/jetHadron.surf.pdf");
	canvas.Clear();*/
	for (unsigned int i = 0; i < nJetPtBins; i++)
	{
		// We don't care about 15-20, so just skip it
		if (i == 1) { i++; }
		for (unsigned int j = 0; j < nParticlePtBins; j++)
		{
			// The canvas is indexed starting at 1
			canvas.cd(index);
			index++;
			jetHadron[i][j]->Draw("colz");
		}
	}
	canvas.Print("output/jetHadron.pdf");
	// Reset the divide canvas;
	canvas.cd();

	// Fit projections and draw
	// 2 particle correlations
	// It has to be a TH1D according to root documentation
	TH1D * phiProjection = dEtaDPhi->ProjectionY("phiProjection");
	minMax nearSideLimits(-TMath::Pi()/4, TMath::Pi()/4);
	minMax awaySideLimits(3.0*TMath::Pi()/4, 5.0*TMath::Pi()/4);

	// Perform fit
	TF1 * nearSideFit = fitPeak(phiProjection, nearSideLimits, "nearSide");
	TF1 * awaySideFit = fitPeak(phiProjection, awaySideLimits, "awaySide");
	
	// Set fit parameters
	nearSideFit->SetLineColor(kRed);
	awaySideFit->SetLineColor(kBlue);

	// Draw projection and fit
	phiProjection->Draw();
	nearSideFit->Draw("same");
	awaySideFit->Draw("same");
	
	// Save projection
	canvas.Print("output/phiProjection.pdf");

	// For eta within +/- 0.5
	TH1D * phiProjectionEtaCut = dEtaDPhi->ProjectionY("phiProjectionEtaCut", 25, 75);
	TF1 * nearSideEtaCutFit = fitPeak(phiProjectionEtaCut, nearSideLimits, "nearSideEtaCut");
	TF1 * awaySideEtaCutFit = fitPeak(phiProjectionEtaCut, awaySideLimits, "awaySideEtaCut");
	
	// Set fit parameters
	nearSideEtaCutFit->SetLineColor(kRed);
	awaySideEtaCutFit->SetLineColor(kBlue);

	// Draw projection and fit
	phiProjectionEtaCut->Draw();
	nearSideEtaCutFit->Draw("same");
	awaySideEtaCutFit->Draw("same");

	// Save projection
	canvas.Print("output/phiProjectionEtaCut.pdf");

	// Jet Hadron
	minMax nearSideJetHLimits(-TMath::Pi()/6, TMath::Pi()/6);
	minMax awaySideJetHLimits(3.0*TMath::Pi()/4, 5.0*TMath::Pi()/4);

	// Setup Canvas
	canvas.Clear();
	// -1 since we are ignoring 15-20 to recreate the STAR plots
	canvas.Divide(nParticlePtBins, nJetPtBins-1);
	index = 1;

	// Define variables
	TH1D * jetHProjection[nJetPtBins-1][nParticlePtBins];
	TF1 * nearSideJetHFit[nJetPtBins-1][nParticlePtBins];
	TF1 * awaySideJetHFit[nJetPtBins-1][nParticlePtBins];
	TString combinedNames;
	// Vectors for widths
	std::vector <double> pt10ASWidths;
	std::vector <double> pt20ASWidths;
	std::vector <double> pt10ASWidthsErr;
	std::vector <double> pt20ASWidthsErr;


	for (unsigned int i = 0; i < nJetPtBins; i++)
	{
		// We don't care about 15-20, so just skip it
		if (i == 1) { i++; }
		for (unsigned int j = 0; j < nParticlePtBins; j++)
		{
			// Set appropraite part of canvas
			canvas.cd(index);
			index++;

			// Perform projection
			combinedNames = jetHadron[i][j]->GetName();
			combinedNames += "_Projection";
			jetHProjection[i][j] = jetHadron[i][j]->ProjectionY(combinedNames.Data());
			
			// Perform fit
			combinedNames = jetHadron[i][j]->GetName();
			combinedNames += "_nsFit";
			nearSideJetHFit[i][j] = fitPeak(jetHProjection[i][j], nearSideJetHLimits, combinedNames.Data());
			combinedNames = jetHadron[i][j]->GetName();
			combinedNames += "_asFit";
			awaySideJetHFit[i][j] = fitPeak(jetHProjection[i][j], awaySideJetHLimits, combinedNames.Data());

			// Set fit parameters
			nearSideJetHFit[i][j]->SetLineColor(kRed);
			awaySideJetHFit[i][j]->SetLineColor(kBlue);

			// Save widths to array
			if (i == 0) { 
        pt10ASWidths.push_back(awaySideJetHFit[i][j]->GetParameter("Sigma")); 
        pt10ASWidthsErr.push_back(awaySideJetHFit[i][j]->GetParError(awaySideJetHFit[i][j]->GetParNumber("Sigma"))); 
      }
			if (i == 2) { 
        pt20ASWidths.push_back(awaySideJetHFit[i][j]->GetParameter("Sigma")); 
        pt20ASWidthsErr.push_back(awaySideJetHFit[i][j]->GetParError(awaySideJetHFit[i][j]->GetParNumber("Sigma"))); 
      }



			// Draw projection and fit
			jetHProjection[i][j]->Draw();
			nearSideJetHFit[i][j]->Draw("same");
			awaySideJetHFit[i][j]->Draw("same");
		}
	}
	// Save projections
	canvas.Print("output/jetHProjection.pdf");

	// Create vector for ptBins for TGraph
	std::vector <double> ptBinsForTGraph;	
	// Take the average of the two values for each bin location
	for (unsigned int i = 0; i < nParticlePtBins; i++) { ptBinsForTGraph.push_back((particlePtBins.at(i) + particlePtBins.at(i+1))/2); }

	// Graph for widths
//	TGraph * pt10Widths = new TGraph(nParticlePtBins, &ptBinsForTGraph[0], &pt10ASWidths[0]);
//	TGraph * pt20Widths = new TGraph(nParticlePtBins, &ptBinsForTGraph[0], &pt20ASWidths[0]);

	TGraphErrors * pt10Widths = new TGraphErrors(nParticlePtBins, &ptBinsForTGraph[0], &pt10ASWidths[0],0,&pt10ASWidthsErr[0]);
	TGraphErrors * pt20Widths = new TGraphErrors(nParticlePtBins, &ptBinsForTGraph[0], &pt20ASWidths[0],0,&pt20ASWidthsErr[0]);

	canvas.Clear();

  pt10Widths->SetMarkerSize(1);
  pt10Widths->SetMarkerStyle(21);
	pt10Widths->SetMarkerColor(kRed);
	pt10Widths->SetLineColor(kRed);
	pt10Widths->GetYaxis()->SetRangeUser(0.09, 1.4);
	pt10Widths->GetYaxis()->SetTitle("#sigma_{as}");
	pt10Widths->GetXaxis()->SetTitle("p_{T}^{associated}");
	pt10Widths->SetTitle("Awayside widths");

	pt10Widths->Draw("ALP");
  pt20Widths->SetMarkerSize(1);
  pt20Widths->SetMarkerStyle(21);
	pt20Widths->Draw("same LP");

	TLegend * leg = new TLegend(0.60, 0.65, 0.80, 0.85);

	leg->AddEntry(pt10Widths, "10 #leq p_{T}^{jet} #leq 15", "p");
	leg->AddEntry(pt20Widths, "20 #leq p_{T}^{jet} #leq 40", "p");

	leg->Draw("same");
	
	canvas.SetLogy(1);
	canvas.Print("output/widths.pdf");
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

	// Parse input
	TString inputFilename = "";
	TString outputFilename = "";
	if (argc == 3)
	{
		inputFilename = argv[1];
    outputFilename = argv[2];
	}
	else
	{
		if (argc != 1)
		{
			Printf("Must specifiy the input filename, or use the default by passing no arguments.");
			Printf("Defaults: inputFilename = \"yajemDefaultProfile\"");
			std::exit(1);
		}
	}
	
	if (inputFilename == "")
	{
		basicAnalysis();
	}
	else
	{
		basicAnalysis(inputFilename,outputFilename);
	}

	#ifdef RE_STYLE
	delete readableStyle;
	#endif

	return 0;
}

