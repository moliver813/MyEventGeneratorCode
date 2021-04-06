// phase1 : take root files, run jetfinding and cuts, make histograms
// 
// Can compile with:
//  make phase1
//
//  or
//
//  g++ -Wl,--no-as-needed -std=c++11 phase1.cc -o phase1 -I`root-config --incdir` `root-config --libs` `fastjet-config --cxxflags --libs --plugins`
//
//  or (if you have readableStyle.h)
//
//  make USER_DEFINED=-DRE_STYLE=1 basicAnalysis

#ifdef RE_STYLE
#include <readableStyle.h>
#endif

#define PI TMath::Pi()

// Set compact root files: use floats instead of doubles, etc.
// 
#define COMPACT_ROOT true

#define YAJEM false

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <iomanip>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>

#include <TPDGCode.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TIterator.h>
#include <TObject.h>
#include <TRandom3.h>

// removing SoftDrop parts
//#include "fastjet/contrib/SoftDrop.hh"

#include "analysis_params.h"

//using namespace fastjet::contrib;

// Functions are defined at the end of file
void printHelp(std::string input = "");
void parseInput(int argc, char * argv[], TString & inputFilename, TString & outputFilename, double & jet_constituent_cut, bool & ptHardBinFlag, TString & responseMatrixFilename, int & iTriggerType);

void AddToyModelParticles(int nGen, TRandom3 * rand, std::vector<fastjet::PseudoJet> &particles) {
	if (!rand) return;

//	float ToyEtaRange = 3; 
	float tau = 2.; // Tau for Pt Spectrum. Mean value? Desired Mean value (+ ptmin)
	float ptMinToy = 0.15; // Minimum pt for generated particles

	// could vary id of particles
	int iPID = 211; // pi+
	double fMass = 0.139; // mass pi+-

	float phi = 0;
	float eta = 0;
	float pt = 0;
	float pz = 0;
	float energy = 0;



  // to do
  // add randomized v3 angle
  // Can pick v3 angle within 1 third of circle (or anywhere)
  //double psi_3 = rand->Rndm() * 2.*TMath::Pi()/3.;
  double psi_3 = rand->Rndm() * 2.*TMath::Pi();


	for (int i = 0; i < nGen * ToyEtaRange; i++) {
		eta = (2.*rand->Rndm()-1)*ToyEtaRange;
		// this expression works for exponential distribution
		pt = ptMinToy + rand->Exp(tau); // sample from exp(-t/tau);

		pz = pt * TMath::SinH(eta);
		energy = TMath::Sqrt(fMass*fMass + TMath::Power(pt*TMath::CosH(eta),2));

		if (iFlowVersion < 0) {
			phi = rand->Rndm()*2*PI;
		} else {


			double v2 = 0;
      double v3 = 0;
			double v4 = 0;
	
			// VN Fit Parameters
			double V2FP_0 = 0;
			double V2FP_1 = 0;
			double V2FP_2 = 0;		
			double V4FP_0 = 0;
			double V4FP_1 = 0;
			double V4FP_2 = 0;		
			switch (iFlowVersion) {
        // Experimental error in comments
				case 3:
					V2FP_0 = 1.51740e+00;//   7.20128e-02
					V2FP_1 = 3.42402e+00;//   1.36342e-01
					V2FP_2 = 1.42951e+00;//   7.27495e-02
          V4FP_0 = 8.96183e-02;//   2.84786e-03
          V4FP_1 = 4.10369e+00;//   2.64685e-01
          V4FP_2 = 1.55910e+00;//   1.06870e-01 
				case 2:
					V2FP_0 = 1.35467e+00;//   4.47105e-02
					V2FP_1 = 3.51387e+00;//   1.64360e-01 
					V2FP_2 = 1.47053e+00;//   8.16601e-02
          V4FP_0 = 2.04272e-01;//   1.01123e-03
          V4FP_1 = 3.78135e+00;//   1.99824e-01
          V4FP_2 = 1.37851e+00;//   8.02073e-02 
				case 1:
					V2FP_0 = 1.04734e+00;//   1.46114e-02
					V2FP_1 = 3.72192e+00;//   1.90239e-01
					V2FP_2 = 1.55223e+00;//   8.94406e-02
          V4FP_0 = 3.26444e-01;//   9.77867e-03
          V4FP_1 = 3.53048e+00;//   1.72260e-01 
          V4FP_2 = 1.28028e+00;//   7.07583e-02
				case 0:
				default:
					V2FP_0 = 5.39182e-01;//   8.36928e-03
					V2FP_1 = 3.70690e+00;//   1.94541e-01
					V2FP_2 = 1.54252e+00;//   9.17354e-02
          V4FP_0 = 3.89071e-01;//   2.13836e-02
          V4FP_1 = 3.43959e+00;//   1.90255e-01
          V4FP_2 = 1.23961e+00;//   8.31609e-02
			}
			v2 = V2FP_0 * TMath::Landau(pt,V2FP_1,V2FP_2,false);
			v4 = V4FP_0 * TMath::Landau(pt,V4FP_1,V4FP_2,false);

      // v3 approximation
      //
      v3 = v2*v2;

			for (int z = 0; z < 30; z++) {
				phi = rand->Rndm()*2*PI;
				double pdf = 1. + 2.*v2*TMath::Cos(2.*phi) + 2. * TMath::Cos(3.*(phi - psi_3)) + 2.*v4*TMath::Cos(4.*phi);
				double testValue = 2.*rand->Rndm();
				if (testValue <= pdf) break;
				// After 30 tries, the phi will just be random
			}
		}

		fastjet::PseudoJet newParticle(pt*cos(phi),pt*sin(phi),pz,energy);

    // todo: randomly assign some toy particles to be pi0s
    // todo: apply a status for toy particles kToyParticleStatus
		newParticle.set_user_index(kToyParticleLabel); // label for Toy model particles

    float fPi0Fraction = 0.2; // total guess for pi0s proportion at highest pt
    if (rand->Rndm() < fPi0Fraction) newParticle.set_user_index(111);


		//printf("DEBUG: Toy adding particle with (px=%f,py=%f,pz=%f,E=%f)\n",newParticle.px(),newParticle.py(),newParticle.pz(),newParticle.e());
		particles.push_back(newParticle);
	}
}

// Function from Marco. Useful for handling pt hard bins with only a few counts
void removeOutliers(TH1* h1, Int_t minBinValue = 2) {
  // set bins with less than minBinValue entries to 0
  Int_t nbins = h1->GetNbinsX()+2;
  TH2* h2 = dynamic_cast<TH2*>(h1);
  if (h2 && h2->GetNbinsY() != 0)
    nbins *= h2->GetNbinsY()+2;
  TH3* h3 = dynamic_cast<TH3*>(h1);
  if (h3 && h3->GetNbinsZ() != 0)
    nbins *= h3->GetNbinsZ()+2;
  std::cout << "hist " << h1->GetName() << " nbins " << nbins << std::endl;
  for (Int_t ibin = 0; ibin < nbins; ibin++) {
    if (h1->GetBinContent(ibin) < minBinValue && h1->GetBinContent(ibin) > 0) {
      //cout << "Histo: " << h1->GetName() << " setting bin " << ibin << " to zero"<< endl;
      h1->SetBinContent(ibin,0);
    }
  }
}

// Accept particles within abs(eta) < 1, pt > 0.2 GeV
bool rejectParticle(const fastjet::PseudoJet & particle)
{
  return !( (particle.pt() > particle_pt_min) && (TMath::Abs(particle.eta()) < eta_cut) );
}


bool rejectJet(const fastjet::PseudoJet & jet)
{
  bool rejectFlag = ! jet.has_constituents();
  // Eta cut on jets. We only want jets within eta-R 
  if (requireHardCore == true) { 
 
    //FIXME FIXME
    // don't remember why I put two of these tags here, but I'd better leave them
    // to be on the safe side 
    rejectFlag = true; 

    std::vector <fastjet::PseudoJet> particles = jet.constituents();
    for (unsigned int i = 0; i < particles.size(); i++)
    {
      if ( particles.at(i).pt() >= globalHardCoreCut) { rejectFlag = false; break; }
      //      if ( particles.at(i).pt() >= jet_constituent_cut) { rejectFlag = false; break; }
    }
  }
  if ((rejectFlag == false) && (TMath::Abs(jet.eta()) >= trigger_eta_cut)) { rejectFlag = true; }

  return rejectFlag;
}

//bool rejectJetHardCore(const fastjet::PseudoJet & jet, double R, double cut)
bool rejectJetHardCore(fastjet::PseudoJet * jet, double R, double cut)
{
  bool rejectFlag = true;
  if (jet->has_constituents() == true)
  {
    std::vector <fastjet::PseudoJet> particles = jet->constituents();
    for (unsigned int i = 0; i < particles.size(); i++)
    {
      if ( particles.at(i).pt() >= cut) { rejectFlag = false; break; }
      //      if ( particles.at(i).pt() >= jet_constituent_cut) { rejectFlag = false; break; }
    }
  }
  // Eta cut on jets. We only want jets within eta-R 
  if ((rejectFlag == false) && (TMath::Abs(jet->eta()) >= trigger_eta_cut)) { rejectFlag = true; }

  return rejectFlag;
}

// method for Dijet.  reject jets whose hard (p_t > 2 GeV/c) constituents
// add up to less than 10 (or 20) GeV/c

//bool rejectJetDijet(const fastjet::PseudoJet & jet, double R, double sum_cut)
bool rejectJetDijet(const fastjet::PseudoJet * jet, double R, double sum_cut)
{
  bool rejectFlag = true;
  double sum_hard = 0;
  if (jet->has_constituents() == true)
  {
    std::vector <fastjet::PseudoJet> particles = jet->constituents();
    for (unsigned int i = 0; i < particles.size(); i++)
    {
      if (particles.at(i).pt() >= 2) sum_hard += particles.at(i).pt();

      //      if ( particles.at(i).pt() >= jet_constituent_cut) { rejectFlag = false; break; }
    }
  }
  if (sum_hard > sum_cut) rejectFlag = false;

  // Eta cut on jets. We only want jets within eta-R 
  if ((rejectFlag == false) && (TMath::Abs(jet->eta()) >= (eta_cut-R))) { rejectFlag = true; }

  return rejectFlag;
}

double getHighestPtConst(const fastjet::PseudoJet & jet) {  
  double highest = 0;
  if (jet.has_constituents())
  {
    std::vector <fastjet::PseudoJet> particles = jet.constituents();
    for (unsigned int i = 0; i < particles.size(); i++)
    {
      highest = TMath::Max(highest,particles.at(i).pt());
    }
  }
  return highest;
}


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

// Take phi of jet, return 1 for in-plane, 2 for mid-plane, 3 for out-of-plane
// Bin 0 is taken to be have all events
int GetEventPlaneBin(const double phi) {
  double fPhi = fmod(phi,PI);
  if ( (fPhi <= PI/6.) || (PI-fPhi <= PI/6.) ) return 1;
  if ( (fPhi <= PI/3.) || (PI-fPhi <= PI/3.) ) return 2;
    return 3;
//  if ( (abs(phi) <= PI/6.) || (abs(PI-phi) <= PI/6.) ) return 1;
//  if ( (abs(phi) <= PI/3.) || (abs(PI-phi) <= PI/3.) ) return 2;
//    return 3;
}


// Take Particle level jet pt, and calculate weights proportional to the probability that it would give a detector level jet
std::vector<double> FindJetReponseWeights(double jetPt, TH2F * responseMatrix) {
  std::vector<double> weights = {};

  if (!responseMatrix) {
    fprintf(stderr,"Missing ResponseMatrix!\n");
    return weights;
  }  
  int particleLevelBin = responseMatrix->GetYaxis()->FindBin(jetPt);

  TH1F * detectorPtRow = (TH1F *) responseMatrix->ProjectionX("detectorPtRow",particleLevelBin,particleLevelBin);
  if (!detectorPtRow) {
    fprintf(stderr,"Something went wrong with projection\n");
    return weights;
  }

  // FIXME normalize detectorPtRow. 
  double integral = detectorPtRow->Integral(1,detectorPtRow->GetNbinsX(),"width");
  // FIXME it is more efficient to scale the final weights than the histogram, right?
  if (integral > 0) detectorPtRow->Scale(1./integral);

  for (int i = 0; i < nJetPtBins; i++) {
    double localWeight = 0.;
    
    detectorPtRow->GetXaxis()->SetRangeUser(jetPtBins[i],jetPtBins[i+1]);
    localWeight = detectorPtRow->Integral("width");

    weights.push_back(localWeight);
  }

  return weights;
}


  // particles
  /** list
   kGamma,             // Photon
   kElectron,          // Electron
   kMuonPlus,          // Muon 
   kPiPlus,            // Pion
   kKPlus,             // Kaon
   kK0Short,           // K0s
   kK0Long,            // K0l
   kProton,            // Proton 
   kNeutron,           // Neutron
   kLambda0,           // Lambda_0
   kSigmaMinus,        // Sigma Minus
   kSigmaPlus,         // Sigma Plus
   3312,               // Xsi Minus 
   3322,               // Xsi 
   3334,               // Omega
   kNuE,               // Electron Neutrino 
   kNuMu,              // Muon Neutrino
   kNuTau              // Tau Neutrino


   pi0??
   kPi0
   */

bool isNeutrino(int pid) {
  pid = TMath::Abs(pid);
  return pid == kNuE || pid == kNuMu || pid == kNuTau;
}


//It would probably be faster to use define statements for these functions,
// to minimize function calls during operation
bool isChargedHadron(int pid) {
    /*  list of stable particles:
      22,             // Photon
      11,          // Electron
      13,          // Muon 
      ,            // Pion
      kKPlus,             // Kaon
      kK0Short,           // K0s
      kK0Long,            // K0l
      kProton,            // Proton 
      kNeutron,           // Neutron
      kLambda0,           // Lambda_0
      kSigmaMinus,        // Sigma Minus
      kSigmaPlus,         // Sigma Plus
      3312,               // Xsi Minus 
      3322,               // Xsi 
      3334,               // Omega
      kNuE,               // Electron Neutrino 
      kNuMu,              // Muon Neutrino
      kNuTau              // Tau Neutrino
      */

  pid = TMath::Abs(pid);
  return pid == kProton || pid == kPiPlus || pid == kKPlus || pid == kSigmaPlus || pid == 3312 || pid == 3334;
}

bool isChargedParticle(int pid) {
  pid = TMath::Abs(pid);
  return (pid == 11 || pid == 13 || pid == kProton || pid == kPiPlus || pid == kKPlus || pid == kSigmaPlus || pid == 3312 || pid == 3334);
}

bool isGamma(int pid) {
  return pid == kGamma;
}

bool isNeutralHadron(int pid) {
  pid = TMath::Abs(pid);
  return pid==kNeutron || pid == kPi0 || pid == kK0Short || pid == kK0Long || pid == kLambda0 || pid == 3322 ;
}

/*
 * Identifies long lived neutral hadrons that ALICE cannot measure
 */
bool isNeutralLongHadron(int pid) {
  pid = TMath::Abs(pid);
  return pid==kNeutron || pid == kK0Long ;
}

// get background from perpendicular cone away from jet.
double getBackground(fastjet::PseudoJet *thisJet, std::vector <fastjet::PseudoJet> *particles) {
  double rho = 0;
  double area = 0;
  double dEtaRange =  0.2;
  double dPhiRange =  PI / 4.0;

  for (int i = 0; i < particles->size(); i++) {
    fastjet::PseudoJet *part = &particles->at(i);
    // printf("particle at (eta,phi) = (%f,%f)\t\t",part->eta(),part->phi());
    //rho+= particles->at(i).pt();
    double dEta = part->eta() - thisJet->eta();
    double dPhi = part->delta_phi_to(*thisJet);
    //  printf(" DEta,DPhi = %f,%f\t",dEta,dPhi);      
    if (TMath::Abs(dEta) > dEtaRange || TMath::Abs(TMath::Abs(dPhi)- PI/2) > dPhiRange || part->pt() < jet_constituent_cut) {
      //   printf("\n");
      continue;
    }
    //  printf("Accepted!\n");
    rho+= part->pt();
  }
  // Divide by area
  rho = rho / (8 * dEtaRange * dPhiRange);
  return rho;
}

void phase1(TString inputFilename, TString outputFilename, double jet_constituent_cut, bool ptHardBinFlag, TString responseMatrixFilename, int iTriggerType)
{
  TRandom3 * rand;

  switch (iTriggerType) {
    case 0: 
      printf("Will build correlations using Jet triggers.\n");
      break;
    case 1:
      printf("Will build correlations using Pi0 triggers.\n");
      break;
    default:
      printf("Unknown Trigger type\n");
  }

  // Load Response Matrix
//    TH2F * responseMatrix = 0;
  TH2F * responseMatrix[nEPBins];
  // the old response matrices
  //const char *responseMatrixNames[4] = {"JESCorrection_Clus6.00_all","JESCorrection_Clus6.00_inPlane","JESCorrection_Clus6.00_midPlane","JESCorrection_Clus6.00_outOfPlane"};
  const char *responseMatrixNames[4] = {"responseMatrix_inclusive","responseMatrix_in_plane","responseMatrix_mid_plane","responseMatrix_out_of_plane"};

  bool usingResponseMatrix = false;

  if (responseMatrixFilename != "") {
    usingResponseMatrix = true;
    printf("Attempting to load RM ...\n");
    printf("Attempting to load response Matrix %s\n",responseMatrixFilename.Data());
    TFile * fResponseMatrix = TFile::Open(responseMatrixFilename.Data());
  
    // Check if file is open
    if (fResponseMatrix->IsOpen() == kFALSE)
    {
      fprintf(stderr,"Response Matrix File \"%s\" is not found!\n", responseMatrixFilename.Data());
      exit(1);
    }

    //responseMatrix = (TH2F *) fResponseMatrix->Get("JESCorrection_Clus6.00_all"); // FIXME update name 
    for (int i = 0; i < nEPBins; i++) {
      responseMatrix[i] = (TH2F *) fResponseMatrix->Get(responseMatrixNames[i]); // FIXME update name 
      if (!responseMatrix[i]) {
        fprintf(stderr,"Could not find %s in file %s\n",responseMatrixNames[i],responseMatrixFilename.Data());
        exit(1);
      }
    }
  } else {
     printf("Not using a response matrix for smearing.\n");
  }
  

  // Open input file
  //  inputFilename = Form("output/%s.root", inputFilename.Data());
  TFile * fIn = TFile::Open(inputFilename.Data());

  // Check if file is open
  if (fIn->IsOpen() == kFALSE)
  {
    Printf("inputFilename \"%s\" is not found!", inputFilename.Data());
    std::exit(1);
  }

	// Initiating random object, using file hash for reproducibility
  if (phiAcceptanceSim || nToyParticles > 0) {
    rand = new TRandom3((UInt_t) fIn->Hash()); // casting ULong_t to UInt_t 
  }

  // The name of the tree is "T"
  TTree * tree = (TTree *) fIn->Get("T");

  TTree * runInfo = (TTree *) fIn->Get("runInfo");

  if (!tree) {
    fprintf(stderr,"Error: Tree \"T\" not found\n");
    std::exit(1);
  }


  // Define our basic variables
  // Basic types don't need to be pointers, but classes do
  //Int_t eventNumber = 0;
  //Int_t nParticles = 0;

  // Need to do something about type?  dynamically chech type in branch?
  // 
  //

#if COMPACT_ROOT 
  std::vector <short> * particleID = 0;
  std::vector <char> * status = 0;
  std::vector <float> * px = 0;
  std::vector <float> * py = 0;
  std::vector <float> * pz = 0;
  std::vector <float> * energy = 0;
  std::vector <float> * mass = 0;

  float vertex_x = 0;
  float vertex_y = 0;
#else
  std::vector <int> * particleID = 0;
  std::vector <int> * status = 0;
  std::vector <double> * px = 0;
  std::vector <double> * py = 0;
  std::vector <double> * pz = 0;
  std::vector <double> * energy = 0;
  std::vector <double> * mass = 0;

  double vertex_x = 0;
  double vertex_y = 0;
#endif


  // for weighted data:
#if YAJEM
  float weight = 0;
#else 
  double weight = 0;
#endif

  double total_weight = 0;

  printf("Beginning to set branch addresses\n");

  // Set Branches to our variables
  //tree->SetBranchAddress("eventNumber", &eventNumber);
  //tree->SetBranchAddress("nParticles", &nParticles);
  tree->SetBranchAddress("particleID", &particleID);
  tree->SetBranchAddress("px", &px);
  tree->SetBranchAddress("py", &py);
  tree->SetBranchAddress("pz", &pz);
  tree->SetBranchAddress("energy", &energy);
  tree->SetBranchAddress("mass", &mass);

  printf("Checking if vertices are on... \n");

  bool vertexOn = 0 != tree->GetBranch("vertex_x");
  if (vertexOn) {
    tree->SetBranchAddress("vertex_x",&vertex_x);
    tree->SetBranchAddress("vertex_y",&vertex_y);
  }
  bool statusOn = 0 != tree->GetBranch("status");
  if (statusOn) tree->SetBranchAddress("status",&status);

/*    printf("Setting up soft drop\n");
    //soft drop analysis
    double z_cut = 0.10;
    double beta  = 0.0;
    SoftDrop sd(beta,z_cut);
*/
    //checking whether weighting is on.  For YaJEM output, weighting might 
  // not be present.
  // TODO: Make an argument
  bool weightingOn = 0 != tree->GetBranch("weight");
  // FIXME just for hard spectrum test
  // weightingOn = 0;


  if (weightingOn) {
    printf("Weighting is on.\n");    
    tree->SetBranchAddress("weight",&weight);
  } else {
    if ( ptHardBinFlag == true )
    {
      printf("Pt Hard Weighting is on. Weighting is set per run.\n");

      // Retreive run info
      // Define temporary variables to calculate weight
      double crossSection = 0;
      Int_t nEvents = 0;
      double totalCrossSection = 0;
      double totalEvents = 0;
      // For this case, it is weighted per run, so the "totalWeight" branch is just the weight for each run.
      // TFileMerger cannot sum up the values itself, so it is calculated here.
      runInfo->SetBranchAddress("crossSection", &crossSection);
      runInfo->SetBranchAddress("nEvents", &nEvents);
      for (Int_t i = 0; i < runInfo->GetEntries(); i++)
      {
        runInfo->GetEntry(i);
        totalCrossSection += crossSection;
        totalEvents += nEvents;
      }
      weight = totalCrossSection / ((Double_t) totalEvents);
    }
    else
    {
      printf("Weighting is off.  All events set to equal weight.\n");
      weight = 1;
    }
  }
  if (iTriggerType == 0) printf("Jet resolution parameter R = %f\n",R);

  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, R);
  double ghost_maxrap= 1.5;
  double ghost_area = 0.01;


  int nTriggerPtBins = nJetPtBins; // Default for Jet-h analysis
  std::vector<double> fTriggerPtBins = jetPtBins;
  if (iTriggerType == 1) {
    printf("Using Bins for Pi0-hadron analysis\n");
    nTriggerPtBins = nPi0PtBins;
    fTriggerPtBins = pi0PtBins;
  }
  printf("Will use trigger pt bins: ");
  for (int i = 0; i < nTriggerPtBins; i++) {
    printf("%.2f, ",fTriggerPtBins.at(i));
  }
  printf("%.2f\n",fTriggerPtBins.back());
	

  // Use either TLorentzVector or fastjet::PseudoJet for the particles
  // Need to use fastjet::PesudoJet for jets regardless
  //std::vector <TLorentzVector> particles;
  //vector for all particles 
  std::vector <fastjet::PseudoJet> particles;
  // subset of particles that are scattering centers from JEWEL
  std::vector <fastjet::PseudoJet> scatCents; 



  // vector for particles to be used for jet reconstruction
  std::vector <fastjet::PseudoJet> constParticles;
  std::vector <fastjet::PseudoJet> jets;
  //std::vector <fastjet::PseudoJet> rec_jets;

  std::vector <fastjet::PseudoJet> hardScatters;

  TString triggerType = "";
  TString jetClass = "";
  TString particleClass = "particlePt_%.2f_%.2f";
  if (iTriggerType == 0) {
    triggerType = "jet";
    jetClass = "jetPt_%.0f_%.0f";
  } else {
    triggerType = "pi0";
    jetClass = "pi0Pt_%.0f_%.0f";
    if (bUseZtBins) particleClass = "zT_%.2f_%.2f";
  }
  TString combinedClass = "";

  // Define histograms
  TH2F * dEtaDPhi = new TH2F("dEtaDPhi", "dEtaDPhi", 100, -2*eta_cut, 2*eta_cut, 100, -TMath::Pi()/2, 3.0*TMath::Pi()/2);
  TH2F * etaPhi = new TH2F("etaPhi", "etaPhi", 100, -5, 5, 100, 0, 2.0*TMath::Pi());

  TH2F * etaPt = new TH2F("etaPt", "etaPt", 100, -5, 5, (int) nTrackPtBins, (double *) &trackPtBins[0]);
  etaPt->Sumw2();
//    TH2F * etaPt = new TH2F("etaPt", "etaPt", 100, -5, 5, 100, 0, 100);
  TH2F * phiPt = new TH2F("phiPt", "p_{T} vs #phi;#phi;p_{T}", 100, 0, 2.0*TMath::Pi(), (int) nTrackPtBins, (double *) &trackPtBins[0]);
//  TH2F * phiPt = new TH2F("phiPt", "p_{T} vs #phi;#phi;p_{T}", 100, 0, 2.0*TMath::Pi(), 100, 0, 100);
  phiPt->Sumw2();

  TH1F * dijetAj = new TH1F("dijetAj","A_{j};A_{j}",100,0,1); dijetAj->Sumw2();
  TH2F * dijetAjJetPt = new TH2F("dijetAjJetPt","A_{j} vs p_{t}^{jet};p_{t}^{jet};A_{j}",(int) nJetPtBins, (double *) &jetPtBins[0],100,0,1); dijetAjJetPt->Sumw2();
  TH1F * dijetXj = new TH1F("dijetXj","X_{j};X_{j}",100,0,1); dijetXj->Sumw2();
  TH2F * dijetXjJetPt = new TH2F("dijetXjJetPt","x_{j} vs p_{t}^{jet};p_{t}^{jet};x_{j}",(int) nJetPtBins, (double *) &jetPtBins[0],100,0,1); dijetXjJetPt->Sumw2();

  TH1F * vertexR = new TH1F("vertexR", "vertexR", 100, 0, 15);
  TH1F * vertexRHardCore = new TH1F("vertexRHardCore", "vertexRHardCore", 100, 0, 15);
  TH1F * vertexRNoHardCore = new TH1F("vertexRNoHardCore", "vertexNoHardCore", 100, 0, 15);
  // testing the leading and subleading cut's effect.
  TH1F * vertexRDijetCut = new TH1F("vertexRDijetCut","vertexRDijetCut;R_{vertex} (fm)",100,0,15);
  TH2F * vertexXY = new TH2F("vertexXY", "vertexXY", 100, -10, 10, 100, -10, 10);
  TH2F * vertexXYHardCore = new TH2F("vertexXYHardCore", "vertexXYHardCore", 100, -10, 10, 100, -10, 10);
  TH2F * vertexXYNoHardCore = new TH2F("vertexXYNoHardCore", "vertexXYNoHardCore", 100, -10, 10, 100, -10, 10);
  //vertex rotated so that the direction of the jet is -x
  // as in the plots made by T. Renk for his jet-h corr theory paper
  TH2F * vertexXYRot = new TH2F("vertexXYRot", "vertexXYRot", 100, -10, 10, 100, -10, 10);
  TH2F * vertexXYRotHardCore = new TH2F("vertexXYRotHardCore", "vertexXYRotHardCore", 100, -10, 10, 100, -10, 10);
  TH2F * vertexXYRotNoHardCore = new TH2F("vertexXYRotNoHardCore", "vertexXYRotNoHardCore", 100, -10, 10, 100, -10, 10);
  // binned by jet bin pt
  //  TH3F * vertexXYRotJetPt = new TH3F("vertexXYRotJetPt", "vertexXYRotJetPt", 100, -10, 10, 100, -10, 10,(int) nJetPtBins, (double *) &jetPtBins[0]);
  //  TH3F * vertexXYRotJetPtHardCore = new TH3F("vertexXYRotJetPtHardCore", "vertexXYRotJetPtHardCore", 100, -10, 10, 100, -10, 10,(int) nJetPtBins, (double *) &jetPtBins[0]);

  // pt and energy loss     
  TH2F * ptLossVertexRLeadingJet = new TH2F("ptLossVertexRLeadingJet","Leading Jet p_{T} Loss vs Hard Scatter R;R_{h.s.} (fm);p_{T}^{parton} - p_{T}^{jet} (GeV/c)",25,0,7.5,50,0,50);
  TH2F * ptLossVertexRSubleadingJet = new TH2F("ptLossVertexRSubleadingJet","Sub-Leading Jet p_{T} Loss vs Hard Scatter R;R_{h.s.} (fm);p_{T}^{parton} - p_{T}^{jet} (GeV/c)",25,0,7.5,50,0,50);


  // Use trigger (jet or pi0) for JetPtBin hists and jet bins for PartonPtBin
  TH1F * ptLossLeadingJetPtBinEP[nEPBins][nTriggerPtBins];
  TH1F * energyLossLeadingJetPtBinEP[nEPBins][nTriggerPtBins];
  TH1F * ptLossSubleadingJetPtBinEP[nEPBins][nTriggerPtBins];
  TH1F * energyLossSubleadingJetPtBinEP[nEPBins][nTriggerPtBins];

  TH1F * ptLossLeadingPartonPtBinEP[nEPBins][nJetPtBins];
  TH1F * energyLossLeadingPartonPtBinEP[nEPBins][nJetPtBins];
  TH1F * ptLossSubleadingPartonPtBinEP[nEPBins][nJetPtBins];
  TH1F * energyLossSubleadingPartonPtBinEP[nEPBins][nJetPtBins];

  //if (iTriggerType == 0) {
    for (int k = 0; k < nEPBins; k++) {
      for (int i = 0; i < nTriggerPtBins; i++) {
        ptLossLeadingJetPtBinEP[k][i] = new TH1F(Form("ptLossLeadingJetPtBinEP_%d_" + jetClass,k,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)),Form("p_{T} Loss (%.1f < p_{T}^{jet} < %.1f);p_{T}^{parton} - p_{T}^{jet} (GeV/c)",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)),125,-5,120);
        energyLossLeadingJetPtBinEP[k][i] = new TH1F(Form("energyLossLeadingJetPtBinEP_%d_" + jetClass,k,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)),Form("E Loss (%.1f < p_{T}^{jet} < %.1f);E^{parton} - E^{jet} (GeV)",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)),125,-5,120);
        ptLossLeadingJetPtBinEP[k][i]->Sumw2();
        energyLossLeadingJetPtBinEP[k][i]->Sumw2();

        ptLossSubleadingJetPtBinEP[k][i] = new TH1F(Form("ptLossSubleadingJetPtBinEP_%d_" + jetClass,k,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)),Form("p_{T} Loss (%.1f < p_{T}^{jet} < %.1f);p_{T}^{parton} - p_{T}^{jet} (GeV/c)",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)),125,-5,120);
        energyLossSubleadingJetPtBinEP[k][i] = new TH1F(Form("energyLossSubleadingJetPtBinEP_%d_" + jetClass,k,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)),Form("E Loss (%.1f < p_{T}^{jet} < %.1f);E^{parton} - E^{jet} (GeV)",fTriggerPtBins.at(i),fTriggerPtBins.at(i+1)),125,-5,120);

        ptLossSubleadingJetPtBinEP[k][i]->Sumw2();
        energyLossSubleadingJetPtBinEP[k][i]->Sumw2();
      }

      for (int i = 0; i < nJetPtBins; i++) {
        ptLossLeadingPartonPtBinEP[k][i] = new TH1F(Form("ptLossLeadingPartonPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1)),Form("p_{T} Loss (%.1f < p_{T}^{Parton} < %.1f);p_{T}^{parton} - p_{T}^{jet} (GeV/c)",jetPtBins.at(i),jetPtBins.at(i+1)),125,-5,120);
        energyLossLeadingPartonPtBinEP[k][i] = new TH1F(Form("energyLossLeadingPartonPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1)),Form("E Loss (%.1f < p_{T}^{Parton} < %.1f);E^{parton} - E^{jet} (GeV)",jetPtBins.at(i),jetPtBins.at(i+1)),125,-5,120);

        ptLossLeadingPartonPtBinEP[k][i]->Sumw2();
        energyLossLeadingPartonPtBinEP[k][i]->Sumw2();

        ptLossSubleadingPartonPtBinEP[k][i] = new TH1F(Form("ptLossSubleadingPartonPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1)),Form("p_{T} Loss (%.1f < p_{T}^{Parton} < %.1f);p_{T}^{parton} - p_{T}^{jet} (GeV/c)",jetPtBins.at(i),jetPtBins.at(i+1)),125,-5,120);
        energyLossSubleadingPartonPtBinEP[k][i] = new TH1F(Form("energyLossSubleadingPartonPtBinEP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1)),Form("E Loss (%.1f < p_{T}^{Parton} < %.1f);E^{parton} - E^{jet} (GeV)",jetPtBins.at(i),jetPtBins.at(i+1)),125,-5,120);

        ptLossSubleadingPartonPtBinEP[k][i]->Sumw2();
        energyLossSubleadingPartonPtBinEP[k][i]->Sumw2();
      }  
    }
  //}

  // Hard Scatter information
  TH1F * PartonPtByJetBin[nEPBins][nTriggerPtBins];
  TH1F * qqbar_PartonPtByJetBin[nEPBins][nTriggerPtBins];
  TH1F * qq_PartonPtByJetBin[nEPBins][nTriggerPtBins];
  // Gluon as leading jet
  TH1F * gq_PartonPtByJetBin[nEPBins][nTriggerPtBins];
  // Quark as leading jet
  TH1F * qg_PartonPtByJetBin[nEPBins][nTriggerPtBins];
  TH1F * gg_PartonPtByJetBin[nEPBins][nTriggerPtBins];

  TH2F * vertexXYJetPt[nTriggerPtBins];
  TH2F * vertexXYRotJetPt[nTriggerPtBins];
  TH2F * vertexXYRotJetPtHardCore[nTriggerPtBins];
  //if (iTriggerType == 0) {
    for (unsigned int i = 0; i < nTriggerPtBins; i++)
    {
      combinedClass = Form("vertexXYJetPt_" + jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
      vertexXYJetPt[i] = new TH2F(combinedClass,combinedClass, 100,-10,10,100,-10,10);
      vertexXYJetPt[i]->GetXaxis()->SetTitle("x (fm)");
      vertexXYJetPt[i]->GetYaxis()->SetTitle("y (fm)");

      // could update combinedClass paradigm to make labeling axes easier
      combinedClass = Form("vertexXYRotJetPt_" + jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
      vertexXYRotJetPt[i] = new TH2F(combinedClass,combinedClass, 100,-10,10,100,-10,10);
      vertexXYRotJetPt[i]->GetXaxis()->SetTitle("x (fm)");
      vertexXYRotJetPt[i]->GetYaxis()->SetTitle("y (fm)");

      if (iTriggerType == 0) {
        combinedClass = Form("vertexXYRotJetPtHardCore_" + jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
        vertexXYRotJetPtHardCore[i] = new TH2F(combinedClass,combinedClass, 100,-10,10,100,-10,10);

        vertexXYRotJetPtHardCore[i]->GetXaxis()->SetTitle("x (fm)");
        vertexXYRotJetPtHardCore[i]->GetYaxis()->SetTitle("y (fm)");
      }
  //PartonPtByJetBin

    }
  //}

  TH2F * leadingJetPtNScatCent[nEPBins];
  //if (iTriggerType == 0) {
    for (unsigned int k = 0; k < nEPBins; k++) {

//      for (unsigned int i = 0; i < nJetPtBins; i++) {
      for (unsigned int i = 0; i < nTriggerPtBins; i++) {

        // Hard Scatter inaformation
  //      combinedClass = Form(jetClass,jetPtBins.at(i),jetPtBins.at(i+1));
        //combinedClass = Form("EP_%d_" + jetClass,k,jetPtBins.at(i),jetPtBins.at(i+1));
        combinedClass = Form("EP_%d_" + jetClass,k,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
        PartonPtByJetBin[k][i] = new TH1F(Form("PartonPtByJetBin_%s",combinedClass.Data()),Form("qqbar_PartonPtByJetBin_%s;p_{T}^{jet} (GeV/c)",combinedClass.Data()), 200,0,200);
        PartonPtByJetBin[k][i]->Sumw2();
        qqbar_PartonPtByJetBin[k][i] = new TH1F(Form("qqbar_PartonPtByJetBin_%s",combinedClass.Data()),Form("qqbar_PartonPtByJetBin_%s;p_{T}^{jet} (GeV/c)",combinedClass.Data()), 200,0,200);
        qqbar_PartonPtByJetBin[k][i]->Sumw2();
        qq_PartonPtByJetBin[k][i]  = new TH1F(Form("qq_PartonPtByJetBin_%s",combinedClass.Data()),Form("qq_PartonPtByJetBin_%s;p_{T}^{jet} (GeV/c)",combinedClass.Data()), 200,0,200);
        qq_PartonPtByJetBin[k][i]->Sumw2();
        // Gluon as leading jet
        gq_PartonPtByJetBin[k][i]  = new TH1F(Form("gq_PartonPtByJetBin_%s",combinedClass.Data()),Form("gq_PartonPtByJetBin_%s;p_{T}^{jet} (GeV/c)",combinedClass.Data()), 200,0,200);
        gq_PartonPtByJetBin[k][i]->Sumw2();
        // Quark as leading jet
        qg_PartonPtByJetBin[k][i]  = new TH1F(Form("qg_PartonPtByJetBin_%s",combinedClass.Data()),Form("qg_PartonPtByJetBin_%s;p_{T}^{jet} (GeV/c)",combinedClass.Data()), 200,0,200);
        qg_PartonPtByJetBin[k][i]->Sumw2();
        gg_PartonPtByJetBin[k][i]  = new TH1F(Form("gg_PartonPtByJetBin_%s",combinedClass.Data()),Form("gg_PartonPtByJetBin_%s;p_{T}^{jet} (GeV/c)",combinedClass.Data()), 200,0,200);
        gg_PartonPtByJetBin[k][i]->Sumw2();
    
      }
    }
    for (int k = 0; k < nEPBins; k++) {

      leadingJetPtNScatCent[k] = new TH2F(Form("leadingJetPtNScatCent_EP_%d",k),"Number of JEWEL Scattering Centers vs Leading Jet p_{T};p_{T}^{jet} (GeV/c);N_{Scat. Cent.}",(int) nTriggerPtBins, (double *) &fTriggerPtBins[0],250,0,1000);
      //leadingJetPtNScatCent[k] = new TH2F(Form("leadingJetPtNScatCent_EP_%d",k),"Number of JEWEL Scattering Centers vs Leading Jet p_{T};p_{T}^{jet} (GeV/c);N_{Scat. Cent.}",(int) nJetPtBins, (double *) &jetPtBins[0],250,0,1000);
      leadingJetPtNScatCent[k]->Sumw2();

    }
  //}

  TH2F * vertexXYRotDijet = new TH2F("vertexXYRotDijet", "vertexXYRotDijet", 100, -10, 10, 100, -10, 10);
  TH2F * vertexXRotJetPt = new TH2F("vertexXRotJetPt","x_{rot} vs p_{t}^{jet};p_{t}^{jet} (GeV/c);x_{rot} (fm)",(int) nTriggerPtBins, (double *) &fTriggerPtBins[0],100,-15,15);
  TH2F * vertexXRotJetPtHardCore = new TH2F("vertexXRotJetPtHardCore","x_{rot} vs p_{t}^{jet} (Hard Core cut);p_{t}^{jet} (GeV/c);x_{rot} (fm)",(int) nTriggerPtBins, (double *) &fTriggerPtBins[0],100,-15,15);
  TH2F * vertexHighestPtXRot = new TH2F("vertexHighestPtXRot","vertexHighestPtXRot;p_{T}^{HC};x (fm)",100,0,100,100,-15,15);  vertexHighestPtXRot->Sumw2();


  // vertexR as function of jetPt:
  TH2F * vertexJetPtR = new TH2F("vertexJetPtR","vertexJetPtR;p_{T}^{jet} (GeV/c);R (fm)",(int) nTriggerPtBins, (double *) &fTriggerPtBins[0],100,0,15);
//  TH2F * vertexJetPtR = new TH2F("vertexJetPtR","vertexJetPtR;p_{T}^{jet} (GeV/c);R (fm)",(int) nTriggerPtBins, (double *) &jetPtBins[0],100,0,15);
  TH2F * vertexJetPtRHardCore = new TH2F("vertexJetPtRHardCore","vertexJetPtRHardCore;p_{t}^{jet} (GeV/c);R (fm)",(int) nTriggerPtBins, (double *) &fTriggerPtBins[0],100,0,15);
  TH2F * vertexJetPtRDijet = new TH2F("vertexJetPtRDijet","vertexJetPtRDijet;p_{t}^{jet} (GeV/c);R (fm)",(int) nTriggerPtBins, (double *) &fTriggerPtBins[0],100,0,15); vertexJetPtRDijet->Sumw2();
  TH2F * vertexDijetAjR = new TH2F("vertexDijetAjR","R_{vertex} vs A_{j}^{Dijet}",100,0,1,25,0,15); vertexDijetAjR->Sumw2();


  TH2F * vertexHighestPtR = new TH2F("vertexHighestPtR","vertexHighestPtR;p_{T}^{HC}",100,0,100,100,0,15);  vertexHighestPtR->Sumw2();
  TH3F * vertexXYRotAjDijet = new TH3F("vertexXYRotAjDijet","vertexXYRotAjDijet;x (fm);y (fm); A_{j}",100,-10,10,100,-10,10,5,0,1);  vertexXYRotAjDijet->Sumw2();

  TH2F * vertexHighestJetZR = new TH2F("vertexHighestJetZR","Vertex R vs Highest z;z;R (fm)",100,0,1,100,0,15);
  TH2F * vertexHighestJetZXRot = new TH2F("vertexHighestJetZXRot","X_{rot} vs Highest z;z;x (fm)",100,0,1,100,0,15);

  // Studies for the Hadron-Lambda/k0 study
  TH2F * vertexXYRot4GeVTrigger = new TH2F("vertexXYRot4GeVTrigger","vertexXYRot4GeVTrigger;x (fm); y (fm)",100,-10,10,100,-10,10);
  TH2F * vertexXYRot5GeVTrigger = new TH2F("vertexXYRot5GeVTrigger","vertexXYRot5GeVTrigger;x (fm); y (fm)",100,-10,10,100,-10,10);
  TH2F * vertexXYRot20GeVTrigger = new TH2F("vertexXYRot20GeVTrigger","vertexXYRot20GeVTrigger;x (fm); y (fm)",100,-10,10,100,-10,10);

  // General Study of Hadron Trigger
  TH2F * vertexXRotTriggerPt = new TH2F("vertexXRotTriggerPt","vertexXRotTriggerPt; p_{T}^{h} (GeV/c); x (fm)",200,0,200,100,-10,10);

  vertexXYRot4GeVTrigger->Sumw2();
  vertexXYRot5GeVTrigger->Sumw2();
  vertexXYRot20GeVTrigger->Sumw2();

  vertexXRotTriggerPt->Sumw2();

  vertexR->Sumw2();
  vertexRHardCore->Sumw2();
  vertexRNoHardCore->Sumw2();
  vertexRDijetCut->Sumw2();
  vertexXY->Sumw2();
  vertexXYHardCore->Sumw2();
  vertexXYNoHardCore->Sumw2();
  vertexXYRot->Sumw2();
  vertexXYRotHardCore->Sumw2();
  vertexXYRotNoHardCore->Sumw2();
  vertexXRotJetPt->Sumw2();
  vertexXRotJetPtHardCore->Sumw2();
  vertexJetPtR->Sumw2();
  vertexJetPtRHardCore->Sumw2();
  vertexHighestJetZR->Sumw2();
  vertexHighestJetZXRot->Sumw2();

  TH1D * jetPt = new TH1D("jetPt", "jetPt;p_{T}^{jet} (GeV/c)", 200, 0, 200);
  TH1D * recJetPt = new TH1D("recJetPt", "recJetPt;p_{T}^{jet} (GeV/c)", 225, -25, 200);
  TH1D * leadingJetPt = new TH1D("leadingJetPt","Leading Jet p_{T};p_{T}^{jet} (GeV/c)",200,0,200);
  TH1D * jetPtUnweighted = new TH1D("jetPtUnweighted", "jetPtUnweighted;p_{T}^{jet} (GeV/c)", 200, 0, 200);
  TH1D * particleEnergy = new TH1D("particleEnergy", "Total Energy of Outgoing Particles; E_{T} (GeV)", 100, 0, 2000);
  //MHO
  TH2F * jetPhiPt = new TH2F("jetPhiPt", "jetPhiPt;#phi;p_{T}^{jet} (GeV/c)", 100, 0, 2.0*TMath::Pi(), 100, 0, 100);
  TH2F * jetEtaPhi = new TH2F("jetEtaPhi", "jetEtaPhi;#Delta#eta;#Delta#phi", 100, -5, 5, 100, 0, 2.0*TMath::Pi());
  TH2F * leadingJetEtaPhi = new TH2F("leadingJetEtaPhi", "leadingJetEtaPhi;#Delta#eta;#Delta#phi", 100, -5, 5, 100, 0, 2.0*TMath::Pi());
  // this histogram basically only exists to keep track of the number of events
  TH1I * hWeight = new TH1I("hWeight","hWeight",nHWeightBins, (double *) &hWeightBins[0]);
//  TH1I * hWeight = new TH1I("hWeight","hWeight",110,0,1.1);
  TH1D * gammaEt = new TH1D("gammaEt","gammaEt",200,0,200);
  TH1D * trackPt = new TH1D("trackPt","trackPt",200,0,200);
  TH1D * chargedTrackPt = new TH1D("chargedTrackPt","chargedTrackPt",200,0,200);
  TH1D * leadingTrackPt = new TH1D("leadingTrackPt","Leading Track Pt",200,0,200);
    

  // Scattering Center Research
  TH1F * scatCentPt = new TH1F("scatCentPt","scatCentPt",100,0,8);
  TH1F * scatCentE = new TH1F("scatCentE","scatCentE",100,0,20);
  TH1F * scatCentEtot = new TH1F("scatCentEtot","scatCentEtot",100,0,1000);

  printf("Initializing some hard scatter histograms\n");

  // Hard Scatter info
  double max_pt_hs_spectrum = 500;
  TH1F * qqbar_pt = new TH1F("qqbar_pt","q#bar{q};p_{t} GeV/c",100,0,max_pt_hs_spectrum);
  TH1F * qq_pt = new TH1F("qq_pt","qq;p_{t} GeV/c",100,0,max_pt_hs_spectrum);
  TH1F * gg_pt = new TH1F("gg_pt","gg;p_{t} GeV/c",100,0,max_pt_hs_spectrum);
  TH1F * qg_pt = new TH1F("qg_pt","qg;p_{t} GeV/c",100,0,max_pt_hs_spectrum);
  TH2F * hsPtLeadingJetPt = new TH2F("hsPtLeadingJetPt","Leading Jet vs Hard Scatter;p_{t}^{hard} GeV/c;p_{t}^{jet} GeV/c",100,0,100,100,0,100);
  TH2F * hsPtSubLeadingJetPt = new TH2F("hsPtSubLeadingJetPt","Sub-Leading Jet vs Hard Scatter;p_{t}^{hard} GeV/c;p_{t}^{jet} GeV/c",100,0,100,100,0,100);
  TH2F * hsDEtaDPhiLeadingJet = new TH2F("hsDEtaDPhiLeadingJet","Hard Scatter - Leading Jet d#eta d#phi;d#eta;d#phi",100,-1,1,100,-1,1);
  TH2F * hsDEtaDPhiSubLeadingJet = new TH2F("hsDEtaDPhiSubLeadingJet","Hard Scatter - Sub-leading Jet d#eta d#phi;d#eta;d#phi",100,-1,1,100,-1,1);

  // Testing hard core extraction from JEWEL:
  TH2F * hs1hs2Pt = new TH2F("hs1hs2Pt","p_{T}^{hs1} vs p_{T}^{hs2};p_{T}^{hs1} (GeV/c); p_{T}^{hs2} (GeV/c)",150,0,150,150,0,150);  
  TH1F * hs1hs2dPhi = new TH1F("hs1hs2dPhi","#Delta#phi HS1 - HS2;#Delta#phi",100, -TMath::Pi()/2, 3.0*TMath::Pi()/2);
  hs1hs2Pt->Sumw2();
  hs1hs2dPhi->Sumw2();

  // hard core versions?

  qqbar_pt->Sumw2();
  qq_pt->Sumw2();
  gg_pt->Sumw2();
  qg_pt->Sumw2();

  hsPtLeadingJetPt->Sumw2();
  hsPtSubLeadingJetPt->Sumw2();
  hsDEtaDPhiLeadingJet->Sumw2();
  hsDEtaDPhiSubLeadingJet->Sumw2();

  leadingJetPt->Sumw2();
  jetPt->Sumw2();
  recJetPt->Sumw2();
  gammaEt->Sumw2();
  trackPt->Sumw2();
  chargedTrackPt->Sumw2();
  leadingTrackPt->Sumw2();
  jetPhiPt->Sumw2();

  // For Event Plane analysis
  TH1D * leadingJetPtEP[nEPBins];
  for (int k = 0; k < nEPBins; k++) {
    combinedClass = Form("leadingJetPtEP_%d", k);
    leadingJetPtEP[k] = new TH1D(combinedClass,Form("Leading Jet p_{T} (EP Bin %d);p_{T}^{jet} (GeV/c)",k),200,0,200);
    leadingJetPtEP[k]->Sumw2();
  }
    


  scatCentPt->Sumw2();
  scatCentE->Sumw2();
  scatCentEtot->Sumw2();

  // Mass analysis
  TH2F *hJetPtMass = new TH2F("hJetPtMass","hJetPtMass",10,0,200,70,0,70);
  TH2F *hSDJetPtMass = new TH2F("hSDJetPtMass","hSDJetPtMass",10,0,200,70,0,70);
  hJetPtMass->Sumw2();
  hSDJetPtMass->Sumw2();


    // For softdrop analysis
    /*TH2F *hJetPtZ  = new TH2F("hJetPtZ",  "hJetPtZ",  nHighPtBins, (double *) &highPtBins[0],50, 0.0, 1.0);
    hJetPtZ->Sumw2();
    TH2F * hJetPtSD2Mass = new TH2F("hJetPtSD2Mass","hJetPtSD2Mass;p_{T}^{jet} (GeV/c);m_{s.d.} (GeV/c^2)", nHighPtBins, (double *) &highPtBins[0],250,0,500);
    hJetPtSD2Mass->Sumw2();

    TH2F * hJetPtDR = new TH2F("hJetPtDR","hJetPtDR;p_{T}^{jet} (GeV/C);#Delta R", nHighPtBins, (double *) &highPtBins[0],50,0,0.4);
    hJetPtDR->Sumw2();
    TH2F *hJetPtZ_DRCut  = new TH2F("hJetPtZ_DRCut",  "hJetPtZ_DRCut;p_{T}^{jet};z",  nHighPtBins, (double *) &highPtBins[0],50, 0.0, 1.0);
    hJetPtZ_DRCut->Sumw2();
*/
  //->Fill(sd_jet.pt(), sd_deltaR, weight);
  //        if (sd_deltaR > 0.1) {
  //          hJetPtZ_DRCut->Fill(sd_jet.pt(), z_g, weight);


  TH1D * hSDJetPt = new TH1D("hSDJetPt","Soft Drop p_{T}^{jet};p_{T}^{jet} (GeV/c)",200,0,200);
  hSDJetPt->Sumw2();
  TH1D * hSDJetPtBkgSub = new TH1D("hSDJetPtBkgSub","Soft Drop p_{T}^{jet};p_{T}^{jet} (GeV/c)",200,0,200);
  hSDJetPtBkgSub->Sumw2();


  // Jet-hadron is done more differentially
  //  enum { nJetPtBins = 3, nAssocParticleBins = 11 };

  

  // what was this for?
  int nParticleBins = nJetAssocPtBins; // either particlePtBins, pi0ParticlePtBins, or zTbins
  
  std::vector<double> assocParticlePtBins = {};
  if (iTriggerType == 0) {
    nAssocParticleBins = nJetAssocPtBins;
    assocParticlePtBins = particlePtBins;
    // trigger_eta_cut = jet_eta_cut by default
  } else { //Pi0
    trigger_eta_cut = pi0_eta_cut;
    if (bUseZtBins) {
      nAssocParticleBins = nZtBins;
      assocParticlePtBins = zTBins;
    }
    else {
      nAssocParticleBins = nPi0ParticlePtBins; 
      assocParticlePtBins = pi0ParticlePtBins;
    }
  }

  printf("Will use associated bins: ");
  for (int i = 0; i < nAssocParticleBins; i++) {
    printf("%.2f, ",assocParticlePtBins.at(i));
  }
  printf("%.2f\n",assocParticlePtBins.back());

  TH2F * jetHadron[nTriggerPtBins][nAssocParticleBins];
//  TH2F * jetHadron[nJetPtBins][nAssocParticleBins];
  TH2F * jetHadronScatCent[nTriggerPtBins][nAssocParticleBins];
  TH2F * jetHadronOnlyJet[nTriggerPtBins][nAssocParticleBins];
  TH2F * jetHadronExcludeJet[nTriggerPtBins][nAssocParticleBins];

  // can I really have non-const number of bins in 2nd axis?
  TH2F * jetHadronEP[nEPBins][nTriggerPtBins][nAssocParticleBins];
  //TH2F * jetHadronEP[nEPBins][nJetPtBins][nAssocParticleBins];
  //
  // for doing correlations with scattering centers
  //TH2F * jetHadronScatCentEP[nEPBins][nTriggerPtBins][nAssocParticleBins];

    // Binning By Event Plane relative angle
  //  const int nEventPlaneBins = 3;
  //  TH2F * jetHadronEventPlane[nEventPlaneBins][nJetPtBins][nAssocParticleBins];

  TH1F * ptBinHadronPt[nTriggerPtBins];
  TH1F * ptBinHadronPtEP[nEPBins][nTriggerPtBins];
  // for Swift background
  TH1F *swiftJetYield[nTriggerPtBins];
  TH1F *swiftParticleYield[nTriggerPtBins][nAssocParticleBins];

  int nBinsDEta = 56;
  // In case of larger eta cut, maintain same bin width as Joel
  if (abs(trigger_eta_cut + eta_cut - 3.25) < 0.001) nBinsDEta = 140;

  int nBinsDPhi = 144;

  for (unsigned int i = 0; i < nTriggerPtBins; i++)
  {
    combinedClass = Form("ptBinHadronPt_" + jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
    //combinedClass = Form("ptBinHadronPt_" + jetClass,jetPtBins.at(i),jetPtBins.at(i+1));
    ptBinHadronPt[i] = new TH1F(combinedClass,combinedClass, 200,0,50);
    ptBinHadronPt[i]->Sumw2();

    //FIXME
    combinedClass = Form("swiftJetYield_" + jetClass,fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
    swiftJetYield[i] = new TH1F(combinedClass,combinedClass,2*nBinsDEta,-trigger_eta_cut,trigger_eta_cut); //x axis is eta
    swiftJetYield[i]->Sumw2();  // useful for fitting.


    for (unsigned int j = 0; j < nAssocParticleBins; j++)
    {
      combinedClass = Form(triggerType+"Hadron_" + jetClass + "_" + particleClass, fTriggerPtBins.at(i), fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
      jetHadron[i][j] = new TH2F(combinedClass, combinedClass, nBinsDEta, -eta_cut - trigger_eta_cut, eta_cut + trigger_eta_cut, nBinsDPhi, -TMath::Pi()/2, 3.0*TMath::Pi()/2);
      combinedClass = Form(triggerType+"HadronScatCent_" + jetClass + "_" + particleClass, fTriggerPtBins.at(i), fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
      jetHadronScatCent[i][j] = new TH2F(combinedClass, combinedClass, nBinsDEta, -eta_cut - trigger_eta_cut, eta_cut + trigger_eta_cut, nBinsDPhi, -TMath::Pi()/2, 3.0*TMath::Pi()/2);
      combinedClass = Form(triggerType+"HadronOnlyJet_" + jetClass + "_" + particleClass, fTriggerPtBins.at(i), fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
      jetHadronOnlyJet[i][j] = new TH2F(combinedClass, combinedClass, nBinsDEta, -eta_cut - trigger_eta_cut, eta_cut + trigger_eta_cut, nBinsDPhi, -TMath::Pi()/2, 3.0*TMath::Pi()/2);
      combinedClass = Form(triggerType+"HadronExcludeJet_" + jetClass + "_" + particleClass, fTriggerPtBins.at(i), fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
      jetHadronExcludeJet[i][j] = new TH2F(combinedClass, combinedClass, nBinsDEta, -eta_cut - trigger_eta_cut, eta_cut + trigger_eta_cut, nBinsDPhi, -TMath::Pi()/2, 3.0*TMath::Pi()/2);

      jetHadron[i][j]->Sumw2();
      jetHadronScatCent[i][j]->Sumw2();
      jetHadronOnlyJet[i][j]->Sumw2();
      jetHadronExcludeJet[i][j]->Sumw2();

      combinedClass = Form("swiftParticleYield_" + jetClass + "_" + particleClass, fTriggerPtBins.at(i),fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
      swiftParticleYield[i][j] = new TH1F(combinedClass,combinedClass, 2*nBinsDEta ,-eta_cut,eta_cut);
      swiftParticleYield[i][j]->Sumw2();

    }
  }
  // Initializing histos by event plane binning
  for (int k = 0; k < nEPBins; k++) {
    for (unsigned int i = 0; i < nTriggerPtBins; i++)
    {

      combinedClass = Form("ptBinHadronPt_EP_%d_" + jetClass, k, fTriggerPtBins.at(i),fTriggerPtBins.at(i+1));
      ptBinHadronPtEP[k][i] = new TH1F(combinedClass,combinedClass, 200,0,50);
      ptBinHadronPtEP[k][i]->Sumw2();

      for (unsigned int j = 0; j < nAssocParticleBins; j++)
      {
        combinedClass = Form(triggerType+"Hadron_EP_%d_" + jetClass + "_" + particleClass, k, fTriggerPtBins.at(i), fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
        jetHadronEP[k][i][j] = new TH2F(combinedClass, combinedClass, nBinsDEta, -eta_cut - trigger_eta_cut, eta_cut + trigger_eta_cut, nBinsDPhi, -TMath::Pi()/2, 3.0*TMath::Pi()/2);
        jetHadronEP[k][i][j]->Sumw2();

        //combinedClass = Form(triggerType+"jetHadronScatCent_EP_%d_" + jetClass + "_" + particleClass, k, fTriggerPtBins.at(i), fTriggerPtBins.at(i+1), assocParticlePtBins.at(j), assocParticlePtBins.at(j+1));
        //jetHadronScatCentEP[k][i][j] = new TH2F(combinedClass, combinedClass, nBinsDEta, -eta_cut - trigger_eta_cut, eta_cut + trigger_eta_cut, nBinsDPhi, -TMath::Pi()/2, 3.0*TMath::Pi()/2);
        //jetHadronScatCentEP[k][i][j]->Sumw2();
      }
    }
  }




  int nGoodEvents = 0;
  
  int nevents = tree->GetEntries();
  // Event loop
  for (Int_t eventNumber = 0; eventNumber < nevents; eventNumber++)
  {
    if (iOddEven > 0) { // Skip Some events
      if ((iOddEven == 1) && (eventNumber % 2 == 1)) continue;
      if ((iOddEven == 2) && (eventNumber % 2 == 0)) continue;
      if (iOddEven > 2) {
        int nMod4 = eventNumber % 4;
        if ((nMod4 == 1) && (iOddEven == 3)) continue;
        if ((nMod4 == 2) && (iOddEven == 4)) continue;
        if ((nMod4 == 3) && (iOddEven == 5)) continue;
        if ((nMod4 == 0) && (iOddEven == 6)) continue;
      }
    }


    tree->GetEntry(eventNumber);
  
    TString hashString = "";

    total_weight += weight * weightingOn;
    // Fill weight histogram
    hWeight->Fill(weight);

    vertexXY->Fill(vertex_x,vertex_y,weight);
    vertexR->Fill(TMath::Sqrt(pow(vertex_x,2) + pow(vertex_y,2)),weight);


    // Fill particles ito array
    for (unsigned int particleNumber = 0; particleNumber < px->size(); particleNumber++)
    {
      // For TLorentzVector
      /*      TLorentzVector particle;
            particle.SetPxPyPzE(px->at(particleNumber), py->at(particleNumber),
            pz->at(particleNumber), energy->at(particleNumber) );
            particles.push_back(particle);*/
      // For Fastjet
      fastjet::PseudoJet jet =  fastjet::PseudoJet( px->at(particleNumber), py->at(particleNumber),
          pz->at(particleNumber), energy->at(particleNumber) );
      hashString += TString::Format(" %e ",energy->at(particleNumber));

      jet.set_user_index(particleID->at(particleNumber)); // storing pid in user index

      if (iTriggerType == 1) {
        if (particleID->at(particleNumber) == 111) {
					if (jet.pt() >= pi0_pt_min && abs(jet.eta()) < pi0_eta_cut) { 
						jets.push_back(jet); // treat the pi0 as a jet in this code
					}
        }
      }

      if (statusOn && status->at(particleNumber) == 3) {
        scatCents.push_back(jet);
        continue; //so that we don't count this in the other particles lists
      }

      if (statusOn && status->at(particleNumber) == 6) { // lucky number
  //        jet.set_user_index(particleID->at(particleNumber)); // storing pid in user index
        hardScatters.push_back(jet);
        continue; // don't use later
      }

      if (isNeutrino(particleID->at(particleNumber))) {
        continue;  // ALICE,STAR don't detect neutrinos.
      }

      // ALICE doesn't detect these
      if (isNeutralLongHadron(particleID->at(particleNumber))) continue;

      if (isGamma(particleID->at(particleNumber))) {
        gammaEt->Fill(jet.Et(),weight);
      }
      // setting a unique id code for each particle
      // to be used later for ID
      jet.set_user_index(particleNumber);       
      // Cut to minimum momentum for reconstruction
      if (jet.pt() >= jet_constituent_cut) { 
        constParticles.push_back(jet);
      } 
      if ((!requireChargedAssoc) || (isChargedParticle(particleID->at(particleNumber)))) {
        particles.push_back(jet);
      }
      //particles.push_back( fastjet::PseudoJet( px->at(particleNumber), py->at(particleNumber),
      //      pz->at(particleNumber), energy->at(particleNumber) ));
      // printf("%d\t%f\t%f\t%f\t%f\n",particleID->at(particleNumber),px->at(particleNumber),py->at(particleNumber),pz->at(particleNumber),energy->at(particleNumber));
    }
		
    // don't add to consituent particles, for convenience
		//int nParticlesBefore = particles.size();
		if (nToyParticles > 0) {
			AddToyModelParticles(nToyParticles,rand,particles);
		}
		//int nParticlesAfter = particles.size();
		//printf("Toy increased n particles from %d to %d\n",nParticlesBefore,nParticlesAfter);
    // Cut particles using erase remove idiom to avoid worrying about removing particles incorrectly
    particles.erase( std::remove_if(particles.begin(), particles.end(), rejectParticle), particles.end() );
    scatCents.erase( std::remove_if(scatCents.begin(),scatCents.end(), rejectParticle), scatCents.end() );

    constParticles.erase( std::remove_if(constParticles.begin(), constParticles.end(), rejectParticle), constParticles.end() );

    fastjet::AreaDefinition areaDef(fastjet::active_area,fastjet::GhostedAreaSpec(ghost_maxrap,1,ghost_area));
    fastjet::ClusterSequenceArea * cs = 0;
    if (iTriggerType == 0) {
      // Refactor from here so that it can be used by everyone
      // Define jet algorithm 
      //fastjet::AreaDefinition areaDef(fastjet::active_area,fastjet::GhostedAreaSpec(ghost_maxrap,1,ghost_area));

      // Perform jet finding
      //  fastjet::ClusterSequenceArea cs(particles, jetDef, areaDef);
      cs = new fastjet::ClusterSequenceArea(constParticles, jetDef, areaDef);
      //fastjet::ClusterSequenceArea cs(constParticles, jetDef, areaDef);
      //    fastjet::ClusterSequence cs(particles, jetDef);

      // Where the jets are extracted from FastJet
      jets = fastjet::sorted_by_pt(cs->inclusive_jets());
      //jets = fastjet::sorted_by_pt(cs.inclusive_jets());

      // removing jets that don't pass eta_cut, global hard core cut (if set)
      jets.erase( std::remove_if(jets.begin(), jets.end(), rejectJet), jets.end() );
    }
    else if (iTriggerType == 1) { // Pi0s
      // FIXME
      jets = fastjet::sorted_by_pt(jets);
    }

    double jet_phi, jet_eta, rot_x, rot_y, vertex_r;

    // Begin the actual analysis 

    // Examine QCD Hard Scatter
    fastjet::PseudoJet hs1,hs2;
    int pid_sum;  // declaring here so I can use it later
    if (hardScatters.size()==2) {
      hs1 = hardScatters.at(0);
      hs2 = hardScatters.at(1);
      nGoodEvents++;
      pid_sum = TMath::Abs(hs1.user_index()) + TMath::Abs(hs2.user_index());

      // check eta cut for each hard scatter
      bool hs1pass = TMath::Abs(hs1.eta()) < eta_cut;
      bool hs2pass = TMath::Abs(hs2.eta()) < eta_cut;
      if (pid_sum < 13) { // 2 quarks/antiquarks
        if (hs1.user_index()*hs2.user_index() < 0 ) { //q qbar
          if (hs1pass) qqbar_pt->Fill(hs1.pt(),weight);
          if (hs2pass) qqbar_pt->Fill(hs2.pt(),weight);
        } else { // qq or qbarqbar
          if (hs1pass) qq_pt->Fill(hs1.pt(),weight);
          if (hs2pass) qq_pt->Fill(hs2.pt(),weight);
        }

      } else if (pid_sum == 42) { // 2 gluons
        if (hs1pass) gg_pt->Fill(hs1.pt(),weight);
        if (hs2pass) gg_pt->Fill(hs2.pt(),weight);
      } else if (pid_sum > 21 && pid_sum < 28) { // 1 gluon, 1 quark
        if (hs1pass) qg_pt->Fill(hs1.pt(),weight);
        if (hs2pass) qg_pt->Fill(hs2.pt(),weight);
      }

      hs1hs2Pt->Fill(hs1.pt(),hs2.pt(),weight);
      hs1hs2dPhi->Fill(calculateDeltaPhi(hs1.phi(),hs2.phi()),weight);
    }

    // Jet-h correlations

    // Alternating hs1,hs2
    fastjet::PseudoJet *hs;
    if (eventNumber % 2 == 0) {hs = &hs1; 
    }
    else {hs = &hs2;   
    }
    
    //leading jet
    if (jets.size()) {
    //  if (jets.size() && (TMath::Abs(hs->eta()) < eta_cut - R)) {
    
      fastjet::PseudoJet * leadingJet = &jets.at(0);
      // Note: redefine leadingJet here to mess with the trigger
      // FIXME remember to remove this
      //fastjet::PseudoJet * leadingJet = hs;

      jet_phi = leadingJet->phi();
      jet_eta = leadingJet->eta();
      rot_x = vertex_x * TMath::Cos(PI - jet_phi) - vertex_y * TMath::Sin(PI - jet_phi);
      rot_y = vertex_x * TMath::Sin(PI - jet_phi) + vertex_y * TMath::Cos(PI - jet_phi);
      vertex_r = TMath::Sqrt(pow(vertex_x,2) + pow(vertex_y,2));

      // FIXME this is not the leading pi0 yet. just the first one
      double leading_pt = leadingJet->pt();
      double true_leading_pt = leading_pt;

      // doing background subtraction, if requested
      if ((iTriggerType == 0) && doJetBkgSub) {
        double rho = 0.0;
        rho = getBackground(leadingJet, &particles);
        double jet_area = leadingJet->area();
        leading_pt = leading_pt - rho * jet_area;
      }

      int iEPBin = 0;
      if (nEPBins > 1) {
        iEPBin = GetEventPlaneBin(leadingJet->phi());
//          leadingJetPtEP[iEPBin]->Fill(leading_pt,weight);
      }

      std::vector<double>::iterator jetPtBinIter;
    
      // For jet-h, same as jetPtBinIter; for pi0, pi0
  //    std::vector<double>::iterator triggerPtBinIter;

      jetPtBinIter = std::lower_bound(fTriggerPtBins.begin(), fTriggerPtBins.end(), leading_pt);
     
  //    if (iTriggerType == 0) {
//        triggerPtBinIter
  //    }
        //jetPtBinIter = std::lower_bound(jetPtBins.begin(), jetPtBins.end(), leading_pt);
      //// jetPtBinIter = std::lower_bound(jetPtBins.begin(), jetPtBins.end(), jets.at(0).pt());
      //} else if (iTriggerType == 1) {
      // // FIXME 
      //  jetPtBinIter = std::lower_bound(pi0PtBins.begin(), pi0PtBins.end(), leading_pt);
      //}
     
      // The -1 is because lower_bound counts below the lowest value as defining the 0 index.
      // We are interested in everything above that index.
      Int_t jetPtBin = jetPtBinIter - fTriggerPtBins.begin() - 1;
      //Int_t jetPtBin = jetPtBinIter - jetPtBins.begin() - 1;


      // Compare to QCD Hard Scatter, if available
      if (hardScatters.size()==2) {
        double hs1Phi = hs1.phi();
        double hs2Phi = hs2.phi();      
        double dPhi2hs1 = calculateDeltaPhi(jet_phi,hs1Phi);
        double dPhi2hs2 = calculateDeltaPhi(jet_phi,hs2Phi);
        double dEta2hs1 = jet_eta - hs1.eta();
        double dEta2hs2 = jet_eta - hs2.eta();

  //        hs1hs2Pt->Fill(hs1.pt(),hs2.pt(),weight);
  //        hs1hs2dPhi->Fill(calculateDeltaPhi(hs1Phi,hs2Phi),weight);


        std::vector<double>::iterator subleadingJetPtBinIter;
        std::vector<double>::iterator partonPtBinIter;

          // safe to assume leading jet is nearly within PI of one hard scatter
  //        if (TMath::Abs(dPhi2hs1) <= TMath::Abs(dPhi2hs2) && TMath::Abs(hs1.eta()) < eta_cut) { // leading jet is hs1.
        if (TMath::Abs(dPhi2hs1) <= TMath::Abs(dPhi2hs2) && TMath::Abs(dEta2hs1) < 1 && TMath::Abs(dPhi2hs1) < 1) { // leading jet is hs1.
            // and hs1 meets cut


          hsPtLeadingJetPt->Fill(hs1.pt(),leading_pt,weight);
          hsDEtaDPhiLeadingJet->Fill(-dEta2hs1,-dPhi2hs1,weight);         
            // Filling type histograms for hs1
          ptLossVertexRLeadingJet->Fill(vertex_r,hs1.pt() - true_leading_pt,weight);

          if ((jetPtBin >= 0) && (jetPtBin < nTriggerPtBins)) {
            PartonPtByJetBin[iEPBin][jetPtBin]->Fill(hs1.pt(),weight);

            ptLossLeadingJetPtBinEP[iEPBin][jetPtBin]->Fill(hs1.pt() - true_leading_pt,weight);
            energyLossLeadingJetPtBinEP[iEPBin][jetPtBin]->Fill(hs1.e() - leadingJet->e(),weight);

            partonPtBinIter = std::lower_bound(jetPtBins.begin(), jetPtBins.end(), hs1.pt());
            Int_t partonPtBin = partonPtBinIter - jetPtBins.begin() - 1;
            if ((partonPtBin >= 0) && (partonPtBin < nJetPtBins)) {
              ptLossLeadingPartonPtBinEP[iEPBin][partonPtBin]->Fill(hs1.pt() - true_leading_pt,weight);
              energyLossLeadingPartonPtBinEP[iEPBin][partonPtBin]->Fill(hs1.e() - leadingJet->e(),weight);
            }

            if (pid_sum < 13) { // 2 quarks/antiquarks
              if (hs1.user_index()*hs2.user_index() < 0 ) { //q qbar
                qqbar_PartonPtByJetBin[iEPBin][jetPtBin]->Fill(hs1.pt(),weight);
              } else { // qq or qbarqbar
                qq_PartonPtByJetBin[iEPBin][jetPtBin]->Fill(hs1.pt(),weight);
              }

            } else if (pid_sum == 42) { // 2 gluons
              gg_PartonPtByJetBin[iEPBin][jetPtBin]->Fill(hs1.pt(),weight);   
            } else if (pid_sum > 21 && pid_sum < 28) { // 1 gluon, 1 quark
              if (hs1.user_index() == 21) { //hs 1 is the gluon
                gq_PartonPtByJetBin[iEPBin][jetPtBin]->Fill(hs1.pt(),weight);
              } else if (hs2.user_index() == 21) { // hs 1 is  the quark
                qg_PartonPtByJetBin[iEPBin][jetPtBin]->Fill(hs1.pt(),weight);
              }
            }
          }

          // check for first subleading jet within PI of hs2 (jets ordered by pt)
          // assume this is hs2. 
          if (TMath::Abs(hs2.eta()) < eta_cut) {
            for (unsigned int i = 1; i < jets.size(); i++)
            {
              if (TMath::Abs(calculateDeltaPhi(hs2Phi,jets.at(i).phi()) < PI/2)) {
                hsPtSubLeadingJetPt->Fill(hs2.pt(),jets.at(i).pt(),weight);
                hsDEtaDPhiSubLeadingJet->Fill(hs2.eta() - jets.at(i).eta(),calculateDeltaPhi(hs2.phi(),jets.at(i).phi()),weight);         
                ptLossVertexRSubleadingJet->Fill(vertex_r,hs2.pt() - jets.at(i).pt(),weight);

                subleadingJetPtBinIter = std::lower_bound(jetPtBins.begin(), jetPtBins.end(), jets.at(i).pt());
                  Int_t subleadingJetPtBin = subleadingJetPtBinIter - jetPtBins.begin() - 1;
                  if ((subleadingJetPtBin >= 0) && (subleadingJetPtBin < nJetPtBins)) {
                    ptLossSubleadingJetPtBinEP[iEPBin][subleadingJetPtBin]->Fill(hs2.pt() - jets.at(i).pt(),weight);
                    energyLossSubleadingJetPtBinEP[iEPBin][subleadingJetPtBin]->Fill(hs2.e() - jets.at(i).e(),weight);
                  }
    
                  partonPtBinIter = std::lower_bound(jetPtBins.begin(), jetPtBins.end(), hs2.pt());
                  Int_t partonPtBin = partonPtBinIter - jetPtBins.begin() - 1;
                  if ((partonPtBin >= 0) && (partonPtBin < nJetPtBins)) {
                    ptLossSubleadingPartonPtBinEP[iEPBin][partonPtBin]->Fill(hs2.pt() - jets.at(i).pt(),weight);
                    energyLossSubleadingPartonPtBinEP[iEPBin][partonPtBin]->Fill(hs2.e() - jets.at(i).e(),weight);
                  }

                  break;
                }
              }
            }
        //  } else if (TMath::Abs(hs2.eta()) < eta_cut) { // leading jet is hs2, hs2 meets cut
          } else if (TMath::Abs(dPhi2hs1) > TMath::Abs(dPhi2hs2) && TMath::Abs(dEta2hs2) < 1 && TMath::Abs(dPhi2hs2) < 1) { // leading jet is hs2, hs2 meets cut


            hsPtLeadingJetPt->Fill(hs2.pt(),leading_pt,weight);
            hsDEtaDPhiLeadingJet->Fill(-dEta2hs2,-dPhi2hs2,weight);         
            ptLossVertexRLeadingJet->Fill(vertex_r,hs2.pt() - true_leading_pt,weight);

            // Filling type histograms for hs1
            if ((jetPtBin >= 0) && (jetPtBin < nTriggerPtBins)) {
              PartonPtByJetBin[iEPBin][jetPtBin]->Fill(hs2.pt(),weight);

              ptLossLeadingJetPtBinEP[iEPBin][jetPtBin]->Fill(hs2.pt() - true_leading_pt,weight);
              energyLossLeadingJetPtBinEP[iEPBin][jetPtBin]->Fill(hs2.e() - leadingJet->e(),weight);

              partonPtBinIter = std::lower_bound(jetPtBins.begin(), jetPtBins.end(), hs2.pt());
              Int_t partonPtBin = partonPtBinIter - jetPtBins.begin() - 1;
              if ((partonPtBin >= 0) && (partonPtBin < nJetPtBins)) {
                ptLossLeadingPartonPtBinEP[iEPBin][partonPtBin]->Fill(hs2.pt() - true_leading_pt,weight);
                energyLossLeadingPartonPtBinEP[iEPBin][partonPtBin]->Fill(hs2.e() - leadingJet->e(),weight);
              }

              if (pid_sum < 13) { // 2 quarks/antiquarks
                if (hs1.user_index()*hs2.user_index() < 0 ) { //q qbar
                  qqbar_PartonPtByJetBin[iEPBin][jetPtBin]->Fill(hs2.pt(),weight);
                } else { // qq or qbarqbar
                  qq_PartonPtByJetBin[iEPBin][jetPtBin]->Fill(hs2.pt(),weight);
                }

              } else if (pid_sum == 42) { // 2 gluons
                gg_PartonPtByJetBin[iEPBin][jetPtBin]->Fill(hs2.pt(),weight);   
              } else if (pid_sum > 21 && pid_sum < 28) { // 1 gluon, 1 quark
                if (hs2.user_index() == 21) { //hs 2 is the gluon
                  gq_PartonPtByJetBin[iEPBin][jetPtBin]->Fill(hs2.pt(),weight);
                } else if (hs1.user_index() == 21) { // hs 2 is  the quark
                  qg_PartonPtByJetBin[iEPBin][jetPtBin]->Fill(hs2.pt(),weight);
                }
              }
            }

            // check for first subleading jet within PI of hs1
            if (TMath::Abs(hs1.eta()) < eta_cut) { // check that hs1 is in eta range
              for (unsigned int i = 1; i < jets.size(); i++)
              {
                if (TMath::Abs(calculateDeltaPhi(hs2Phi,jets.at(i).phi()) < PI/2)) {
                  hsPtSubLeadingJetPt->Fill(hs1.pt(),jets.at(i).pt(),weight);
                  hsDEtaDPhiSubLeadingJet->Fill(hs1.eta() - jets.at(i).eta(),calculateDeltaPhi(hs1.phi(),jets.at(i).phi()),weight);         
                  ptLossVertexRSubleadingJet->Fill(vertex_r,hs1.pt() - jets.at(i).pt(),weight);

                  subleadingJetPtBinIter = std::lower_bound(jetPtBins.begin(), jetPtBins.end(), jets.at(i).pt());
                  Int_t subleadingJetPtBin = subleadingJetPtBinIter - jetPtBins.begin() - 1;
                  if ((subleadingJetPtBin >= 0) && (subleadingJetPtBin < nJetPtBins)) {
                    ptLossSubleadingJetPtBinEP[iEPBin][subleadingJetPtBin]->Fill(hs1.pt() - jets.at(i).pt(),weight);
                    energyLossSubleadingJetPtBinEP[iEPBin][subleadingJetPtBin]->Fill(hs1.e() - jets.at(i).e(),weight);
                  }

                  partonPtBinIter = std::lower_bound(jetPtBins.begin(), jetPtBins.end(), hs1.pt());
                  Int_t partonPtBin = partonPtBinIter - jetPtBins.begin() - 1;
                  if ((partonPtBin >= 0) && (partonPtBin < nJetPtBins)) {
                    ptLossSubleadingPartonPtBinEP[iEPBin][partonPtBin]->Fill(hs1.pt() - jets.at(i).pt(),weight);
                    energyLossSubleadingPartonPtBinEP[iEPBin][partonPtBin]->Fill(hs1.e() - jets.at(i).e(),weight);
                  }

                  break;
                }
              }
            }
          }
        }
        // End of Hard Scatter section

        // Phi Acceptance Simulation
        // Simulates the limited acceptance of the ALICE EMCAL by looking at the leading jet in each of three 120 degree sections.
        // The starting phi of these sections must be randomized, to avoid correlation with the event plane.

        int phiSteps = 1;
        float initPhi = 0;  // initial phi for the rotating phi acceptance cut

        if (phiAcceptanceSim) {
          phiSteps = 3;
          initPhi = rand->Rndm()*2*PI;

        }

        // place a loop here for the number of jets to analyze (1 for leading jet only, jets.size() for all)
        int nTriggerJets = jets.size();
        if (bLeadingJetOnly) nTriggerJets = 1;
        for (int iTrigger = 0; iTrigger < nTriggerJets; iTrigger++) { 

					fastjet::PseudoJet * trigger = &jets.at(iTrigger);
					double trigger_pt = trigger->pt();

					//iterate once or three times, depending on whether the phi acceptance simulation is on
					for (int z = 0; z < phiSteps; z++) {
						if (phiAcceptanceSim && bLeadingJetOnly) {
							float phiStep = fmod(initPhi + 2. * z * PI / 3, 2*PI);
							// If we do the phi acceptance sim, we re-determine the leading jet in the acceptance region
							leadingJet = 0;
							float phiAcceptanceRange = 107. * PI / 180.; // phi range of EMCal main
							for (int i = 0; i < jets.size(); i++) { //Looking for leading jet within acceptance
								if (abs(calculateDeltaPhi(jets.at(i).phi(),phiStep)) < phiAcceptanceRange / 2.) {
									leadingJet = &jets.at(i); // this works because jets are sorted by pt
									break;
								}
							}
							if (!leadingJet) { // No jet found within this phi acceptance
								continue;
							}
							trigger = leadingJet;
						}

						jet_phi = trigger->phi();
						jet_eta = trigger->eta();
						rot_x = vertex_x * TMath::Cos(PI - jet_phi) - vertex_y * TMath::Sin(PI - jet_phi);
						rot_y = vertex_x * TMath::Sin(PI - jet_phi) + vertex_y * TMath::Cos(PI - jet_phi);
						vertex_r = TMath::Sqrt(pow(vertex_x,2) + pow(vertex_y,2));

						trigger_pt = trigger->pt();
						double true_trigger_pt = trigger_pt;;
						// doing background subtraction, if requested
						if ((iTriggerType == 0) && doJetBkgSub) {
							double rho = 0.0;
							rho = getBackground(trigger, &particles);
							double jet_area = trigger->area();
							trigger_pt = trigger_pt - rho * jet_area;
						}

						int nScatCent = scatCents.size();
						if (nScatCent != 0) leadingJetPtNScatCent[0]->Fill(trigger_pt,nScatCent,weight);
						iEPBin = 0;
						if (nEPBins > 1) {
							iEPBin = GetEventPlaneBin(trigger->phi());
			
							if (nScatCent != 0) leadingJetPtNScatCent[iEPBin]->Fill(trigger_pt,nScatCent,weight);
						}

						jetPtBinIter = std::lower_bound(fTriggerPtBins.begin(), fTriggerPtBins.end(), trigger_pt);
						// The -1 is because lower_bound counts below the lowest value as defining the 0 index.
						// We are interested in everything above that index.
						jetPtBin = jetPtBinIter - fTriggerPtBins.begin() - 1;

						leadingJetEtaPhi->Fill(jet_eta,jet_phi,weight);

						vertexXRotJetPt->Fill(trigger_pt,rot_x,weight);
						if (!rejectJetHardCore(trigger, R, jet_hard_core_cut)) { 
							vertexXRotJetPtHardCore->Fill(trigger_pt,rot_x,weight);
						}

						// jet energy resolution 
						// check if res matrix exists.  if so, loop over jetptbins with an array of weights (0. - 1.)
						// function: get weight --> Return the array

						std::vector<double> responseWeights = {};

						int nTimes = 1;  
//          if (responseMatrix) {
//            responseWeights = FindJetReponseWeights(jets.at(0).pt(),responseMatrix);
//            nTimes = nJetPtBins;
//          }
						if (usingResponseMatrix) {
							responseWeights = FindJetReponseWeights(trigger->pt(),responseMatrix[iEPBin]);
							nTimes = nJetPtBins;
						}
						for (int m = 0; m < nTimes; m++) {
							double localWeight = 1.;
							if (usingResponseMatrix) {
								jetPtBin = m;
								localWeight = responseWeights[m];

								// hack to make sure the leading jet pt histos are useful for normalization
								trigger_pt  = (fTriggerPtBins.at(m+1) + fTriggerPtBins.at(m)) / 2.;
								//leading_pt = (jetPtBins.at(m+1) + jetPtBins.at(m)) / 2.; 

							}

							leadingJetPt->Fill(trigger_pt,weight*localWeight);
							leadingJetPtEP[0]->Fill(trigger_pt,weight*localWeight);
							leadingJetPtEP[iEPBin]->Fill(trigger_pt,weight*localWeight);

							if ((jetPtBin >= 0) && (jetPtBin < nTriggerPtBins))
							{
								vertexXYJetPt[jetPtBin]->Fill(vertex_x,vertex_y,weight);
								vertexXYRotJetPt[jetPtBin]->Fill(rot_x,rot_y,weight);
								if ((iTriggerType == 0) && !rejectJetHardCore(trigger, R, jet_hard_core_cut)) { 
									vertexXYRotJetPtHardCore[jetPtBin]->Fill(rot_x,rot_y,weight);
								}
	//              swiftJetYield[jetPtBin]->Fill(jets.at(0).eta(),weight);
								swiftJetYield[jetPtBin]->Fill(trigger->eta(),weight);
								// for jets we are using: assigning each particle to its jet using
								// the user_index feature
								// int   user_index () const 
								// void   set_user_index (const int index)
								/*      if (jets.at(i).has_constituents() == true)
											{
											std::vector <fastjet::PseudoJet> const_particles = jets.at(i).constituents();
											for (unsigned int j = 0; j < const_particles.size(); j++) 
											const_particles.at(j).set_user_index(i);
											}
											*/
								//Printf("jetPt: %f, Bin: %d", jets.at(i).pt(), jetPtBin);
								std::vector <double>::iterator particlePtBinIter;
								Int_t particlePtBin = -1;
								for (unsigned int j = 0; j < particles.size(); j++)
								{
									double particleWeight = 1; // Additional weight factor for associated particles
									// Could simulate flow with weighting here
									if(particles.at(j).user_index() == kToyParticleLabel) particleWeight = kToyParticleWeightScale;
									double particleObservable = particles.at(j).pt();
									if (bUseZtBins) {
										particleObservable = particleObservable / trigger_pt;
									}

									particlePtBinIter = std::lower_bound(assocParticlePtBins.begin(), assocParticlePtBins.end(), particleObservable);
									// The -1 is because lower_bound counts below the lowest value as defining the 0 index.
									// We are interested in everything above that index.
									particlePtBin = particlePtBinIter - assocParticlePtBins.begin() - 1;
									if ((particlePtBin >= 0) && (particlePtBin < nAssocParticleBins))
									{
										bool jetContainsParticle = false;
										if (iTriggerType == 0) {
											unsigned int particle_user_index = particles.at(j).user_index();
											if (trigger->has_constituents()) {
												std::vector <fastjet::PseudoJet> const_particles = trigger->constituents();
												for (unsigned int k = 0; k < const_particles.size(); k++) {
													if (particle_user_index == const_particles.at(k).user_index()) {
														jetContainsParticle = true;
														break;
													}
												}
											}
										}
										ptBinHadronPt[jetPtBin]->Fill(particles.at(j).pt(),weight*particleWeight);
										//jetHadron[jetPtBin][particlePtBin]->Fill(leadingJet->eta()-particles.at(j).eta(), calculateDeltaPhi(leadingJet->phi(), particles.at(j).phi()),weight);
										jetHadron[jetPtBin][particlePtBin]->Fill(trigger->eta()-particles.at(j).eta(), calculateDeltaPhi(trigger->phi(), particles.at(j).phi()),weight*localWeight*particleWeight);
										swiftParticleYield[jetPtBin][particlePtBin]->Fill(particles.at(j).eta(),weight*localWeight*particleWeight);

										jetHadronEP[0][jetPtBin][particlePtBin]->Fill(trigger->eta() - particles.at(j).eta(), calculateDeltaPhi(trigger->phi(),particles.at(j).phi()),weight*localWeight*particleWeight);
										ptBinHadronPtEP[0][jetPtBin]->Fill(particles.at(j).pt(),weight*localWeight*particleWeight);
										if (iEPBin) {
											jetHadronEP[iEPBin][jetPtBin][particlePtBin]->Fill(trigger->eta() - particles.at(j).eta(), calculateDeltaPhi(trigger->phi(),particles.at(j).phi()),weight*localWeight*particleWeight);
											ptBinHadronPtEP[iEPBin][jetPtBin]->Fill(particles.at(j).pt(),weight*localWeight*particleWeight);
										}
										if (iTriggerType == 0) {
											if (jetContainsParticle) {
												jetHadronOnlyJet[jetPtBin][particlePtBin]->Fill(trigger->eta()-particles.at(j).eta(), calculateDeltaPhi(trigger->phi(), particles.at(j).phi()),weight*localWeight*particleWeight);
											} else {
												jetHadronExcludeJet[jetPtBin][particlePtBin]->Fill(trigger->eta()-particles.at(j).eta(), calculateDeltaPhi(trigger->phi(), particles.at(j).phi()),weight*localWeight*particleWeight);
											}
										}
									}
								}
								// Filling ScatCent Only histograms
								for (unsigned int j = 0; j < scatCents.size(); j++)
								{ 
									particlePtBinIter = std::lower_bound(particlePtBins.begin(), particlePtBins.end(), scatCents.at(j).pt());
									// The -1 is because lower_bound counts below the lowest value as defining the 0 index.
									// We are interested in everything above that index.
									particlePtBin = particlePtBinIter - particlePtBins.begin() - 1;
									if ((particlePtBin >= 0) && (particlePtBin < nAssocParticleBins))
									{
										jetHadronScatCent[jetPtBin][particlePtBin]->Fill(trigger->eta()-scatCents.at(j).eta(), calculateDeltaPhi(trigger->phi(), scatCents.at(j).phi()),weight*localWeight);

										//jetHadronScatCentEP[0][jetPtBin][particlePtBin]->Fill(trigger->eta() - scatCents.at(j).eta(), calculateDeltaPhi(trigger->phi(),scatCents.at(j).phi()),weight*localWeight);
										//if (iEPBin) {
											//jetHadronScatCentEP[iEPBin][jetPtBin][particlePtBin]->Fill(trigger->eta() - scatCents.at(j).eta(), calculateDeltaPhi(trigger->phi(),scatCents.at(j).phi()),weight*localWeight);
										//}
									}
								}
							}
						}
          }
        }
      }

      // Study of geo bias for hadron-Lambda/k0 paper
      for (unsigned int i = 0; i < particles.size(); i++) {
        double part_pt = particles.at(i).pt();

        if ( part_pt > 4. ) {
          double hadron_phi = particles.at(i).phi();
          rot_x = vertex_x * TMath::Cos(PI - hadron_phi) - vertex_y * TMath::Sin(PI - hadron_phi);
          rot_y = vertex_x * TMath::Sin(PI - hadron_phi) + vertex_y * TMath::Cos(PI - hadron_phi);
          vertexXYRot4GeVTrigger->Fill(rot_x,rot_y,weight);
          if ( part_pt > 5. ) {
            vertexXYRot5GeVTrigger->Fill(rot_x,rot_y,weight);
            if ( part_pt > 20. ) {
              vertexXYRot20GeVTrigger->Fill(rot_x,rot_y,weight);
            }
          }
        }
        vertexXRotTriggerPt->Fill(part_pt,rot_x,weight);
      }

      // More jet-related measurements
      for (unsigned int i = 0; i < jets.size(); i++)
      {
  //      if (rejectJet(jets.at(i), R)) continue;
        
        double jet_pt = jets.at(i).pt();


        if ((iTriggerType == 0) && doJetBkgSub) {
          double rho = 0.0;
          rho = getBackground(&jets.at(i), &particles);
          double jet_area = jets.at(i) .area();
          jet_pt = jet_pt - rho * jet_area;
        }



        // Fill jet spectra
        jetPt->Fill(jet_pt,weight);
        jetPtUnweighted->Fill(jet_pt);
        jetEtaPhi->Fill(jets.at(i).eta(),jets.at(i).phi(),weight);
        jetPhiPt->Fill(jets.at(i).phi(),jet_pt,weight);
       // jetPhiPt->Fill(jets.at(i).phi(),jets.at(i).pt(),weight);

        // Mass analysis
        hJetPtMass->Fill(jet_pt,sqrt(jets.at(i).m2()),weight);
        
        //Soft Drop analysis  
        /*
        fastjet::PseudoJet sd_jet = sd(jets.at(i));
        assert(sd_jet != 0); //because soft drop is a groomer (not a tagger), it should always return a soft-dropped jet
   
        std::vector<fastjet::PseudoJet> piecesSoft = sd_jet.pieces();
        if ( (sd_jet != 0) && (piecesSoft.size()) >= 2 ) {
          fastjet::PseudoJet sd_jet_0 = sd(piecesSoft[0]);
          fastjet::PseudoJet sd_jet_1 = sd(piecesSoft[1]);
  //        float pt1 = sd_jet_0.pt();
  //        float pt2 = sd_jet_1.pt();
          float pt1 = piecesSoft[0].pt();
          float pt2 = piecesSoft[1].pt();
          float z_g = 0.0;
          if (pt1 <= pt2) z_g = pt1/(pt1+pt2);
          else            z_g = pt2/(pt1+pt2);

          fastjet::PseudoJet sd_jet_sum = sd_jet_0 + sd_jet_1;

          hJetPtSD2Mass->Fill(sd_jet.pt(),sqrt(sd_jet_sum.m2()),weight);
          hJetPtZ->Fill(sd_jet.pt(), z_g, weight);
  //        hjetptz_CMS->Fill(sd_jet.pt(), z_g, tweight);

  //        float sd_deltaR = sd_jet_0.delta_R(sd_jet_1);
          float sd_deltaR = piecesSoft[0].delta_R(piecesSoft[1]);
          hJetPtDR->Fill(sd_jet.pt(), sd_deltaR, weight);
          if (sd_deltaR > 0.1) {
            hJetPtZ_DRCut->Fill(sd_jet.pt(), z_g, weight);
  //          hJetPtZ_CMS_DRCut->Fill(sd_jet.pt(), z_g, tweight);
          }
        }
        else {
          hJetPtZ->Fill(sd_jet.pt(), 0.0, weight);
        }

         hJetPtZ->Fill(sd_jet.pt(), sd_jet.structure_of<SoftDrop>().symmetry(),weight);
        hSDJetPtMass->Fill(sd_jet.pt(),sqrt(sd_jet.m2()),weight);   
        */

        fastjet::PseudoJet thisJet = jets.at(i);

        // Try some background subtraction here
        double rho = 0.0;
        rho = getBackground(&thisJet, &particles);
        double jet_area = 0;
        if (iTriggerType == 0) jet_area = thisJet.area();
       
//        hSDJetPt->Fill(sd_jet.pt());
        // experiment
//        hSDJetPtBkgSub->Fill(sd_jet.pt() - rho * jet_area); // this makes no sense

        // After filling basic spectra, require hard core (ie particle > 6 GeV)
        // Removing the jets would break the loop above
        //jets.erase( std::remove_if(jets.begin(), jets.end(), rejectJet), jets.end() );
        jet_phi = jets.at(i).phi();
        rot_x = vertex_x * TMath::Cos(PI - jet_phi) - vertex_y * TMath::Sin(PI - jet_phi);
        rot_y = vertex_x * TMath::Sin(PI - jet_phi) + vertex_y * TMath::Cos(PI - jet_phi);
        vertex_r = TMath::Sqrt(pow(vertex_x,2) + pow(vertex_y,2));
        vertexXYRot->Fill(rot_x,rot_y,weight);

        //finding the highest pt constituent of the jet
        double highestPtConst = getHighestPtConst(jets.at(i));

        vertexHighestPtXRot->Fill(highestPtConst,rot_x,weight);
        vertexHighestPtR->Fill(highestPtConst,vertex_r,weight);
        vertexHighestJetZR->Fill(highestPtConst/jet_pt,vertex_r,weight);
        vertexHighestJetZXRot->Fill(highestPtConst/jet_pt,rot_x,weight);


        //Looking at surface bias with Dijet method:
        //looking for leading jets, then subleading jets with delta_phi > 2 pi / 3
        if ((iTriggerType == 0) && !rejectJetDijet(&jets.at(i), R, 20.)) {
          for ( int j = i+1; j < jets.size(); j++ ) {
            // to save time, take advantage of ordering:
            double subleading_pt = jets.at(j).pt();
            if (subleading_pt < 10.) break;

            if (!rejectJetDijet(&jets.at(j), R, 10.) && TMath::Abs(calculateDeltaPhi(jet_phi,jets.at(j).phi())) > 2. * PI/3.) {
              // Then these jets have passed the dijet cut      
              vertexRDijetCut->Fill(vertex_r,weight);
              vertexXYRotDijet->Fill(rot_x,rot_y,weight);
              vertexJetPtRDijet->Fill(jet_pt,vertex_r,weight);

              double a_j = (jet_pt - subleading_pt) / (jet_pt + subleading_pt);
              double x_j = subleading_pt / jet_pt;

              dijetAj->Fill(a_j,weight);
              dijetAjJetPt->Fill(jet_pt,a_j,weight);
              dijetXj->Fill(x_j,weight);
              dijetXjJetPt->Fill(jet_pt,x_j,weight);
              vertexDijetAjR->Fill(a_j,vertex_r,weight);
              vertexXYRotAjDijet->Fill(rot_x,rot_y,a_j,weight);
            }
          }
        }

        vertexJetPtR->Fill(jet_pt, TMath::Sqrt(pow(vertex_x,2) + pow(vertex_y,2)),weight);

        if (iTriggerType == 0) {
          //if (rejectJetHardCore(jets.at(i), R, jet_hard_core_cut) == true) {
          if (rejectJetHardCore(&jets.at(i), R, jet_hard_core_cut) == true) {
            vertexXYNoHardCore->Fill(vertex_x,vertex_y,weight);
            vertexXYRotNoHardCore->Fill(rot_x,rot_y,weight);
            vertexRNoHardCore->Fill(vertex_r,weight);
          } else {
            vertexXYHardCore->Fill(vertex_x,vertex_y,weight);
            vertexXYRotHardCore->Fill(rot_x,rot_y,weight);
            vertexRHardCore->Fill(vertex_r,weight);
            vertexJetPtRHardCore->Fill(jet_pt, TMath::Sqrt(pow(vertex_x,2) + pow(vertex_y,2)),weight);
          }
    //      if (rejectJet(jets.at(i), R)) continue;
      
          //std::vector <fastjet::PseudoJet> rec_jets;
            if (doJetBkgSub) recJetPt->Fill(jet_pt,weight);
            else recJetPt->Fill(thisJet.pt() - jet_area*rho,weight);
          }
        }
      // Do angular correlations for all particles
      double particleEnergyTemp = 0;
      double fLeadingTrackPt = 0;
      for (unsigned int i = 0; i < particles.size(); i++)
      {
        // Eta-Phi of each particle
        etaPhi->Fill(particles.at(i).eta(), particles.at(i).phi(),weight);

        phiPt->Fill(particles.at(i).phi(), particles.at(i).pt(),weight);
        // Eta-Pt of each particle
        // the purpose of this plot is to investigate the double ridge 
        // in all particles in AA
        etaPt->Fill(particles.at(i).eta(), particles.at(i).pt(),weight);

        trackPt->Fill(particles.at(i).pt(),weight);
        if (isChargedParticle(particles.at(i).user_index())) chargedTrackPt->Fill(particles.at(i).pt(),weight);

        fLeadingTrackPt = TMath::Max(fLeadingTrackPt,particles.at(i).pt());

        //MHO fix: we need to double count here 
        // to avoid biases in the ordering of particles
        //      for (unsigned int j = i+1; j < particles.size(); j++)
        for (unsigned int j = 0; j < particles.size(); j++)
        {
          // MHO
          if (i == j) continue;
          // deltaEta-deltaPhi or each particle pair
          // For TLorentzVector
          //dEtaDPhi->Fill(particles.at(i).Eta()-particles.at(j).Eta(), calculateDeltaPhi(particles.at(i).Phi(), particles.at(j).Phi()),weight);
          // For Fastjet
          dEtaDPhi->Fill(particles.at(i).eta()-particles.at(j).eta(), calculateDeltaPhi(particles.at(i).phi(), particles.at(j).phi()) ,weight);
        }
        particleEnergyTemp += particles.at(i).E();
      }
      leadingTrackPt->Fill(fLeadingTrackPt,weight);
      // Fill total particle energy
      particleEnergy->Fill(particleEnergyTemp,weight);

      // Scattering Center research
      double scatCentTotalE = 0;

      for (unsigned int i = 0; i < scatCents.size(); i++) {
        scatCentPt->Fill(scatCents.at(i).pt(),weight); // note: no rejection criteria
        scatCentE->Fill(scatCents.at(i).E(),weight);
        scatCentTotalE += scatCents.at(i).E();
      }
      scatCentEtot->Fill(scatCentTotalE,weight);

      // Clear vector of PseudoJets after each event
      particles.clear();
      scatCents.clear();
      hardScatters.clear();
      constParticles.clear();
      jets.clear();

      //if ((eventNumber % 1000 == 0) && (eventNumber != 0)) { printf("eventNumber = %d\n", eventNumber); }
      if ((eventNumber % 10 == 0) && (eventNumber != 0)) { printf("eventNumber = %d\n", eventNumber); }

    //hashString
//     printf("H %d\nW %e\n",(int) hashString.Hash(),weight);

    // For testing purposes
    //if (eventNumber > 99) {break;}
  }
  

  // If the event plane analysis is running, add up the ep bins to get the inclusive bin
  for (int k = 1; k < nEPBins; k++) { // this actually checks that nEPBins > 1 !
    for (int i = 0; i < nTriggerPtBins; i++) {
      ptLossLeadingJetPtBinEP[0][i]->Add(ptLossLeadingJetPtBinEP[k][i]);
      energyLossLeadingJetPtBinEP[0][i]->Add(energyLossLeadingJetPtBinEP[k][i]);

      ptLossSubleadingJetPtBinEP[0][i]->Add(ptLossSubleadingJetPtBinEP[k][i]);
      energyLossSubleadingJetPtBinEP[0][i]->Add(energyLossSubleadingJetPtBinEP[k][i]);

      PartonPtByJetBin[0][i]->Add(PartonPtByJetBin[k][i]);
      qqbar_PartonPtByJetBin[0][i]->Add(qqbar_PartonPtByJetBin[k][i]);
      qq_PartonPtByJetBin[0][i]->Add(qq_PartonPtByJetBin[k][i]);
      gq_PartonPtByJetBin[0][i]->Add(gq_PartonPtByJetBin[k][i]);
      qg_PartonPtByJetBin[0][i]->Add(qg_PartonPtByJetBin[k][i]);
      gg_PartonPtByJetBin[0][i]->Add(gg_PartonPtByJetBin[k][i]);
    }
    for (int i = 0; i < nJetPtBins; i++) {
      ptLossLeadingPartonPtBinEP[0][i]->Add(ptLossLeadingPartonPtBinEP[k][i]);
      energyLossLeadingPartonPtBinEP[0][i]->Add(energyLossLeadingPartonPtBinEP[k][i]);
      ptLossSubleadingPartonPtBinEP[0][i]->Add(ptLossSubleadingPartonPtBinEP[k][i]);
      energyLossSubleadingPartonPtBinEP[0][i]->Add(energyLossSubleadingPartonPtBinEP[k][i]);
    }
  }
  
  
  //  outputFilename = Form("output/%s.root", outputFilename.Data());
  TFile * fOut = TFile::Open(outputFilename.Data(),"RECREATE");
  if (fOut->IsOpen() == kFALSE)
  {
    Printf("outputFilename \"%s\" failed to open!", outputFilename.Data());
    std::exit(1);
  }
  
  fOut->Add(dEtaDPhi);
  fOut->Add(etaPhi);
  fOut->Add(phiPt);
  fOut->Add(etaPt);
  fOut->Add(dijetAj);
  fOut->Add(dijetAjJetPt);
  fOut->Add(dijetXj);
  fOut->Add(dijetXjJetPt);

  if (leadingJetPtNScatCent[0]->GetEntries() > 0) {
    for (int k = 0; k < nEPBins; k++) {
      fOut->Add(leadingJetPtNScatCent[k]);
    }
  }

  fOut->Add(vertexR);
  fOut->Add(vertexRHardCore);
  fOut->Add(vertexRNoHardCore);
  fOut->Add(vertexXY);
  fOut->Add(vertexXYHardCore);
  fOut->Add(vertexXYNoHardCore);
  fOut->Add(vertexXYRot);
  fOut->Add(vertexXYRotHardCore);
  fOut->Add(vertexXYRotNoHardCore);
  fOut->Add(vertexXYRotDijet);
  fOut->Add(vertexXYRotAjDijet);
  fOut->Add(vertexXRotJetPt);
  fOut->Add(vertexXRotJetPtHardCore);
  fOut->Add(vertexHighestPtXRot);
  
  fOut->Add(ptLossVertexRLeadingJet);
  fOut->Add(ptLossVertexRSubleadingJet);


  for (int k = 0; k < nEPBins; k++) {
    for (int i = 0; i < nJetPtBins; i++) {
      fOut->Add(ptLossLeadingJetPtBinEP[k][i]);
      fOut->Add(energyLossLeadingJetPtBinEP[k][i]);
      fOut->Add(ptLossLeadingPartonPtBinEP[k][i]);
      fOut->Add(energyLossLeadingPartonPtBinEP[k][i]);

      fOut->Add(ptLossSubleadingJetPtBinEP[k][i]);
      fOut->Add(energyLossSubleadingJetPtBinEP[k][i]);
      fOut->Add(ptLossSubleadingPartonPtBinEP[k][i]);
      fOut->Add(energyLossSubleadingPartonPtBinEP[k][i]);
    }
  }

  fOut->Add(vertexJetPtR);
  fOut->Add(vertexJetPtRHardCore);
  fOut->Add(vertexJetPtRDijet);
  fOut->Add(vertexRDijetCut);
  fOut->Add(vertexDijetAjR);
  fOut->Add(vertexHighestPtR);
  fOut->Add(vertexHighestJetZR);
  fOut->Add(vertexHighestJetZXRot);

  fOut->Add(jetPt);
  fOut->Add(leadingJetPt);

  for (int k = 0; k < nEPBins; k++) {
    fOut->Add(leadingJetPtEP[k]);
  }

  fOut->Add(recJetPt);
  fOut->Add(gammaEt);
  fOut->Add(trackPt);
  fOut->Add(chargedTrackPt);
  fOut->Add(leadingTrackPt);
  fOut->Add(jetPtUnweighted);
  fOut->Add(particleEnergy);
  fOut->Add(jetEtaPhi);
  fOut->Add(leadingJetEtaPhi);
  fOut->Add(jetPhiPt);
  fOut->Add(hWeight);
  
  fOut->Add(scatCentPt);
  fOut->Add(scatCentE);
  fOut->Add(scatCentEtot);
  
  fOut->Add(qqbar_pt);
  fOut->Add(qq_pt);
  fOut->Add(gg_pt);
  fOut->Add(qg_pt);
  fOut->Add(hsPtLeadingJetPt);
  fOut->Add(hsPtSubLeadingJetPt);
  fOut->Add(hsDEtaDPhiLeadingJet);
  fOut->Add(hsDEtaDPhiSubLeadingJet);
  
  fOut->Add(vertexXYRot4GeVTrigger);
  fOut->Add(vertexXYRot5GeVTrigger);
  fOut->Add(vertexXYRot20GeVTrigger);
  
  fOut->Add(vertexXRotTriggerPt);

/*  fOut->Add(hJetPtZ);
  fOut->Add(hJetPtSD2Mass);
  fOut->Add(hJetPtDR);
  fOut->Add(hJetPtZ_DRCut);

  fOut->Add(hJetPtMass);  
  fOut->Add(hSDJetPtMass);  
  fOut->Add(hSDJetPt);
  fOut->Add(hSDJetPtBkgSub);
*/
  fOut->Add(hs1hs2Pt);
  fOut->Add(hs1hs2dPhi);

  for (unsigned int i = 0; i < nTriggerPtBins; i++)
  {
    fOut->Add(ptBinHadronPt[i]);
    fOut->Add(swiftJetYield[i]);
    fOut->Add(vertexXYJetPt[i]);
    fOut->Add(vertexXYRotJetPt[i]);
    if (iTriggerType == 0) fOut->Add(vertexXYRotJetPtHardCore[i]);
    
    for (unsigned int j = 0; j < nAssocParticleBins; j++)
    {
      fOut->Add(jetHadron[i][j]);
      fOut->Add(jetHadronScatCent[i][j]);
      fOut->Add(jetHadronOnlyJet[i][j]);
      fOut->Add(jetHadronExcludeJet[i][j]);
      fOut->Add(swiftParticleYield[i][j]);
    }
  }

  for (unsigned int k = 0; k < nEPBins; k++)
  {
    for (unsigned int i = 0; i < nJetPtBins; i++)
    {
      fOut->Add(PartonPtByJetBin[k][i]);
      fOut->Add(qqbar_PartonPtByJetBin[k][i]);
      fOut->Add(qq_PartonPtByJetBin[k][i]);
      fOut->Add(gq_PartonPtByJetBin[k][i]);
      fOut->Add(qg_PartonPtByJetBin[k][i]);
      fOut->Add(gg_PartonPtByJetBin[k][i]);
    }
  }
  for (unsigned int k = 0; k < nEPBins; k++)
  {
    for (unsigned int i = 0; i < nTriggerPtBins; i++)
    {

      fOut->Add(ptBinHadronPtEP[k][i]);
      for (unsigned int j = 0; j < nAssocParticleBins; j++)
      {
        fOut->Add(jetHadronEP[k][i][j]);  
        //fOut->Add(jetHadronScatCentEP[k][i][j]);  
      }
    }
  }


  // useful code from rene
  // https://root.cern.ch/root/roottalk/roottalk02/3719.html
  
  // July 21: 
  // new idea: store sum of weights in a single entry in a branch of a TTree
  // then when merging the files, the number of entries will be the number of 
  // files.  Scale everything by 1 / n_files.
  // Can also add cross section later.  This requires reading the log file.
  // Depending on where I store the log files, it may be necessary to add another 
  // argument to phase1, indicating the path to the log.
  //
  // We could also add a branch for the energy, and for the type of event
  // That information is also in the logfile
  
  // August 19:
  // This will now be created first in the HepMC2Root stage.
  // For compatibility, it will also be created here if it does not exist
  // FIXME
  
    

  if (!runInfo ) { //for compatibility 
 // if (!runInfo) delete runInfo; // Making a new runInfo in case the previous one was missing something
    printf("TTree runInfo not found.  Creating new one ... \n");
    runInfo = new TTree("runInfo","runInfo");
    TBranch *bTotalWeight = runInfo->Branch("totalWeight", &total_weight, "totalWeight/D");  
    TBranch *bNEvents = runInfo->Branch("nEvents", &nevents, "nEvents/I");  
//    TBranch *bCrossSection = runInfo->Branch("crossSection", &storedCrossSection, "crossSection/D");

    // to add: cross secion, sqrt(s), type
    // pt range?
    runInfo->Fill();
  }
  else {
    printf("TTree runInfo found!\n");
    // check if total weight, cross section info filled.
    // Not running if multiple entries in runInfo (ptHardBins on )
    if (runInfo->GetEntries() == 1) {
      double vTotalWeight = 0;
      int vNEvents = 0;
      double vNEventsDouble = 0;
      double vCrossSection = 0;
      runInfo->SetBranchAddress("totalWeight",&vTotalWeight);
      if (!strcmp(runInfo->GetBranch("nEvents")->GetTitle(),"nEvents/I")) runInfo->SetBranchAddress("nEvents",&vNEvents);
      else runInfo->SetBranchAddress("nEvents",&vNEventsDouble);
      runInfo->SetBranchAddress("crossSection",&vCrossSection);
      runInfo->GetEntry(0);
      printf("MHO: got nEvents = %d, totalWeight = %f, crossSection = %f\n",vNEvents,vTotalWeight,vCrossSection);

      if (vTotalWeight == 0 || (vNEvents + (int) vNEventsDouble == 0)) {  
        printf("Correcting total weight or nEvents branch.\n");
    //    vTotalWeight = total_weight;
    //    vNEvents = nevents;  
        // create a new runinfo tree cause this one is broken.
        delete runInfo;
        runInfo = new TTree("runInfo","runInfo");
        TBranch *bTotalWeight = runInfo->Branch("totalWeight", &total_weight, "totalWeight/D");  
        TBranch *bNEvents = runInfo->Branch("nEvents", &nevents, "nEvents/I");  
        TBranch *bCrossSection = runInfo->Branch("crossSection", &vCrossSection, "crossSection/D");  
        runInfo->Fill();
  
      } 
      else { 
        // Everything is fine.
        runInfo->CloneTree(); 
      }
    } 
    else if (runInfo->GetEntries() == 0) { 
      printf("TTree runInfo has zero entries\n");
    } else {
      runInfo->CloneTree(); 
    }
  }
  
  // Creating TTree for parameter information
  TTree * parameters = new TTree("parameters","parameters");
  TBranch * bJetR = parameters->Branch("jetR",&R,"jetR/D");
  TBranch * bConstCut = parameters->Branch("constCut",&jet_constituent_cut,"constCut/D");
  TBranch * bEtaCut = parameters->Branch("etaCut",&eta_cut,"etaCut/D");
  TBranch * bParticlePtCut = parameters->Branch("particlePtCut",&particle_pt_min,"particlePtCut/D");
  
  parameters->Fill();
  
  
  
  
  
  
  
  
  if (false && weightingOn) {
  //if (weightingOn) {
  
    // I have no idea what i should be doing here
  
    //total_weight = nevents / total_weight;
    //FIXME  
    //double rescale = 1 / total_weight;
    double rescale = nevents / total_weight;
    printf("Rescaling histograms...\n");
  
    TObjLink *lnk = gDirectory->GetList()->FirstLink();
    TObject *obj;
    TH1 *hobj;
    while (lnk) {
      obj = gDirectory->FindObject(lnk->GetObject()->GetName());    
      obj = gDirectory->Get(lnk->GetObject()->GetName());    
      //   printf("name: %s named %s\n",obj->ClassName(),obj->GetName()); 
      lnk = lnk->Next();
      if (obj->InheritsFrom("TH1") && strcmp(obj->GetName(),"hWeight") && strcmp(obj->GetName(),"jetPtUnweighted")) {
        //      printf("\tscaling %s by %f \n",obj->GetName(),total_weight);
        hobj = (TH1 *) obj;
        hobj->Scale(rescale);
      }    
    }
  } 
  
  printf("Analyzed %d good events out of %d events!\n",nGoodEvents,nevents); 
  
  
  fOut->Write();
  fOut->Close();
  
  
  
  
  return;
  
  //CUT HERE-ISH

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
  TString inputFilename = "root/pp.root";
  TString outputFilename = "root/pp.hist.root";
  TString responseMatrixFilename = ""; 
  // Defined in analysis_params.h
  //double jet_constituent_cut = 0;
  bool ptHardBinFlag = false;

  int iTriggerType = 0; // 0 for jet, 1 for pi0

  // Parse input and set the appropriate values
  parseInput(argc, argv, inputFilename, outputFilename, jet_constituent_cut, ptHardBinFlag, responseMatrixFilename,iTriggerType);

  if (!inputFilename.CompareTo(outputFilename)) {
    fprintf(stderr,"Error: output file is same as input file.\n");
    exit(1);
  }

  phase1(inputFilename, outputFilename, jet_constituent_cut, ptHardBinFlag, responseMatrixFilename, iTriggerType);

  #ifdef RE_STYLE
  delete readableStyle;
  #endif

  return 0;
}

// Process command line inputs. This is not nearly as robust as using something like Boost::program_options,
// but that seemed a bit heavy for this instance, making program_options less useful. The code here is simple,
// but it worked well and allowed for quick development, while getting the job done.
void parseInput(int argc, char * argv[], TString & inputFilename, TString & outputFilename, double & jet_constituent_cut, bool & ptHardBinFlag, TString & responseMatrixFilename, int & iTriggerType)
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
    if ((argv[i] == std::string("--ptHardBin")) || (argv[i] == std::string("-p")))
    {
      // Use pt hard bin 
      ptHardBinFlag = true;
      continue;
    }
    if ((argv[i] == std::string("--eventPlane")) || (argv[i] == std::string("-e")))
    {
      // Bin Jet-Hadron correlations by angle relative to event plane
//      nEPBins_dyn = 4;
      nEPBins = 4;
      continue;
    }
    if ((argv[i] == std::string("--phi")) || (argv[i] == std::string("-a")))
    {
      phiAcceptanceSim = true;
      continue;
    }
    if ((argv[i] == std::string("--charged")) || (argv[i] == std::string("-ch")))
    {
      requireChargedAssoc = true;
      continue;
    }
    if ((argv[i] == std::string("--triggerType")) || (argv[i] == std::string("-t")) ){
      if (argc > i+1)
      {
        if (argv[i+1] == std::string("jet")) {
          iTriggerType = 0; // jet
          i++;
          continue;
        }
        else if (argv[i+1] == std::string("pi0"))
        {
          iTriggerType = 1; //pi0
          i++;
          continue;
        } else {
          std::cout << "Error: unexpected trigger type " << argv[i+1] << std::endl;
        }
      }
      else
      {
        std::cout <<"A choice of trigger must follow the " << argv[i] << " flag." << std::endl;
      }
    }
    if ((argv[i] == std::string("--Zt")) || (argv[i] == std::string("-z")))
    {
      bUseZtBins = true;
      continue;
    }


    if ((argv[i] == std::string("--sign")) || (argv[i] == std::string("-s")) )
    {
      // Check for input value
      if (argc > i+1)
      {
        iOddEven = std::stoi(argv[i+1]);
        i++;
        continue;
      }
      else
      {
        std::cout << "An value for the sign (or modulus) of events to use!" << std::endl;
        printHelp();
      }
    }


    if ((argv[i] == std::string("--numToyParticles")) || (argv[i] == std::string("-n")) )
    {
      // Check for input value
      if (argc > i+1)
      {
        nToyParticles = std::stoi(argv[i+1]);
        i++;
        continue;
      }
      else
      {
        std::cout << "An value for the number of toy particles must be passed!" << std::endl;
        printHelp();
      }
    }
		if ((argv[i] == std::string("--flow")) || (argv[i] == std::string("-v")) ){
		// Check for input value
			if (argc > i+1) {
				iFlowVersion = std::stoi(argv[i+1]);
				i++;
				continue;
			}	
			else 
			{
        std::cout << "A choice of flow parameters must be given!" << std::endl;
        printHelp();
      }
		}
    if ((argv[i] == std::string("--constituentCut")) || (argv[i] == std::string("-c")) )
    {
      // Check for consituent cut value
      if (argc > i+1)
      {
        jet_constituent_cut = std::stod(argv[i+1]);
        i++;
        continue;
      }
      else
      {
        std::cout << "An value for the constituent cut must be passed!" << std::endl;
        printHelp();
      }
    }
    if ((argv[i] == std::string("--jetResolutionParameter")) || (argv[i] == std::string("-r")) )
    {
      if (argc > i+1)
      {
        R = std::stod(argv[i+1]);
        i++;
        continue;
      }
      else
      {
        std::cout << "An value for the jet resolution parameter must be passed!" << std::endl;
        printHelp();
      }
    }
//    << std::setw(5) << std::left << "\t-m" << "\t--responseMatrix <filename>" 
    if ((argv[i] == std::string("--responseMatrix")) || (argv[i] == std::string("-m")) )
    {
      // Sets response matrix filename
      if (argc > i+1)
      {
        responseMatrixFilename = argv[i+1];
        i++;
        continue;
      }
      else
      {
        std::cout << "An reponse matrix filename must be passed!" << std::endl;
        printHelp();
      }
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
    if ((argv[i] == std::string("--hardCore")) || (argv[i] == std::string("-k")) )
    {
      // Check for hard core cut value
      if (argc > i+1)
      {
        globalHardCoreCut = std::stod(argv[i+1]);
        if (globalHardCoreCut > 0) requireHardCore = true; 
        i++;
        continue;
      }
      else
      {
        std::cout << "A value for the global hard core cut must be passed!" << std::endl;
        printHelp();
      }
    }
    if ((argv[i] == std::string("--doJetBkgSub")) || (argv[i] == std::string("-b")) ) 
    {
      doJetBkgSub = true; 
    }
    else
    {
      printHelp(argv[i]);
    }
  }
}

// Print help for using pionDecay and exit afterwards
// Default argument defined at the top of the file
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
    << "\t-> Sets the output filename. Default: \"root/pp.hist.root\"" << std::endl      
    << std::setw(5) << std::left << "\t-t" << "\t--triggerType [pi0 or jet]" 
    << "\t-> Sets whether pi0s or jets or used as the trigger. Default: jet." << std::endl
    << std::setw(5) << std::left << "\t-c" << "\t--constituentCut <value>" 
    << "\t-> Sets the constituent cut for jet finding. Default: 0." << std::endl
    << std::setw(5) << std::left << "\t-r" << "\t--jetResolutionParameter <value>" 
    << "\t-> Sets the jet resolution parameter R for jet finding. Default: 0." << std::endl
    << std::setw(5) << std::left << "\t-k" << "\t--hardCore <value>" 
    << "\t-> Sets the requirement of a hard core to be true, and sets global hard core cut. Default: false." << std::endl
    << std::setw(5) << std::left << "\t-z" << "\t--Zt"
    << "\t\t\t-> Uses z_t = (p_t^assoc) / (p_t^trigger) for associated particle binning. Default: false." << std::endl 
    << std::setw(5) << std::left << "\t-p" << "\t--ptHardBin"
    << "\t\t\t-> Enables pt hard bin weighting. Default: false." << std::endl 
    << std::setw(5) << std::left << "\t-e" << "\t--eventPlane"
    << "\t\t\t-> Enables event plane binning. Default: false." << std::endl 
    << std::setw(5) << std::left << "\t-ch" << "\t--charged"
    << "\t\t\t-> Requires associated particles be charged. Default: false." << std::endl 
    << std::setw(5) << std::left << "\t-m" << "\t--responseMatrix <filename>" 
    << "\t-> The path to the root file with the jet energy response matrix." << std::endl      
    << std::setw(5) << std::left << "\t-a" << "\t--phi"
    << "\t\t\t-> Applies a phi acceptance simulation. Default: false." << std::endl 
    << std::setw(5) << std::left << "\t-b" << "\t--doJetBkgSub"
    << "\t\t\t-> Uses an out of cone background subtraction for jet pT for jetPt bins. Default: false." << std::endl
    
    << std::setw(5) << std::left << "\t-s" << "\t--sign"
    << "\t\t\t\t\t-> Sets the usage of only even or odd numbered events"<<std::endl<<std::endl

    << std::setw(5) << std::left << "\t-n" << "\t--numToyParticles"
    << "\t\t\t-> Sets the number of particles to produce from a toy model. Set to 0 to disable toy model. Default: 100." << std::endl << std::endl
    << std::setw(5) << std::left << "\t-v" << "\t--flow [version]"
    << "\t\t\t-> Enables parametrized flow simulation. Version 0,1,2,3 correspond to 5.02 PbPb centralities 0-10,10-30,30-50,50-80. Set to -1 to disable. Default: -1" << std::endl << std::endl;

  std::exit(-1);
}
