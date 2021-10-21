#ifndef ANALYSIS_PARAMS_H
#define ANALYSIS_PARAMS_H

// include file for pt bins, and cut information.

  const int cTriggerColorList[5] = {kViolet+10,kMagenta-3,kPink,kOrange+7,kSpring-1};

  // Choices for 2D background subtraction
  const int Nbkg2DMethods = 7;  
  const char *bkg2DMethodArray[Nbkg2DMethods] = {"NoBkgSub","2DBkgSub","dEtaBkgSub","dEtaFarSub","dEtaGenGausFitSub","dEtaMixedGausFitSub"};
/*
  const int aNoBkgSub = 0;
  const int a2DFitSub = 1;
  const int aDEtaFitSub = 2;
  const int aDEtaFarSub = 3;
	const int aDEtaGenGausFitSub = 4;
*/
  // Choices for dPhi Projection Fit Function
  const int NdPhiFitMethods = 7;  
  const char *dPhiFitMethodArray[NdPhiFitMethods] = {"1G1G","1M1M","2G1M","2G1GG","1GG1GG","2G1SH","2G2G"};
  // Format a[NSFit][ASFit]
  // G = Gaussian
  // GG = Generalized Gaussian 
  // MG = Modified Gaussian
  // SH = sech(|x|^beta)

  // Both Fwhm and Sigma versions are used

/*  const int a1M1M = 0; 
  const int a2G1M = 1;
  const int a2G1GG = 2;
  const int a2G1SH = 3;
  // const int a1G1G = 4; // Does this still work?
*/
  //Choices for dPhi Projection Fit Background
  const int NdPhiBkgMethods = 3;  
  const char *dPhiBkgMethodArray[NdPhiBkgMethods] = {"NoDPhiBkg","FreeDPhiBkg","ZYAM"};

  // Choices for the delta eta range for the far eta range
  const int NdEtaCutModes = 3;
  const char *dEtaCutModesArray[NdEtaCutModes] = {"Default","3Sigma","FarOut"};


//  int aNoDPhiBkg = 0; // B = 0
//  int aFreeDPhiBkg = 1; // B is a free parameter
//  int aZYAM = 2;  // B is fixed by ZYAM

	int nHWeightBins = 55;
	std::vector <double> hWeightBins = {0,1e-6,2e-6,3e-6,4e-6,5e-6,6e-6,7e-6,8e-6,9e-6,1e-5,2e-5,3e-5,4e-5,5e-5,6e-5,7e-5,8e-5,9e-5,1e-4,2e-4,3e-4,4e-4,5e-4,6e-4,7e-4,8e-4,9e-4,1e-3,2e-3,3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3,1e-2,2e-2,3e-2,4e-2,5e-2,6e-2,7e-2,8e-2,9e-2,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};


  // Only used for Phi vs Pt 
	int nTrackPtBins = 33;
	std::vector <double> trackPtBins = {0.15,0.20,0.25,0.30,0.35,0.5,0.75,1.0,
		1.25,1.5,2,2.5,3.0,3.5,4,5,6,7,8,9,10,12.5, 
		15,20,25,30,35,40,45,50,60,80,100,200};


 // enum { nJetPtBins = 4, nParticlePtBins = 11 };
    int nEPBins = 1;  // 1 For no EP bins, 4 for 4 bins. Duh.

    int nJetPtBins = 4;
    int nJetAssocPtBins = 11;
    
    int nPi0PtBins = 5;
    int nPi0ParticlePtBins = 8;//11
    int nZtBins    = 7 ; 


    int nGammaPtBins = 7;


    int nAssocParticleBins = nJetAssocPtBins; // either particlePtBins, pi0ParticlePtBins, or zTbins

    // oh my god this is poorly written 

    // Define the pt bins
    // Size of the bins is from STAR Jet-h paper
    // Desired bins are 0 - 2 based on the convention below.
  //std::vector <double> jetPtBins = {10,15,20,40};
  std::vector <double> jetPtBins = {10,15,20,40,60};
    // Desired bins are 0 - 10 based on the convention below.
  std::vector <double> particlePtBins = {0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 10, 15};
  //std::vector <double> particlePtBins = {0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 17};


    // Bins for Pi0-h analysis
    std::vector <double> pi0PtBins = {5,7,9,11,14,17};
    std::vector <double> gammaPtBins = {11,14,17,20,25,30,40,50};
    //std::vector <double> pi0ParticlePtBins = {0.15,0.4,0.8,1.45,2.5,4.2,6.95,11.4,18.6};
    // Cleaned up pi0 associated particle pt bins
    std::vector <double> pi0ParticlePtBins = {0.2,0.4,0.8,1.5,2.5,4,7,11,17};
    //{0.15, 0.25, 0.5,   1, 1.5, 2,   3, 4, 5,   6, 8, 10};


    //std::vector <double> zTBins    = {0,1./6, 2./6,3./6,4./6,5./6,1.};
    //std::vector<double> zTBins    = {0,1./6, 2./6,3./6,4./6,5./6,1.,7./6,8./6,9./6,2.};
    std::vector <double> zTBins = {0,1./6,2./6, 3./6,4./6,5./6,  6./6,7./6};

  // Bins for HighPt analyses: eg. soft drop
    int nHighPtBins = 14;
  std::vector <double> highPtBins = {20,40,60,80,100,120,140,160,180,200,250,300,400,500,1000};


    // Minimum pt for all particles
    double particle_pt_min = 0.15;
    // Minimum pt for particles included in Jet Reconstruction
    // Warning: this can be changed in run-time by phase1.cc 
    double jet_constituent_cut = 3.0;
    // Minimum pt for at least one constituent of a jet with a hard core
    // used for specifc analyses.  For a global hard core cut on jets, see
    // globalHardCoreCut
    double jet_hard_core_cut = 6.0;
    double eta_cut = 3.0; // 0.9; FIXME trying large eta again
  //	double eta_cut = 7; // 0.9; FIXME trying large eta again
    double jet_eta_cut = 0.5;	

    double pi0ThetaCutHighPt = 0.017; // High pT cut
    double pi0ThetaCutLowPt = 0.023; // Low pT cut


  double pi0_pt_min = 5;
  double pi0_eta_cut = 0.7;

  // FIXME should be <= 1.8 (avoid influence of non-QGP region
	double delta_eta_cut = 1.8; // 1.2

  double trigger_eta_cut = jet_eta_cut;
	double rmsRange = 1.047;
//	double rmsRange = 3.14159/2;

// jet cone radius / resolution parameter
	double R = 0.2;

//	double eta_range = 2. * (eta_cut  - R);
//	double eta_range = 2. * (eta_cut  - R);
	double eta_range = 2 * jet_eta_cut; // total eta range of jets

	double delta_eta_range = jet_eta_cut + eta_cut; // delta eta is in +- delta_eta_range

  // Delta eta cut for the NS peak
  // FIXME I think I overrode the original version of this
//  double fDeltaEtaPeakCut = 0.6;
  double fFarOutDeltaEtaMin = 1.7;
  double fFarOutDeltaEtaMax = 2.2;

  double fDeltaEtaPeakCut = 0.8;


  bool requireHardCore = false;
  double globalHardCoreCut = 6.0;

  double c_width = 600;
  double c_height = 600;

	double c_width_large = 1200;
	double c_height_large = 900;


  bool use2DFit = false;
  bool useSecHFit = false;
  bool useModGausFit = true;
  bool useZYAM = true;


  // default settings
  int bkg2DMethod = 0;
  int dPhiFitMethod = 3;
  int dPhiBkgMethod = 1;
  int dEtaCutMode = 0;

  bool doJetBkgSub = false;  // use background subtraction for jetPt bins? (affects surf bias, jet-h correlations)

	bool noAwaySide = false; // Switch just for comparing the nearside peak from single peak simulations.
	bool noNearSide = false; // Switch just for comparing the awayside peak from single peak simulations.
	bool noFitting = false; // Don't do final fits, just show pre-fits.

	bool plotSigmaFits = true; // Use the sigma fits when ploting azimuthal projections with fits

  // what should default be ?
	bool acceptanceCorrection = false; 
	bool fastMixed = true;
	//int fastMixingMode = 2; // By default, use both jet and associate particle yields
  int fastMixingMode = 1; // Convolve the jet yield with a constant yield in associated particles

  int iOddEven = 0; // 0:all,1:odd,2:even; 3 = 1/4, 4 = 2/4, 5 = 3/4 6 = 4/4

  bool bEnableToyModel = true;
  float ToyEtaRange = 5; //2 //5
  int nToyParticles = 200; // how many to generate per event
  int iFlowVersion = -1; // -1 = disable, 0 = cent0, 1 = cent1, 2 = cent2, 3 = cent3
  const int kToyParticleStatus = 8; // another way of labelling toy particles
  const int kToyParticleLabel = -7; // label for toy particles. technically this is the PDG bbar-prime
  const int kToyPi0Label = -111; // For toy particles sampled to be treated as pi0 triggers.
  const int kToyPhotonLabel = -22; // Toy photons
  const float kToyParticleWeightScale = 1.0; // Downscale of toy particles to keep reasonable S/B and avoid toy model particles uncertainties overwhelming signal

  bool bLeadingJetOnly = false;
	bool phiAcceptanceSim = false; // only useful for leading jet requirement
	bool requireChargedAssoc = false;
  bool bUseZtBins = false; // Switch for using Zt bins instead of assoc pt bins

  bool ApplySystUncertBeforeRatio = false; // Whether to apply sys uncert before taking EP ratios

#endif
