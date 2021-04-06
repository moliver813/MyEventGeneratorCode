# Detect OS
# From: http://stackoverflow.com/a/12099167
ifeq ($(shell uname -s),Darwin)
  extraLinkOption=
else
  extraLinkOption=-Wl,--no-as-needed
endif

default:  phase1 phase2 phase3 compareConstCut compare calcSystUncert

phase1: phase1.cc analysis_params.h
	g++ $(USER_DEFINED) -O2 -std=c++11 $(extraLinkOption) phase1.cc -o $@ -I`root-config --incdir` -I$(MYINSTALL)/include `root-config --libs` `fastjet-config --cxxflags --libs --plugins` -lRecursiveTools

phase2: phase2.cc fitAlgos.o analysis_params.h phase2.h
	g++ $(USER_DEFINED) -std=c++11 $(extraLinkOption) phase2.cc fitAlgos.o -o $@ -I`root-config --incdir` -I$(MYINSTALL)/include `root-config --libs` `fastjet-config --cxxflags --libs --plugins`

phase3: phase3.cc analysis_params.h
	g++ $(USER_DEFINED) -std=c++11 $(extraLinkOption) phase3.cc -o $@ -I`root-config --incdir` -I$(MYINSTALL)/include `root-config --libs` 

fitAlgos.o: fitAlgos.cc analysis_params.h phase2.h
	g++ $(USER_DEFINED) -std=c++11 $(extraLinkOption) -c fitAlgos.cc -o $@ -I`root-config --incdir` -I$(MYINSTALL)/include `root-config --libs` `fastjet-config --cxxflags --libs --plugins` 

compareConstCut: compareConstCut.cc analysis_params.h
	g++ $(USER_DEFINED) -std=c++11 $(extraLinkOption) compareConstCut.cc -o $@ -I`root-config --incdir` -I$(MYINSTALL)/include `root-config --libs` 

compare: compare.cc
	g++ $(USER_DEFINED) -std=c++11 $(extraLinkOption) compare.cc -o $@ -I`root-config --incdir` -I$(MYINSTALL)/include `root-config --libs` 

calcSystUncert: calcSystUncert.cc
	g++ $(USER_DEFINED) -std=c++11 $(extraLinkOption) calcSystUncert.cc -o $@ -I`root-config --incdir` -I$(MYINSTALL)/include `root-config --libs` 
