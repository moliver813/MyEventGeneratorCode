#!/bin/bash

# pp Section

#INPUTDIR=eventPlane/eta_09/phiAcceptanceSim/chargeReq
INPUTDIR=eventPlane/eta_09/latest/chargeReq

INFILE=$INPUTDIR/pp_2760GeV_all_R_0.2_K_0_CCUT_3_.hist.root


./analysis/phase2 -i $INFILE -o eventPlane/pp_2760GeV_FastMixed_ZYAM.final.root -b 0 -B ZYAM

##./analysis/phase2 -i $INFILE -o eventPlane/pp_2760GeV_FastMixed_dEtaGenGausFit.final.root -b dEtaGenGausFitSub -B 0
##./analysis/phase2 -i $INFILE -o eventPlane/pp_2760GeV_FastMixed_dEta2GausFit.final.root -b dEtaBkgSub -B 0

##./analysis/phase2 -i $INFILE -o eventPlane/pp_2760GeV_FastMixed_FreeB_f_0.final.root -b 0 -B FreeDPhiBkg -f 0
##./analysis/phase2 -i $INFILE -o eventPlane/pp_2760GeV_FastMixed_FreeB_f_1.final.root -b 0 -B FreeDPhiBkg -f 1

./analysis/phase2 -i $INFILE -o eventPlane/pp_2760GeV_FastMixed_FreeB_f_2.final.root -b 0 -B FreeDPhiBkg -f 2
./analysis/phase2 -i $INFILE -o eventPlane/pp_2760GeV_FastMixed_FreeB_f_3.final.root -b 0 -B FreeDPhiBkg -f 3
./analysis/phase2 -i $INFILE -o eventPlane/pp_2760GeV_FastMixed_FreeB_f_4.final.root -b 0 -B FreeDPhiBkg -f 4


./analysis/calcSystUncert eventPlane/SysUncert_pp_2760GeV_R_0.2_K_0_CCUT_3.root eventPlane/pp_2760GeV_FastMixed_ZYAM.final.root "ZYAM" eventPlane/pp_2760GeV_FastMixed_FreeB_f_2.final.root "FreeB, F 2" eventPlane/pp_2760GeV_FastMixed_FreeB_f_3.final.root "FreeB, F 3" eventPlane/pp_2760GeV_FastMixed_FreeB_f_4.final.root "FreeB, F 4"
#./analysis/calcSystUncert eventPlane/SysUncert_pp_2760GeV_R_0.2_K_0_CCUT_3.root eventPlane/pp_2760GeV_FastMixed_ZYAM.final.root "ZYAM" eventPlane/pp_2760GeV_FastMixed_FreeB_f_2.final.root "FreeB, F 2" eventPlane/pp_2760GeV_FastMixed_FreeB_f_3.final.root "FreeB, F 3" eventPlane/pp_2760GeV_FastMixed_FreeB_f_4.final.root "FreeB, F 4" eventPlane/pp_2760GeV_FastMixed_dEtaGenGausFit.final.root "dEta GenGaus Fit" eventPlane/pp_2760GeV_FastMixed_dEta2GausFit.final.root "dEta 2Gaus Fit"

mkdir -p eventPlane/pp_2760GeV_MethodComparison 

./analysis/compare eventPlane/pp_2760GeV_MethodComparison eventPlane/pp_2760GeV_FastMixed_ZYAM.final.root "ZYAM" eventPlane/pp_2760GeV_FastMixed_FreeB_f_2.final.root "FreeB, F 2" eventPlane/pp_2760GeV_FastMixed_FreeB_f_3.final.root "FreeB, F 3" eventPlane/pp_2760GeV_FastMixed_FreeB_f_4.final.root "FreeB, F 4"
#./analysis/compare eventPlane/pp_2760GeV_MethodComparison eventPlane/pp_2760GeV_FastMixed_ZYAM.final.root "ZYAM" eventPlane/pp_2760GeV_FastMixed_FreeB_f_2.final.root "FreeB, F 2" eventPlane/pp_2760GeV_FastMixed_FreeB_f_3.final.root "FreeB, F 3" eventPlane/pp_2760GeV_FastMixed_FreeB_f_4.final.root "FreeB, F 4" eventPlane/pp_2760GeV_FastMixed_dEtaGenGausFit.final.root "dEta GenGaus Fit" eventPlane/pp_2760GeV_FastMixed_dEta2GausFit.final.root "dEta 2Gaus Fit"

#./analysis/compare eventPlane/pp_2760GeV_MethodComparison eventPlane/pp_2760GeV_FastMixed_ZYAM.final.root "FastMixed, ZYAM" eventPlane/pp_2760GeV_FastMixed_FreeB_f_0.final.root "FastMixed, FreeB, F 0" eventPlane/pp_2760GeV_FastMixed_FreeB_f_1.final.root "FastMixed, FreeB, F 1" eventPlane/pp_2760GeV_FastMixed_FreeB_f_2.final.root "FastMixed, FreeB, F 2" eventPlane/pp_2760GeV_FastMixed_FreeB_f_3.final.root "FastMixed, FreeB, F 3" eventPlane/pp_2760GeV_FastMixed_FreeB_f_4.final.root "FastMixed, FreeB, F 4"



#LIST="eventPlane/pp_2760GeV_FastMixed_ZYAM.final.root \"FastMixed, ZYAM\" eventPlane/pp_2760GeV_FastMixed_FreeB_f_0.final.root "FastMixed, FreeB, F 0" eventPlane/pp_2760GeV_FastMixed_FreeB_f_1.final.root "FastMixed, FreeB, F 1" eventPlane/pp_2760GeV_FastMixed_FreeB_f_2.final.root "FastMixed, FreeB, F 2" eventPlane/pp_2760GeV_FastMixed_FreeB_f_3.final.root "FastMixed, FreeB, F 3" eventPlane/pp_2760GeV_FastMixed_FreeB_f_4.final.root "FastMixed, FreeB, F 4"
#LIST='eventPlane/pp_2760GeV_FastMixed_ZYAM.final.root "FastMixed, ZYAM" eventPlane/pp_2760GeV_FastMixed_FreeB_f_0.final.root "FastMixed, FreeB, F 0" eventPlane/pp_2760GeV_FastMixed_FreeB_f_1.final.root "FastMixed, FreeB, F 1" eventPlane/pp_2760GeV_FastMixed_FreeB_f_2.final.root "FastMixed, FreeB, F 2" eventPlane/pp_2760GeV_FastMixed_FreeB_f_3.final.root "FastMixed, FreeB, F 3" eventPlane/pp_2760GeV_FastMixed_FreeB_f_4.final.root "FastMixed, FreeB, F 4"'
