#!/bin/bash


#INPUTDIR=eventPlane/eta_09/phiAcceptanceSim/chargeReq
INPUTDIR=eventPlane/eta_09/latest/chargeReq

INFILE=$INPUTDIR/PbPb_NR_T_485_2760GeV_CENT_30_50_all_R_0.2_K_0_CCUT_3_.hist.root


./analysis/phase2 -i $INFILE -o eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_ZYAM.final.root -b 0 -B ZYAM

##./analysis/phase2 -i $INFILE -o eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_dEtaGenGausFit.final.root -b dEtaGenGausFitSub -B 0
##./analysis/phase2 -i $INFILE -o eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_dEta2GausFit.final.root -b dEtaBkgSub -B 0

##./analysis/phase2 -i $INFILE -o eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_0.final.root -b 0 -B FreeDPhiBkg -f 0
##./analysis/phase2 -i $INFILE -o eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_1.final.root -b 0 -B FreeDPhiBkg -f 1
./analysis/phase2 -i $INFILE -o eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_2.final.root -b 0 -B FreeDPhiBkg -f 2
./analysis/phase2 -i $INFILE -o eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_3.final.root -b 0 -B FreeDPhiBkg -f 3
./analysis/phase2 -i $INFILE -o eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_4.final.root -b 0 -B FreeDPhiBkg -f 4


./analysis/calcSystUncert eventPlane/SysUncert_PbPb_NR_T_485_2760GeV_CENT_30_50_R_0.2_K_0_CCUT_3.root eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_ZYAM.final.root "ZYAM" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_2.final.root "FreeB, F 2" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_3.final.root "FreeB, F 3" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_4.final.root "FreeB, F 4"

mkdir -p eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_MethodComparison 

./analysis/compare eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_MethodComparison eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_ZYAM.final.root "ZYAM" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_2.final.root "FreeB, F 2" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_3.final.root "FreeB, F 3" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_4.final.root "FreeB, F 4"

#./analysis/compare eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_MethodComparison eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_ZYAM.final.root "FastMixed, ZYAM" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_0.final.root "FastMixed, FreeB, F 0" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_1.final.root "FastMixed, FreeB, F 1" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_2.final.root "FastMixed, FreeB, F 2" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_3.final.root "FastMixed, FreeB, F 3" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_4.final.root "FastMixed, FreeB, F 4"



#LIST="eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_ZYAM.final.root \"FastMixed, ZYAM\" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_0.final.root "FastMixed, FreeB, F 0" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_1.final.root "FastMixed, FreeB, F 1" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_2.final.root "FastMixed, FreeB, F 2" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_3.final.root "FastMixed, FreeB, F 3" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_4.final.root "FastMixed, FreeB, F 4"
#LIST='eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_ZYAM.final.root "FastMixed, ZYAM" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_0.final.root "FastMixed, FreeB, F 0" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_1.final.root "FastMixed, FreeB, F 1" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_2.final.root "FastMixed, FreeB, F 2" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_3.final.root "FastMixed, FreeB, F 3" eventPlane/PbPb_NR_T_485_2760GeV_CENT_30_50_FastMixed_FreeB_f_4.final.root "FastMixed, FreeB, F 4"'
