#!/bin/bash

IND=(0 1 2 3 4 5)
#IND=(0 5)
TAUIS=(02 04 06 08 10 12)
TAUIVALUES=(0.2 0.4 0.6 0.8 1.0 1.2)

LISTOFOBJECTS="trackPt particleEnergy jetPt leadingJetPt \
OutOverIn_AS_pi0Pt_5_7 \
MidOverIn_AS_pi0Pt_5_7 \
OutOverIn_NS_pi0Pt_5_7 \
MidOverIn_NS_pi0Pt_5_7 \
OutMinusIn_AS_pi0Pt_5_7 \
MidMinusIn_AS_pi0Pt_5_7 \
OutMinusIn_NS_pi0Pt_5_7 \
MidMinusIn_NS_pi0Pt_5_7 \
OutOverIn_AS_pi0Pt_11_14 \
MidOverIn_AS_pi0Pt_11_14 \
OutOverIn_NS_pi0Pt_11_14 \
MidOverIn_NS_pi0Pt_11_14 \
OutMinusIn_AS_pi0Pt_11_14 \
MidMinusIn_AS_pi0Pt_11_14 \
OutMinusIn_NS_pi0Pt_11_14 \
MidMinusIn_NS_pi0Pt_11_14 \
OutOverIn_AS_pi0Pt_14_17 \
MidOverIn_AS_pi0Pt_14_17 \
OutOverIn_NS_pi0Pt_14_17 \
MidOverIn_NS_pi0Pt_14_17 \
OutMinusIn_AS_pi0Pt_14_17 \
MidMinusIn_AS_pi0Pt_14_17 \
OutMinusIn_NS_pi0Pt_14_17 \
MidMinusIn_NS_pi0Pt_14_17 \
Trigger_V2 \
Trigger_V4 \
Track_V2 \
Track_V4 \
Track_V6
"


#ptBinHadronPt_pi0Pt_5_7 \
#energyLossLeadingJetPtBinEP_0_pi0Pt_5_7"

ARGSTRING=""
FILELIST=""
#for TAUI in ${TAUIS[@]} 
for INDEX in ${IND[@]}
#for TAUI in ${TAUIS[@]} 
do
	TAUI=${TAUIS[$INDEX]}
	TAUIVALUE=${TAUIVALUES[$INDEX]}
	OUTFILE=pi0Output/TauScan_NR/NoToy/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set3_Pt_NoToy_FarAve/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set3_Pt_NoToy_FarAve.pi0final.root
	ARGSTRING="$ARGSTRING $OUTFILE \"#tau_{i} = $TAUIVALUE\""
	FILELIST="$FILELIST $OUTFILE"

done

CompareDir=pi0Output/TauScan_NR/NoToy/Cent_30_50/Cmp_FarAve
mkdir -p $CompareDir
echo ./analysis/compare $CompareDir $ARGSTRING

echo ~/cern/gammaHadron/wrk/analysis/sysCompare.py -l "$LISTOFOBJECTS" -d $CompareDir -o  CmpOutput.root -f $FILELIST -t "${TAUIVALUES[0]}" "${TAUIVALUES[1]}" "${TAUIVALUES[2]}" "${TAUIVALUES[3]}" "${TAUIVALUES[4]}" "${TAUIVALUES[5]}"
~/cern/gammaHadron/wrk/analysis/sysCompare.py -l "$LISTOFOBJECTS" -d $CompareDir -o  CmpOutput.root -f $FILELIST -t "${TAUIVALUES[0]}" "${TAUIVALUES[1]}" "${TAUIVALUES[2]}" "${TAUIVALUES[3]}" "${TAUIVALUES[4]}" "${TAUIVALUES[5]}"

#~/cern/gammaHadron/wrk/analysis/sysCompare.py -l "$LISTOFOBJECTS" -d $CompareDir -o  CmpOutput.root -f $FILELIST -t "${TAUIVALUES[0]}" "${TAUIVALUES[2]}" "${TAUIVALUES[5]}"

#~/cern/gammaHadron/wrk/analysis/sysCompare.py -l "$LISTOFOBJECTS" -d $CompareDir -o  CmpOutput.root -f $FILELIST -t "${TAUIVALUES[0]}" "${TAUIVALUES[1]}" "${TAUIVALUES[2]}" "${TAUIVALUES[3]}" "${TAUIVALUES[4]}" "${TAUIVALUES[5]}"


#~/cern/gammaHadron/wrk/analysis/sysCompare.py -l "$LISTOFOBJECTS" -d $CompareDir -o  CmpOutput.root -f $FILELIST -t "${TAUIVALUES[0]}" "${TAUIVALUES[5]}"


