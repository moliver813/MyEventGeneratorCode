#!/bin/bash

IND=(0 1 2 3 4 5)
TAUIS=(02 04 06 08 10 12)
TAUIVALUES=(0.2 0.4 0.6 0.8 1.0 1.2)

ARGSTRING=""
#for TAUI in ${TAUIS[@]} 
for INDEX in ${IND[@]}
#for TAUI in ${TAUIS[@]} 
do
	TAUI=${TAUIS[$INDEX]}
	TAUIVALUE=${TAUIVALUES[$INDEX]}
	OUTFILE=pi0Output/TauScan_NR/NoFlow/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set1_Pt_NoToy_ZYAM/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set1_Pt_NoToy_ZYAM.pi0final.root
	ARGSTRING="$ARGSTRING $OUTFILE \"#tau_{i} = $TAUIVALUE\""
	#OUTFILE=pi0Output/TauScan_NR/NoFlow/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set1_Pt_NoToy_ZYAM.pi0final.root
#	OUTDIR=pi0Output/TauScan_NR/NoFlow/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set1_Pt_NoToy_ZYAM/
#	./analysis/phase2 -i $INFILE -o $OUTFILE -od $OUTDIR -B ZYAM
done

#echo "ARGSTRING = $ARGSTRING"
CompareDir=pi0Output/TauScan_NR/NoFlow/Cmp_ZYAM
mkdir -p $CompareDir
echo ./analysis/compare $CompareDir $ARGSTRING



