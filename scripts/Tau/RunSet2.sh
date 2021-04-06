#!/bin/bash

TAUIS=(02 04 06 08 10 12)

for TAUI in ${TAUIS[@]} 
do
	INFILE=pi0Root/TauScan_NR/NoFlow/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set2_Pt_NoFlow_.pi0hist.root
	OUTFILE=TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set2_Pt_NoFlow_ZYAM.pi0final.root
	OUTDIR=pi0Output/TauScan_NR/NoFlow/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set2_Pt_NoFlow_ZYAM/
	./analysis/phase2 -i $INFILE -o $OUTFILE -od $OUTDIR -B ZYAM
done

exit 0

for TAUI in ${TAUIS[@]} 
do
	INFILE=pi0Root/TauScan_NR/NoFlow/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set2_Pt_NoFlow_.pi0hist.root
	OUTFILE=TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set2_Pt_NoFlow_ZYAM.pi0final.root
	OUTDIR=pi0Output/TauScan_NR/NoFlow/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set2_Pt_NoFlow_ZYAM/
	./analysis/phase2 -i $INFILE -o $OUTFILE -od $OUTDIR -B ZYAM
done



