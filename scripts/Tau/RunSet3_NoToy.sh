#!/bin/bash

TAUIS=(02 04 06 08 10 12)


# Far Eta Average

for TAUI in ${TAUIS[@]} 
do
	INFILE=pi0Root/TauScan_NR/NoToy/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set3_Pt_NoToy_.pi0hist.root
	OUTFILE=TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set3_Pt_NoToy_FarAve.pi0final.root
	OUTDIR=pi0Output/TauScan_NR/NoToy/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set3_Pt_NoToy_FarAve/
	./analysis/phase2 -i $INFILE -o $OUTFILE -od $OUTDIR --background2D 3 -B 0
done

# exit 0

# ZYAM

for TAUI in ${TAUIS[@]} 
do
	INFILE=pi0Root/TauScan_NR/NoToy/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set3_Pt_NoToy_.pi0hist.root
	OUTFILE=TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set3_Pt_NoToy_ZYAM.pi0final.root
	OUTDIR=pi0Output/TauScan_NR/NoToy/TauScan_PbPb_NR_Tau_$TAUI''_Cent_30_50_S_5023_Set3_Pt_NoToy_ZYAM/
	./analysis/phase2 -i $INFILE -o $OUTFILE -od $OUTDIR -B ZYAM
done

