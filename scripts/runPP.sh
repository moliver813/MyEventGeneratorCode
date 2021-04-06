#!/bin/bash

#TAUIS=(02 04 06 08 10 12)

# Far Eta Average

INFILE=pi0Root/pp/NoToy/Prod_07_11_pp_S_5023_Pt_NoToy_.pi0hist.root
OUTFILE=pp_S_5023_Set4_Pt_NoToy_FarAve.pi0final.root
OUTDIR=pi0Output/pp/NoToy/pp_S_5023_Set4_Pt_FarAve/
./analysis/phase2 -i $INFILE -o $OUTFILE -od $OUTDIR --background2D 3 -B 0

# ZYAM

INFILE=pi0Root/pp/NoToy/Prod_07_11_pp_S_5023_Pt_NoToy_.pi0hist.root
OUTFILE=pp_S_5023_Set4_Pt_NoToy_ZYAM.pi0final.root
OUTDIR=pi0Output/pp/NoToy/pp_S_5023_Set4_Pt_ZYAM/
./analysis/phase2 -i $INFILE -o $OUTFILE -od $OUTDIR -B ZYAM


