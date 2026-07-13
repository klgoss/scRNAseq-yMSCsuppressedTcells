# scRNAseq-yMSCsuppressedTcells
This repository contains code used to analyze data in "Interferon-γ activated mesenchymal stromal cell prophylaxis of GVHD by a T cell Nrf2-dependent mechanism" by Foppiani et al, 2026.

There are three scRNAseq datasets included in the paper. Each dataset has it's own script of code used to analyze and generate figures:

96-hour, bead activated T cells (activated and yMSC suppressed) by cell cycle: 
- Script: 96hour_BeadActivated_Tcells_Analysis.R

96-hour, transwell (activated and bystander T cells, control and yMSC suppressed):
- Script: 96hour_Transwell_Analysis.R

24-hour, bead activated T cells (control and yMSC suppressed):
- Script: 24h_BeadActivated_Analysis.R


For all analyses, filtered feature BC matrices were used to generate Seurat objects. These folders are also included in this repository.
