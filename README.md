# scRNAseq-yMSCsuppressedTcells
This repository contains code used to analyze data in "Interferon-γ activated mesenchymal stromal cell prophylaxis of GVHD by a T cell Nrf2-dependent mechanism" by Foppiani et al, 2026.

There are three scRNAseq datasets included in the paper. Each dataset has its own code used to analyze and generate figures:

96-hour, bead activated T cells (activated and yMSC suppressed) by cell cycle: 96hour_BeadActivated_Tcells_Analysis.R

96-hour, transwell (activated and bystander T cells, control and yMSC suppressed): 96hour_Transwell_Analysis.R

24-hour, bead activated T cells (control and yMSC suppressed): 24h_BeadActivated_Analysis.R

Note: For all analyses, filtered feature BC matrices were used to generate Seurat objects. These folders are also included in this repository.
