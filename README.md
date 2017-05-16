# compare_surf_tools
A BrainHack Project to compare thickness outputs from different pipelines run on ABIDE I

# To prepare the data and run the example analyses:

$ cd compare-surf-tools  
# Start R  
$ R  
R> source("analysis/data_prep.R")  
R> source("analysis/pairs_plot.R")  
R> source("analysis/heatmap_plot.R")  


The main data file is ~analysis/abide_ct.RData, which is a row-wise collection of cortical thicknesses for each study ID by method (ANTS, Freesurfer v5.1, and Freesurfer v5.3) by cortical unit. (i.e. a "long" data format). For example, to see the first 6 rows of data for the "left cuneus" label:

R> head(abide[abide[, "Label Name"] %in% "ctx-lh-cuneus",c(1:3, 10)])
 tex table generated in R 3.3.3 by xtable 1.8-2 package

% Tue May 16 16:44:06 2017
\begin{table}[ht]
\centering
\begin{tabular}{rlrll}
  \hline
 & SubjID & thickness & method & label\_abbrev \\ 
  \hline
6 & Caltech\_0051456 & 1.77 & FS53 & L.CN \\ 
  17 & Caltech\_0051456 & 0.88 & ANTS & L.CN \\ 
  20 & Caltech\_0051456 & 1.94 & FS51 & L.CN \\ 
  219 & Caltech\_0051457 & 1.49 & ANTS & L.CN \\ 
  306 & Caltech\_0051457 & 2.06 & FS51 & L.CN \\ 
  364 & Caltech\_0051457 & 2.36 & FS53 & L.CN \\ 
   \hline
\end{tabular}
\end{table}

