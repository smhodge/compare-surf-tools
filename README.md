# compare_surf_tools
A BrainHack Project to compare thickness outputs from different pipelines run on ABIDEI

# To prepare the data and run the example analyses:

$ cd compare-surf-tools
# Start R
$ R
R> source("analysis/data_prep.R")
R> source("analysis/pairs_plot.R")
R> source("analysis/heatmap_plot.R")


The main data file is ~analysis/abide_ct.RData, which is a row-wise collection of cortical thicknesses for each study ID by method (ANTS, Freesurfer v5.1, and Freesurfer v5.3) by cortical unit. (i.e. a "long" data format). For example, to see the first 6 rows of data for the "left cuneus" label:

R> head(abide[abide[, "Label Name"] %in% "ctx-lh-cuneus",c(1:3, 10)])
             SubjID thickness method label_abbrev
6   Caltech_0051456   1.76900   FS53         L.CN
17  Caltech_0051456   0.88290   ANTS         L.CN
20  Caltech_0051456   1.94100   FS51         L.CN
219 Caltech_0051457   1.49429   ANTS         L.CN
306 Caltech_0051457   2.05900   FS51         L.CN
364 Caltech_0051457   2.35700   FS53         L.CN
 
