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

<!-- html table generated in R 3.3.3 by xtable 1.8-2 package -->
<!-- Tue May 16 16:48:27 2017 -->
<table border=1>
<tr> <th>  </th> <th> SubjID </th> <th> thickness </th> <th> method </th> <th> label_abbrev </th>  </tr>
  <tr> <td align="right"> 6 </td> <td> Caltech_0051456 </td> <td align="right"> 1.77 </td> <td> FS53 </td> <td> L.CN </td> </tr>
  <tr> <td align="right"> 17 </td> <td> Caltech_0051456 </td> <td align="right"> 0.88 </td> <td> ANTS </td> <td> L.CN </td> </tr>
  <tr> <td align="right"> 20 </td> <td> Caltech_0051456 </td> <td align="right"> 1.94 </td> <td> FS51 </td> <td> L.CN </td> </tr>
  <tr> <td align="right"> 219 </td> <td> Caltech_0051457 </td> <td align="right"> 1.49 </td> <td> ANTS </td> <td> L.CN </td> </tr>
  <tr> <td align="right"> 306 </td> <td> Caltech_0051457 </td> <td align="right"> 2.06 </td> <td> FS51 </td> <td> L.CN </td> </tr>
  <tr> <td align="right"> 364 </td> <td> Caltech_0051457 </td> <td align="right"> 2.36 </td> <td> FS53 </td> <td> L.CN </td> </tr>
   </table>

