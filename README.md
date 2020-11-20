# Temperature-DEBIPM
R code for generating the size- and temperature depenendent DEB-IPM for co-authors. 

** Info repository content **

https://www.slu.se/fishinfoodwebs

* To run and plot the model, download the content of the repository (green button, right and up). Unzip and store the two R-scripts (Thunell_DEBIPM_model_20201120 and Thunell_DEBIPM_plot_20201120) in your working directory. When scripts are open in R, first change the path in the beginning of each script to your working directory.

* The Thunell_DEBIPM_model_20201120 contains all functions and parameter values for running the model while Thunell_DEBIPM_plot_20201120 contains code for plotting the DEB, the IPM vital rate functions and the DEBIPM output (Kappa and temperature dependent long term population growth rate lambda, the population stable structure and reproductive output).

* Run both scripts (first model and then plot) or just the plot script which sources Thunell_DEBIPM_model_20201120.

* Running the model for many Kappa and temperature values takes a while so I included a somewhat more high resolution plot of Lambda vs Kappa for many temperature (kappaTlamdba_20201120.pdf). 

* Anna and Magnus, you dont have access to the Windermere data, you'll find plots  of how the DEB performs given current parameter values in DEB-Windplots_20201120.



