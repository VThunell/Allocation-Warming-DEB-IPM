# Temperature-DEBIPM
R code for generating the size- and temperature depenendent DEB-IPM for co-authors. 

## Info on current repository content

* To run and plot the model, download the content of the repository (green button that says Code, right and up). Unzip and store the two R-scripts (Thunell_DEBIPM_model_20201120 and Thunell_DEBIPM_plot_20201120) in your working directory. When scripts are open in R, first change the path in the beginning of each script to your working directory.

* The Thunell_DEBIPM_model_20201120 contains all functions and parameter values for running the model while Thunell_DEBIPM_plot_20201120 contains code for plotting the DEB, the IPM vital rate functions and the DEBIPM output (Kappa and temperature dependent long term population growth rate lambda, the population stable structure and reproductive output).

* Run both scripts (first model and then plot) or just the plot script which sources Thunell_DEBIPM_model_20201120.

* Running the model for many Kappa and temperature values takes a while so I included a somewhat more high resolution plot of Lambda vs Kappa for many temperature (kappaTlamdba_20201120.pdf). The temepratue dependence of growth There is still some work left on 

* Anna and Magnus, you dont have access to the Windermere data, you'll find plots of how the DEB performs in relaion to data on Lake Windemere pike given current parameter values in DEB-Windplots_20201120.pdf.

* Im working on a document to give a bit of background on my assumptions regarding vital rate functions and parameter values.

* You might consider the plot of fecundity of the DEB and Windermere pike. Given the kappa-scaling that we use, this underestimates reproductive output of large but not small individuals. The only way that I can think of to fix this is to include an allometric relationship of Kappa as we need to increase allocable energy for reproduction but not growth for large individuals.

* I have used the survival of eggs to scale the long-term population growth rate to reasonable values (to slightly above one).