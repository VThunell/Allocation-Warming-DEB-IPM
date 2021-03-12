# Temperature-DEBIPM
R code for generating the size- and temperature depenendent DEB-IPM for co-authors. 

## Info on current repository content

Despite the Branch name Cost-of-reproduction, this branch contains the new code for including a Size dependent kappa into the DEBIPM. For switching between size dependent and independent kappa, see function kap_fun(). When using the size indepedent kappa, the size scaling paramters of intake, should be reset (to eps1=0.64 and eps2=0.55).

To run and plot the model, download the content of the repository (green button that says Code, right and up). 

The Thunell_DEBIPM_model_2021xxxx contains all functions and parameter values for running the model while Thunell_DEBIPM_plot_2021xxxx contains code for plotting the DEB, the IPM vital rate functions and the DEBIPM output (Kappa and temperature dependent long term population growth rate lambda, the population stable structure and reproductive output).
