This is supplemental material for the paper entitled 'Bayesian
Functional Principal Components Analysis using Relaxed Mutually Orthogonal Processes'.
In the 3 directories are code used to replicate the results presented in the paper. 
Below, we describe the contents of the directories and discuss replication workflow. 

/cebu_data_analysis - 

 file descriptions:

  cebu_mcmc.R - This script imports data from the /data directory and runs the MCMC algorithm to 
replicate the results in Section 5 of the main paper.

  cebu_mcmc_vis.R - This script loads MCMC samples from the /results directory to replicate
all figures in section 5 of the main paper and appendix F.2, F.4, and F.5. 

  cebu_wb.R - This script imports data from the /data directory and runs the MCMC algorithm to 
implement the CG model in Appendix F.3.

  cebu_wb_model.txt - This is a WinBUGS model file that is called within cebu_wb.R

  cebu_wb_vis.R - This script loads MCMC samples from the /results directory to replicate
all figures in Appendix F.3.

  create_tent_basis.R - This script creates a tent basis used in the historical functional
regression model described in Section 5 of the main paper.

  /data, /figures, /results - These are all subdirectories that are called within the 
R files in the /cebu_data_analysis directory.

 workflow: To succesfully run the code in this directory, the directory will need to be 
the current working directory in the user's R session. 

Since the authors do not have explicit permission to share the data from the Cebu 
longitudinal health and nutrition survey, the user must download the 'mlong.dta', 
'mmort.dta', and 'mbirth2.dta' files from the Cebu Mother Baseline Survey, 
1983-1986 directory in the following url: https://dataverse.unc.edu/dataset.xhtml?persistentId=hdl:1902.29/11680.
Store these files in the /data directory.

With the data files in place, the algorithm to run the model described in Section 5
of the paper can be ran using cebu_wb.R. The results of the script will be stored in the 
/results directory. If the user would like a file that contains the authors MCMC samples,
this is available upon request. The file was too large to share on this platform. 

With saved MCMC samples in the /resuts directory, a user can remake all figures related 
to this example using the cebu_mcmc_vis.R script. The script also computes WAIC from the 
model. The seed is set for the MCMC algorithm for reproducibility. 

If a user would like to reproduce the results from the CG model (Appendix F.3), the user must
download the data and run the cebu_wb.R script. This script calls the WinBUGS model file,
cebu_wb_model.txt, and assumes that WinBUGS is installed in the directory: "c:/Program Files/WinBUGS14/".
If this is not the case, then the user must specify the correct location of WinBUGS. 
The results of the script will be stored in the /results directory. If the user would 
like a file that contains the authors MCMC samples, this is available upon request. 
The file was too large to share on this platform. 

With saved MCMC samples in the /resuts directory, a user can remake all figures related 
to this CG example (Appendix F.3) using the cebu_wb_vis.R script. The script also computes WAIC from the 
model. The seed is set for the MCMC algorithm for reproducibility. 

/simulations - 

 file descriptions:

  bfpca_simulation.R - This script replicates the simulation experiment performed in Section 
4.2 of the main paper.

  fpca_growth_simulation.R - This script replicates the simulation experiment performed in Appendix E.  

  fpca_reg_simulation.R - This script replicates the simulation experiment performed in Section 
4.3 of the main paper. 

  fpca_simulation.R - This script replicates the simulation experiment performed in Section 
4.1 of the main paper. 

  lambda_prior_study.R - This script replicates the prior illustration in Section 2.2.

  nu_sensitivity_study.R - This script replicates the simulation experiment performed in 
Appendix C. 

  wb_growth_model.txt - This is a WinBUGS model file that is called within the fpca_growth_simulation.R script.

  wb_model.txt - This is a WinBUGS model file that is called within the fpca_simulation.R script.

  wb_reg_model.txt - This is a WinBUGS model file that is called within the fpca_reg_simulation.R script.
 
 workflow: To succesfully run the code in this directory, the directory will need to be 
the current working directory in the user's R session. 

The scripts bfpca_simulation.R, fpca_reg_simulation.R, and fpca_simulation.R can all be ran to 
replicate the results in section 4.2, 4.3, and 4.1 of the main paper respecively. The script 
fpca_growth_simulation.R can be ran to replicate the results in Appendix E. Each script
begins with generating replicates of simulated datasets, fitting competing models on the same 
generated dataset, and concluding with illustrating the performance othe competing methods. 
A user can select whether they want to replicate the simulation experiment in a 
'low' or 'high' variability setting, which is discussed in Section 4 of the paper. The seed is 
set for the replication experiments for reproducibility. 


The fpca_growth_simulation.R, fpca_reg_simulation.R and fpca_simulation.R scripts call the WinBUGS 
model files: wb_growth_model.txt wb_reg_model.txt and wb_model.txt. It is assumed that WinBUGS is 
installed in the directory: "c:/Program Files/WinBUGS14/".If this is not the case, then the user 
must specify the correct location of WinBUGS. 

The lambda_prior_study.R script can be ran to recreate all figures from Section 2.2 in the 
main paper. The seed is set for the replication experiments for reproducibility. 

The nu_sensitivity_study.R script can be ran to recreate all figures from Appendix C. 
The seed is set for the replication experiments for reproducibility. 

/src - 

 file descriptions:

  msf.cpp - This c++ file has functions used to resolve label switching and sign ambiguity
in MCMC samples.

  remo_fpca.cpp - This c++ file has functions used in steps of MCMC algorithms presented in the paper. 

 workflow: These files contain functions that are used to speed up the MCMC implementation in R. 