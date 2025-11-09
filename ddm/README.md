To run the scripts contained in these folders, a couple of preparations are necessary and some information might be helpful.

1) R (https://www.r-project.org/) and R Studio (https://www.rstudio.com/) need to be preinstalled to run the scripts. 
   
2) The C++ boost libraries (https://www.boost.org/) need to be downloaded and saved under "C:Program Files" (we used the version 1_76_0) to run the scripts containing the code for the Drift Diffusion Models.

3) The data file "ddmdata.csv" contains the raw data and can be found in the "Data and Analyses" folder. 
   The raw data contains the subject id ("id"), the trial number ("trialNo"), whether it is the first or second decision in the respective trial ("decNo"), the reaction time in seconds ("rt"), the choice outcome ("choice"; 0 = avoid, 1 = engage), the trial type ("TT3"; 1 = self trials, 2 = partner trials, 3 = group trials), the self performance ("S"), the partner performance ("P"), the performance of the relevant opponent ("Or"), the performance of the irrelevant performance ("Oi"), the performance of the first opponent ("O1"), the performance of the second opponent ("O2"), and the bonus ("bonus").

4) The scripts containing the code to run the time-varying Drift Diffusion Model (the model first presented by Maier et al. (2020), adapted for the purpose of our study) can be found in the folder "tDDM".
   
5) The parameter estimates and analyses scripts on which the results reported in the manuscript are based can be found in the folder "fits". The individual scripts should save the modelfits in there as well. If you have all the right packages, and the paths are in order, simply running mainscript, and then waiting for a long time, should get you there.