
set scheme cleanplots
graph set window fontface "Helvetica"
* Regression on likelihood of initial prescription
clear

cd "/Users/blocke/Box Sync/Residency Personal Files/Stats/RCT Duplicate"

program define datetime 
end
capture mkdir "Results and Figures"
capture mkdir "Results and Figures/$S_DATE/" //make new folder for figure output if needed
capture mkdir "Results and Figures/$S_DATE/Logs/" //new folder for stata logs
local a1=substr(c(current_time),1,2)
local a2=substr(c(current_time),4,2)
local a3=substr(c(current_time),7,2)
local b = "RCT Duplicate Analysis.do" // do file name
copy "`b'" "Results and Figures/$S_DATE/Logs/(`a1'_`a2'_`a3')`b'"


import excel "table 1 abstraction.xlsx", sheet("Sheet1") firstrow case(lower)

gen calc_se = (log(upper)-log(lower))/(3.92)
gen calc_est = log(point)
gen z = abs((calc_est)/calc_se)
gen p = exp((-0.717*z) - (0.416 * z^2) )
recode p (0/0.05 = 1) (0.05/max=0) (.=0), gen(orig_sig)

//was the replication p<0.05?
gen rep_sig = .
replace rep_sig = 1 if (orig_sig == 1 & replicated == 1)
replace rep_sig = 0 if (orig_sig == 1 & replicated == 0)
replace rep_sig = 0 if (orig_sig == 0 & replicated == 1)
replace rep_sig = 1 if (orig_sig == 0 & replicated == 0)

//what was the predicted likelihood of replication?
gen pred_rep = .
replace pred_rep = pred_p05 if orig_sig == 1
replace pred_rep = 1-pred_p05 if orig_sig == 0



//Visualizations

scatter pred_rep p if orig_sig == 1, mlabel(trialname) mlabsize(5pt) graphregion(margin(small)) mlabgap(0) mlabcolor(%75) jitter(10) mlabangle(-45) xscale(log extend) yscale(range(0.4, 1) extend)
scatter pred_rep p if orig_sig == 0, mlabel(trialname) mlabsize(5pt) graphregion(margin(small)) mlabgap(0) mlabcolor(%75) jitter(10) mlabangle(-45) yscale(range(0.4, 1))


logistic replicated pred_rep
estat gof, group(10)

//does not seem like the ones that were more likely not to replicate didn't?


//method - simulation
