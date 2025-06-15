* Steroid PNa
* Meta-analysis

cd "/Users/blocke/Box Sync/Residency Personal Files/Scholarly Work/Locke Research Projects/statistics_sandbox/Steroids PNa/"
set scheme cleanplots

clear

import excel "Steroids PNa MA.xls", sheet("Sheet1") firstrow case(lower)
list study year chest_ma int_death int_alive pla_death pla_alive
label variable chest_ma "Old Trials"
label define chest_ma_label 0 "Old Trials" 1 "NEW"
label values chest_ma chest_ma_label
label variable study "Study"
//encode study, gen(study_n)
//label variable study_n "Study"

meta esize int_death int_alive pla_death pla_alive, random(reml) esize(lnor) studylabel(study year) 
meta summarize, eform predinterval(95)
meta forestplot, eform  nullrefline(lcolor(gs3)) esrefline(lcolor(gs3) lpattern(dash_dot)) subgroup (chest_ma) title("Steroids OR for Death, Random-Effects Meta-analysis")  xsize(11) ysize(7) xscale(log range(0.01 64)) xlabel(0.03125 "1/32" 0.0625 "1/16" 0.125 "1/8" 0.25 "1/4" 0.5 "1/2" 1 "1" 2 "2" 4 "4" 8 "8" 16 "16" 32 "32") scheme(tab2)
//95% prediction interval for exp(theta): [ 0.54, 1.320]



