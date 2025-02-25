* Meta-analysis
* O2 Target Data

//Mac Version
//cd "/Users/reblocke/Box Sync/Residency Personal Files/MSCI/MDCRC 6200 SRMA/O2 TargetMeta-anlysis"

//PC Version
cd "/Users/blocke/Box Sync/Residency Personal Files/Stats/ASA Baseline Use MA/"

//set scheme s1mono
set scheme cleanplots

clear
import excel "ASA BL MA.xls", sheet("Sheet1") firstrow case(lower)

meta esize int_out int_no_out pla_out pla_no_out, random(reml) studylabel(study) esize(lnrratio)

meta summarize, eform predinterval(95)
meta forestplot, random(reml) eform nullrefline(lcolor(gs3) favorsright("Favors Placebo") favorsleft("Favors ASA")) note("") esrefline(lcolor(gs3) lpattern(dash_dot)) title("Incidence Ratio of Primary Outcmoe") 
graph export "ASA BL SRMA.png", as(png) name("Graph") replace

