* Meta-analysis
* O2 Target Data

//Mac Version
cd "/Users/blocke/Box Sync/Residency Personal Files/Scholarly Work/Locke Research Projects/statistics_sandbox/Obesity NIV Reintubation"

//PC Version
//cd "C:\Users\reblo\Box\Residency Personal Files\Stats\Heterotopic Ossification"



//set scheme s1mono
set scheme cleanplots

clear
import excel "ObesNIV MA.xls", sheet("Sheet1") firstrow case(lower)

//drop if study == "Hernandez"

meta esize int_ho int_no_ho pla_ho pla_no_ho, random(reml) studylabel(study) esize(lnrratio)

meta summarize, eform predinterval(95)
meta forestplot, random(reml) eform nullrefline(lcolor(gs3) favorsright("Favors Placebo") favorsleft("Favors Intervention")) note("") esrefline(lcolor(gs3) lpattern(dash_dot)) title("Reintubation Rate") 
graph export "ObesNIV SRMA.png", as(png) name("Graph") replace

