* Meta-analysis
* O2 Target Data

//Mac Version
//cd "/Users/reblocke/Box Sync/Residency Personal Files/MSCI/MDCRC 6200 SRMA/O2 TargetMeta-anlysis"

//PC Version
cd "C:\Users\reblo\Box\Residency Personal Files\Stats\Heterotopic Ossification"

//set scheme s1mono
set scheme cleanplots

clear
import excel "HO MA.xls", sheet("Sheet1") firstrow case(lower)

meta esize int_ho int_no_ho pla_ho pla_no_ho, random(reml) studylabel(study) esize(lnrratio)

meta summarize, eform predinterval(95)
meta forestplot, random(reml) eform nullrefline(lcolor(gs3) favorsright("Favors Placebo") favorsleft("Favors Intervention")) note("") esrefline(lcolor(gs3) lpattern(dash_dot)) title("Incidence Ratio of Heterotopic Ossification") 
graph export "HO SRMA.png", as(png) name("Graph") replace

