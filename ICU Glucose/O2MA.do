* Meta-analysis
* O2 Target Data

//Mac Version
//cd "/Users/reblocke/Box Sync/Residency Personal Files/MSCI/MDCRC 6200 SRMA/O2 TargetMeta-anlysis"

//PC Version
cd "C:\Users\reblo\Box\Residency Personal Files\MSCI\MDCRC 6200 SRMA\O2 TargetMeta-anlysis\"

//set scheme s1mono
set scheme cleanplots

clear
import excel "O2 icu target MA.xls", sheet("Sheet1") firstrow case(lower)
keep if !missing(name) // drop lines that don't correspond to studies


*For subgroups
/*
gen scr_contrast=0
replace scr_contrast=1 if baseline_scr=="lessthan2"
replace scr_contrast =2 if baseline_scr=="morethan2" & contrast_volume=="morethan100"
replace scr_contrast =3 if baseline_scr=="morethan2" & contrast_volume=="lessthan100"
label define scr_contrast 1 "scr<2" 2 "scr>2_contrast>100" 3 "scr>2_contrast<100"
label values scr_contrast scr_contrast
*/

* For meta-regression 
/*
gen level_Scr = 0
replace level_Scr = 1 if baseline_scr =="lessthan2"
replace level_Scr = 2 if baseline_scr =="morethan2"
*/

meta esize lowo2died lowo2alive higho2died higho2alive, random(reml) studylabel(name) 

* [ ] need to fix color -> dashes of lines. 
* random(sjonkman) option is good if sample size is very small.; random(dl) otherwise ok. reml is default which is generally ok
meta summarize, eform subgroup (outcome) predinterval(95)
meta forestplot, random(reml) eform nullrefline(lcolor(gs3) favorsright("Favors High O2") favorsleft("Favors Low O2")) note("") esrefline(lcolor(gs3) lpattern(dash_dot)) subgroup ( outcome ) title("Odds Ratio of Mortality, Lower vs Higher  O2 Targets") 
graph export "O2 Target SRMA.png", as(png) name("Graph") replace


//currently not working because ldbounds is not imported as "bounds" in the R environment... seems like two options are two=

//play around with metacumbounds for TSA 
//windows version
metacumbounds lowo2died lowo2alive higho2died higho2alive, data(count) rrr(0.05) alpha(.05) beta(.20) effect(r) spending(1) is(ais) id(author) stat(rr) listRout listRin rdir("C:\Program Files\R\R-4.3.1\bin\") shwRRR graph xtitle(Information Size) ytitle(z-score)

/*
* Other subgruops
* baseline S.Cr >2 & used contrast volume >100
meta summarize, eform subgroup ( scr_contrast )
meta forestplot, eform nullrefline(lcolor(gs3)) esrefline(lcolor(gs3) lpattern(dash_dot)) subgroup ( scr_contrast ) title("Risk Ratio, Random-Effects Meta-analysis")
*/

//left over code
/*
* risk of bias
meta summarize, eform subgroup ( rob )
meta forestplot, eform  nullrefline(lcolor(gs3)) esrefline(lcolor(gs3) lpattern(dash_dot)) subgroup ( rob ) title("Risk Ratio, Random-Effects Meta-analysis")

* meta-regression 
meta regress level_Scr
* Key factors: I2 (amount of variation) and P-value (is this a predictor)
* These are more important than beta-coefficients per se

* [ ] need to fix colors
meta funnelplot, random(dl) contours(1 5 10)


meta trimfill, estimator(run) funnel (contours(1 5 10)) eform
meta bias, egger

* log-linear dose-response regression models: glst command

*L'Abbe plot
labbe nac_cin nac_nocin pla_cin pla_nocin, null rr(0.72) id(author)

* OR
clear
import excel "NAC_CIN Binary_Data for MA_workshop.xls", sheet("Sheet1") firstrow case(lower)
list age author baseline_scr contrast_volume hydration_type nac_cin nac_nocin pla_cin pla_nocin rob year

gen odds_nac= nac_cin/nac_nocin
gen odds_pla= pla_cin /pla_nocin
gen or= odds_nac/ odds_pla

meta esize nac_cin nac_nocin pla_cin pla_nocin, random(dl) esize(lnor) studylabel(author year)
meta summarize, eform subgroup ( baseline_scr)
meta forestplot, eform  nullrefline(lcolor(gs3)) esrefline(lcolor(gs3) lpattern(dash_dot))  subgroup ( baseline_scr) title("Odds Ratio, Random-Effects Meta-analysis")

* RD
clear
import excel "NAC_CIN Binary_Data for MA_workshop.xls", sheet("Sheet1") firstrow case(lower)
list age author baseline_scr contrast_volume hydration_type nac_cin nac_nocin pla_cin pla_nocin rob year

meta esize nac_cin nac_nocin pla_cin pla_nocin, random(dl) esize(rdiff) studylabel(author year)
meta summarize, subgroup ( baseline_scr)
meta forestplot,  nullrefline(lcolor(gs3)) esrefline(lcolor(gs3) lpattern(dash_dot))  subgroup ( baseline_scr) title("Risk Difference Random-Effects Meta-analysis")



*meta esize = uses 4 variables  t_event, t_noevent, c_event, c_noevent
*meta set = 3 variables RR and 95 CI
*meta set = 3 set lnrr, 


* Continuous data demo

clear
import excel "Pharmacist_continous data_MA_workshop.xlsx", sheet("Continous") firstrow case(lower)
meta esize n1 m1 sd1 n2 m2 sd2, eslabel(MD) random(dl) studylabel (authoryear) esize( mdiff )
meta summarize, predinterval(95)
meta forestplot

* HKSJ method - most accurate in small sample size (because normality assumption of DL is problematic)


*/

*MD if have same variables

*SMD if different variable scales
