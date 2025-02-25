/* PFT Simulator */ 

/* AIM: to identify the extent to which PRISM "evoluation" is just a function of measurement error */ 


// Clear the existing dataset
clear

// Create variable 'a' with values from 0 to 100
set obs 125
gen fvc = _n
label variable fvc "FVC Percent Predicted"

// Expand the dataset for each value of 'b'
expand 125

// Generate variable 'b'
bysort fvc: gen fev1 = _n
label variable fev1 "FEV1 Percent Predicted"

// Sort the dataset by 'a' and 'b'
sort fvc fev1


gen interp = ""
replace interp = "normal" if (fev1/fvc >= 0.7) & (fev1 >= 80) & (fvc >= 80)
replace interp = "prism" if (fev1/fvc >= 0.7) & (fev1 < 80)
replace interp = "restriction suggested" if fvc < 80
replace interp = "obstruction" if fev1/fvc < 0.7
replace interp = "" if (fev1/fvc >= 1.4)
//replace interp = "" if (fvc < 0.3) //also probably 
label define interp_lab 0 "normal" 1 "prism" 2 "restriction suggested" 3 "obstruction"
encode interp, generate(interpretation) label(interp_lab)
label variable interpretation "Spirometry Interpretation"
drop interp

heatplot interpretation fev1 fvc, xbins(125) ybins(125) legend(off) color(HCL qualitative) ///
    text(75 45 "Restriction" "Suggested") ///
	text(30 80 "Obstruction") ///
	text(100 100 "Normal") ///
	text(75 90 "PRISM")

	//Ref std's
	
	
//TODO: could you simulate test-retest reliability here? to quantify stability. 

//TODO: why not model fev1, fvc, and an interaction term (to represent that a small fev1 in the context of a large fvc - ie obstruction- is something different than a small fev1 in the context of a small fvc - ie restriction)

/*
ATS/ERS guidance says within-day differences
in FEV1  5%, week-to-week differences 11%, and
year-to-year differences  15% for normal subjects
should not be interpreted as clinically meaningful;
greater differences apply to chronic obstructive
pulmonary disease (or COPD) patients (within-day
 13%; week-to-week  20%) (Pellegrino, Viegi
et al., 2005).
*/ 
