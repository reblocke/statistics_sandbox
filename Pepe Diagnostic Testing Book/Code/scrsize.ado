*! 1.0.4 GML 15 Jan 2003
*
* 1.0.4 : eliminate dependence on user program -propci-;
*         and changes for vers 6 to 7 (forvalues)
* 1.0.3 GML 13 Dec 2000
pro def scrsize, rclass
    version 7
    preserve
    drop _all
    gettoken tpr_1 0 : 0, parse(" ,")
    gettoken fpr_1 0 : 0, parse(" ,")
    local dalpha = 1 - $S_level/100
    syntax , ND(int) NDB(int) TPNull(real) FPNull(real) [Alpha(real `dalpha') /*
      */     NSim(int 500) Dots noDISplay ]
    cap confirm number `tpr_1'
     if _rc~=0 {
         di in red "tp_rate argument needs to be a number"
         error _rc
     }
    cap confirm number `fpr_1'
     if _rc~=0 {
         di in red "fp_rate argument needs to be a number"
         error _rc
     }
    if `tpr_1' > 1 | `tpr_1' < 0 {
        di in red "tp_rate must be between 0 & 1"
        exit 198
    }
    if `fpr_1' > 1 | `fpr_1' < 0 {
        di in red "fp_rate must be between 0 & 1"
        exit 198
    }
    if `tpnull' > 1 | `tpnull' < 0 {
        di in red "argument to tpnull(#) must be between 0 & 1"
        exit 198
    }
    if `fpnull' > 1 | `fpnull' < 0 {
        di in red "argument to fpnull(#) must be between 0 & 1"
        exit 198
    }
    if `alpha' > 1 | `alpha' < 0 {
        di in red "argument to alpha(#) must be between 0 & 1"
        exit 198
    }

    tempfile resfile
    simres `nsim' `resfile' `alpha' `nd' `ndb' `tpr_1' `fpr_1' `dots'

    use `resfile',clear
    qui gen byte reject = (tprlb > `tpnull') & (fprub < `fpnull')
    qui sum reject,meanonly
    local power = r(sum)/_N

    if `"`display'"' == `""' {
       di
       di
       di in gr "         Null   True "
       di in gr "         ----   ---- " /*
         */    in gr "                 nd:" in y %5.0f `nd'
       di in gr "   TPR   " in y %4.2f `tpnull' "   " %4.2f `tpr_1' /*
         */    in gr "                 ndb:" in y %5.0f `ndb'
       di in gr "   FPR   " in y %4.2f `fpnull' "   " %4.2f `fpr_1' /*
         */    in gr "       # simulations:" in y %5.0f `nsim'
       di
       di in gr "         alpha: " in y %5.2g `alpha'  in g "      1-sided test"
       di
       di in gr "   ==>   power: " in y %5.3g `power'
       di
    }
    return scalar power = `power'
end
*************************************************************
program define simres
    args nsim resfile alpha nd ndb tpr fpr dots
    tempname alpha_2 mysim

    qui postfile `mysim' tprlb tprub fprlb fprub using `resfile', replace
    local n = `nd' + `ndb'
    qui set obs `n'
    qui gen byte d = (_n<=`nd')
    qui gen byte y = .

    /* mult x 2 since using 1-sided test */
    scalar `alpha_2' = 2 * (1 - sqrt(1-`alpha'))

    local level = 100.*(1.0 - `alpha_2')

    local dots = cond("`dots'"=="", "*", "noisily")
    forvalues i = 1/`nsim' {
       `dots' di in gr "." _c
       doasim `tpr' `fpr' `level'
       post `mysim' (r(tpr_lb)) (r(tpr_ub)) (r(fpr_lb)) (r(fpr_ub))
    }
    postclose `mysim'
end
********************************************************
pro def doasim, rclass
    args tpr fpr level

    qui replace y = (uniform() < `tpr') if d==1
    qui replace y = (uniform() < `fpr') if d==0

    ciexact y if d==1, level(`level')
    return scalar tpr_lb = r(lb)
    return scalar tpr_ub = r(ub)

    ciexact y if d==0, level(`level')
    return scalar fpr_lb = r(lb)
    return scalar fpr_ub = r(ub)
end
*********************************************************
pro def ciexact, rclass
    * calculate and return exact confidence limits for a binomial proportion
    * note Stata's -ci- will not permit non-integer 'level'
    syntax varname if, Level(real)
    assert inlist(`varlist',0,1,.) /* remove after trial run */
    qui sum `varlist' `if', meanonly
    local n = r(N)
    local k = int(r(mean)*`n' + .5)
    return scalar lb = invbinomial(`n', `k',(100-`level')/200)
    return scalar ub = invbinomial(`n', `k',1-(100-`level')/200)
end
