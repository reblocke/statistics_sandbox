*! 1.3.3 GML Mar 2 2001
* supercedes ccaucbs
* Empirical AUC estimate and bootstrapped se
* v 1.1 adds cohort sampling (default, case/control sampling now
*   specified as option), cluster sampling option, point estimate
*   for empirical ROC, (version 1.0 was ccaucbs)
* v 1.2 addition of option for 2nd test variable and difference statistics
* v 1.3 ROC argument scale change (% from 1-100) for consistency with partial()
* v 1.3.2 error in subprog roc fixed: d to `d' in -count- cmd
* v 1.3.3 adjustment for partial AUC calculation added for deviation of
*         specified t from corresponding discrete quantile of y_db
prog def aucbs, rclass
    version 7.0
    syntax varlist(min=2 max=3) [if] [in], [ PARtial(real 0) ROC_e(real 0) /*
    */ CCsamp  Level(integer $S_level) Nsamp(integer 50) CLuster(varlist) /*
    */ RESfile(string) REPLACE  * ]
    local ntest: word count `varlist'
    local ntest = `ntest' - 1
    tokenize `varlist'
    local y1 "`1'"
    if `ntest'==1 {
        local d "`2'"
    }
    else {
        local y2 "`2'"
        local d "`3'"
    }
    marksample touse
    markout `touse' `y1' `y2' `d'

    qui ta `d' if `touse'
    if r(r)~=2 {
        di in red "`d' must take on two values"
        exit 198
    }
    qui sum `d' if `touse',meanonly
    if r(min)~=0 | r(max)~=1 {
        di in red "`d' must be 0/1"
        exit 198
    }
    if `level'<10 | `level'>99 {
        local level 95
    }
    if `partial'~=0 {
        cap assert `partial' >=1.0 & `partial' <=100.0
        if _rc~=0 {
            di in red "argument for partial() option must be between 1 & 100"
            exit
        }
    }
    /* so that roc() and partial() arguments are on same scale */
    if `roc_e'~=0 {
        cap assert `roc_e' >=1.0 & `roc_e' < 100.0
        if _rc~=0 {
            di in red "fp% argument for roc_e() option must be between 1 & 100"
            exit
        }
    }

    if "`cluster'" ~="" {
        confirm variable `cluster'
        local clustop `",cl(`cluster')"'
    }
    tempname pfile roc roc1 roc2 auc1 auc2

    if "`resfile'"=="" {tempfile resfile}
    else {local ressave "yes" }

    if "`replace'"!="" { local replacm ",`replace'" }

    if `ntest'==1 {
        if `roc_e'~=0 {postfile `pfile' auc roc_t using `resfile' `replacm'}
        else          {postfile `pfile' auc using `resfile' `replacm'}
        local testvar : var label `y1'
        if `"`testvar'"'==`""' { local testvar `"`y1'"' }
    }
    else {
        if `roc_e'~=0 {postfile `pfile' auc1 auc2 roc1_t roc2_t using `resfile' `replacm'}
        else          {postfile `pfile' auc1 auc2 using `resfile' `replacm'}
        local tv1 : var label `y1'
        if `"`tv1'"'==`""' { local tv1 `"`y1'"' }
        local tv2 : var label `y2'
        if `"`tv2'"'==`""' { local tv2 `"`y2'"' }
        local testvar `"`tv1' and `tv2'"'
        if `"`testvar'"'==`""' { local testvar `"`y1' and `y2'"' }
    }
    preserve
    keep `y1' `y2' `d' `touse' `cluster'

    if `roc_e'~=0 {
        roc `y1' `d' `touse' `roc_e'
        if `ntest'==1 { scalar `roc' = r(roc) }
        else {
            scalar `roc1'= r(roc)
            roc `y2' `d' `touse' `roc_e'
            scalar `roc2'= r(roc)
        }
    }

    auc `y1' `d' `touse' `partial'
    if `ntest'==1 {
        if `roc_e'~=0 {post `pfile' (r(auc)) (`roc')}
        else          {post `pfile' (r(auc)) }
    }
    else {
        scalar `auc1' = r(auc)
        auc `y2' `d' `touse' `partial'
        scalar `auc2' = r(auc)
        if `roc_e'~=0 {post `pfile' (`auc1') (`auc2') (`roc1') (`roc2')}
        else          {post `pfile' (`auc1') (`auc2') }
    }

    if "`ccsamp'" == "" {
        tempfile obsfile
        qui save `obsfile'
    }
    else {
        tempfile caseobs cassamp contobs
        qui keep if `d'==1 & `touse'
        qui save `caseobs'
        restore, preserve
        qui keep if `d'==0 & `touse'
        qui save `contobs'
    }

    local i = 1
    while `i' <= `nsamp' {
        qui {
            if "`ccsamp'" == "" {
                use `obsfile',clear
                bsample `clustop'
            }
            else {
                use `caseobs',clear
                bsample `clustop'
                save `cassamp'
                use `contobs',clear
                bsample `clustop'
                append using `cassamp'
                erase `cassamp'
            }

            if `roc_e'~=0 {
                roc `y1' `d' `touse' `roc_e'
                if `ntest'==1 { scalar `roc' = r(roc) }
                else {
                    scalar `roc1'= r(roc)
                    roc `y2' `d' `touse' `roc_e'
                    scalar `roc2'= r(roc)
                }
            }

            auc `y1' `d' `touse' `partial'
            if `ntest'==1 {
                if `roc_e'~=0 {post `pfile' (r(auc)) (`roc')}
                else          {post `pfile' (r(auc)) }
            }
            else {
                scalar `auc1' = r(auc)
                auc `y2' `d' `touse' `partial'
                scalar `auc2' = r(auc)
                if `roc_e'~=0 {post `pfile' (`auc1') (`auc2') (`roc1') (`roc2')}
                else          {post `pfile' (`auc1') (`auc2') }
            }

        }
        local i = `i' + 1
    }
    postclose `pfile'
    use `resfile',clear

    if `ntest'==1 {
        local x = auc[1]       /* first obs is from observed data */
        char auc[bstrap] `x'
        if `roc_e' ~= 0 {
            local x = roc_t[1]
            char roc_t[bstrap] `x'
        }
    }
    else {
        gen aucdelta = auc2 - auc1
        local x1 = auc1[1]
        local x2 = auc2[1]
        local x3 = aucdelta[1]
        char auc1[bstrap] `x1'
        char auc2[bstrap] `x2'
        char aucdelta[bstrap] `x3'
        if `roc_e'~=0 {
            gen rocdelta = roc2_t - roc1_t
            local x1 = roc1_t[1]
            local x2 = roc2_t[1]
            local x3 = rocdelta[1]
            char roc1_t[bstrap] `x1'
            char roc2_t[bstrap] `x2'
            char rocdelta[bstrap] `x3'
        }
    }
    qui drop in 1
    if "`ressave'" ~= ""{
        erase `resfile'.dta
        qui save `resfile'
    }

    di
    if `ntest'==1 {
        di in g "Empirical AUC estimate for " in yel "`testvar'"
        if `roc_e' ~= 0 {
            di in g "  and empirical ROC(t) at t = " in yel `roc_e' "%"
        }
        di
    }
    else {
        di in g "Empirical AUC estimates for"
        di in g "    test1: " in yel "`tv1'"
        di in g "    test2: " in yel "`tv2'"
        if `roc_e' ~= 0 {
            di in g "  empirical ROC(t) at t = " in yel `roc_e' "%"
            di
            di in g "  AUC and ROC(t) differences for test 2 - test 1 (aucdelta & rocdelta)"
            di
        }
        else {
            di
            di in g "  and AUC difference for test 2 - test 1 (aucdelta)"
            di in g
        }
    }

***********************************************

    if "`ccsamp'" == "" {
        di in g "  with bootstrap standard error estimates
        di in g "   (sampling w/o respect to case/control status)"
    }
    else {
        di in g "  bootstrap standard error estimates based on
        di in g "   sampling separately from cases and controls"
    }
    if `partial'~=0 {
        di
        di in g "  auc is Partial AUC for FPR < " in yel "`partial'%"
    }
    bstat, level(`level')
*****************************************
* returned results:

    if `ntest'==1 {
        qui bstat auc
        if `partial'==0 {
            return scalar se_auc = r(se)
            return scalar auc = r(stat)
        }
        else {
            return scalar pauc_t = `partial'
            return scalar se_pauc = r(se)
            return scalar pauc = r(stat)
        }
        if `roc_e'~=0 {
            qui bstat roc_t
            return scalar t = `roc_e'/100.0
            return scalar se_roct = r(se)
            return scalar roct = r(stat)
        }
    }
    else /* ** 2 tests ** */ {
        qui bstat auc1
        if `partial'==0 {
            return scalar se_auc1 = r(se)
            return scalar auc1 = r(stat)
        }
        else {
            return scalar pauc_t = `partial'
            return scalar se_pauc1 = r(se)
            return scalar pauc1 = r(stat)
        }
        qui bstat auc2
        if `partial'==0 {
            return scalar se_auc2 = r(se)
            return scalar auc2 = r(stat)
        }
        else {
            return scalar se_pauc2 = r(se)
            return scalar pauc2 = r(stat)
        }
        qui bstat aucdelta
            return scalar se_aucdl = r(se)
            return scalar aucdelta = r(stat)
        if `roc_e'~=0 {
            qui bstat roc1_t
            return scalar t = `roc_e'/100.0
            return scalar se_roc1t = r(se)
            return scalar roc1t = r(stat)
            qui bstat roc2_t
            return scalar se_roc2t = r(se)
            return scalar roc2t = r(stat)
            qui bstat rocdelta
            return scalar se_rocdl = r(se)
            return scalar rocdelta = r(stat)
        }
    }
end
*********************************************************
prog def auc,rclass
    /* calculate auc or partial auc */
    args y d touse partial
    tempvar yt
    tempname ydmax nd ndb tpr_t fpr_t y_perc auc_adj
    qui gen `yt' = `y'
    if `partial'~=0 {
        local perc = 100.0 - `partial'
        qui sum `yt' if `d'==1 & `touse',meanonly
        sca `ydmax' = r(max)
        sca `nd' = r(N)
        qui count if `d'==0 & `touse'
        scalar `ndb' = r(N)
        qui _pctile `yt' if `d'==0 & `touse',p(`perc')
        scalar `y_perc' = r(r1)
        qui count if `d'==1 & `touse' & `yt' >= `y_perc'
        sca `tpr_t' = r(N)/`nd'
        qui count if `d'==0 & `touse' & `yt' >= `y_perc'
        sca `fpr_t' = r(N)/`ndb'
        /* partial area adjustment - for deviation of discrete non-dis test quantile
           from the specified partial t : */
        sca `auc_adj' = `tpr_t' * (`fpr_t' - (`partial'/100.0))
        qui replace `yt' = `ydmax' + 1.0 if `d'==0 & `yt'< `y_perc'
    }
    qui ranksum `yt' if `touse',by(`d')
    return scalar auc = 1.0 - (r(sum_obs) - r(N_1)*(r(N_1) + 1)/2.0) / (r(N_1) * r(N_2))
    if `partial'~=0 {
        return scalar auc = return(auc) - `auc_adj'
    }
end
***********************************************************
prog def roc,rclass
    /* calculate empirical roc for specified t
        add check on range of t (0<t<1)
    */
    args y d touse t
    local perc = 100 - `t'
    qui _pctile `y' if `d'==0 & `touse', p(`perc')
    qui count if `y' > r(r1) & `d'==1 & `touse'
    local num = r(N)
    qui count if `d'==1 & `touse'
    return scalar roc = `num'/r(N)
end
