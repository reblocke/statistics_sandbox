*! 1.1.2 GML 07 Sept 2005
* originally written for MP Diagnostic Methods course
* 1.0.1 03 March 2000
* 1.0.2:  fixed d -> `d' references
* 1.0.3: trap for logit error with zero fpr
* 1.0.4: 29 June 2001
*        add'l returned results, display format option,
*        macro length error fix pre stata 7 for DLR+ joint interval (with DLR-)
* 1.1.0: 07 Oct 2003
*        remove dependence on -propci-
*        replace Wald interval (for proportions) with Wilson interval
* 1.1.1: GML 21 Jan 2004
*        include conf. bounds for tpr, fpr, npv, fpv in returned results
* 1.1.2: GML 07 Sept 2005
*        help converted to smcl
*        returned joint intervals for PPV, NPV renamed for consistency
********
* consider including option for Agresti Coull Interval?
prog def binscrn1, rclass
    version 7.0
    syntax varlist(min=2 max=2) [if] [in], [ Level(integer $S_level)   /*
    */     Format(string) ]
    tempvar tpr fpr dn dbn dnsum dbnsum
    tokenize `varlist'
    local y "`1'"
    local d "`2'"
    marksample touse
    markout `touse' `y' `d'

    /*  check that disease and test vars are 0/1 only: */

    qui ta `d' if `touse'
    if r(r)~=2 {
        di in red "disease indicator `d' must take on two values"
        exit 198
    }
    qui sum `d' if `touse',meanonly
    if r(min)~=0 | r(max)~=1 {
        di in red "disease indicator `d' must be 0/1"
        exit 198
    }
    /*  error check for y: 0/1 only */
    qui ta `y' if `touse'
    if r(r)~=2 {
        di in red "test result variable `y' must take on two values"
        exit 198
    }
    qui sum `y' if `touse',meanonly
    if r(min)~=0 | r(max)~=1 {
        di in red "test result variable `y' must be 0/1"
        exit 198
    }
    if `level'<10 | `level'>99 {
        local level 95
    }
    if `"`format'"'==`""' {
        local format "%4.2g"
    }
    else {
        cap qui di `format' _pi
        if _rc~=0 | index("`format'","%")~=1 {
            di in red "format(#) option is invalid"
            exit 120
        }
    }

    di
    di in g "Summary measures for binary screening tests
    di
    di in g "Disease status by test outcome:
    di in g "-------------------------------

    ta `y' `d' if `touse' ,

    ******************************************************************
    *  calculate confidence intervals:

    local ylab : var label `y'
    tempvar not_d
    qui gen byte `not_d' = (`d'==0)
    di
    di in g "screening measures:
    di in g "-------------------
    di
    di in g "  test measure `y': `ylab'"
    di
    di in g "                                |------------    `level'% CI    -------------|
    di in g "                                     exact       Wilson        logit
    di
    cicalc "True Positive rate:   " `y' "`d'==1" `level' `touse' `format'

    return sca tp_lb_ex = r(lb_ex)
    return sca tp_ub_ex = r(ub_ex)
    return sca tp_lb_wi = r(lb_wi)
    return sca tp_ub_wi = r(ub_wi)
    return sca tp_lb_lg = r(lci)
    return sca tp_ub_lg = r(uci)

    cicalc "False Positive rate:  " `y' "`d'==0" `level' `touse' `format'
    return sca fp_lb_ex = r(lb_ex)
    return sca fp_ub_ex = r(ub_ex)
    return sca fp_lb_wi = r(lb_wi)
    return sca fp_ub_wi = r(ub_wi)
    return sca fp_lb_lg = r(lci)
    return sca fp_ub_lg = r(uci)

    di
    cicalc "Disease Prevalence:   " `d' "`d'~=." `level' `touse' `format'
    return sca disprev = r(p)

    cicalc "Pos. Predictive Value:" `d' "`y'==1" `level' `touse' `format'
    return sca ppv = r(p)
    return sca ppv_lb_ex = r(lb_ex)
    return sca ppv_ub_ex = r(ub_ex)
    return sca ppv_lb_wi = r(lb_wi)
    return sca ppv_ub_wi = r(ub_wi)
    return sca ppv_lb_lg = r(lci)
    return sca ppv_ub_lg = r(uci)

    cicalc "Neg. Predictive Value:" `not_d' "`y'==0" `level' `touse' `format'
    return sca npv = r(p)
    return sca npv_lb_ex = r(lb_ex)
    return sca npv_ub_ex = r(ub_ex)
    return sca npv_lb_wi = r(lb_wi)
    return sca npv_ub_wi = r(ub_wi)
    return sca npv_lb_lg = r(lci)
    return sca npv_ub_lg = r(uci)

    ************************************************************
    *  calculate various accuracy measures:

    tempname tp fp dlrpos dlrneg z se
    tempname dlrp_lb dlrp_ub dlrn_lb dlrn_ub

    qui count if `d'==1 & `touse'
    local nd = r(N)              /* number with disease     */
    qui count if `d'==0 & `touse'
    local ndb = r(N)             /* number without disease  */

    qui count if `d'==1 & `y'==1 & `touse'
    scalar  `tp' = r(N)/`nd'     /* true positive proportion */

    qui count if `d'==0 & `y'==1 & `touse'
    scalar `fp' = r(N)/`ndb'     /* false positive proportion */

    scalar `dlrpos' = `tp'/`fp'
    scalar `dlrneg' = (1.0 - `tp')/(1.0 - `fp')

    ***  CI's for DLR+ and DLR-  *******************
    ***     (using Simel's formula)

    local alpha = 1.0 - (`level'/100.)
    sca `z' = invnorm(1.0 - (`alpha'/2))

    * for DLR+ :
    scalar `se' = sqrt((1 - `tp')/(`tp'*`nd') + (1-`fp')/(`fp'*`ndb'))
    scalar `dlrp_lb' = exp(ln(`dlrpos') - `z'*`se')
    scalar `dlrp_ub' = exp(ln(`dlrpos') + `z'*`se')
    return sca selndlrp = `se'

    * for DLR- :
    scalar `se' = sqrt((`tp')/((1-`tp')*`nd') + (`fp')/((1-`fp')*`ndb'))
    scalar `dlrn_lb' = exp(ln(`dlrneg') - `z'*`se')
    scalar `dlrn_ub' = exp(ln(`dlrneg') + `z'*`se')
    return sca selndlrn = `se'

    /* display CI's for DLR+ and DLR-  */
    di
    di
    di in g "Diagnostic Likelihood Ratios (`level'% Simel CI's) : "
    di in g "-----------------------------------------------
    di
    di in g "    DLR+:  " in y `format' `dlrpos' "      (" `format' `dlrp_lb' "," `format' `dlrp_ub' ")  "
    di in g "    DLR-:  " in y `format' `dlrneg' "      (" `format' `dlrn_lb' "," `format' `dlrn_ub' ")  "
    di

    **************************************************************************
    /* iv) Joint confidence intervals:

      rectangles based on sqrt(1-alpha) level logit based confidence limits for
      the component measures (TP, FP, PPV, NPV)

      Confidence regions for DLR+ and DLR- taken directly from the above calculated
        95% confidence limits
    */

    tempname tp_lb tp_ub fp_lb fp_ub
    tempname ppv_lb ppv_ub npv_lb npv_ub alpha_2
    tempname dlrp_l2 dlrp_u2 dlrn_l2 dlrn_u2

    scalar `alpha_2' = 1 - sqrt(1-`alpha')
    local level_2 = 100.*(1.0 - `alpha_2')

    /* first get logit based limits calculated by cicalc */
    qui cicalc "TP" `y' "`d'==1" `level_2' `touse'
    scalar `tp_lb' = r(lci)
    scalar `tp_ub' = r(uci)
    qui cicalc "FP" `y' "`d'==0" `level_2' `touse'
    scalar `fp_lb' = r(lci)
    scalar `fp_ub' = r(uci)

    return sca tpj_lb = `tp_lb'
    return sca tpj_ub = `tp_ub'
    return sca fpj_lb = `fp_lb'
    return sca fpj_ub = `fp_ub'

    qui cicalc "PPV" `d' "`y'==1" `level_2' `touse'
    scalar `ppv_lb' = r(lci)
    scalar `ppv_ub' = r(uci)
    qui cicalc "NPV" `not_d' "`y'==0" `level_2' `touse'
    scalar `npv_lb' = r(lci)
    scalar `npv_ub' = r(uci)

    return sca ppvj_lb = `ppv_lb'
    return sca ppvj_ub = `ppv_ub'
    return sca npvj_lb = `npv_lb'
    return sca npvj_ub = `npv_ub'

    *** calculate joint DLR limits
    *    using CI's for TP and FP  (sqrt(1-alpha) level CI's rather than 1-alpha):

    * for DLR+ :
    scalar `dlrp_l2' = `tp_lb'/`fp_ub'
    scalar `dlrp_u2' = `tp_ub'/`fp_lb'
    return sca dlrpj_ub = `dlrp_u2'
    return sca dlrpj_lb = `dlrp_l2'

    * for DLR- :
    scalar `dlrn_l2' = (1-`tp_ub')/(1-`fp_lb')
    scalar `dlrn_u2' = (1-`tp_lb')/(1-`fp_ub')
    return sca dlrnj_ub = `dlrn_u2'
    return sca dlrnj_lb = `dlrn_l2'

   di
   di in g "Joint Confidence Regions (`level'%): "
   di in g "-------------------------------
   di
   di in g "   (TP, FP)   =>   TP(l,u) x FP(l,u)  : " /*
      */ in y "(" `format' `tp_lb' "," `format' `tp_ub'                /*
      */  ")" in g" x " in y "(" `format' `fp_lb' "," `format' `fp_ub' ")
   di in g "   (PPV,NPV)  =>  PPV(l,u) x NPV(l,u) : " /*
      */ in y  "(" `format' `ppv_lb' "," `format' `ppv_ub'                /*
      */  ")" in g" x " in y "(" `format' `npv_lb' "," `format' `npv_ub' ")
   di in g "   (DLR+,DLR-)                        : " /*
      */ in y "(" `format' `dlrp_l2' "," `format' `dlrp_u2'                /*
      */  ")" in g" x " in y "(" `format' `dlrn_l2' "," `format' `dlrn_u2' ")
   di

   return sca tp = `tp'
   return sca fp = `fp'
   return sca dlrneg = `dlrneg'
   return sca dlrpos = `dlrpos'
   return sca nd = `nd'
   return sca ndb = `ndb'
end
***************************************************
program define cicalc, rclass
    /* calculate and return logit based confidence limits,
       obtain Wald and Exact limits, and display intervals for
       all three methods.
       (for accuracy measures which are simple proportions)

       arguments:
        mname -  string with label for accuracy measure (proportion)
        x -  binomial indicator variable defining the proportion
        include - condition defining denominator subset for proportion
        * alpha - for 1-alpha level confidence bounds
        level
        touse
    */
    args mname x include level touse format
    tempname z lci uci lb_exct ub_exct lb_wilson ub_wilson
    local alpha = 1.0 - (`level'/100.)
    scalar `z' = invnorm(1.0 - (`alpha'/2))

    /* CI based on logit model */

    cap logit `x' if `include' & `touse'
    if _rc~=0 {
        sca `lci' = .
        sca `uci' = .
        return sca lci = .
        return sca uci = .
        }
    else {
        sca `lci' = _b[_cons] - `z'* _se[_cons]
        sca `lci' = exp(`lci')/(exp(`lci') + 1.0)
        sca `uci' = _b[_cons] + `z'* _se[_cons]
        sca `uci' = exp(`uci')/(exp(`uci') + 1.0)
        return sca lci = `lci'  /* allow program to return confidence bounds for later use */
        return sca uci = `uci'
    }
    /* obtain Exact (Clopper-Pearson) and Wilson CI's */

    ciexact `x' if `include' & `touse', level(`level')
    sca `lb_exct'    = r(lb)
    return sca lb_ex = r(lb)
    sca `ub_exct'    = r(ub)
    return sca ub_ex = r(ub)

    ciwilson `x' if `include' & `touse', level(`level')
    sca `lb_wilson'  = r(lb)
    return sca lb_wi = r(lb)
    sca `ub_wilson'  = r(ub)
    return sca ub_wi = r(ub)

    di in g "  `mname'  " in y `format' r(p_hat) "    ("`format' `lb_exct' ","`format' `ub_exct' ")" /*
    */                           "  ("`format' `lb_wilson' ","`format' `ub_wilson'  ")" /*
    */                           "  ("`format' `lci' "," `format' `uci'  ")"
    return sca p = r(p_hat)
end
*********************************************************
pro def ciexact, rclass
    * calculate and return exact confidence limits for a binomial proportion
    * note vers. 7 Stata's -ci- will not permit non-integer 'level'
    syntax varname if, Level(real)
    assert inlist(`varlist',0,1,.) /* remove after trial run */
    qui sum `varlist' `if', meanonly
    local n = r(N)
    local k = int(r(mean)*`n' + .5)
    return scalar lb = invbinomial(`n', `k',(100-`level')/200)
    return scalar ub = invbinomial(`n', `k',1-(100-`level')/200)
end
*********************************************************
pro def ciwilson, rclass
    * calculate and return Wilson confidence limits for a binomial proportion
    * note version 7 Stata's -ci- will not permit non-integer 'level'
    syntax varname if, Level(real)
    tempname kappa phat
    local alpha = 1.0 - (`level'/100.)
    sca `kappa' = invnorm(1 - (`alpha'/2))
    qui sum `varlist' `if', meanonly
    local n = r(N)
    local k = int(r(mean)*`n' + .5)
    sca `phat' = `k'/`n'
    return scalar lb = (`k' + ((`kappa'^2)/2))/(`n' + `kappa'^2) /*
        */              - ((`kappa'*sqrt(`n'))/(`n' + `kappa'^2)) * sqrt(`phat'*(1-`phat') + (`kappa'^2)/(4*`n'))
    return scalar ub = (`k' + ((`kappa'^2)/2))/(`n' + `kappa'^2) /*
        */              + ((`kappa'*sqrt(`n'))/(`n' + `kappa'^2)) * sqrt(`phat'*(1-`phat') + (`kappa'^2)/(4*`n'))
    return scalar p_hat = `phat'
end
