*! 1.3.0 GML Sep 20 2000
*! Semi-parametric ROC estimator
*!   with bootstrapped se's for binormal parameters and auc
* 1.1: addition of partial() option
* 1.2: addition of ccsamp and roc(t) options.
* 1.3: a) change in Fdb(y) calculation: now simply #Yi_db >= y / n_db
*      b) addition of NObstrap option  (don't do bootstrap estimation of se's)
prog def dfroc, rclass
    version 6.0
    syntax varlist(min=2 max=2) [if] [in], [ noGraph noBStrap Nsamp(integer 50)  /*
    */  CCsamp Level(integer $S_level) Link(string) PARtial(real 0) ROCt(real 0)  /*
    */  RESfile(string) REPLACE  * ]
    tokenize `varlist'
    local y "`1'"
    local d "`2'"
    marksample  touse
    markout `touse' `y' `d'
    local link = lower("`link'")
    if "`link'" ~= "" {
        if "`link'"~="logit" & "`link'"~="probit" {
            di in red "link option must be logit or probit"
            exit 198
        }
    }
    else { local link "probit" }
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
    if `partial'~=0 {
        cap assert `partial' >=1.0 & `partial' <=100.0
        if _rc~=0 {
            di in red "fp% argument for partial() option must be between 1 & 100"
            exit
        }
    }
    /* so that roc() and partial() arguments are on same scale */
    if `roct'~=0 {
        cap assert `roct' >=1.0 & `roct' < 100.0
        if _rc~=0 {
            di in red "fp% argument for roc() option must be between 1 & 100"
            exit
        }
        local rocp = `roct'
        local roct = `roct'/100.0
    }
    if `level'<10 | `level'>99 {
        local level 95
    }
    if "`resfile'"=="" {tempfile resfile}
    else {local ressave "yes" }

    tempname Cov a_obs b_obs auc_obs ase bse pfile

    local roctitl : var label `y'
    if `"`roctitl'"'==`""' { local roctitl `"`y'"' }

    preserve
    qui keep if `touse'

    if `"`bstrap'"'=="" {
        if "`replace'"!="" { local replacm ",`replace'" }
        postfile `pfile' a b using `resfile' `replacm'

        /* first call to mpest prog and post with observed data */

        tempfile obsfile
        qui save `obsfile'

        if "`ccsamp'" ~= "" {
            tempfile caseobs cassamp contobs
            qui {
                keep if `d'==1 & `touse'
                save `caseobs'
                use `obsfile'
                keep if `d'==0 & `touse'
                save `contobs'
                use `obsfile'
            }
        }
    }

    /* estimate and post results for observed data */

    mpest `y' `d' `link' `partial'

    if `"`bstrap'"'~= ""  {    /* i.e. no bs se's */
        scalar `a_obs' = r(a)
        scalar `b_obs' = r(b)
        if `roct' ~= 0 {
            tempname roctobs
            scalar `roctobs' = normprob(`a_obs' + `b_obs'*(invnorm(`roct')))
        }
        if `partial'==0 {
            scalar `auc_obs' = normprob(`a_obs'/sqrt(1+`b_obs'^2))
            local aucmac = `auc_obs'
            return scalar auc = `auc_obs'
        }
    }
    else  {    /* i.e. if getting boostrapped se's */

        post `pfile' (r(a)) (r(b))

        local i = 1
        while `i' <= `nsamp' {
            qui {
                if "`ccsamp'" == "" {
                    use `obsfile',clear
                    bsample `clustop'
                }
                else {
                    use `caseobs',clear
                    bsample
                    save `cassamp'
                    use `contobs',clear
                    bsample
                    append using `cassamp'
                    erase `cassamp'
                }
                mpest `y' `d' `link' `partial'
                post `pfile' (r(a)) (r(b))
            }
            local i = `i' + 1
        }
        postclose `pfile'
        qui use `resfile',clear
        local ax = a[1]
        local bx = b[1]
        char a[bstrap] `ax'
        char b[bstrap] `bx'
        scalar `a_obs' = a[1]
        scalar `b_obs' = b[1]
        if `roct' ~= 0 {
            gen roc_t = normprob(a + b*(invnorm(`roct')))
            local rx = roc_t[1]
            char roc_t[bstrap] `rx'
            tempname roctobs
            scalar `roctobs' = roc_t[1]
        }
        if `partial'==0 {
            scalar `auc_obs' = normprob(a[1]/sqrt(1+b[1]^2))
            local aucmac = `auc_obs'
            return scalar auc = `auc_obs'
        }

        qui drop in 1   /* first record contains results from observed data */

        if `partial'==0 {
            qui gen double auc = normprob(a/sqrt(1+b^2))
            char auc[bstrap] `aucmac'
            qui sum auc
            return scalar auc_se = r(sd)
        }

        qui matrix accum `Cov' = a b, deviations noconstant
        matrix `Cov' = `Cov'/(r(N)-1)
        qui sum a
        scalar `ase' = r(sd)
        qui sum b
        scalar `bse' = r(sd)
        local nb = r(N)
        if `roct' ~= 0 {
            tempname  rocse
            qui sum roc_t
            scalar `rocse' = r(sd)
        }
        return scalar b_se = `bse'
        return scalar a_se = `ase'
    }  /* end of bootstrap specific section */

    di
    di in gr "Distribution-free ROC estimator for " in yel "`roctitl', `link'" in gr " model"

    if `partial'==0 {
        di in gr "binormal ROC parameters and corresponding AUC estimate"
        di
    }
    else {
        di in gr "  binormal ROC parameters, estimated for fpr < " in yel `partial' "%"
        di
    }
    if `roct' ~= 0 {
        di in gr "  parametric (binormal) ROC estimate at t = " in yel `rocp' "%"
        di
    }
    if `"`bstrap'"' == "" {
        if "`ccsamp'" == "" {
            di in g "  with bootstrap standard error estimates
            di in g "   (sampling w/o respect to case/control status)"
            di
        }
        else {
            di in g "  bootstrap standard error estimates based on
            di in g "   sampling separately from cases and controls"
            di
        }

        bstat, level(`level')

        di
        if `nb'~=`nsamp' {
            di "Var estimates based on " `nb' " of " `nsamp' " samples with valid estimates for b"
            di
        }
        di
        di in gr "         Var(a,b)"
        mat l `Cov',nohead
        di
        return matrix Cov `Cov'
        if "`ressave'" ~= ""{
            erase `resfile'.dta
            qui save `resfile'
        }
    }
    else {
        di in gr "      a: " in yel %6.0g `a_obs'
        di in gr "      b: " in yel %6.0g `b_obs'
        if `partial'==0 { di in gr "    AUC: " in yel %6.0g `auc_obs' }
        if `roct' ~= 0  { di in gr " ROC(t): " in yel %6.0g `roctobs' }
    }
    if `roct' ~= 0 {
        return scalar t = `roct'
        return scalar roc = `roctobs'
        if `"`bstrap'"' == "" { return scalar roc_se = `rocse' }
    }
    return scalar b = `b_obs'
    return scalar a = `a_obs'

    if `"`graph'"' == `""' {
        local a = `a_obs'
        local b = `b_obs'

        if `"`bstrap'"' ~= "" {
            local se_auc = -1
            local nsamp = -1
        }

        if `partial'==0 {
            local auc = `auc_obs'
            if `"`bstrap'"' == "" { local se_auc = return(auc_se) }
        }
        else {
            local auc = -1
            local se_auc = -1
        }
        local partprp = `partial'/100.

        Graphit , a(`a') b(`b') auc(`auc') se_auc(`se_auc') /*
           */   nsamp(`nsamp') tstlb(`roctitl') link(`link')/*
           */   partial(`partprp') `options'
    }
    restore
end
****************************************************************
program def mpest, rclass
    /* get moderately parametric ROC estimates
       estimates are returned in
          r(a)
          r(b)
    */
    version 6.0
    args y d link partial
    keep `y' `d'

    tempname nd ndb
    tempvar tdbnp uid idd uddb b1np yd sumdb

    qui sum `y' if `d'==1, meanonly
    sca `nd' = r(N)            /* nd: # with disease */
    qui sum `y' if `d'==0, meanonly
    sca `ndb' = r(N)         /* ndb: # without disease */

    /* STEP #1: CALCULATE nonparametric FDB(Y) */

    gsort -`y'
    qui gen int `sumdb' = sum(`d'==0)
    sort `y' `sumdb'
    qui by `y': gen `tdbnp' = `sumdb'[_N]
    qui replace `tdbnp' = `tdbnp'/`ndb'
    drop `sumdb'

    if `partial'~=0 {
        qui drop if (`tdbnp'> (`partial'/100.)) & d==0
        * and recalc ndb:
        qui sum `y' if `d'==0, meanonly
        sca `ndb' = r(N)
    }

    /* STEP #2 : RESTRUCTURE DATA TO RECORDS WITH PAIRS (YD,YDB) */

    qui gen `uid' = _n  /* unique id for each record */
    sort `d' `y'
    qui by `d' :gen `idd' = _n if `d'==1  /* idd unique id each `d'  */

    qui expand `nd' if `d'==0
    sort `uid'
    qui by `uid':replace `idd'= _n if `d'==0
     /* assigns a set of d's to each db */

    local ydb "`y'"
    gsort `idd' -`d'
    qui by `idd': gen `yd' = `ydb'[1]
    qui drop if `d'==1
    drop `uid' `idd' `d'

    /* STEP #3 : GENERATE BINARY INDCATORS ETC FROM PAIRS */
    qui gen byte `uddb' = (`yd' >= `ydb')
    qui gen `b1np'=invnorm(`tdbnp') if `tdbnp'<1

    /* STEP 4 : CALULATE THE Moderately PARAMETRIC ROC */

    qui `link' `uddb' `b1np'
    return scalar a = _b[_cons]
    cap return scalar b = _b[`b1np']
    if _rc~=0 {return scalar b = .}
end
******************************************************************
prog def Graphit
    version 6.0
    syntax , A(real) B(real) AUC(real) SE_auc(real) NSamp(int) /*
      */  Tstlb(string) Link(string) Partial(real) [*]
    drop _all
    tempvar t roct forty5
    qui set obs 200
    qui range `t' 0 1
    qui gen `forty5' = `t' if _n==1 | _n==_N
    qui gen `roct' = normprob(`a'+`b'*invnorm(`t'))
    qui replace `roct' = 0 if `t'==0  /* otherwise missing/undefined */
    qui replace `roct' = 1 if `t'==1

    local gridpts `".20,.40,.60,.80"'
    local gridlb   `"0,.20,.40,.60,.80,1"'

    local 0 `",`options'"'
    syntax [, Symbol(string) XLAbel(string) YLAbel(string) XLIne(string) /*
        */    YLIne(string) L2title(string) B2title(string) L1title(string) /*
        */    T1title(string) T2title(string) Pen(string) Connect(string) *]

    if `partial'~=0 {
        tempvar expi prtline
        qui gen byte `expi' = `t'==0
        qui expand 3 if `expi'
        sort `expi'
        qui by `expi': replace `t' = `partial' if `expi'==1 & _n>1
        qui by `expi': gen     `prtline' = 0 if `expi'==1 & _n==2
        qui by `expi': replace `prtline' = 1 if `expi'==1 & _n==3
        qui by `expi': replace `roct'   = . if `expi'==1 & _n>1
        qui by `expi': replace `forty5' = . if `expi'==1 & _n>1

        local yvars "`roct' `forty5' `prtline'"
        if `"`symbol'"' == `""' { local symbol `"iii"' }
        if `"`connect'"' == `""' { local connect `"sll"' }
        if `"`pen'"' == `""' { local pen `"239"' }
        if `"`t2title'"'  == `""' {  /*
        */ local t2title: di "based on model estimated for t < " %3.2f `partial' " only"
        }
    }
    else {
        local yvars "`roct' `forty5'"
        if `"`symbol'"' == `""' { local symbol `"ii"' }
        if `"`connect'"' == `""' { local connect `"sl"' }
        if `"`pen'"' == `""' { local pen `"23"' }
        if `"`t2title'"'  == `""' {
            if `se_auc' == -1 {
                local t2title: di "AUC = " %4.3f `auc' "
            }
            else {
                local t2title: di "AUC = " %4.3f `auc' ";  boostrap se(AUC) = " %4.3f `se_auc'  /*
                    */   "  (`nsamp'  samples)"
            }
        }
    }
    if `"`xlabel'"' == `""' { local xlabel `"`gridlb'"' }
    if `"`ylabel'"' == `""' { local ylabel `"`gridlb'"' }
    if `"`xtick'"'  == `""' { local xtick  `"`gridlb'"' }
    if `"`ytick'"'  == `""' { local ytick  `"`gridlb'"' }
    if `"`xline'"'  == `""' { local xline  `"`gridpts'"' }
    if `"`yline'"'  == `""' { local yline  `"`gridpts'"' }
    if `"`l2title'"'  == `""' { local l2title  `"ROC(t)"' }
    if `"`l1title'"'  == `""' { local l1title  `" "'}
    if `"`b2title'"'  == `""' { local b2title  `"t"' }
    if `"`t1title'"'  == `""' { /*
    */ local t1title  `"Distribution-free ROC estimator for `tstlb';  `link' model"' }

    sort `t'
    #delimit ;
    graph `yvars' `t',
        c(`connect') s(`symbol') border
        t1(`"`t1title'"') t2(`"`t2title'"')
        xlabel(`xlabel') ylabel(`ylabel')
        xline(`xline') yline(`yline')
        pen(`pen')
        l2(`l2title') b2(`b2title') l1(`"`l1title'"')
        `options' ;
    #delimit cr
end
