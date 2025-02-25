*! 1.3.5 GML 25 Feb 2003
* ROC confidence bands and Var(AUC) added in version 1.1
* error in sample used for kdensity estimation fixed in vers 1.2
* optional return of kdensity estimates added in vers 1.2
* partial AUC option added in 1.3
* absence of `touse' in partial auc calc remedied in 1.3.1
* correction of se(ROC) in 1.3.2 (used for confidence bands):
*     i.e.  L = fdc/fdbc replaced with L^2 in subprog cbgen
* adjustment for partial AUC calculation added in 1.3.3
*    for deviation of specified t from corresponding discrete
*    quantile of y_db
* 1.3.4: optional return of variable with roc se based on kdensity
*    estimates included in 1.3.4
*    and slight modification to varauc subprog.
* 1.3.5: trap added for fix of graph title error when test variable label exceeds 33 characters
*        (t(1) limit is 69 characters.  graphit conscructs t(1) from test var label and 36
*        other characters)
* 1.0.0 GML 12 November 1999
* 1.1.0 GML 17 December 1999
* 1.2.0 GML 16 February 2000
* 1.3.1 GML 05 April 2000
* 1.3.2 GML 14 December 2000
* 1.3.3 GML 4 May 2001
* 1.3.4 GML 24 May 2002
*
prog def emroc, rclass
    version 7.0
    syntax varlist(min=2 max=2) [if] [in], [ noGraph Level(integer $S_level)   /*
    */  GENSPec(string) GENSEns(string) SE_roc(string) PARtial(real 0)  /*
    */  KERNal(string) KWIDth(real -1.0) CONFband GENKd(string) REPLACE  * ]
    tempvar tpr fpr dn dbn dnsum dbnsum
    tokenize `varlist'
    local y "`1'"
    local d "`2'"
    marksample touse
    markout `touse' `y' `d'

    if "`replace'"!=`""' {
        cap drop `genspec'
        cap drop `gensens'
        cap drop `genkd'
        cap drop `se_roc'
    }

    /*  need to add error check for d: 0/1 only */
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

    if `"`se_roc'"' ~= "" {
        if `"`confband'"' == "" {
            di in red " confband option must be included if se_roc() option is specified"
            exit 198
        }
    }

    if `"`genkd'"' ~= "" {
        if `"`confband'"' == "" {
            di in red " confband option must be included if genkd() option is specified"
            exit 198
        }
    }

    /* get AUC */
    auc `y' `d' `touse' 0
    return scalar auc = r(auc)

    if `partial'~=0 {
        auc `y' `d' `touse' `partial'
        return scalar pauc = r(auc)
        return scalar pauc_t = `partial'
        }

    qui count if `d'==1 & `touse'
    local dtot = r(N)
    qui count if `d'==0 & `touse'
    local dbtot = r(N)

    qui egen `dn'  = sum(`d'==1) if `touse',by(`y')
    qui egen `dbn' = sum(`d'==0) if `touse',by(`y')

    /* get variance(AUC) - old & new at this point, for comparison */
    varauc `y' `d' `touse' `dtot' `dbtot'
    return scalar aucvar = r(vauc)

    /* `touse2' is restricted, beyond `touse', to observations corresponding
       to distinct test values */

    tempvar touse2
    gen byte `touse2' = `touse'
    sort `touse2' `y'
    qui by `touse2' `y': replace `touse2' = 0 if _n>1

    /* nsum = no. >= y */

    gsort -`y'

    qui {
    gen `dnsum'= sum(`dn') if `touse2'
    gen `dbnsum' = sum(`dbn') if `touse2'

    gen `tpr' = `dnsum'/`dtot'
    gen `fpr' = `dbnsum'/`dbtot'

    sort `y' `touse2'
    by `y': replace `tpr' = `tpr'[_N]
    by `y': replace `fpr' = `fpr'[_N]
    }

    /* confidence bands */

    if `"`confband'"'~="" {
        tempvar ub lb kdd kddb se
        qui gen double `ub'=.
        qui gen double `lb'=.
        qui gen `kdd'=.
        qui gen `kddb'=.
        qui gen `se'=.
        /* kwidth: kernal half-width argument for kdensity program */

        if `kwidth' == -1 { local kwidth "" }
        else {local kwidth `"width(`kwidth')"'}
        cbgen `tpr' `fpr' `d' `level' `dtot' `dbtot' `y' `touse'  /*
           */ `se' `lb' `ub' `kdd' `kddb' `kernal' `kwidth'
        qui replace `ub' = 1.0 if `ub'>1 & `ub' ~=.
        qui replace `lb' = 0 if `lb'<0 & `lb' ~=.
        local kwidth = r(width)

        if "`genkd'"!=`""' {
            local wc : word count `genkd'
            cap assert `wc' == 2
            if _rc~=0 {
                di in red "must have 2 newvar names in genkd( ) option"
                error _rc
            }
            tokenize `genkd'
            cap confirm new variable `1'
            if _rc~=0 {
                di in red "variable `1' already exists,"
                di in red `"  use "replace" option to replace `1' "'
                error _rc
            }
            cap confirm new variable `2'
            if _rc~=0 {
                di in red "variable `2' already exists,"
                di in red `"  use "replace" option to replace `2' "'
                error _rc
            }
            rename `kdd' `1'
            rename `kddb' `2'
        }
        if "`se_roc'"!=`""' {
            local gvarnum : word count `se_roc'
            cap assert `gvarnum' == 1
            if _rc~=0 {
                di in red "more than one newvar name in se_roc( ) option"
                error _rc
            }
            cap confirm new variable `se_roc'
            if _rc~=0 {
                di in red "variable `se_roc' already exists,"
                di in red `"  use "replace" option to replace `se_roc'"'
                error _rc
            }
            qui gen `se_roc' = `se'
        }
    }

    /* restrict to unique test values for graph */

    if `"`graph'"' == `""' {
        local aucarg = return(auc)
        Graphit `tpr' `fpr' `y' `lb' `ub' if `touse2', auc(`aucarg') l(`level') kwd(`kwidth') `options'
        }

/* turn this into a subprogram later? : */

    if "`gensens'"!=`""' {
        local gvarnum : word count `gensens'
        cap assert `gvarnum' == 1
        if _rc~=0 {
            di in red "more than one newvar name in gensens( ) option"
            error _rc
        }
        cap confirm new variable `gensens'
        if _rc~=0 {
            di in red "variable `gensens' already exists,"
            di in red `"  use "replace" option to replace `gensens'"'
            error _rc
        }
        rename `tpr' `gensens'
    }
    if "`genspec'"!=`""' {
        local gvarnum : word count `genspec'
        cap assert `gvarnum' == 1
        if _rc~=0 {
            di in red "more than one newvar name in genspec( ) option"
            error _rc
        }
        cap confirm new variable `genspec'
        if _rc~=0 {
            di in red "variable `genspec' already exists,"
            di in red `"  use "replace" option to replace `genspec'"'
            error _rc
        }
        qui gen `genspec' = 1.0 - `fpr'
    }
end
*********************************************************************
prog def Graphit
    version 6.0
    syntax varlist(min=3 max=5) [if], Level(integer) AUC(real) kwd(string) [ * ]
    tempvar origin forty5
    local varnum: word count `varlist'
    tokenize `varlist'
    local tpr "`1'"
    local fpr "`2'"
    local y "`3'"
    local lb "`4'"
    local ub "`5'"

    marksample touse

    local aucstr : display %4.3f `auc'
    if `varnum' ~=5 { local kwdt ""}
    else { local kwdt : display %5.0g `kwd' }
    preserve

    /* need to create (tpr,fpr) = (0,0) for ROC curve  */
    qui {
    count
    local newobs = r(N) + 1
    set obs `newobs'
    gen byte `origin' = ( _n==_N )
    replace `tpr' = 0 if `origin'
    replace `fpr' = 0 if `origin'
    }

    local roctitl : var label `y'
    if `"`roctitl'"'==`""' |  length(`"`roctitl'"') > 32 { local roctitl `"`y'"' }

    local gridpts `".20,.40,.60,.80"'
    local gridlb   `"0,.20,.40,.60,.80,1"'

    local 0 `",`options'"'
    syntax [, Symbol(string) XLAbel(string) YLAbel(string) XLIne(string) /*
        */    YLIne(string) L2title(string) B2title(string) L1title(string) /*
        */    T1title(string) T2title(string) *]

    if `"`symbol'"' == `""' { local symbol `"o"' }
    if `"`xlabel'"' == `""' { local xlabel `"`gridlb'"' }
    if `"`ylabel'"' == `""' { local ylabel `"`gridlb'"' }
    if `"`xline'"'  == `""' { local xline  `"`gridpts'"' }
    if `"`yline'"'  == `""' { local yline  `"`gridpts'"' }

    if `"`l2title'"'  == `""' { local l2title  `"True Positive rate"' }
    if `"`l1title'"'  == `""' { local l1title  `" "'}
    if `"`b2title'"'  == `""' { local b2title  `"False Positive rate"' }
    if `"`t1title'"'  == `""' { local t1title  `"empirical ROC curve for `roctitl'; AUC = `aucstr'"' }

    local symbol `"`symbol'i"'
    local connect `"ll"'
    if `varnum' == 5 {
        local symbol `"`symbol'ii"'
        local connect `"`connect'll"'
        if `"`t2title'"'  == `""' { /*
         */  local t2title `"`level'% confidence bands     kernal half-width = `kwdt' "'
        }
    }
    qui gen `forty5' = .
    qui replace `forty5' = cond(`fpr'==0, 0,cond(`fpr'==1, 1, .))
    format `tpr' `forty5' `fpr' %4.1f
    sort `y'
    #delimit ;
    graph `tpr' `forty5' `lb' `ub' `fpr' if `touse',
        c(`connect') s(`symbol') border
        t1(`"`t1title'"') t2(`"`t2title'"')
        xlabel(`xlabel') ylabel(`ylabel')
        xline(`xline') yline(`yline')
        l2(`l2title') b2(`b2title') l1(`"`l1title'"')
        `options' ;
    #delimit cr
end
********************************************************************************
prog def cbgen
    args roct t d level nd ndb y touse se lb ub kdd kddb kernal kwidth
    tempname z
    tempvar kdx fdc kdx2 fdbc var
    scalar `z' = invnorm(.5 + `level'/200.)
    kdensity `y' if `d'==0 & `touse',gen(`kdx' `fdbc') at(`y') nogr `kernal' `kwidth'
    kdensity `y' if `d'==1 & `touse',gen(`kdx2' `fdc') at(`y') nogr `kernal' `kwidth'
    qui {
    replace `fdbc' = 1e-30 if `fdbc'==0
    gen double `var' = `roct'*(1.0 - `roct')/`nd' + ((`fdc'/`fdbc')^2)*(`t'*(1.0-`t')/`ndb')
    replace `se' = sqrt(`var')
    replace `lb' = `roct' - `z'*`se'
    replace `ub' = `roct' + `z'*`se'
    replace `kdd' = `fdc'
    replace `kddb' = `fdbc'
    }
    /* save as permanent vbls for troubleshooting (for now)
        cap drop se-lb
        gen se = `se'
        gen fdc = `fdc'
        gen fdbc = `fdbc'
        gen t = `t'
        gen roct = `roct'
        gen ratio = `fdc'/`fdbc'
        gen ub = `ub'
        gen lb = `lb'
    */
end
******************************************************************
prog def varauc,rclass
    /* calculate and return variance(AUC)
      per DeLong, DeLong, and Clarke-Pearson (1988)
        but w/o correction for ties
      (as described in Pepe, OUP book, Ch.5
    */
    args y d touse nd ndb
    tempname varauc
    tempvar vd v
    gsort `y' -`d'
    qui gen double `vd' = 1 - sum(`d'==0 & `touse')/`ndb'
    gsort -`y' -`d'
    qui gen double `v' = sum(`d'==1 & `touse')/`nd'
    qui sum `vd' if `d'==1
    sca `varauc' = r(Var)
    qui sum `v' if `d'==0
    return scalar vauc = `varauc'/`nd' + r(Var)/`ndb'
end
******************************************************************
prog def auc,rclass
    /* calculate auc or partial auc */
    args y d touse partial
    tempvar yt
    tempname ydmax nd ndb tpr_t fpr_t y_perc auc_adj
    qui gen `yt' = `y'
    if `partial'~=0 {
        cap assert `partial' >=1.0 & `partial' <=100.0
        if _rc~=0 {
            di in red "argument for partial() option must be between 1 & 100"
            exit
            }
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
