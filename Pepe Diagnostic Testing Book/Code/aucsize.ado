*! 1.1.2 GML 5 Feb 2003
*  vers 1.1.0 : added logit option (includes power calculations using logit(auc)-based
*             confidence bound.
*  vers 1.1.1 logit based power calc included as returned result
*  vers 1.1.2: eliminate call to user written program -univar-
*
* 1.0.0 GML 5 Sep 2001
* 1.1.1 GML 7 Nov 2001
*
* for ROC book ch. 9
prog def aucsize, rclass
    version 7.0
    preserve
    drop _all
    local dalpha = 1 - $S_level/100
    syntax,  MDB(real) SDB(real) MD(real) SD(real)  /*
       */  ND(int) NDB(int)                         /*
       */  [AUCnull(real .5) Alpha(real `dalpha') DELong MP LOgit /*
       */   NSim(int 500) Dots noDISplay RESfile(string) REPLACE ]

    if `aucnull' > 1 | `aucnull' < 0 {
        di in red "argument to aucnull(#) must be between 0 & 1"
        exit 198
    }
    if `alpha' > 1 | `alpha' < 0 {
        di in red "argument to alpha(#) must be between 0 & 1"
        exit 198
    }
    if "`delong'" !="" & "`mp'"!="" {
        di in red "can currently specify at most ONE of the alternative AUC variance"
        di in red "  estimator options (DeLong or mp)"
        exit 198
    }

    if "`delong'"  != "" {
        local varmeth = "`delong'"
        local methstr = "DeLong, DeLong, and Clarke-Pearson"
    }
    else if "`mp'" != "" {
        local varmeth = "`mp'"
        local methstr = "MP version of DeLong"
    }
    else {
        local varmeth = "hanley"
        local methstr = "Hanley & McNeil"
    }

    if "`resfile'"=="" {tempfile resfile}
    else {local ressave "yes" }

    if "`replace'"!="" { local replacm ",`replace'" }

    local dots = cond("`dots'"=="", "*", "noisily")

    tempname mysim
    if `"`logit'"' ~= "" {local aucvname "lauc_lb"}
    qui postfile `mysim' auc auc_lb `aucvname' using `resfile' `replacm'

    simres `nsim' `resfile' `alpha' `nd' `ndb' `md' `mdb' `sd' `sdb' /*
       */  `aucnull' `mysim' `varmeth' `dots' `logit'

    use `resfile',clear
    qui gen byte reject1 = (auc_lb > `aucnull') if auc_lb~=.
    qui sum reject1,meanonly
    local n1 = r(N)
    local power1 = r(sum)/r(N)
    local nbad1 = `nsim' - r(N)

    if `"`logit'"' ~= `""' {
        qui gen byte reject2 = (lauc_lb > `aucnull')    if lauc_lb~=.
        qui sum reject2,meanonly
        local n2 = r(N)
        local power2 = r(sum)/r(N)
        local nbad2 = `nsim' - r(N)
    }

    tempname a b bnauc
    sca `a' = (`md' - `mdb')/`sd'
    sca `b' = `sdb'/`sd'
    sca `bnauc' = normprob(`a'/sqrt(1+`b'^2))

    if `"`display'"' == `""' {
       di
       di as text "Power to reject Ho: AUC = " as result `aucnull'
       di
       if `"`logit'"' == "" {
           di as text _skip(8) "==> " as result %5.3g `power1'
           di
       }
       else {
           di as text _skip(8) "==> " as result %5.3g `power1' as text "  (Wald interval rejection)"
           di as text _skip(8) "==> " as result %5.3g `power2' as text "  (logit(auc)-based interval rejection)"
           di
       }
       di as text " based on empirical (non-parametric) AUC estimate
       di
       di as text " for sample size  n_d  = " as result `nd'
       di as text "                  n_db = " as result `ndb'
       di
       di as text "  and test value distributions,"
       di as text "       Y_db  ~ N(" as result `mdb' as text ", " as result `sdb' in g ")"
       di as text "       Y_d   ~ N(" as result `md' as text ", " as result `sd' in g ")"
       di
       di as text "          # simulations:" as result %5.0f `nsim'
       di
       di as text "    # simulations with invalid conf. bound  (Wald):" as result %4.0f `nbad1'
       if `"`logit'"' ~= "" {
           di as text "    # simulations with invalid conf. bound (logit):" as result %4.0f `nbad2'
       }
       di as text "     (these are disregarded in power calculations)"
       di
       di as text "  power based on lower 1-sided " as result (1-`alpha')*100 as text "% confidence bound"
       di as text "    and " as result "`methstr'" in g " variance estimator
       di
       di as text "  AUC based on binormal ROC for the specified distributions: " as result %5.3g `bnauc'
       di
       di
       di in text " AUC distribution for simulated samples: "
       di
       di in text " {ralign 23:#}{ralign 55:{hline 4} quantiles {hline 4}}
       di in text "                    samples   mean     sd      min    max     .25    mdn    .75
       di in text " {hline 17}{c +}{hline 60}
       qui sum auc,d
       di in text " {ralign 16: non-param. AUC} {c |}" as result %6.0f `r(N)' %9.2f `r(mean)' /*
        */  %8.3f `r(sd)' %8.2f `r(min)' %7.2f `r(max)' %8.2f `r(p25)' %7.2f `r(p50)' %7.2f `r(p75)'
       qui sum auc_lb,d
       di in text " {ralign 16:lower conf bound} {c |}" as result %6.0f `r(N)' %9.2f `r(mean)' /*
        */  %8.3f `r(sd)' %8.2f `r(min)' %7.2f `r(max)' %8.2f `r(p25)' %7.2f `r(p50)' %7.2f `r(p75)'
       di
    }
    return scalar power_W = `power1'
    if "`logit'" ~= "" {return scalar power_L = `power2'}
end
*************************************************************
program define simres
    version 7
    /* drive simulation loop with calls to doasim,
       post simulations results,
       and close result file
    */
    * set trace on
    args nsim resfile alpha nd ndb md mdb sd sdb aucnull mysim varmeth dots logit
    tempname alpha_2 mu_1

    local n = `nd' + `ndb'
    qui set obs `n'
    qui gen byte d = (_n<=`nd')
    qui gen y = .
    local level = 100. * (1.0 - `alpha')
    local level2s =  100. * (1.0 - 2*`alpha')
    * local dots = cond("`dots'"=="", "*", "noisily")
    if "`logit'" ~= "" {local logitvarnm "(r(lauc_lb))"}
    forvalues i = 1/`nsim' {
       `dots' di as text "." _c
       * set trace on
       doasim `md' `mdb' `sd' `sdb' `aucnull'  `nd' `ndb' `level' `level2s' `varmeth' `logit'
       post `mysim' (r(auc)) (r(auc_lb)) `logitvarnm'
    }
    postclose `mysim'
end
********************************************************
pro def doasim, rclass
    version 7.0
    args md mdb sd sdb aucnull nd ndb level level2s varmeth logit
    tempname auc auc_1 se logitauc lb
    qui replace y = `md'*(d==1) + `mdb'*(d==0) + /*
                 */   (`sd'*(d==1) + `sdb'*(d==0)) * invnorm(uniform())
    if "`varmeth'"=="mp" {
        auc y d
        scalar `auc' = r(auc)
        return scalar auc = `auc'
        varauc y d `nd' `ndb'
        sca `se' = sqrt(r(vauc))
        return scalar auc_lb =  `auc' - ( (invnorm(`level'/100.)) * `se' )
    }
    else if "`varmeth'" == "delong" {
        qui roctab d y, lev(`level2s')
        return scalar auc = r(area)
        return scalar auc_lb = r(lb)
        sca `se' = r(se)
    }
    else {
        qui roctab d y, lev(`level2s') hanley
        return scalar auc = r(area)
        return scalar auc_lb = r(lb)
        sca `se' = r(se)
    }
    if `"`logit'"' ~= "" {
        if return(auc) == 0 {
            scalar `auc_1' = .0001
        }
        else if return(auc) == 1 {
            scalar `auc_1' = .9999
        }
        else {
            scalar `auc_1' = return(auc)
        }
        scalar `logitauc' = log( (`auc_1')/(1.0 - `auc_1'))
        scalar `lb' = `logitauc' - (invnorm(`level'/100.)) * (`se'/ (`auc_1'*(1 - `auc_1')))
        return scalar lauc_lb = exp(`lb')/(1.0 + exp(`lb'))
    }
end
******************************************************************
prog def varauc,rclass
    /* calculate and return variance(AUC) */
    args y d nd ndb
    tempname varauc
    tempvar vd v
    gsort `y' -`d'
    qui gen double `vd' = sum(`d'==0)/`ndb'
    gsort -`y' `d'
    qui gen double `v' = sum(`d'==1)/`nd'
    qui replace `v' = `vd' if `d'==1
    qui sum `v' if `d'==1
    sca `varauc' = r(Var)
    qui sum `v' if `d'==0
    return scalar vauc = `varauc'/`nd' + r(Var)/`ndb'
end
*********************************************************
prog def auc,rclass
    /* calculate auc */
    args y d
    tempvar yt
    tempname ydmax nd ndb tpr_t fpr_t y_perc auc_adj
    qui gen `yt' = `y'
    qui ranksum `yt',by(`d')
    return scalar auc = 1.0 - (r(sum_obs) - r(N_1)*(r(N_1) + 1)/2.0) / (r(N_1) * r(N_2))
end
