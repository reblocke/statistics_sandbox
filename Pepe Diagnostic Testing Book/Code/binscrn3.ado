*! 1.1.0 01 Sept 2005
* 1.0.0 09 March 2000
* comparison of binary screening tests, paired data
* originally for MP Diagnostic Methods course
* 7/2005 modified for SEMT section 3.3
prog def binscrn3, rclass
    version 9
    syntax varlist(min=3 max=3) [if] [in], [ Level(integer $S_level) Format(string) NIter(integer 20)]
    tokenize `varlist'
    local y1 "`1'"
    local y2 "`2'"
    local d "`3'"
    marksample touse
    markout `touse' `y1' `y2' `d'

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
    qui ta `y1' if `touse'
    if r(r)~=2 {
        di in red "test result variable `y1' must take on two values"
        exit 198
    }
    qui sum `y1' if `touse',meanonly
    if r(min)~=0 | r(max)~=1 {
        di in red "test result variable `y1' must be 0/1"
        exit 198
    }
    qui ta `y2' if `touse'
    if r(r)~=2 {
        di in red "test result variable `y2' must take on two values"
        exit 198
    }
    qui sum `y2' if `touse',meanonly
    if r(min)~=0 | r(max)~=1 {
        di in red "test result variable `y2' must be 0/1"
        exit 198
    }
    if `level'<10 | `level'>99 {
        local level 95
    }
    if `"`format'"'==`""' {
        local format "%5.3g"
    }
    else {
        cap qui di `format' _pi
        if _rc~=0 | index("`format'","%")~=1 {
            di in red "format(#) option is invalid"
            exit 120
        }
    }
    local alpha = (100. - `level')/100.0
    tempname alpha_jnt
    sca `alpha_jnt' = 1 - sqrt(1-`alpha')

    local t1lb: var label `y1'
    local t2lb: var label `y2'
    if "`t1lb'"=="" {
        local t1lb "`y1'"
    }
    if "`t2lb'"=="" {
        local t2lb "`y2'"
    }
    local t1par " (`t1lb')"
    local t2par " (`t2lb')"

    di
    di "{txt} Comparison of binary tests, paired data"
    di
    di "{txt}   Test 1: `t1lb'"
    di "{txt}   Test 2: `t2lb'"
    di
    di "{txt}************************************************************ "
    di
    di "{txt}True positive fraction:"
    tempname M1 tpf1 tpf2 tpr varlogr z tpr_lb tpr_ub p_mcc z p_z
    di
    di "{txt} Test result concordance among those with disease:
    ta `y1' `y2' if `d'==1 & `touse',matcell(`M1')
    di

    local nd = r(N)
    cheng `=`M1'[1,1]' `=`M1'[1,2]' `=`M1'[2,1]' `=`M1'[2,2]' `r(N)' `alpha'
    sca `z' = r(z)
    sca `p_z' = 2*(normal(-abs(`r(z)')))
    return scalar rtpf = r(pratio)
    return scalar rtpf_lb = r(lci)
    return scalar rtpf_ub = r(uci)
    sca `tpf1' = r(p1)
    sca `tpf2' = r(p2)

    di "{txt}       TPF1: " as res `format' `tpf1'
    di "{txt}       TPF2: " as res `format' `tpf2'
    di
    di "{txt}   rTPF = TPF1/TPF2 (`level'% CI) : " as res `format' `return(rtpf)' ///
       "  ("`format' `return(rtpf_lb)' ","`format' `return(rtpf_ub)' ")"
    di
    qui mcc `y1' `y2' if `d'==1 & `touse'  /* will obtain McNemar's test statistic from this
                                    (see [R] epitab or help for -epitab- ) */
    sca `p_mcc' = chiprob(1,r(chi2))
    di "{txt}    test of Ho: rTPF = 1.0"
    di
    di "{txt}        McNemar's statistic"
    di "{txt}            chi2(1) = " as res `format' r(chi2) as txt ",  p = " as res %6.5g `p_mcc'
    di
    di "{txt}        test based on asymptotic normality of ln(rTPF) "
    di "{txt}                  Z = " as res `format' `z' as txt ",  p = " as res %6.5g `p_z'
    di
    * get CI's for adjusted alpha_jnt, for joint conf region
    cheng `=`M1'[1,1]' `=`M1'[1,2]' `=`M1'[2,1]' `=`M1'[2,2]' `nd' `alpha_jnt'
    return scalar rtpfj_lb = r(lci)
    return scalar rtpfj_ub = r(uci)

    **** Compare FPF1 vs FPF2   *****************
    tempname fpf1 fpf2 fpr fpr_lb fpr_ub
    di "{txt} *********************************************************
    di
    di "{txt}False positive fraction:
    di
    di "{txt} Test result concordance among those without disease:
    ta `y1' `y2' if `d'==0 & `touse',matcell(`M1')
    di

    local ndb = r(N)
    cheng `=`M1'[1,1]' `=`M1'[1,2]' `=`M1'[2,1]' `=`M1'[2,2]' `r(N)' `alpha'
    sca `z' = r(z)
    sca `p_z' = 2*(normal(-abs(`r(z)')))
    sca `fpf1' = r(p1)
    sca `fpf2' = r(p2)
    return scalar rfpf = r(pratio)
    return scalar rfpf_lb = r(lci)
    return scalar rfpf_ub = r(uci)

    di "{txt}       FPF1: " as res `format' `fpf1'
    di "{txt}       FPF2: " as res `format' `fpf2'
    di
    di "{txt}   rFPF = FPF1/FPF2 (`level'% CI) : " as res `format' `return(rfpf)' ///
       "  ("`format' `return(rfpf_lb)' ","`format' `return(rfpf_ub)' ")"
    di
    qui mcc `y1' `y2' if `d'== 0 & `touse'
    sca `p_mcc' = chiprob(1,r(chi2))
    di "{txt}    test of Ho: rFPF = 1.0"
    di
    di "{txt}        McNemar's statistic"
    di "{txt}            chi2(1) = " as res `format' r(chi2) as txt ",  p = " as res %6.5g `p_mcc'
    di
    di "{txt}        test based on aymptotic normality of ln(rFPF) "
    di "{txt}                  Z = " as res `format' `z' as txt ",  p = " as res %6.5g `p_z'
    di
    * get CI's for adjusted alpha_jnt, for joint conf region
    cheng `=`M1'[1,1]' `=`M1'[1,2]' `=`M1'[2,1]' `=`M1'[2,2]' `ndb' `alpha_jnt'
    return scalar rfpfj_lb = r(lci)
    return scalar rfpfj_ub = r(uci)

    **********************************************************
    * joint CI's for rFPF, rTPF

    di
    di "{txt}****************"
    di "{txt} `level'% joint rectangular confidence region for (rFPF, rTPF):"
    di
    di _col(10)"{txt}rFPF`jointpar' = {res}" `format' `return(rfpf)' " (" `format' `return(rfpfj_lb)' ///
                "," `format' `return(rfpfj_ub)' ")"
    di _col(10)"{txt}rTPF`jointpar' = {res}" `format' `return(rtpf)' " (" `format' `return(rtpfj_lb)' ///
                "," `format' `return(rtpfj_ub)' ")"
    di

    /* get GLM based confidence limits for rPPV and rNPV */
    glm_pv `y1' `y2' `d' `touse' `level' `niter'
    return scalar rnpv_lb = r(rnpv_lb)
    return scalar rnpv_ub = r(rnpv_ub)
    return scalar rppv_lb = r(rppv_lb)
    return scalar rppv_ub = r(rppv_ub)
    local glm_converge = r(glm_converge)

    **** Compare PPV1 vs PPV2   *******************
    tempname M2
    di "{txt} *********************************************************
    di
    di "{txt} Positive Predictive values:
    di
    qui ta `d' if `y1' & `touse',matcell(`M2')
    return sca ppv1 = `M2'[2,1]/r(N)
    qui ta `d' if `y2' & `touse',matcell(`M2')
    return sca ppv2 = `M2'[2,1]/r(N)
    return sca rppv = `return(ppv1)'/`return(ppv2)'

    di "{txt}        PPV1: " as res `format' `return(ppv1)'
    di "{txt}        PPV2: " as res `format' `return(ppv2)'
    di
    di "{txt}   rPPV = PPV1/PPV2 (`level'% CI*) : " as res `format' `return(rppv)' ///
         "  ("`format' `return(rppv_lb)' ","`format' `return(rppv_ub)' ")"
    di

    **** Compare NPV1 vs NPV2   ************************
    tempname npv1 npv2 npvr
    di "{txt} *********************************************************
    di
    di "{txt} Negative Predictive values:
    di
    qui ta `d' if `y1'==0 & `touse',matcell(`M2')
    return sca npv1 = `M2'[1,1]/r(N)
    qui ta `d' if `y2'==0 & `touse',matcell(`M2')
    return sca npv2 = `M2'[1,1]/r(N)
    return sca rnpv = `return(npv1)'/`return(npv2)'

    di "{txt}        NPV1: " as res `format' `return(npv1)'
    di "{txt}        NPV2: " as res `format' `return(npv2)'
    di
    di "{txt}   rNPV = NPV1/NPV2 (`level'% CI*) : " as res `format' `return(rnpv)' ///
         "  ("`format' `return(rnpv_lb)' ","`format' `return(rnpv_ub)' ")"
    di
    di
    di "{txt} *CI's for rPPV and rNPV based on GLM agreement model "
    di

    **** Compare DLR+_1 vs DLR+_2 ************************
    tempname dlrp1 dlrp2 dlrprat
    di "{txt} ********************************************************* "
    di
    di "{txt} Diagnostic Likelihood Ratios (positive):
    di

    sca `dlrp1' = `tpf1'/`fpf1'
    sca `dlrp2' = `tpf2'/`fpf2'
    sca `dlrprat' = `dlrp1'/`dlrp2'

    di "{txt}          DLR+_1: " as res `format' `dlrp1'
    di "{txt}          DLR+_2: " as res `format' `dlrp2'
    di "{txt}   DLR+_1/DLR+_2: " as res `format' `dlrprat'
    di

    **** Compare DLR-_1 vs DLR-_2 ************************
    tempname dlrn1 dlrn2 dlrnrat
    di "{txt} *********************************************************
    di
    di "{txt} Diagnostic Likelihood Ratios (negative):
    di

    sca `dlrn1' = (1-`tpf1')/(1-`fpf1')
    sca `dlrn2' = (1-`tpf2')/(1-`fpf2')
    sca `dlrnrat' = `dlrn1'/`dlrn2'

    di "{txt}          DLR-_1: " as res `format' `dlrn1'
    di "{txt}          DLR-_2: " as res `format' `dlrn2'
    di "{txt}   DLR-_1/DLR-_2: " as res `format' `dlrnrat'
    di
end
****************************
prog def cheng, rclass
    /* program to calculate and return confidence limits for ratio of
        independent proportions; per 2x2 table for paired data,
        and statistic Z = ln(p1/p2)/se
        SEMT sec 3.3.4   Cheng and Macaluso (1997)
       arguments (table cells):
           nn np pn pp
           alpha - for 1-alpha level confidence limits
    */
    args nn np pn pp n alpha
    tempname p1 p2 z se pratio
    sca `z' = invnorm(1.0 - (`alpha'/2.0))
    return sca p1 = (`pn' + `pp')/`n'
    return sca p2 = (`np' + `pp')/`n'

    sca `se' = sqrt( (`np' + `pn')/((`np' + `pp')*(`pn' + `pp')) )

    return sca pratio = return(p1)/return(p2)
    return sca lci = exp( ln(return(pratio)) - `z'*`se')
    return sca uci = exp( ln(return(pratio)) + `z'*`se')
    return sca z = ln(return(pratio))/`se'
end
****************************
prog def glm_pv, rclass
    version 9
    /* obtain rPPV and rNPV and confidence bounds
       via GLM agreement model (SEMT sec 3.6.2)   */
    args y1 y2 d touse level niter
    preserve
    tempvar subj test y ytest w
    tempname z
    local alpha = (100. - `level')/100.0
    sca `z' = invnorm(1.0 - (`alpha'/2.0))
    qui {
        gen int `subj' = _n
        expand 2 if `touse'
        sort `subj'
        by `subj': gen int `test' = (_n==1)
        by `subj': gen int    `y' = `y1' if _n==1
        by `subj': replace    `y' = `y2' if _n==2
        gen int `w' = (`y'==`d')
        gen int `ytest' = `y'*`test'
    }
    qui cap glm `w' `y' `test' `ytest' if `touse', ///
         family(binomial) link(log) cluster(`subj') iterate(`niter') level(`level')
    if e(converged) == 0 {
        qui glm `w' `y' `test' `ytest' if `touse', ///
             family(binomial) link(log) cluster(`subj') difficult iterate(`niter') level(`level')
    }
    if e(converged) == 0 {
        return sca glm_converge = 0
        exit
    }
    return sca rnpv = exp(_b[`test'])
    return sca rnpv_lb = exp(_b[`test'] - `z'*_se[`test'])
    return sca rnpv_ub = exp(_b[`test'] + `z'*_se[`test'])
    qui lincom `test' + `ytest'
    return sca rppv = exp(r(estimate))
    return sca rppv_lb = exp( r(estimate) - `z'*r(se) )
    return sca rppv_ub = exp( r(estimate) + `z'*r(se) )
    return sca glm_converge = 1
end
