*! 1.1.0 20 July 2005
* 1.0.2 09 March 2000
* comparison of binary screening tests, unpaired data
* originally written for MP Diagnostic Methods course
* version 1.1.0 updated for MP SEMT book, section 3.2
*   additions:  - returned values for relative measures and conf. bounds
*               - calculation of joint confidence regions
*               - return of Cov matrix for (lnDLR+,lnDLR-)
prog def binscrn2, rclass
    version 8.0
    syntax varlist(min=3 max=3) [if] [in], [ Level(integer $S_level)  Format(string) ]
    tokenize `varlist'
    local y "`1'"
    local d "`2'"
    local gp "`3'"
    pause on
    marksample touse
    markout `touse' `y' `d' `gp'
    /*  assure that disease and test vars are 0/1 only: */
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
    qui ta `gp' if `touse'
    if r(r)~=2 {
        di in red "test type variable `gp' must take on two values"
        exit 198
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
    qui sum `gp' if `touse', meanonly
    local t1val = r(min)
    local t2val = r(max)

    local t1lb: label (`gp') `t1val'
    local t2lb: label (`gp') `t2val'
    if "`t1val'"~="" & "`t2val"~="" {
        local t1par " (`t1lb')"
        local t2par " (`t2lb')"
        local jointpar "(`t1lb',`t2lb')"
    }
    local nolab: value label `gp'
    if "`nolab'" == "" {
        local t1par ""
        local t2par ""
        local jointpar ""
        }

    if `level'<10 | `level'>99 {
        local level 95
    }
    local alpha = (100. - `level')/100.0
    di
    di in g "Comparison of binary tests, unpaired data
    di
    tempname M1 tpf1 tpf2 fpf1 fpf2 nd1 nd2 ndb1 ndb2
    di in g "************************************************************
    di
    di in g "True Postive fraction: test `t1val'`t1par' vs test `t2val'`t2par'"
    di
    di in g " Test outcome by test type among diseased:"

    ta `gp' `y' if `d'==1 & `touse',r chi2 exact matcell(`M1')

    ***  CI's for TPF_1 / TPF_2 using Simel's formula ***********

    sca `nd1' = `M1'[1,1] + `M1'[1,2]
    sca `nd2' = `M1'[2,1] + `M1'[2,2]
    sca `tpf1' = `M1'[1,2]/`nd1'
    sca `tpf2' = `M1'[2,2]/`nd2'
    simel `tpf1' `tpf2' `nd1' `nd2' `alpha'

    di
    di in g "                        TPF1:  " in y `format' `tpf1'
    di in g "                        TPF2:  " in y `format' `tpf2'
    di in g "    rTPF = TPF1/TPF2 (`level'% CI): " in y `format' r(ratio) " (" `format' r(lci) /*
        */                 "," `format' r(uci) ")"
    di
    return scalar rtpf = r(ratio)
    return scalar rtpf_lb = r(lci)
    return scalar rtpf_ub = r(uci)
    **********************************************************
    * compare FPF_1 vs FPF_2:

    di in g "***************************************************"
    di
    di in g "False Postive fraction: test `t1val'`t1par' vs test `t2val'`t2par'"
    di
    di in g " test outcome by test type among non-diseased:"

    ta `gp' `y' if `d'==0 & `touse',r chi2 exact matcell(`M1') nokey

    ***  CI's for TPF_1 / TPF_2 using Simel's formula ***********

    sca `ndb1' = `M1'[1,1] + `M1'[1,2]
    sca `ndb2' = `M1'[2,1] + `M1'[2,2]
    sca `fpf1' = `M1'[1,2]/`ndb1'
    sca `fpf2' = `M1'[2,2]/`ndb2'
    simel `fpf1' `fpf2' `ndb1' `ndb2' `alpha'

    di
    di in g "                        FPF1:  " in y `format' `fpf1'
    di in g "                        FPF2:  " in y `format' `fpf2'
    di in g "   rFPF = FPF1/FPF2 (`level'% CI): " in y `format' r(ratio) " (" `format' r(lci) /*
        */                 "," `format' r(uci) ")"
    di
    return scalar rfpf = r(ratio)
    return scalar rfpf_lb = r(lci)
    return scalar rfpf_ub = r(uci)
    **********************************************************
    * joint CI's for rFPF, rTPF

    tempname alpha_jnt
    sca `alpha_jnt' = 1 - sqrt(1-`alpha')
    di
    di "{txt}****************"
    di "{txt} `level'% joint rectangular confidence region for (rFPF, rTPF):"
    di
    simel `fpf1' `fpf2' `ndb1' `ndb2' `alpha_jnt'
    di _col(10)"{txt}rFPF`jointpar' = {res}" `format' r(ratio) " (" `format' r(lci) "," `format' r(uci) ")"
    return scalar rfpfj = r(ratio)
    return scalar rfpfj_lb = r(lci)
    return scalar rfpfj_ub = r(uci)

    simel `tpf1' `tpf2' `nd1' `nd2' `alpha_jnt'
    di _col(10)"{txt}rTPF`jointpar' = {res}" `format' r(ratio) " (" `format' r(lci) "," `format' r(uci) ")"
    di
    return scalar rtpfj = r(ratio)
    return scalar rtpfj_lb = r(lci)
    return scalar rtpfj_ub = r(uci)
    **********************************************************
    * compare PPV_1 vs PPV_2:

    tempname npos1 npos2 ppv1 ppv2 alpha_jnt
    di in g "************************************************************
    di
    di in g "Positive Predictive Values: test `t1val'`t1par' vs test `t2val'`t2par'"
    di
    di in g " disease status by test type among test positive individuals:"

    ta `gp' `d' if `y'==1 & `touse',r chi2 exact matcell(`M1') nokey

    sca `npos1' = `M1'[1,1] + `M1'[1,2]
    sca `npos2' = `M1'[2,1] + `M1'[2,2]
    sca `ppv1' = `M1'[1,2]/`npos1'
    sca `ppv2' = `M1'[2,2]/`npos2'
    simel `ppv1' `ppv2' `npos1' `npos2' `alpha'

    di
    di in g "                      PPV1:  " in y `format' `ppv1'
    di in g "                      PPV2:  " in y `format' `ppv2'
    di in g " rPPV = PPV1/PPV2 (`level'% CI): " in y `format' r(ratio) " (" `format' r(lci) /*
        */                 "," `format' r(uci) ")"
    di
    return scalar rppv = r(ratio)
    return scalar rppv_lb = r(lci)
    return scalar rppv_ub = r(uci)
    **********************************************************
    * compare NPV_1 vs NPV_2:

    tempname nneg1 nneg2 npv1 npv2
    di in g "***************************************************"
    di
    di in g "Negative Predictive Values: test `t1val'`t1par' vs test `t2val'`t2par'"
    di
    di in g " disease status by test type among test negative individuals:"

    ta `gp' `d' if `y'==0 & `touse',r chi2 exact matcell(`M1') nokey

    sca `nneg1' = `M1'[1,1] + `M1'[1,2]
    sca `nneg2' = `M1'[2,1] + `M1'[2,2]
    sca `npv1' = `M1'[1,1]/`nneg1'
    sca `npv2' = `M1'[2,1]/`nneg2'
    simel `npv1' `npv2' `nneg1' `nneg2' `alpha'

    di
    di in g "                      NPV1:  " in y `format' `npv1'
    di in g "                      NPV2:  " in y `format' `npv2'
    di in g " rNPV = NPV1/NPV2 (`level'% CI): " in y `format' r(ratio) " (" `format' r(lci) /*
        */                 "," `format' r(uci) ")"
    di
    return scalar rnpv = r(ratio)
    return scalar rnpv_lb = r(lci)
    return scalar rnpv_ub = r(uci)
    **********************************************************
    * joint CI's for rNPV, rPPV

    tempname alpha_jnt
    sca `alpha_jnt' = 1 - sqrt(1-`alpha')
    di
    di "{txt}****************"
    di "{txt} `level'% joint rectangular confidence region for (rNPV, rPPV):"
    di
    simel `npv1' `npv2' `nneg1' `nneg2' `alpha_jnt'
    di _col(10)"{txt}rNPV`jointpar' = {res}" `format' r(ratio) " (" `format' r(lci) "," `format' r(uci) ")"
    return scalar rnpvj = r(ratio)
    return scalar rnpvj_lb = r(lci)
    return scalar rnpvj_ub = r(uci)

    simel `ppv1' `ppv2' `npos1' `npos2' `alpha_jnt'
    di _col(10)"{txt}rPPV`jointpar' = {res}" `format' r(ratio) " (" `format' r(lci) "," `format' r(uci) ")"
    di
    return scalar rppvj = r(ratio)
    return scalar rppvj_lb = r(lci)
    return scalar rppvj_ub = r(uci)
    **********************************************************
    * compare DLR+_1 vs DLR+_2:

    tempname z Z dlrp1 dlrp2 dlrprat se_dlrprat dlrplci dlrpuci pval
    sca `dlrp1' = `tpf1'/`fpf1'
    sca `dlrp2' = `tpf2'/`fpf2'
    sca `dlrprat' = `dlrp1'/`dlrp2'
    sca `se_dlrprat' = sqrt(  (1-`tpf1')/(`tpf1'*`nd1') + (1-`fpf1')/(`fpf1'*`ndb1')     /*
        */        + (1-`tpf2')/(`tpf2'*`nd2') + (1-`fpf2')/(`fpf2'*`ndb2') )
    sca `z' = invnorm(1.0 - (`alpha'/2.0))
    sca `dlrplci' = exp( (ln(`dlrp1')-ln(`dlrp2')) - `z'*`se_dlrprat' )
    sca `dlrpuci' = exp( (ln(`dlrp1')-ln(`dlrp2')) + `z'*`se_dlrprat' )

    * test of Ho: DLR1+ = DLR2+:

    sca `Z' =  ( ln(`dlrp1')-ln(`dlrp2') )/`se_dlrprat'
    sca `pval' = 2*(1.0-normprob(abs(`Z')))
    di in g "*******************************************************
    di
    di in g "Diagnostic Likelihood Ratios (positive)
    di
    di in g  "               DLR1+: " in y `format' `dlrp1'
    di in g  "               DLR2+: " in y `format' `dlrp2'
    di
    di in g " rDLR+ = DLR1+/DLR2+ (`level'% CI) : " in y `format' `dlrprat' " (" `format' `dlrplci' "," `format' `dlrpuci' ") "
    di
    di in g " Test of Ho: DLR1+ = DLR2+ :   Z = " in y %5.3f `Z' in g ",  p = " in y `format' `pval'
    di
    return scalar rdlrp = `dlrprat'
    return scalar rdlrp_lb = `dlrplci'
    return scalar rdlrp_ub = `dlrpuci'

    di in g "***************************************************"
    di
    di in g "Diagnostic Likelihood Ratios (negative)

    * compare DLR-_1 vs DLR-_2:

    tempname dlrn1 dlrn2 dlrnrat se_dlrnrat dlrnlci dlrnuci r2 Cov
    sca `dlrn1' = (1-`tpf1')/(1-`fpf1')
    sca `dlrn2' = (1-`tpf2')/(1-`fpf2')
    sca `dlrnrat' = `dlrn1'/`dlrn2'
    sca `se_dlrnrat' = sqrt(  `tpf1'/((1-`tpf1')*`nd1') + `fpf1'/((1-`fpf1')*`ndb1')     /*
        */        + `tpf2'/((1-`tpf2')*`nd2') + `fpf2'/((1-`fpf2')*`ndb2') )
    sca `z' = invnorm(1.0 - (`alpha'/2.0))
    sca `dlrnlci' = exp( (ln(`dlrn1')-ln(`dlrn2')) - `z'*`se_dlrnrat' )
    sca `dlrnuci' = exp( (ln(`dlrn1')-ln(`dlrn2')) + `z'*`se_dlrnrat' )

    * test of Ho: DLR1- = DLR2-:

    sca `Z' =  ( ln(`dlrn1')-ln(`dlrn2') )/`se_dlrnrat'
    sca `pval' = 2*(1.0-normprob(abs(`Z')))

    di
    di in g "                DLR1-: " in y `format' `dlrn1'
    di in g "                DLR2-: " in y `format' `dlrn2'
    di
    di in g " rDLR- = DLR1-/DLR2- (`level'% CI) : " in y `format' `dlrnrat' " (" `format' `dlrnlci' "," `format' `dlrnuci' ") "
    di
    di in g " Test of Ho: DLR1- = DLR2- :   Z = " in y %5.3f `Z' in g ",   p = " in y `format' `pval'
    di
    return scalar rdlrn = `dlrnrat'
    return scalar rdlrn_lb = `dlrnlci'
    return scalar rdlrn_ub = `dlrnuci'

    * these are not independent, so cannot construct joint rectangular conf region.
    *  rather, return covariance matrix, and log(DLR) estimate vector - can be used to
    *  generate elliptical conf region for log(rDLR+),log(rDLR-)

    sca `r2' = -( 1/`nd1' + 1/`ndb1' + 1/`nd2' + 1/`ndb2' )

    matrix `Cov' = (`se_dlrprat'^2,`r2' \ `r2', `se_dlrnrat'^2)
    matrix rownames `Cov' = log(rDLR+) log(rDLR-)
    matrix colnames `Cov' = log(rDLR+) log(rDLR-)
    matrix rDLR = (ln(`dlrprat') \ ln(`dlrnrat'))
    matrix rownames rDLR = log(rDLR+) log(rDLR-)
    return matrix Cov_ln_rDLR = `Cov'
    return matrix ln_rDLR = rDLR
end
*****************************
prog def simel, rclass
    /* program to calculate and return Simel confidence limits for ratio of
        independent proportions,  p1/p2
       arguments:
           p1 - proportion 1
           p2 - proportion 2
           n1 - denominator for p1
           n2 - denominator for p2
           alpha - for 1-alpha level confidence limits
    */
    args p1 p2 n1 n2 alpha
    tempname z se pratio
    sca `z' = invnorm(1.0 - (`alpha'/2.0))
    sca `se' = sqrt((1-`p1')/(`p1'*`n1') + (1-`p2')/(`p2'*`n2'))
    sca `pratio' = `p1'/`p2'
    return sca lci = exp(ln(`pratio') - `z'*`se')
    return sca uci = exp(ln(`pratio') + `z'*`se')
    return sca ratio = `pratio'
end
