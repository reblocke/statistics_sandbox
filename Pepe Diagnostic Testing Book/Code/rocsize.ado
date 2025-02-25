*! 1.1.1 GML Nov 11 2001
*
* 1.0.0 GML Jan 11 2001
*  additional calculations using confidence bounds based on
*    asymptotic normalitity of logit(ROC(t)) added in vers 1.1
*
*  1.1.1 from 1.1.0 minor display changes
*        name changes for returned results
*
pro def rocsize, rclass
    version 7
    preserve
    drop _all
    gettoken tpr_1 0 : 0, parse(" ,")
    gettoken fpr_1 0 : 0, parse(" ,")
    local dalpha = 1 - $S_level/100
    syntax , ND(int) NDB(int) TPNull(real) FPNull(real) [Alpha(real `dalpha') /*
      */     NSim(int 500) KERNal(string) KWIDth(real -1.0) Dots noDISplay    /*
      */     RESfile(string) REPLACE ]
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

    if "`resfile'"=="" {tempfile resfile}
    else {local ressave "yes" }

    if "`replace'"!="" { local replacm ",`replace'" }

    tempname mysim
    qui postfile `mysim' rocf rocf_lb lrocf_lb rocinvt rocinvub lrocinub using `resfile' `replacm'

    simres `nsim' `resfile' `alpha' `nd' `ndb' `tpr_1' `fpr_1' /*
       */  `fpnull' `tpnull' `mysim' `kwidth' `dots' `kernal'

    use `resfile',clear
    qui gen byte reject1 = (rocf_lb > `tpnull')
    qui replace reject1 = . if rocf_lb == .
    qui sum reject1,meanonly
    local power1 = r(sum)/r(N)
    qui gen byte reject2 = (rocinvub < `fpnull')
    qui replace reject2 = . if rocinvub == .
    qui sum reject2,meanonly
    local power2 = r(sum)/r(N)
    qui gen byte reject3 = reject1 | reject2
    qui replace reject3 = . if reject1==. | reject2==.
    qui sum reject3,meanonly
    local power3 = r(sum)/r(N)

    qui gen byte reject4 = (lrocf_lb > `tpnull')
    qui replace reject4 = . if lrocf_lb == .
    qui sum reject4,meanonly
    local power4 = r(sum)/r(N)
    qui gen byte reject5 = (lrocinub < `fpnull')
    qui replace reject5 = . if lrocinub == .
    qui sum reject5,meanonly
    local power5 = r(sum)/r(N)
    qui gen byte reject6 = reject4 | reject5
    qui replace reject6 = . if reject4==. | reject5==.
    qui sum reject6,meanonly
    local power6 = r(sum)/r(N)

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
       di "                                                    Power
       di "                                       based on asymptotic normality of:
       di "
       di "                                             roc(t)   logit(roc(t))
       di
       di in gr "    1) reject if      ROC(fp_0)_lb > tp_0    " in y %5.3g `power1' "       " %5.3g `power4'
       di in gr "    2) reject if ROC^(-1)(tp_0)_ub < fp_0    " in y %5.3g `power2' "       " %5.3g `power5'
       di in gr "    3) reject if (1) or (2)                  " in y %5.3g `power3' "       " %5.3g `power6'
       di
       * univar rocf rocf_lb lrocf_lb rocinvt rocinvub lrocinub, dec(3)
    }
    return scalar pow_l3 = `power6'
    return scalar pow_l2 = `power5'
    return scalar pow_l1 = `power4'
    return scalar pow_w3 = `power3'
    return scalar pow_w2 = `power2'
    return scalar pow_w1 = `power1'
end
*************************************************************
program define simres
    *version 6
    version 7
    * set trace on
    args nsim resfile alpha nd ndb tpr_1 fpr_1 fpr_0 tpr_0 mysim kwidth /* dots kernal */
    macro shift 10
    local 0 ",`*'"
    syntax [, Dots *]
    tempname alpha_2 mu_1

    local n = `nd' + `ndb'
    qui set obs `n'
    qui gen byte d = (_n<=`nd')
    qui gen byte d_inv = 1-d
    qui gen y = .
    qui gen ycut = .
    qui gen byte at = (_n==1) /* at arg (ycut) to kdensity must be a var,
                               but only need density for a single val */
    scalar `mu_1' = invnorm(1.0 - `fpr_1') + invnorm(`tpr_1')
    * set trace on

    /* still true?? */

    /* mult x 2 since using 1-sided test */
    * scalar `alpha_2' = 2 * (1 - sqrt(1-`alpha'))
    scalar `alpha_2' = 2 * `alpha'

    local level = 100.*(1.0 - `alpha_2')
    local dots = cond("`dots'"=="", "*", "noisily")
    local i = 1
    while `i' <= `nsim' {
       `dots' di in gr "." _c
       * set trace on
       doasim `tpr_1' `fpr_1' `fpr_0' `tpr_0' `mu_1' `nd' `ndb' `level' `kwidth' `options'
       post `mysim' (r(roc_f)) (r(rocf_lb)) (r(lrocf_lb)) (r(rocinvt)) (r(rocinvub)) (r(lrocinvub))
       local i = `i' + 1
    }
    postclose `mysim'
end
********************************************************
pro def doasim, rclass
    *version 6.0
    version 7
    args tpr1 fpr1 fpr0 tpr0 mu_1 nd ndb level kwidth /* kernal */
    macro shift 9
    tempname rocf rocinvt
    qui replace y = invnorm(uniform()) + d*`mu_1'

    roc y d `fpr0'
    scalar `rocf' = r(roc)
    qui replace ycut = r(ycut) if at
    cbgen `rocf' `fpr0' d `level' `nd' `ndb' y ycut `kwidth' `*'
    return scalar roc_f = `rocf'
    return scalar rocf_lb = r(roc_lb)
    return scalar lrocf_lb = r(lroc_lb)

    roc y d_inv `tpr0'
    scalar `rocinvt' = r(roc)
    qui replace ycut = r(ycut) if at
    cbgen `rocinvt' `tpr0' d_inv `level' `ndb' `nd' y ycut `kwidth' `*'
    return scalar rocinvt = `rocinvt'
    return scalar rocinvub = r(roc_ub)
    return scalar lrocinvub = r(lroc_ub)
end
***********************************************************
prog def roc,rclass
    /* calculate empirical roc for specified t
        add check on range of t (0<t<1)
    */
    args y d t
    local perc = 100 - (100*`t')
    qui _pctile `y' if `d'==0, p(`perc')
    return scalar ycut = r(r1)
    qui count if `y' > r(r1) & `d'==1
    local num = r(N)
    qui count if `d'==1
    return scalar roc = `num'/r(N)
end
********************************************************************************
prog def cbgen, rclass
    args roct t d level nd ndb y atvar kwidth kernal
    tempname z s
    tempvar kdx fdc kdx2 fdbc se
    scalar `z' = invnorm(.5 + `level'/200.)

    if `kwidth' == -1 { local kwidth "" }
    else {local kwidth `"width(`kwidth')"'}
    kdensity `y' if `d'==0, gen(`kdx' `fdbc') at(`atvar') nogr `kernal' `kwidth'
    kdensity `y' if `d'==1, gen(`kdx2' `fdc') at(`atvar') nogr `kernal' `kwidth'
    qui {
    replace `fdbc' = 1e-30 if `fdbc'==0
    replace `fdc' = 1e-30 if `fdc'==0
    gen double `se' = `roct'*(1.0 - `roct')/`nd' + ((`fdc'/`fdbc')^2)*(`t'*(1.0-`t')/`ndb')
    replace `se' = sqrt(`se')
    sum `se',meanonly
    assert r(N)==1
    sca `s' = r(min)
    return scalar roc_lb = `roct' - (`z' * `s')
    return scalar roc_ub = `roct' + (`z' * `s')

    * CI's based on logit ROC
    * set roct to values close to 0 and 1 when 0 and 1, respectively:
    *   (otherwise logit undefined)
    tempname logitroc lb ub roct_1
    if `roct' == 0 {
        scalar `roct_1' = .0001
    }
    else if `roct' == 1 {
        scalar `roct_1' = .9999
    }
    else {
        scalar `roct_1' = `roct'
    }
    scalar `logitroc' = log( (`roct_1')/(1.0 - `roct_1'))
    scalar `lb' = `logitroc' - `z' * (`s'/ (`roct_1'*(1 - `roct_1')))
    return scalar lroc_lb = exp(`lb')/(1.0 + exp(`lb'))
    scalar `ub' = `logitroc' + `z' *(`s'/ (`roct_1'*(1 - `roct_1')))
    return scalar lroc_ub = exp(`ub')/(1.0 + exp(`ub'))
    }
    /*  troubleshoot missing vals for logit-based bounds :
    if return(lroc_lb)==. {
        di "roct:     "`roct'
        di "logitroc: "`logitroc'  "           lroc_lb: " return(lroc_lb)
        di "s:       "`s'    "                 lroc_ub: " return(lroc_ub)
        di "lb:      "`lb'   "                  roc_lb: " return(roc_lb)
        di "exp(lb): "exp(`lb')  "              roc_ub: " return(roc_ub)
        di "fdc:  "`fdc'
        di "fdbc: "`fdbc'
        pause
    }
    if return(lroc_ub)==. {
        di "roct:     "`roct'
        di "logitroc: "`logitroc'  "           lroc_lb: " return(lroc_lb)
        di "s:       "`s'    "                 lroc_ub: " return(lroc_ub)
        di "ub:      "`ub'   "                  roc_lb: " return(roc_lb)
        di "exp(ub): "exp(`ub')  "              roc_ub: " return(roc_ub)
        di "fdc:  "`fdc'
        di "fdbc: "`fdbc'
        pause
    }
    */
end
