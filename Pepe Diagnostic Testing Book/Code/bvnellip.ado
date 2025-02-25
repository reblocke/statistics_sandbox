*! 1.1.0 GML May 24 2000
*! Joint Confidence ellipse for 2 parameters with bivariate normal distrib.
prog def bvnellip,rclass
    version 6
    syntax [varlist(default=none max=2)] [if] [in], PLotvrs(string) [SORtvar(string) /*
        */ Level(integer $S_level)                            /*
        */ COVar(string) Beta(string) DF2(integer 0) NPoints(integer 200) /*
        */ REPLACE ]

    if "`varlist'"=="" & ("`covar'"=="" | "`beta'"=="" | `df2'==0) {
        di in red "either varlist OR beta(), covar() and df2() options must be specified"
        exit 198
    }
    tempname F beta1 sbeta1 beta2 sbeta2 rho
    tempvar theta new1 new2

    ***** 1) varlist option:

    if "`varlist'"~="" {
        if "`beta'" ~= "" | "`covar'" ~= "" {
            di in red "beta() and covar() options cannot be specified if a varlist is specified"
            exit 198
        }
        local varct: word count `varlist'
        if `varct' < 2 {
            di in red "TWO variables must be included in varlist"
            exit 198
        }
        tokenize `varlist'
        local x1 "`1'"
        local x2 "`2'"
        marksample touse
        markout `touse' `x1' `x2'
        qui sum `x1' if `touse'
        local nobs = r(N)
        scalar `beta1' = r(mean)
        scalar `sbeta1'  = r(sd)/sqrt(`nobs')
        qui sum `x2' if `touse'
        scalar `beta2' = r(mean)
        scalar `sbeta2'  = r(sd)/sqrt(`nobs')
        qui cor `x1' `x2' if `touse'
        scalar `rho' = r(rho)
        local df2 = `nobs' - 2
    }
    ***** 2) beta,covar arg options
    else {
        tempname CorM
        if  ("`covar'"=="" | "`beta'"=="" | `df2'==0) {
            di in red "beta(), covar() and df2() options must all be specified
            di in red "  when varlist is not."
            exit 198
        }
        if (rowsof(`covar')~=2 | colsof(`covar')~=2) {
            di in red "Covariance matrix must be 2x2"
            exit 198
        }
        if (rowsof(`beta')~=2 | colsof(`beta')~=1) {
            di in red "beta matrix (vector) must be 2x1"
            exit 198
        }
        scalar `beta1' = `beta'[1,1]
        scalar `beta2' = `beta'[2,1]
        scalar `sbeta1' = sqrt(`covar'[1,1])
        scalar `sbeta2' = sqrt(`covar'[2,2])
        matr `CorM' = corr(`covar')
        scalar `rho' = `CorM'[2,1]
    }
    *********
    scalar `F' = invfprob(2,`df2',(100-`level')/100)
    qui range `theta' 0 2*_pi `npoints'
    qui gen `new1' = `beta1' + `sbeta1'*sqrt(2*`F')*cos(`theta'+acos(`rho'))
    qui gen `new2' = `beta2' + `sbeta2'*sqrt(2*`F')*cos(`theta')

    if "`replace'"!=`""' {
        cap drop `plotvrs'
        cap drop `sortvar'
    }

    local wc : word count `plotvrs'
    cap assert `wc' == 2
    if _rc~=0 {
        di in red "must have 2 newvar names in plotvrs( )"
        error _rc
    }
    tokenize `plotvrs'
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
    rename `new1' `1'
    rename `new2' `2'

    if "`sortvar'" ~= "" {
        local wc : word count `sortvar'
        cap assert `wc' == 1
            if _rc~=0 {
                di in red "must have 1 newvar name in sortvar( )"
                error _rc
            }
        cap confirm new variable `sortvar'
        if _rc~=0 {
            di in red "variable `sortvar' already exists,"
            di in red `"  use "replace" option to replace `sortvar' "'
            error _rc
        }
        rename `theta' `sortvar'
    }
    return scalar beta2 = `beta2'
    return scalar beta1 = `beta1'
end
