*! version 1.0.2    29apr1999
* 1.0.1 did not work with version 6
*   changing sens and spec from tempvars to scalars
*    with the temporary names appeared to fix this.
program define lrreg
    version 5.0
    pause on
    local options "Level(integer $S_level) Eform"
    if substr("`1'",1,1)==","| "`*'"=="" {
        if "$S_E_cmd"~="lrreg" {
            error 301
        }
        parse "`*'"
        if "`eform'" ~= "" { local eform "eform(DLR ratio)" }
        ml mlout lrreg, level(`level') `eform'
        exit
    }
    parse "`*'", parse(" ,")

    capture confirm variable `1'
    if _rc~=0 {
        di in red "`1' found where binary test outcome variable expected"
        di in red "`1' is not an existing variable"
        exit
    }
    capture assert `1'==0 | `1'==1 | `1'==.
    if _rc~=0 {
        di in r "test outcome variable " in g "`1'" in r " is not 0/1"
        exit
    }
    capture confirm variable `2'
    if _rc~=0 {
        di in red "`2' found where binary disease variable expected"
        di in red "`2' is not an existing variable"
        exit
    }
    capture assert `2'==0 | `2'==1 | `2'==.
    if _rc~=0 {
        di in r "disease variable " in g "`2'" in r " is not 0/1"
        exit
    }

    local y "`1'"
    local d "`2'"

    capture confirm existence `4'
    if _rc~=0 | substr("`4'",1,1)=="," {
        eq ? `3'
        eq LRpos: $S_1
        eq ? `3'
        eq LRneg: $S_1
        local lrp LRpos
        local lrn LRneg
        macro shift 3
    }
    else {
        local lrp "`3'"
        local lrn "`4'"
        mac shift 4
    }
    eq ? `lrp'
    local varlrp $S_1
    eq ? `lrn'
    local varlrn $S_1

    local if "opt"
    local in "opt"
    local options "`options' Robust CLuster(string) noCONs *"

    parse "`*'"

    if "`eform'" ~= "" { local eform "eform(DLR Ratio)" }

    tempvar touse mltouse
    tempname B0 ll D sens spec

/* Determine estimation sample. */
    mark `touse' `if' `in'
    markout `touse' `y' `d' `varlrp' `varlrn'
    if "`cluster'"~="" {
        confirm variable `cluster'
        markout `touse' `cluster', strok
        local clopt "cluster(`cluster')"
    }
/* Get starting values */

    matrix `B0' = J(1,2,0)
    matrix coleq `B0' = `lrp' `lrn'

    qui summ `y' if `d'==1 & `touse'
    sca `sens' = _result(3)
    qui summ `y' if `d'==0 & `touse'
    sca `spec' = 1-_result(3)

/*  a) for model with constant  */

    if "`cons'" == "" {
        matrix colnames `B0' = _cons _cons
        matrix `B0'[1,1] = ln(`sens'/(1-`spec'))
        matrix `B0'[1,2] = ln((1-`sens')/`spec')
    }
/*  b) for model without constant */

    else {
        local c1var : word 1 of `varlrp'
        local c2var : word 1 of `varlrn'
        matrix colnames `B0' = `c1var' `c2var'
        sum `c1var',meanonly
        mat `B0'[1,1] = ln(`sens'/(1-`spec')) /* /_result(3) */
        sum `c2var',meanonly
        mat `B0'[1,2] = ln((1-`sens')/`spec')  /* /_result(3) */
    }
    * di "for cons = `cons'"
    * mat list `B0'
    * pause
    * display in white "checkpoint 2"

    * di in white "beginning of ml commands"

/* Do ml */
    ml begin
    ml func lrreg_ll
    ml method lf

    * di in white "checkpoint 2b"

    ml model b = `lrp' `lrn', depv(00) from(`B0') `cons'

    * di in white "checkpoint 3"
    * set trace on

    ml sample `mltouse' if `touse'
    ml depn `y' `d'

*   set more on
*   set trace on

    ml max `ll' `D', `options'

    * di in white "checkpoint 4"

/* reset depname in order to trick mlout into displaying
   the dep var. name (will do so if only a single name
   exists in $S_depnam */

    global S_mldepn : word 2 of $S_mldepn

    ml post lrreg, title(Likelihood Ratio Regression)
    if "`robust'"~="" | "`cluster'"~="" {
        tempvar s1 s2
        lrr_scor `y' `d' `lrp' `lrn' `s1' `s2'
    _robust `s1' `s2' if `touse', `clopt'
    }

    * di in wh "checkpoint 5 "

    ml mlout lrreg, level(`level') `eform'
    * cap eq drop LRpos LRneg
end

program define lrr_scor
    version 5.0

    * di in white "made it to prog scores"
    local y "`1'"
    local d "`2'"
    local lrp "`3'"
    local lrn "`4'"
    local s1 "`5'"
    local s2 "`6'"
    tempvar x1b x2b g1 g2 dg1 dg2

        quietly predict double `x1b', equation(`lrp')
        quietly predict double `x2b', equation(`lrn')
        quietly gen double `g1' = exp(`x1b')
        quietly gen double `g2' = exp(`x2b')
        quietly gen double `dg1' = `g1'
        quietly gen double `dg2' = `g2'

        #delimit ;
        quietly gen double `s1' = cond(`y',
            cond(`d', (`dg1'*`g2')/(`g1'*(`g2'-`g1')),
                (`dg1')/(`g2'-`g1') ) ,
                 ((`dg1')*(1-`g2'))/((1-`g1')*(`g2'-`g1')) );
        quietly gen double `s2' = cond(`y',
            (`dg2'*(1-`g1'))/((`g2'-1)*(`g2'-`g1')) ,
                cond(`d', (-`dg2'*`g1')/(`g2'*(`g2'-`g1')),
                (-`dg2')/(`g2'-`g1') ) );
        #delimit cr
end
