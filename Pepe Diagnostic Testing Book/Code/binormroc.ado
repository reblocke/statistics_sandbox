*! 1.0.2 GML 20 Nov 2001
prog def binormroc
    version 7.0
    syntax,  mu_db(real) s_db(real) mu_d(real) s_d(real) [ NPoints(int 500) /*
       */  GENRoc(string) GENt(string) REPLACE]

    if `"`genroc'"'=="" & `"`gent'"'=="" { preserve }

    tempfile g1 g2 gnull

    cap assert `s_db' > 0 & `s_d' > 0
     if _rc~=0 {
         di in red "numeric arguments for s_db() and s_d() must be positive"
         exit
     }

    if "`replace'"!=`""' {
         cap drop `genroc'
         cap drop `gent'
    }

    tempvar y1 y2 t roc x forty5

    local xmin = `mu_db' - 3.5 * `s_db'
    local xmax = `mu_d' + 3.5 * `s_d'

    qui range `x' `xmin' `xmax' `npoints'

    qui normd `x', d(`y1') m(`mu_db') s(`s_db') replace
    qui normd `x', d(`y2') m(`mu_d') s(`s_d') replace

    set textsize 175
    #delimit ;
    qui gr `y1' `y2' `x', s(..) c(..)
       xlabel xtick
       t1("mu_db = `mu_db',    s_db = `s_db'")
       t2(" mu_d = `mu_d',     s_d = `s_d'")
       saving(`g1')
       ;

    #delimit cr
    * after density curves:

    tempname a b
    sca `a' = (`mu_d' - `mu_db')/`s_d'
    sca `b' = `s_db'/`s_d'

    qui {
       range `t' 0 1 `npoints'
       gen `roc' = norm(`a' + `b'*invnorm(`t'))
       replace `roc' = `t' if `t'==0 | `t'==1
       sort `t'
       gen `forty5' = `t' if `t'==0 | `t'==1
    }

    local bstr = string(round(`b',.01))
    local astr = string(round(`a',.01))

    qui gr using , saving(`gnull')
    #delimit ;
    qui gr `roc' `forty5' `t',
      c(sl) s(ii)
      xscale(0,1) yscale(0,1)
      xlabel(0,.2,.4,.6,.8,1) ylabel(0,.2,.4,.6,.8,1)
      xtick(.1,.3,.5,.7,.9) ytick(.1,.3,.5,.7,.9)
      /* xline(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)
         yline(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1) */
      xline(0,1) yline(0,1) noaxis
      b1("FP rate") b2(" ")
      l1("TP rate")
      t1("binormal ROC; a = `astr', b = `bstr' ")
      saving(`g2')
    ;
    #delimit cr
    set textsize 100
    gr using `g2' `gnull' `g1', margin(5)

    if "`genroc'"!=`""' {
        local gvarnum : word count `genroc'
        cap assert `gvarnum' == 1
        if _rc~=0 {
            di in red "more than one newvar name in genroc( ) option"
            error _rc
        }
        cap confirm new variable `genroc'
        if _rc~=0 {
            di in red "variable `genroc' already exists,"
            di in red `"  use "replace" option to replace `genroc'"'
            error _rc
        }
        rename `roc' `genroc'
    }
    if "`gent'"!=`""' {
        local gvarnum : word count `gent'
        cap assert `gvarnum' == 1
        if _rc~=0 {
            di in red "more than one newvar name in gent( ) option"
            error _rc
        }
        cap confirm new variable `gent'
        if _rc~=0 {
            di in red "variable `gent' already exists,"
            di in red `"  use "replace" option to replace `gent'"'
            error _rc
        }
        rename `t' `gent'
    }
    if `"`genroc'"'=="" & `"`gent'"'=="" { restore }
end
pro def normd
    version 7
    /* this was already written as a stand alone, hence the syntax overkill */
    syntax varlist(min=1 max=1), Dens(string) [ Mu(real 0) Sigma(real 1) REPLACE ]
    local x "`varlist'"
    cap assert `sigma'>0
    if _rc~=0 {
        di in red "arg to sigma() must be > 0"
        error _rc
    }
    if `"`replace'"'!=`""' {
        cap drop `dens'
    }
    else {
        cap confirm new variable `dens'
        if _rc~=0 {
            di in red "variable `dens' already exists,"
            di in red `"  use "replace" option to replace `dens'"'
            error _rc
        }
    }
    gen `dens' = (1.0/(sqrt(2.0*_pi)*`sigma')) * exp(-((`x'-`mu')^2)/(2.0*(`sigma'^2)))
end
/*
 future modifications:

  - need to allow npoints less than than no. obs in current dataset -  use -mark- cmds to restrict?

  - alternative specification of binormal parameters a,b ?

*/


