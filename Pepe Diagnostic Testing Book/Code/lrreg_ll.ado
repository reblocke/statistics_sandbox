*! version 1.0.0	08apr1997
*! likelihood function called by lrreg
program define lrreg_ll
	version 5.0
	
	* di in white "in lrr_mll"

	tempvar g1 g2  
	local lnf "`1'"
	local X1b1 "`2'"
	local X2b2 "`3'"

	quietly gen double `g1' = exp(`X1b1')
	quietly gen double `g2' = exp(`X2b2')

	local y : word 1 of $S_mldepn
	local D : word 2 of $S_mldepn

*	set more on
*	sort `y' `D' `g1' `g2'
*	by `y' `D': l `X1b1' `X2b2' `g1' `g2' 

	#delimit ;
   	quietly replace `lnf' =  cond( `y', 
	  cond(`D', ln( ((`g2'- 1)*(`g1'))/(`g2'-`g1')),
   	        ln((`g2'- 1)/(`g2'-`g1')) ), 
   	  cond(`D' ,ln( (`g2'*(1-`g1'))/(`g2'-`g1') ), 
	         ln( (1-`g1')/(`g2'-`g1') ) )    ) ;
	#delimit cr

end
