*! version 1.0.1  15jul2022
cap program drop tot_tut
program tot_tut, eclass
	version 17.0
	
	syntax varlist(min=3 max=3) [if] [in] [, vce(string)] 
	
	
	gettoken var rest : varlist
	gettoken Z choose : rest
	
	qui replace `choose' = 0 if missing(`choose')
	marksample touse
	
	*Check randomization range 0-2
	
	*Check choice range 0-1

	
	tempname Y X1 X0 W WPY WPX1i WPX0i theta1 theta0 theta1_0 theta xbhat1 xbhat0 U V Suu Svv Suv Svu WPU WPV cov1 cov0 cov1_0 vartottut cov
	tempvar x0 x1 z0_ z0 z1 clustervar
	
	
	*Cluster - robust
    if `"`vce'"' != "" {
        my_vce_parse , vce(`vce') 
        local vcetype     "robust"
        local clustervar  "`r(clustervar)'"
        if "`clustervar'" != "" {
            markout `touse' `clustervar'
        }
    }
	else {
		gen `clustervar' = _n if `touse'
		local robust = "robust"
	}
	
	sort `clustervar' 
	
	gen `x0' = -(`Z'==2)*(`choose'==0)
	gen `x1' = (`Z'==2)*(`choose'==1)
	gen `z0_' = -(`Z'==0)
	gen `z0' = (`Z'==0)
	gen `z1' = (`Z'==1)

	mata: mata clear
	qui putmata `Y' = `var' `X1' = (`x1' `z1' 1) `X0' = (`x0' `z0_' 1) `W' = (`z1' `z0' 1) if `touse'
	
	*'Projection' matrices
	mata : `WPY' = quadcross(`W',`Y')
	mata : `WPX1i' = pinv(quadcross(`W',`X1'))
	mata : `WPX0i' = pinv(quadcross(`W',`X0'))

	
		*ToT
	mata : `theta1' = (`WPX1i'*`WPY')
		*TuT
	mata : `theta0' = (`WPX0i'*`WPY')	
		*ToT-TuT
	mata : `theta1_0' = (1,0,0)*(`theta1' - `theta0')
	
	
	*Residuals
	mata : `U' = `Y' - `X1'*`theta1'
	mata : `V' = `Y' - `X0'*`theta0'	

	
	*Weighted Sum of Residuals
	mata : `Suu' = J(3,3,0)
	mata : `Svv' = J(3,3,0)
	mata : `Suv' = J(3,3,0)

	if "`robust'"=="robust" {
		qui count if `touse'
		local Nc = `r(N)'
		forvalues i = 1/`r(N)' {
			mata : `WPU' = `W'[`i',1..3]'*`U'[`i',1]
			mata : `WPV' = `W'[`i',1..3]'*`V'[`i',1]
			mata : `Suu' = `Suu' + `WPU'*`WPU''
			mata : `Svv' = `Svv' + `WPV'*`WPV''
			mata : `Suv' = `Suv' + `WPU'*`WPV''
		}
	}
	else {
		local k = 1
		qui levelsof `clustervar', local(levels) 
		local Nc = `r(r)'
		foreach l of local levels {
			*Number of obs in cluster
			qui count if `clustervar' == `l' & `touse'
			
			mata : `WPU' = J(3,1,0)
			mata : `WPV' = J(3,1,0)
			forvalues i = `k'/`=`k'+`r(N)'-1' {
				mata : `WPU' = `WPU' + `W'[`i',1..3]'*`U'[`i',1]
				mata : `WPV' = `WPV' + `W'[`i',1..3]'*`V'[`i',1]
			}
			mata : `Suu' = `Suu' + `WPU'*`WPU''
			mata : `Svv' = `Svv' + `WPV'*`WPV''
			mata : `Suv' = `Suv' + `WPU'*`WPV''
			
			local k = `k' + `r(N)'
		}
	}	

	*Covariance matrix
	mata : `cov0' = `WPX0i'*`Svv'*`WPX0i''
	mata : `cov1' = `WPX1i'*`Suu'*`WPX1i''
	mata : `cov1_0' = `WPX1i'*`Suv'*`WPX0i''
	mata : `cov'  = (`cov1' , `cov1_0' \ `cov1_0'', `cov0')
	
	mata : `vartottut' = (`WPX1i', -`WPX0i')*(`Suu', `Suv' \ `Suv'', `Svv')*(`WPX1i'' \ -`WPX0i'')
	mata : `vartottut' = (1,0,0)*`vartottut'*(1\0\0)
	
	*Rearranging/stacking
	mata : `cov' = (`cov'[2,2], `cov'[2,1], `cov'[2,4], `cov'[2,6], `cov'[2,3], -1	 \ ///
					`cov'[1,2], `cov'[1,1], `cov'[1,4], `cov'[1,6], `cov'[1,3], -1	 \ ///
					`cov'[4,2], `cov'[4,1], `cov'[4,4], `cov'[4,6], `cov'[4,3], -1	 \ ///
					`cov'[6,2], `cov'[6,1], `cov'[6,4], `cov'[6,6], `cov'[6,3], -1	 \ ///
					`cov'[3,2], `cov'[3,1], `cov'[3,4], `cov'[3,6], `cov'[3,3], -1	 \ ///
						-1	  , 	-1	  , 	-1	   , 	-1	    , 		-1	, `vartottut')
							
	mata : `theta' = (`theta1'[2,1], `theta1'[1,1], `theta0'[1,1], `theta1'[3,1], `theta0'[3,1], `theta1_0'[1,1])
	
	mata : st_matrix("`theta'", `theta')
	mata : st_matrix("`cov'", `cov')

	matrix colnames `theta' = ATE ToT TuT E[Y1] E[Y0] ToT-TuT
	matrix colnames `cov' = ATE ToT TuT E[Y1] E[Y0] ToT-TuT	
	matrix rownames `cov' = ATE ToT TuT E[Y1] E[Y0] ToT-TuT	


	*Display
	qui count if `touse'
	local N = `r(N)'
	di " "
	di "------------------------------------------------------------------------------"
	di "ToT & TuT 2sls stacked regression            	Number of obs   =    `r(N)'"	
	di " "
	if "`robust'"=="robust" {
	di "								(Robust std. err.)"
	}
	else {
	di "				(Std. err. adjusted for `Nc' clusters in suc_x_dia)"
	}	
	di " "
	ereturn post `theta' `cov', esample(`touse') buildfvinfo obs(`N') 
	ereturn scalar df_r = `Nc' - 1
    ereturn local vce "`vce'"
    ereturn local vcetype "`vcetype'"
    ereturn local clustvar "`clustvar'"
    ereturn local cmd "tot_tut"	
	ereturn display
end

cap program drop my_vce_parse
program define my_vce_parse, rclass
    syntax  [, vce(string) ]
 
    local case : word count `vce'
     
    if `case' > 2 {
        my_vce_error , typed(`vce')
    }
 
    local 0 `", `vce'"' 
    syntax  [, Robust CLuster * ]
 
    if `case' == 2 {
        if "`robust'" == "robust" | "`cluster'" == "" {
            my_vce_error , typed(`vce')
        }
 
        capture confirm numeric variable `options'
        if _rc {
            my_vce_error , typed(`vce')
        }
 
        local clustervar "`options'" 
    }
    else {    // case = 1
        if "`robust'" == "" {
            my_vce_error , typed(`vce')
        }
 
    }
 
    return clear    
    return local clustervar "`clustervar'" 
end
 
cap program drop my_vce_error 
program define my_vce_error
    syntax , typed(string)
 
    display `"{red}{bf:vce(`typed')} invalid"'
    error 498
end
