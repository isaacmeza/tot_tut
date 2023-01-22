*! version 1.0.3  19jan2023
cap program drop tot_tut
program tot_tut, eclass
	version 17.0
	
	syntax varlist(min=3 max=3) [if] [in] [, vce(string)] 
	
	
	gettoken var rest : varlist
	gettoken Z choose : rest

	
	*Check randomization range 0-2
	qui levelsof `Z'
	if r(levels)!="0 1 2" {
		display "Randomization range outside of 0-1-2"
		exit
	}
	*Check choice range 0-1
	qui levelsof `choose'
	if r(levels)!="0 1" {
		display "Not a binary choice variable"
		exit
	}
	
	tempname Y Ybias Ylevel X1 X0 Xb1 Xb0 Xl1 Xl0 W Wb1 Wb0 Wl1 Wl0 WPY WPX1i WPX0i Wb1PYb Wb1PXb1i Wb0PYb Wb0PXb0i Wl1PYl Wl1PXl1i Wl0PYl Wl0PXl0i theta1 theta0 theta1_0 theta asb asl U V Ub1 Ub0 Ul1 Ul0 Suu Svv Suv Svu Sb1 Sb0 Sub Sl1 Sl0 Sul WPU WPV Wb1PUb1 Wb0PUb0 Wl1PUl1 Wl0PUl0 cov1 cov0 cov1_0 vartottut varasb varasl cov
	tempvar choose_ ddd x0 x0_ x1 x1_ z0 z0_ z1 z2 uno_z1 uno_z0 yb yl clustervar
	
	qui gen `choose_' = `choose'
	qui replace `choose' = 0 if missing(`choose')
	marksample touse

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
		qui gen `clustervar' = _n if `touse'
		local robust = "robust"
	}
	
	sort `clustervar' 
	
	gen `ddd' = (`Z'==1) | (`choose'==1)
	gen `x0' = -(`Z'==2)*(`choose'==0)
	gen `x0_' = -`x0'
	gen `x1' = (`Z'==2)*(`choose'==1)
	gen `x1_' = -`x1'
	gen `z0' = (`Z'==0)
	gen `z0_' = -`z0'
	gen `z1' = (`Z'==1)
	gen `z2' = (`Z'==2)
	gen `uno_z1' = (1-`z1')
	gen `uno_z0' = (1-`z0')


	*Convert dep var for selection outcomes
	gen `yb' = `var'*(1-`z1')*(1-`ddd')
	gen `yl' = `var'*(1-`z0')*(`ddd')


	mata: mata clear
	qui putmata `Y' = `var' `X1' = (`x1' `z1' 1) `X0' = (`x0' `z0_' 1) `W' = (`z1' `z0' 1) `Ybias' = `yb' `Xb1' = (`x1_' `uno_z1') `Xb0' = (`x0_' `z0') `Wb1' = (`z2' `uno_z1') `Wb0' = (`z2' `z0') `Ylevel' = `yl' `Xl1' = (`x1' `z1') `Xl0' = (`x0' `uno_z0') `Wl1' = (`z2' `z1') `Wl0' = (`z2' `uno_z0') if `touse'
	

	*'Projection' matrices
	mata : `WPY' = quadcross(`W',`Y')
	mata : `WPX1i' = pinv(quadcross(`W',`X1'))
	mata : `WPX0i' = pinv(quadcross(`W',`X0'))

	mata : `Wb1PYb' = quadcross(`Wb1',`Ybias')
	mata : `Wb1PXb1i' = pinv(quadcross(`Wb1',`Xb1'))

	mata : `Wb0PYb' = quadcross(`Wb0',`Ybias')
	mata : `Wb0PXb0i' = pinv(quadcross(`Wb0',`Xb0'))

	mata : `Wl1PYl' = quadcross(`Wl1',`Ylevel')
	mata : `Wl1PXl1i' = pinv(quadcross(`Wl1',`Xl1'))

	mata : `Wl0PYl' = quadcross(`Wl0',`Ylevel')
	mata : `Wl0PXl0i' = pinv(quadcross(`Wl0',`Xl0'))


	
		*ToT
	mata : `theta1' = (`WPX1i'*`WPY')
		*TuT
	mata : `theta0' = (`WPX0i'*`WPY')	
		*ToT-TuT
	mata : `theta1_0' = (1,0,0)*(`theta1' - `theta0')
	
		*ASB
	mata : `asb' = (1,0)*(`Wb1PXb1i'*`Wb1PYb'-`Wb0PXb0i'*`Wb0PYb')
		*ASB
	mata : `asl' = (1,0)*(`Wl1PXl1i'*`Wl1PYl'-`Wl0PXl0i'*`Wl0PYl')


	*Residuals
	mata : `U' = `Y' - `X1'*`theta1'
	mata : `V' = `Y' - `X0'*`theta0'	

	mata : `Ub1' = `Ybias' - `Xb1'*`Wb1PXb1i'*`Wb1PYb'
	mata : `Ub0' = `Ybias' - `Xb0'*`Wb0PXb0i'*`Wb0PYb'
	mata : `Ul1' = `Ylevel' - `Xl1'*`Wl1PXl1i'*`Wl1PYl'
	mata : `Ul0' = `Ylevel' - `Xl0'*`Wl0PXl0i'*`Wl0PYl'
	
	*Weighted Sum of Residuals
	mata : `Suu' = J(3,3,0)
	mata : `Svv' = J(3,3,0)
	mata : `Suv' = J(3,3,0)

	mata : `Sb1' = J(2,2,0)
	mata : `Sb0' = J(2,2,0)	
	mata : `Sub' = J(2,2,0)
	mata : `Sl1' = J(2,2,0)
	mata : `Sl0' = J(2,2,0)		
	mata : `Sul' = J(2,2,0)


	if "`robust'"=="robust" {
		qui count if `touse'
		local Nc = `r(N)'
		forvalues i = 1/`r(N)' {
			mata : `WPU' = `W'[`i',1..3]'*`U'[`i',1]
			mata : `WPV' = `W'[`i',1..3]'*`V'[`i',1]
			mata : `Suu' = `Suu' + `WPU'*`WPU''
			mata : `Svv' = `Svv' + `WPV'*`WPV''
			mata : `Suv' = `Suv' + `WPU'*`WPV''

			mata : `Wb1PUb1' = `Wb1'[`i',1..2]'*`Ub1'[`i',1]
			mata : `Wb0PUb0' = `Wb0'[`i',1..2]'*`Ub0'[`i',1]
			mata : `Sb1' = `Sb1' + `Wb1PUb1'*`Wb1PUb1''
			mata : `Sb0' = `Sb0' + `Wb0PUb0'*`Wb0PUb0''			
			mata : `Sub' = `Sub' + `Wb1PUb1'*`Wb0PUb0''

			mata : `Wl1PUl1' = `Wl1'[`i',1..2]'*`Ul1'[`i',1]
			mata : `Wl0PUl0' = `Wl0'[`i',1..2]'*`Ul0'[`i',1]
			mata : `Sl1' = `Sl1' + `Wl1PUl1'*`Wl1PUl1''
			mata : `Sl0' = `Sl0' + `Wl0PUl0'*`Wl0PUl0''	
			mata : `Sul' = `Sul' + `Wl1PUl1'*`Wl0PUl0''			
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

			mata : `Wb1PUb1' = J(2,1,0)
			mata : `Wb0PUb0' = J(2,1,0)
			mata : `Wl1PUl1' = J(2,1,0)
			mata : `Wl0PUl0' = J(2,1,0)			

			forvalues i = `k'/`=`k'+`r(N)'-1' {
				mata : `WPU' = `WPU' + `W'[`i',1..3]'*`U'[`i',1]
				mata : `WPV' = `WPV' + `W'[`i',1..3]'*`V'[`i',1]

				mata : `Wb1PUb1' = `Wb1PUb1' + `Wb1'[`i',1..2]'*`Ub1'[`i',1]
				mata : `Wb0PUb0' = `Wb0PUb0' + `Wb0'[`i',1..2]'*`Ub0'[`i',1]

				mata : `Wl1PUl1' = `Wl1PUl1' + `Wl1'[`i',1..2]'*`Ul1'[`i',1]
				mata : `Wl0PUl0' = `Wl0PUl0' + `Wl0'[`i',1..2]'*`Ul0'[`i',1]
			}

			mata : `Suu' = `Suu' + `WPU'*`WPU''
			mata : `Svv' = `Svv' + `WPV'*`WPV''
			mata : `Suv' = `Suv' + `WPU'*`WPV''

			mata : `Sb1' = `Sb1' + `Wb1PUb1'*`Wb1PUb1''
			mata : `Sb0' = `Sb0' + `Wb0PUb0'*`Wb0PUb0''		
			mata : `Sub' = `Sub' + `Wb1PUb1'*`Wb0PUb0''

			mata : `Sl1' = `Sl1' + `Wl1PUl1'*`Wl1PUl1''
			mata : `Sl0' = `Sl0' + `Wl0PUl0'*`Wl0PUl0''				
			mata : `Sul' = `Sul' + `Wl1PUl1'*`Wl0PUl0''			
			
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

	mata : `varasb' = (`Wb1PXb1i', -`Wb0PXb0i')*(`Sb1', `Sub' \ `Sub'', `Sb0')*(`Wb1PXb1i'' \ -`Wb0PXb0i'')
	mata : `varasb' = (1,0)*`varasb'*(1\0)
	mata : `varasl' = (`Wl1PXl1i', -`Wl0PXl0i')*(`Sl1', `Sul' \ `Sul'', `Sl0')*(`Wl1PXl1i'' \ -`Wl0PXl0i'')
	mata : `varasl' = (1,0)*`varasl'*(1\0)
	
	*Rearranging/stacking
	mata : `cov' = (`cov'[2,2], `cov'[2,1], `cov'[2,4], `cov'[2,6], `cov'[2,3], 0 , 0 , 0	 \ ///
					`cov'[1,2], `cov'[1,1], `cov'[1,4], `cov'[1,6], `cov'[1,3], 0 , 0 , 0	 \ ///
					`cov'[4,2], `cov'[4,1], `cov'[4,4], `cov'[4,6], `cov'[4,3], 0 , 0 , 0	 \ ///
					`cov'[6,2], `cov'[6,1], `cov'[6,4], `cov'[6,6], `cov'[6,3], 0 , 0 , 0	 \ ///
					`cov'[3,2], `cov'[3,1], `cov'[3,4], `cov'[3,6], `cov'[3,3], 0 , 0 , 0	 \ ///
						0	  , 	0	  , 	0	   , 	0	    , 		0	, `vartottut' , 0 , 0 \ ///
						0	  , 	0	  , 	0	   , 	0	    , 		0	, 0 , `varasb' , 0 \ ///
						0	  , 	0	  , 	0	   , 	0	    , 		0	, 0 , 0 , `varasl')
							
	mata : `theta' = (`theta1'[2,1], `theta1'[1,1], `theta0'[1,1], `theta1'[3,1], `theta0'[3,1], `theta1_0'[1,1], `asb', `asl')
	
	mata : st_matrix("`theta'", `theta')
	mata : st_matrix("`cov'", `cov')
	
	matrix colnames `theta' = ATE ToT TuT E[Y1] E[Y0] ToT-TuT ASB ASL
	matrix colnames `cov' = ATE ToT TuT E[Y1] E[Y0] ToT-TuT	ASB ASL
	matrix rownames `cov' = ATE ToT TuT E[Y1] E[Y0] ToT-TuT	ASB ASL

	*Renormalization
	qui replace `choose' = `choose_'

	*Display
	qui count if `touse'
	local N = `r(N)'
	di " "
	di "------------------------------------------------------------------------------"
	di "ToT & TuT 2sls stacked regression            	Number of obs   =    `r(N)'"	
	di " "
	if "`robust'"=="robust" {
	di "						(Robust std. err.)"
	}
	else {
	di "						(Std. err. adjusted for `Nc' clusters)"
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
