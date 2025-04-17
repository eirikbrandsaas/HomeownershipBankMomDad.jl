clear all
cap graph drop _all

************* 
** Share of households who remain as owners
************* 

use "$PSIDpath/regression_sample.dta", clear
tempfile leadmeans
tempname  post

postfile `post' int fv                       /// lead: 0,2,?,16
                     double m_first          /// mean for age==first_own
                     double m_owner          /// mean for current owners
                     using `leadmeans', replace
		     
qui forvalues fv = 0(2)16 {
	sum F(`fv').owner if age == first_own & first_own != firstage & inrange(age,25,44) [aw = famwgt]
	local mean_first = `r(mean)'
	sum F(`fv').owner if owner == 1 & first_own != firstage & inrange(age,25,44) [aw = famwgt]
	local mean_owner = `r(mean)'
	post `post' (`fv') (`mean_first') (`mean_owner')
}
postclose `post'      // finish writing

use `leadmeans', clear
twoway  (line  m_first m_owner  fv, lcolor(black orange) lpattern(solid solid) msymbol(square circle)) , ///
	 legend(order(1 "First-Time Owners" 2 "All Owners") ring(0) pos(1) cols(1))     ///
	graphregion(margin(1 1 0 0))                 /// shrink outer frame
         xlabel(0(4)16) xtitle("Years (leads)") ytitle("Share") name(hazard)

graph display, xsize(2) ysize(2) scale(1.8) 
graph export "tabfig/descr/PSID_ownerexit.pdf", replace


************* 
** Impact of wealth and parental wealth on ownership outcomes
************* 

use "$PSIDpath/regression_sample.dta", clear
tempfile results
tempname  post

postfile `post' int fv ///
        double b_hh seb_hh b_par seb_par ///
	long Nb ///
	double c_hh sec_hh c_par sec_par ///
	long Nc ///
        using `results', replace

forvalues fv = 0(2)16 {
    
    
	local controls "c.wealth c.income i.hs i.coll i.white i.married i.year i.state c.age##c.age c.famsize c.age_prnt##c.age_prnt"
	local lag_main "cashonhand_prnt "

	qui regr F(`fv').owner `lag_main' `controls' if owner == 1 & inrange(age,25,44)
    * write fv, coeffs and s.e.'s of B and C
	tempname bB sebB bparB sebparB NB                // stash the numbers
	
        scalar   `bB' = _b[wealth]
        scalar   `sebB' = _se[wealth]
	scalar   `bparB' = _b[cashonhand_prnt]
        scalar   `sebparB' = _se[cashonhand_prnt]
	scalar  `NB' =  `e(N)'
	
	qui regr F(`fv').owner `lag_main' `controls' if age == first_own & first_own != firstage & inrange(age,25,44)
	post `post' (`fv')  /// 
	(scalar(`bB')) (scalar(`sebB')) (scalar(`bparB')) (scalar(`sebparB')) (scalar(`NB')) ///
	(_b[wealth]) (_se[wealth]) (_b[cashonhand_prnt]) (_se[cashonhand_prnt]) (`e(N)')
}

postclose `post'

*------------------------------------------------------------
* 3.  Bring the new dataset into memory
*------------------------------------------------------------
use `results', clear
list, clean

foreach mod in "b" "c" {
	gen l`mod'_hh = `mod'_hh - 1.68*se`mod'_hh
	gen u`mod'_hh = `mod'_hh + 1.68*se`mod'_hh
	gen l`mod'_par = `mod'_par- 1.68*se`mod'_par
	gen u`mod'_par = `mod'_par + 1.68*se`mod'_par

}
cap graph drop _all
twoway  (line b_hh ub_hh  lb_hh fv, lcolor(black black black) lpattern(solid dash  dash) msymbol(square none none)) ///
		(line b_par ub_par  lb_par fv, lcolor(orange orange orange) lpattern(solid dash dash)  msymbol(circle none none)), ///
		 graphregion(margin(1 1 0 0))                 /// shrink outer frame
		 legend(order(1 "Household Wealth" 4 "Parent Wealth")    ///
           ring(0) pos(5) cols(1))  xlabel(0(4)16) ylabel(-0.1(0.05)0.1) xtitle("Lead Length") ytitle("Coefficient") name(all)
	   
graph display, xsize(2) ysize(2) scale(1.8) 
graph export "tabfig/descr/PSID_coefowner.pdf", replace

	
twoway  (line c_hh uc_hh  lc_hh fv, lcolor(black black black) lpattern(solid dash  dash) msymbol(square none none)) ///
		(line c_par uc_par  lc_par fv, lcolor(orange orange orange) lpattern(solid dash dash)  msymbol(circle none none)), ///
		 graphregion(margin(1 1 0 0))                 /// shrink outer frame
		 legend(order(1 "Household Wealth" 4 "Parent Wealth")    ///
           ring(0) pos(5) cols(1)) xlabel(0(4)16)  ylabel(-0.1(0.05)0.1)  xtitle("Lead Length") ytitle("Coefficient") name(ft)
graph display, xsize(2) ysize(2) scale(1.8) 
graph export "tabfig/descr/PSID_coefftowner.pdf", replace


   graph combine all ft, ycommon

