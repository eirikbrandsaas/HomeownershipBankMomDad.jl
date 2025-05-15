clear all
cap graph drop _all

************* 
** Share of households who remain as owners
************* 

use "$PSIDpath/regression_sample.dta", clear
tempfile leadmeans
tempname  post

postfile `post' int fv                       /// lead: 0,2,?,16
                     double m_first se_first         /// mean for age==first_own
                     double m_owner se_owner         /// mean for current owners
                     using `leadmeans', replace
		     
qui forvalues fv = 0(2)16 {
	sum F(`fv').owner if age == first_own & first_own != firstage & inrange(age,25,44) [aw = famwgt]
	local mean_first = `r(mean)'
	local se_first   = r(sd)/sqrt(r(N))      
	sum F(`fv').owner if owner == 1 & first_own != firstage & inrange(age,25,44) [aw = famwgt]
	local mean_owner = `r(mean)'
	local se_owner   = r(sd)/sqrt(r(N))      
	post `post' (`fv') (`mean_first') (`se_first') (`mean_owner') (`se_owner')
}
postclose `post'      // finish writing

use `leadmeans', clear
foreach mod in "owner" "first" {
	gen m_`mod'_l = m_`mod' - 1.68*se_`mod'
	gen m_`mod'_u = m_`mod' + 1.68*se_`mod'


}

twoway  (line  m_first m_first_u m_first_l m_owner m_owner_u m_owner_l  fv, lcolor(black black black orange orange orange ) lpattern(solid dash dash solid dash dash ) msymbol(square circle)) , ///
	 legend(order(1 "First-Time Owners" 4 "All Owners") ring(0) pos(1) cols(1))     ///
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
        double all_k allse_k all_p allse_p ///
	long Nb ///
	double first_k firstse_k first_p firstse_p ///
	long Nc ///
        using `results', replace

forvalues fv = 2(2)16 {
    
    
 	local controls "c.wealth c.income i.hs i.coll i.married i.white c.famsize i.famchange i.year i.state c.age##c.age c.age_prnt##c.age_prnt"
	local lag_main "cashonhand_prnt "

	qui regr F(`fv').owner `lag_main' `controls' if owner == 1 & inrange(age,25,44) & first_own != firstage, vce(cluster famid )
	noi disp `e(N)'

	tempname bB sebB bparB sebparB NB                // stash the numbers
	
        scalar   `bB' = _b[wealth]
        scalar   `sebB' = _se[wealth]
	scalar   `bparB' = _b[cashonhand_prnt]
        scalar   `sebparB' = _se[cashonhand_prnt]
	scalar  `NB' =  `e(N)'
	
	qui regr F(`fv').owner `lag_main' `controls' if age == first_own & first_own != firstage & inrange(age,25,44), vce(cluster famid )
	post `post' (`fv')  /// 
	(scalar(`bB')) (scalar(`sebB')) (scalar(`bparB')) (scalar(`sebparB')) (scalar(`NB')) ///
	(_b[wealth]) (_se[wealth]) (_b[cashonhand_prnt]) (_se[cashonhand_prnt]) (`e(N)')
}

postclose `post'


// Plot 

use `results', clear
list, clean

foreach mod in "all" "first" {
	gen `mod'_l_k = `mod'_k - 1.68*`mod'se_k
	gen `mod'_u_k = `mod'_k + 1.68*`mod'se_k
	gen `mod'_l_p = `mod'_p- 1.68*`mod'se_p
	gen `mod'_u_p = `mod'_p + 1.68*`mod'se_p

}
cap graph drop _all
twoway  (line all_k all_u_k  all_l_k fv, lcolor(black black black) lpattern(solid dash  dash) msymbol(square none none)) ///
		(line all_p all_u_p  all_l_p fv, lcolor(orange orange orange) lpattern(solid dash dash)  msymbol(circle none none)), ///
		 graphregion(margin(1 1 0 0))                 /// shrink outer frame
		 legend(order(1 "Household Net Worth" 4 "Parent Wealth")    ///
           ring(0) pos(5) cols(1))  xlabel(0(4)16) ylabel(-0.1(0.05)0.1) xtitle("Lead Length") ytitle("Coefficient") name(all)
	   
graph display, xsize(2) ysize(2) scale(1.8) 
graph export "tabfig/descr/PSID_coefowner.pdf", replace

	
twoway  (line first_k first_u_k  first_l_k fv, lcolor(black black black) lpattern(solid dash  dash) msymbol(square none none)) ///
		(line first_p first_u_p  first_l_p fv, lcolor(orange orange orange) lpattern(solid dash dash)  msymbol(circle none none)), ///
		 graphregion(margin(1 1 0 0))                 /// shrink outer frame
		 legend(order(1 "Household Net Worth" 4 "Parent Wealth")    ///
           ring(0) pos(5) cols(1)) xlabel(0(4)16)  ylabel(-0.1(0.05)0.1)  xtitle("Lead Length") ytitle("Coefficient") name(ft)
graph display, xsize(2) ysize(2) scale(1.8) 
graph export "tabfig/descr/PSID_coefftowner.pdf", replace


graph combine all ft, ycommon
clear

************
** combine files to use with julia later possibly
use `results', clear
merge 1:1 fv using `leadmeans' 
sort fv

export delimited _all using "data/PSID/maintain.csv", replace // Creates a file with all the data that you can use in Julia
