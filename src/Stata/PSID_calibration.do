cap scalar drop _all
cap graph drop _all
********************************
// Initial wealth
********************************
local basepath $basepath
use  "`basepath'/data/PSID/panel_winsor", clear
// keep if year == $minyear | year == $minyear  + 2
keep if $firstage - 2 <=  age  & age <= $firstage + 2

keep age age_prnt agediff wealth wealth_prnt income_prnt famwgt income

keep if $minagediff <= agediff  & agediff <= $maxagediff
corr wealth*
corr income*

xtile xk = wealth 	[pw = famwgt], nq($Nxk)
xtile xp = wealth_prnt  [pw = famwgt], nq($Nxp)
xtile vp = income_prnt  [pw = famwgt], nq($nv)
xtile vk = income 	[pw = famwgt], nq($nv)


replace wealth = $wealth_threshold if wealth < $wealth_threshold // No notion of negative net worth in the model...

gen count = 1
tab xk vk if xp ==1 & vp == 1
collapse (count) count (median) wealth (max) wealth_prnt , by(xk xp vp vk)
order xk vk xp vp

bys xp vp: egen prob = total(count)
replace prob = count/prob
drop count
drop if xp == .
sort vp xp vk xk

// Deall with what happens if there are no observations within the bin
** Fill the indicators with true values, set probability to zero and set the wealth's to zero 
** (value doesnt matter since prob = 0)!*/
forv ivp = 1/$nvp {
	forv ivk = 1/$nv {
		forv ixk = 1/$Nxk {
			forv ixp = 1/$Nxp {
				count if xk == `ixk' & xp == `ixp' & vp == `ivp'  & vk == `ivk'
				if `r(N)' == 0 { 
				  count
				  expand 2 in `r(N)'
				  replace xk = `ixk' if _n==`r(N)'
				  replace xp = `ixp' if _n==`r(N)'
				  replace vk = `ivk' if _n==`r(N)'
				  replace vp = `ivp' if _n==`r(N)'
				  replace wealth = 0 if _n==`r(N)'
				  replace wealth_prnt = 0 if _n==`r(N)'
				  replace prob = 0 if _n==`r(N)'
				}	
			}
		}
	}
}
sort vp xp vk xk

preserve // Quick detour to plot the distrubtions
	egen group = group(vk xk)
	forv xp = 1/$Nxp {
		sum wealth_prnt  if xp == `xp'
		if `xp' < $Nxp {
			local xp`xp' = round(`r(mean)',1)
			local xp`xp' = "`xp`xp''"
		}
		else {
			local xp$Nxp = "Inf."
		}
		disp "`xp`xp''"
		
	}
	
	
	sort xp vp group
	bys xp vp: gen CDF = sum(prob)
	sort xp vp group
	
	forv ivp = 1/$nvp {
		local vp = `ivp'
		twoway (connect CDF group if xp == 1 & vp ==`vp') (connect CDF group if xp == 2 & vp ==`vp') (connect CDF group if xp == 3 & vp ==`vp') (connect CDF group if xp == 4 & vp ==`vp') , xlabel(1 5 9 12) ylabel(0(0.25)1) xline(5) xline(9) ///
			legend(label(1 "1st Quartile (x<`xp1')") label(2 "2nd Quartile (x<`xp2')") label(3 "3rd Quartile (x<`xp3')") label(4 "4th Quartile (x<`xp4')") rows(1)) title("F(:,`vp')") xtitle("Wealth x & Income y Decile") ///
			name(plt`ivp') nodraw
	}
	grc1leg plt1 plt2 plt3, row(1)
	graph display, xsize(12) ysize(4) scale(2.2)
	graph export "tabfig/empirical/inidistr.eps", replace
	graph export "tabfig/empirical/inidistr.pdf", replace
restore

export delimited xp vp xk vk prob wealth using "data/model/calibration/theta_$nv.csv", replace

collapse (first) wealth_prnt, by(xp)
keep if _n <$Nxp // Drop final value (e.g., for quartiles there are only three limits)
drop xp
export delimited _all using "data/model/calibration/thetaxplims_$nv.csv", replace novarnames
**************************************
// Pattern of deterministic income
**************************************
local basepath $basepath
use  "`basepath'/data/PSID/panel_winsor", clear
keep if age >=$firstage & age<=$lastage

collapse (mean) income [pw = famwgt], by(age year)
collapse (mean) income , by(age)
forv power = 1/$polyorder {
	gen age`power' = age^`power'/10^`power' if age<=$retage-1
}
regr income c.age?
scalar cons = _b[_cons]
forv power = 1/$polyorder {
	scalar incage`power'=_b[age`power']
}

predict yhat, xb
sum yhat if age == $retage-1 // Income at 66
local lastworkingincome = `r(mean)'
sum income if age >= $retage // Income for all households aged 67 or more
local retirementincome = `r(mean)'
scalar retinc = `retirementincome'/`lastworkingincome'
replace yhat =  `retirementincome' if age >= $retage
replace yhat = yhat*2
replace income = income*2

label var yhat "Predicted (Polynomial) Income"
label var income "Average Income"
twoway line yhat income age, xtitle("Age") ytitle("Income $1000") legend(off)
graph display, xsize(6) ysize(3) scale(2.2)
graph export "tabfig/empirical/lifecycleincome.eps", replace
graph export "tabfig/empirical/lifecycleincome.pdf", replace


clear
local obs =  $polyorder + 1
set obs `obs'
gen betas = .
replace betas = cons if _n==1
forv power = 1/$polyorder {
	replace betas = incage`power' if _n == `power'+1
}

export delimited betas using "data/model/calibration/ageincome_betas.csv", novarnames replace

*************************************
** Racial income gaps
**************************************
use  "$basepath/data/PSID/panel_winsor", clear
replace black = . if black != 1
replace black = 0 if white == 1

qui{
	preserve
	keep if age >=$firstage & age<=$lastage
	collapse (mean) income [pw = famwgt], by(age year)
	collapse (mean) income , by(age)
	noi sum income
	local meanall = `r(mean)'
	restore

	preserve
	keep if black==1
	keep if age >=$firstage & age<=$lastage
	collapse (mean) income [pw = famwgt], by(age year)
	collapse (mean) income , by(age)
	noi sum income
	local blackgap= `r(mean)'/`meanall'
	restore

	preserve
	keep if white==1
	keep if age >=$firstage & age<=$lastage
	collapse (mean) income [pw = famwgt], by(age year)
	collapse (mean) income , by(age)
	noi sum income
	local whitegap= `r(mean)'/`meanall'
	restore
}

clear
set obs 2
gen race = "White" if _n==1
gen gap = `whitegap' if _n==1
replace race = "Black" if _n==2
replace gap = `blackgap' if _n==2
export delimited _all  using "data/model/calibration/racialincomegaps.csv", replace

**************************
// Income shocks - Old
**************************
local basepath $basepath
use  "`basepath'/data/PSID/panel_winsor", clear
// keep if year == $minyear | year == $minyear  + 2
keep if age >=$minagepar & age<=$lastage
bys id_hd: keep if _N >= 2
keep id_hd year age income famwgt expend_healthcare
egen agecl = cut(age), at(25(2)$lastage)
xtset id_hd year
keep if income > 10
replace income = income - expend_healthcare
replace income = ff.income / income

gen median = .
gen quant=.

drop if income > 2
drop if income < -1
levelsof agecl, local(ages)
foreach age of local ages {
	capture drop xq
	xtile xq=income [pw = famwgt] if agecl==`age', nq($nvp)
	replace quant=xq if agecl==`age' & missing(quant)
	summarize income [aw = famwgt] if agecl==`age', det
	replace median = `r(p50)'  if agecl == `age'
}

collapse (median) income median (count) obs = income [pw = famwgt], by(agecl quant year)
replace income = income/median
collapse (mean) income (sum) obs, by(agecl quant) // Give each year equal weight
rename income  vp_grd
drop if agecl == .
drop if quant == .
reshape wide vp_grd obs , i(agecl) j(quant)
gen agecl2 = agecl^2
forv ivp = 1/$nvp {
	regr vp_grd`ivp' agecl agecl2
	predict vpval`ivp', xb
}

twoway (line vp_grd* agecl, lpattern(dash dash dash dash))  (line vpval* agecl, lpattern(solid solid solid solid) color(black gs10 sky turqoise)) , xtitle("Age") ytitle("Productivity") legend(off) ylabel(0(0.5)2.5)
graph display, xsize(4) ysize(3) scale(1.5)
graph export "tabfig/empirical/lifecycleproductivity_old_$nvp.eps", replace
graph export "tabfig/empirical/lifecycleproductivity_old_$nvp.pdf", replace
keep agecl vpval*
export delimited _all using "data/model/calibration/vpgrid_$nvp.csv", replace

************************
// Income shocks Young
************************
local basepath $basepath
use  "`basepath'/data/PSID/panel_winsor", clear
// keep if year == $minyear | year == $minyear  + 2
keep if age >=$firstage & age<=$lastage
bys id_hd: keep if _N >= 2
keep id_hd year age income famwgt
xtset id_hd year
local maxage = $minagepar + 2
egen agecl = cut(age), at(25(2)`maxage')
drop if missing(agecl)
drop age
gen quant=.


levelsof agecl, local(ages)

foreach age of local ages {
	capture drop xq
	xtile xq=income [pw = famwgt] if agecl==`age', nq($nv)
	replace quant=xq if agecl==`age' & missing(quant)
}

preserve // Find the values for the v-grid
gen median = .
levelsof agecl, local(ages)
foreach age of local ages {
	summarize income [aw = famwgt] if agecl==`age', det
	replace median = `r(p50)'  if agecl == `age'
}

collapse (p50) income (first) median (count) obs = income [pw = famwgt], by(agecl quant year)
bys agecl: gen vgrid = income / median

collapse (mean) vgrid (count) obs, by(quant agecl)
reshape wide vgrid obs, i(agecl) j(quant)
gen agecl2 = agecl^2
forv iv = 1/$nv {
	regr vgrid`iv' agecl agecl2
	predict vval`iv', xb
}
twoway (line vgrid* agecl, lpattern(dash dash dash dash))  (line vval* agecl, lpattern(solid solid solid solid) color(black gs10 sky turqoise)) , xtitle("Age") ytitle("Productivity") legend(off) ylabel(0(0.5)2.5)
graph display, xsize(6) ysize(3) scale(2.2)
graph export "tabfig/empirical/lifecycleproductivity_$nv.eps", replace
graph export "tabfig/empirical/lifecycleproductivity_$nv.pdf", replace

keep agecl vval*
export delimited agecl vval* using "data/model/calibration/vgrid_k_$nv.csv", replace
restore

// Find the transition matrices
levelsof agecl, local(ages)
forv iv =1/$nv {
	forv ivp = 1/$nv {
		gen v`iv'to`ivp' = 0
		foreach age of local ages {
			replace v`iv'to`ivp' = v`iv'to`ivp' + 1 if quant==`iv' & ff.quant == `ivp' & agecl == `age'
		}
	}
}
sum(year)
keep if year <= `r(max)' - 2 // Cannot transition in last period
forv iv =1/$nv {
	gen rowsum`iv' = 0
	forv ivp = 1/$nv {
		replace rowsum`iv' = rowsum`iv' + v`iv'to`ivp]'
	}
}

collapse (mean) v* (mean) rowsum*, by(agecl)
forv iv =1/$nv {
	forv ivp = 1/$nv {
		replace v`iv'to`ivp' = v`iv'to`ivp'/rowsum`iv'
	}
}

gen agecl2 = agecl^2

forv iv = 1/$nv {
	forv ivn = 1/$nv {
		regr v`iv'to`ivn' agecl agecl2
		predict Piv`iv'to`ivn', xb
	}
	egen Piv`iv'_tot = rowtotal(Piv`iv'*)
	forv ivn = 1/$nv {
		replace Piv`iv'to`ivn' =  Piv`iv'to`ivn'/ Piv`iv'_tot 
	}
	drop Piv`iv'_tot
}
cap graph drop _all
local graphs ""
forv iv = 1/$nv {
	twoway (line v`iv'to* agecl, lpattern(dash dash dash dash))  (line Piv`iv'* agecl, lpattern(solid solid solid solid) color(black gs10 sky turquoise)) , xtitle("Age") ytitle("Probability") legend(off) name(gr`iv') nodraw title("Pi(:,`iv')") ylabel(0(0.25)1)
	local graphs = "`graphs' gr`iv'"
}
graph combine `graphs', rows(1)
graph display, xsize(6) ysize(2) scale(1.6)

graph export "tabfig/empirical/Pi_vk_$nv.eps", replace
graph export "tabfig/empirical/Pi_vk_$nv.pdf", replace
keep if agecl < $minagepar
keep agecl Piv*

export delimited agecl Piv* using "data/model/calibration/Piv_k_$nv.csv", replace
************************************
// Find market value of houses
************************************
local basepath $basepath
use  "`basepath'/data/PSID/panel_winsor", clear
// keep if year == $minyear | year == $minyear  + 2
keep if age >=25  & age<=44

svy: mean valhouse if owner == 1
mat b = e(b)
local valhouse = b[1,1]

// File with all single parameters
clear
set obs 1
gen lambda = retinc
gen retage = $retage
gen valhouse = `valhouse'

export delimited _all using "data/model/calibration/parameters.csv", replace
