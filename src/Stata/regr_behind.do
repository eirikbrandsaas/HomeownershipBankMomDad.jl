qui {

cap scalar drop _all
cap graph drop _all
eststo clear

use "$PSIDpath/regression_sample.dta", clear

local controls "c.wealth c.income i.hs i.coll i.married i.white c.famsize i.famchange i.year i.state c.age##c.age c.age_prnt##c.age_prnt "
local lag_main "cashonhand_prnt "

eststo r0: regr Fbehind `lag_main'  if Fage==first_own & first_own != firstage & inrange(age,25,44), vce(cluster famid ) 
summarize Fbehind if e(sample)
estadd scalar meanoutcome   = `r(mean)': r0
estadd local controls = "N": r0

eststo r1: regr Fbehind `lag_main' `controls' if Fage==first_own & first_own != firstage & inrange(age,25,44), vce(cluster famid ) 
summarize Fbehind if e(sample)
estadd scalar meanoutcome   = `r(mean)': r1
estadd local controls = "Y": r1

eststo r2: reg Fbehind  `lag_main'  if Fowner == 1 &  behind == 0 &  inrange(age,25,44), vce(cluster famid )
summarize Fbehind if e(sample)
estadd scalar meanoutcome   = `r(mean)': r2
estadd local controls = "N": r2

eststo r3: reg Fbehind  `lag_main' `controls' if Fowner == 1 &  behind == 0 &  inrange(age,25,44), vce(cluster famid )
summarize Fbehind if e(sample)
estadd scalar meanoutcome   = `r(mean)': r3
estadd local controls = "Y": r3

eststo  r4: xtreg Fbehind  `lag_main' `controls' if Fowner==1  & behind == 0 & inrange(age,25,44), fe vce(cluster famid )
summarize Fbehind if e(sample)
estadd scalar meanoutcome   = `r(mean)': r4
estadd local controls = "Y": r4



estadd local firsttime = "Y":r0
estadd local firsttime = "Y":r1
estadd local firsttime = "N":r2
estadd local firsttime = "N":r3
estadd local firsttime = "N":r4

label var cashonhand_prnt "\;Parental Wealth"
label var wealth "\;Net Worth "
label var income "\;Income "
label var hs "\;High School"
label var coll "\;College"
label var white "\;White"
label var famsize "\;Family Size"
label var married "\;Married"
label var owner "\;Owner"

esttab using "tabfig/regr/hypoII.tex", replace drop(age* *.age* *state *.year 0.* _cons *famchange*)  ///
	stats(N meanoutcome firsttime controls, label("N" "Outcome (mean)" "First-Time Buyers Only" "Other Controls") fmt(%9.0fc 3 3 3) ) compress star(+ 0.1 * 0.05 ** 0.01 *** 0.001) se booktabs ///
	 mtitle("OLS" "OLS" "OLS" "OLS" "FE") ///
	 label cells(b(fmt(3) star)  se(par fmt(3)))  collabels(none) 	 
	 
 esttab using "tabfig/regr/hypoII_short.tex", replace drop(*age* *income* *wealth* *.*age* *state *.*year 0.* _cons *hs* *coll* *white* *famsize* *married*  *famchange*)  ///
	stats(N meanoutcome firsttime controls, label("N" "Outcome (mean)" "First-Time Buyers Only" "Other Controls") fmt(%9.0fc 3 3 3) ) compress star(+ 0.1 * 0.05 ** 0.01 *** 0.001) se booktabs ///
	 mtitle("OLS" "OLS" "OLS" "OLS" "FE") ///
	 label cells(b(fmt(3) star)  se(par fmt(3))) collabels(none)
}
 esttab, drop(age* *.age* *state *.year 0.* _cons *famchange*)  ///
	stats(N meanoutcome firsttime controls, label("N" "Outcome (mean)" "First-Time Buyers Only" "Other Controls") fmt(%9.0fc 3 3 3) ) compress star(+ 0.1 * 0.05 ** 0.01 *** 0.001) se  ///
	 mtitle("OLS" "OLS" "OLS" "OLS" "FE") ///
	 label cells(b(fmt(3) star)  se(par fmt(3))) collabels(none)

