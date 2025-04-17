qui {

cap scalar drop _all
cap graph drop _all
eststo clear

use "$PSIDpath/regression_sample.dta", clear
gen FFbehind = FFFF.behind
gen Fbehind2 = (FFbehind == 1 | Fbehind == 1)


local controls "c.wealth c.income i.hs i.coll i.white i.married i.year i.state c.age##c.age c.famsize c.age_prnt##c.age_prnt"
local lag_main "cashonhand_prnt "

eststo r1: regr Fbehind `lag_main' `controls' if Fage==first_own & first_own != firstage & inrange(age,25,44)
summarize Fbehind if e(sample)
estadd scalar meanoutcome   = `r(mean)': r1

eststo r2: regr Feverbehind `lag_main' `controls' if Fage==first_own & first_own != firstage  & inrange(age,25,44)
summarize Feverbehind if e(sample)
estadd scalar meanoutcome   = `r(mean)': r2

eststo r3: reg Fbehind  `lag_main' `controls' if Fowner == 1 &  behind == 0 &  inrange(age,25,44)
summarize Fbehind if e(sample)
estadd scalar meanoutcome   = `r(mean)': r3

eststo  r4: xtreg Fbehind  `lag_main' `controls' if Fowner==1  & behind == 0 & inrange(age,25,44), fe
summarize Fbehind if e(sample)
estadd scalar meanoutcome   = `r(mean)': r4



estadd local firsttime = "Y":r1
estadd local firsttime = "Y":r2
estadd local firsttime = "N":r3
estadd local firsttime = "N":r4

label var cashonhand_prnt "\;Par. Wealth"
label var wealth "\;Net Worth "
label var income "\;Income "
label var hs "\;High School"
label var coll "\;College"
label var white "\;White"
label var famsize "\;Family Size"
label var married "\;Married"
label var owner "\;Owner"

esttab using "tabfig/regr/hypoII.tex", replace drop(age* *.age* *state *.year 0.* _cons)  ///
	stats(N meanoutcome firsttime, label("N" "Outcome (mean)" "First-Time Buyers Only") fmt(%9.0fc 3 3) ) compress star(+ 0.1 * 0.05 ** 0.01 *** 0.001) se booktabs ///
	 mtitle("Behind" "Ever Behind" "Behind " "Behind FE") ///
	 label cells(b(fmt(3) star)  se(par fmt(3)))  /// 
	 refcat(wealth "\textit{Child}" cashonhand_prnt "\textit{Parent}", nolabel) collabels(none) 
	 
 esttab using "tabfig/regr/hypoII_short.tex", replace drop(*age* *.*age* *state *.*year 0.* _cons *hs* *coll* *white* *famsize* *married*)  ///
	stats(N meanoutcome firsttime, label("N" "Outcome (mean)" "First-Time Buyers Only") fmt(%9.0fc 3 3) ) compress star(+ 0.1 * 0.05 ** 0.01 *** 0.001) se booktabs ///
	 mtitle("Behind" "Ever Behind" "Behind " "Behind FE") ///
	 label cells(b(fmt(3) star)  se(par fmt(3))) 
}
esttab , replace drop(age* *.age* *state *.year 0.* _cons)  ///
	stats(N meanoutcome firsttime, label("N" "Outcome (mean)" "First-Time Buyers Only") fmt(%9.0fc 3 3) ) compress star(+ 0.1 * 0.05 ** 0.01 *** 0.001) se ///
	 mtitle("Behind" "Ever Behind" "Behind " "Behind FE") ///
	 label cells(b(fmt(3) star)  se(par fmt(3)))  /// 
	 refcat(wealth "\textit{Child}" cashonhand_prnt "\textit{Parent}", nolabel) collabels(none) 

