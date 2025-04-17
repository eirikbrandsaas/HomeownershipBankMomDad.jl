qui {

cap scalar drop _all
cap graph drop _all
use  "$basepath/data/PSID/panel_ind", clear

by id_hd: gen firsthouse = owner[1]
by id_hd: gen firstage = age[1]
by id_hd: egen first_own = min(age / (owner == 1))

gen mortg = valhouse - homeequity
replace mortg = 0 if mortg < 0
replace mortg = . if owner == 0 
replace mortg = . if mortg > 1000

gen LTV = mortg/valhouse
gen LTVatpurchase = LTV
replace LTVatpurchase = . if age != first_own
replace LTVatpurchase = . if first_own == firstage
replace LTVatpurchase = . if LTVatpurchase > 1 & LTVatpurchase < .
replace LTVatpurchase = . if LTVatpurchase == 0 

gen cashonhand_prnt = wealth_prnt + income_prnt
foreach var in income wealth wealth_prnt mortg valhouse income_prnt cashonhand_prnt famsize {
	replace `var' = log(`var')
}

eststo clear

replace behind = . if owner!= 1 // Can't be behind if not owner
bys id_hd: egen everbehind = max(behind)

xtset id_hd year
foreach var in behind everbehind age owner {
	gen F`var' = ff.`var'
}


local controls "c.wealth c.income i.hs i.coll i.white i.married i.year i.state i.owner c.age##c.age c.famsize c.age_prnt##c.age_prnt"
local lag_main "cashonhand_prnt "

eststo r1: regr Fbehind `lag_main' `controls' if Fage==first_own & first_own != firstage & inrange(age,25,44)
summarize Fbehind if e(sample)
estadd scalar meanoutcome   = `r(mean)': r1

eststo r2: regr Feverbehind `lag_main' `controls' if Fage==first_own & first_own != firstage  & inrange(age,25,44)
summarize everbehind if e(sample)
estadd scalar meanoutcome   = `r(mean)': r2

eststo r3: reg Fbehind  `lag_main' `controls' if Fowner == 1 & inrange(age,25,44)
summarize behind if e(sample)
estadd scalar meanoutcome   = `r(mean)': r3

eststo  r4: xtreg Fbehind  `lag_main' `controls' if Fowner==1  & inrange(age,25,44), fe
summarize behind if e(sample)
estadd scalar meanoutcome   = `r(mean)': r4



label var cashonhand_prnt "\;Wealth(t-2)"
label var wealth "\;Net Worth(t-2)"
label var income "\;Income(t-2)"
label var hs "\;High School"
label var coll "\;College"
label var white "\;White"
label var famsize "\;Family Size"
label var married "\;Married"
label var owner "\;Owner"

esttab using "tabfig/regr/hypoII.tex", replace drop(age* *.age* *state *.year 0.* _cons)  ///
	stats(N meanoutcome, label("N" "Outcome (mean)") fmt(%9.0fc 3) ) compress star(+ 0.1 * 0.05 ** 0.01 *** 0.001) se booktabs ///
	 mtitle("Behind First" "Ever Behind" "Behind RE" "Behind FE") ///
	 label cells(b(fmt(3) star)  se(par fmt(3)))  /// 
	 refcat(wealth "\textit{Child}" cashonhand_prnt "\textit{Parent}", nolabel) collabels(none) 
}
esttab , replace drop(age* *.age* *state *.year 0.* _cons)  ///
	stats(N meanoutcome, label("N" "Outcome (mean)") fmt(%9.0fc 3) ) compress star(+ 0.1 * 0.05 ** 0.01 *** 0.001) se ///
	 mtitle("Behind First" "Ever Behind" "Behind RE" "Behind FE") ///
	 label cells(b(fmt(3) star)  se(par fmt(3)))  /// 
	 refcat(wealth "\textit{Child}" cashonhand_prnt "\textit{Parent}", nolabel) collabels(none) 
asd
// Useful to get the numbers reported in the text.
sum behind [aw = famwgt]
noi disp as text "Share of households behind: `r(mean)'"
sum everbehind [aw = famwgt]
noi disp as text "Share of households ever behind: `r(mean)'"
