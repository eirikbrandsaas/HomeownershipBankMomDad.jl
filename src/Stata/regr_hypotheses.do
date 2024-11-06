local datapath $datapath
local basepath $basepath 
local period $period

cap scalar drop _all
cap graph drop _all
use  "`basepath'/data/PSID/panel_winsor", clear

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
replace first_own = . if first_own == firstage
local controls "i.state i.hs i.coll i.white i.year c.age##c.age c.famsize c.age_prnt##c.age_prnt"
foreach var in cashonhand_prnt wealth income {
	gen L`var' = ll.`var'
}
local lag_main "Lcashonhand_prnt Lwealth Lincome"
local behindcontrols = "" // No extra controls for these regressions
keep if age <= 44 & age >=25

eststo: regr valhouse `lag_main' `controls' if age==first_own 
eststo: regr behind `lag_main' `controls' `behindcontrols' if age==first_own 
eststo: regr everbehind `lag_main' `controls' `behindcontrols' if age==first_own 
eststo: xtreg behind  `lag_main' `controls' `behindcontrols' c.age##c.age  if first_own != . & owner==1
eststo: xtreg behind  `lag_main' `controls' `behindcontrols' c.age##c.age if first_own != . & owner==1 , fe

label var Lcashonhand_prnt "\;Wealth(t-2)"
label var Lwealth "\;Net Worth(t-2)"
label var Lincome "\;Income(t-2)"
label var hs "\;High School"
label var coll "\;College"
label var white "\;White"
label var famsize "\;Family Size"

esttab using "tabfig/regr/hypoII.tex", replace drop(age* *.age* *state *.year 0.* _cons)  ///
	stats(N, fmt(%9.0fc)) compress star(+ 0.1 * 0.05 ** 0.01 *** 0.001) se booktabs ///
	 mtitle("House Value" "Behind First" "Ever Behind" "Behind RE" "Behind FE") ///
	 label cells(b(fmt(3) star)  se(par fmt(3)))  /// 
	 refcat(Lwealth "\textit{Child}" Lcashonhand_prnt "\textit{Parent}", nolabel) collabels(none) 


// Useful to get the numbers reported in the text.
sum behind [aw = famwgt]
noi disp as text "Share of households behind: `r(mean)'"
sum everbehind [aw = famwgt]
noi disp as text "Share of households ever behind: `r(mean)'"