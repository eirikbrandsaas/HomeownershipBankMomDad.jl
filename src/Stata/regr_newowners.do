qui {
cap scalar drop _all
cap graph drop _all

**************
** clean the data
***************

use  "$basepath/data/PSID/panel_ind", clear

by id_hd: gen firsthouse = owner[1]
by id_hd: gen firstage = age[1]
by id_hd: egen first_own = min(age / (owner == 1))
replace first_own = . if first_own == firstage // Cant be first-time ownership if it is the first observation

gen cashonhand_prnt = wealth_prnt + income_prnt
foreach var in income wealth wealth_prnt valhouse income_prnt cashonhand_prnt famsize {
	replace `var' = log(`var')
}

replace behind = . if owner!= 1 // Can't be behind if not owner
bys id_hd: egen everbehind = max(behind)


gen totrcved_fromparent1 = totrcved_fromparent
replace rcvr = 0 if totrcved_fromparent< $mingift
replace totrcved_fromparent1 = 0 if totrcved_fromparent< $mingift
replace totrcved_fromparent1 = . if totrcved_fromparent1==0


replace transf2 = 1 if (totrcved_fromparent > 10) & year == 2013 & !missing(totrcved_fromparent) // Include 2013 transfers
gen totrcved_fromparent3 = transf2val 
replace totrcved_fromparent3 = totrcved_fromparent3 + totrcved_fromparent if !missing(totrcved_fromparent) & (totrcved_fromparent > 10)
replace totrcved_fromparent3 = . if totrcved_fromparent3 == 0 

xtset id_hd year
foreach var in owner transf2 rcvr transf2val totrcved_fromparent totrcved_fromparent1 totrcved_fromparent3 age behind everbehind {
	gen F`var' = ff.`var'
}

save "$PSIDpath/regression_sample.dta", replace

*************
** run the regressions
***************
use "$PSIDpath/regression_sample.dta", clear

local lag_main "cashonhand_prnt "
local controls "c.wealth c.income i.hs i.coll i.married i.white c.famsize i.famchange i.year i.state c.age##c.age c.age_prnt##c.age_prnt"

eststo clear
eststo r0: reg Fowner i.Frcvr  if owner == 0 & inrange(age,25,44) & age <= first_own & !missing(Frcvr), vce(cluster famid )
summarize Frcvr if e(sample)
estadd scalar rcvr   = `r(mean)': r0
summarize Fowner if e(sample)
estadd scalar Fown   = `r(mean)': r0
summarize Ftotrcved_fromparent1 if e(sample), det
estadd scalar trsize = `r(p50)': r0


eststo r1: reg Fowner `lag_main' i.Frcvr  `controls' if owner == 0 & inrange(age,25,44) & age <= first_own & !missing(Frcvr), vce(cluster famid )
summarize Frcvr if e(sample)
estadd scalar rcvr   = `r(mean)': r1
summarize Fowner if e(sample)
estadd scalar Fown   = `r(mean)': r1
summarize Ftotrcved_fromparent1 if e(sample), det
estadd scalar trsize = `r(p50)': r1

eststo r3: reg Fowner  i.Ftransf2 if owner == 0 & inrange(age,25,44) & age < first_own , vce(cluster famid )
summarize Ftransf2 if e(sample)
estadd scalar rcvr   = `r(mean)'
summarize Fowner if e(sample)
estadd scalar Fown = `r(mean)': r3
summarize Ftotrcved_fromparent3 if e(sample), det
estadd scalar trsize= `r(p50)': r3


eststo r4: reg Fowner `lag_main' i.Ftransf2 `controls' if owner == 0 & inrange(age,25,44) & age < first_own , vce(cluster famid )
summarize Ftransf2 if e(sample)
estadd scalar rcvr   = `r(mean)'
summarize Fowner if e(sample)
estadd scalar Fown = `r(mean)': r4
summarize Ftotrcved_fromparent3 if e(sample), det
estadd scalar trsize= `r(p50)': r4


eststo r5: xtreg Fowner `lag_main' i.Ftransf2 `controls' if owner == 0 & inrange(age,25,44) & age <= first_own , fe vce(cluster famid )
summarize Ftransf2 if e(sample)
estadd scalar rcvr   = `r(mean)'
summarize Fowner if e(sample)
estadd scalar Fown = `r(mean)': r5
summarize Ftotrcved_fromparent3 if e(sample), det
estadd scalar trsize= `r(p50)': r5



estadd local controls   = "N": r0
estadd local controls   = "Y": r1
estadd local controls   = "N": r3
estadd local controls   = "Y": r4
estadd local controls   = "Y": r5


label var cashonhand_prnt "\;Par. Wealth"
label var wealth "\;Net Worth"
label var income "\;Income"
label var hs "\;High School"
label var coll "\;College"
label var white "\;White"
label var famsize "\;Family Size"
label var married "\;Married"

local giftthreshold = "0$mingift"
display "`giftthreshold'"

noi esttab r0 r1 r3 r4 r5 , replace drop(*age* *.*age* *state *.*year 0.* *famchange* _cons)  ///
	stats(N rcvr Fown trsize controls, label("N" "Receipt Rate" "Rent-to-Own Rate" "Median Transfer (1000 USD)" "Other Controls") fmt(%9.0fc 3 3 1 1)) compress star(+ 0.1 * 0.05 ** 0.01 *** 0.001) se  ///
	 mtitle("OLS" "OLS" "OLS" "OLS" "FE") ///
	 label cells(b(fmt(3) star)  se(par fmt(3)))  order(1.Frcvr  1.Ftransf2 ) /// 
	 varlabels(1.Ftransf2 "\;Any Transfer ($>$10k)" 1.Frcvr "\;Parent Transfer ($>$`giftthreshold'k)" ) ///
	 refcat(1.Frcvr "\textit{Transfer}" cashonhand_prnt "\textit{Other Controls}", nolabel) collabels(none) 

noi esttab r0 r1 r3 r4 r5  using "tabfig/regr/newowners_short.tex", replace drop(*age* *.*age* *state *.*year 0.* _cons *hs* *coll* *white* *famsize* *married* *income* *wealth* *cashonhand* *famchange*  )  ///
	stats(N rcvr Fown trsize controls, label("N" "Receipt Rate" "Rent-to-Own Rate" "Median Transfer (1000 USD)" "Other Controls") fmt(%9.0fc 3 3 1 1)) compress star(+ 0.1 * 0.05 ** 0.01 *** 0.001) se booktabs ///
	 mtitle("OLS" "OLS" "OLS" "OLS" "FE") ///
	 label cells(b(fmt(3) star)  se(par fmt(3)))  order(1.Frcvr  1.Ftransf2 ) /// 
	 varlabels(1.Ftransf2 "\;Any Transfer ($>$10k)" 1.Frcvr "\;Parent Transfer ($>$`giftthreshold'k)" ) ///
	 collabels(none) 

noi esttab r0 r1 r3 r4 r5  using "tabfig/regr/newowners.tex", replace drop(*age* *.*age* *state *.*year 0.* _cons *famchange* )  ///
	stats(N rcvr Fown trsize controls, label("N" "Receipt Rate" "Rent-to-Own Rate" "Median Transfer (1000 USD)" "Other Controls") fmt(%9.0fc 3 3 1 1))  compress star(+ 0.1 * 0.05 ** 0.01 *** 0.001) se booktabs ///
	 mtitle("OLS" "OLS" "OLS" "OLS" "FE") ///
	 label cells(b(fmt(3) star)  se(par fmt(3)))  order(1.Frcvr 1.Ftransf2 ) /// 
	 varlabels(1.Ftransf2 "\;Any Transfer ($>$10k)" 1.Frcvr "\;Parent Transfer ($>$`giftthreshold'k)" ) ///
	 refcat(1.Frcvr "\textit{Transfer}" cashonhand_prnt "\textit{Other Controls}", nolabel) collabels(none) 

***********************
** Plot the two transfer distributions
***********************

use "$PSIDpath/regression_sample.dta", clear
cap graph drop _all



sum Ftransf2val if  Ftransf2val>0 &  owner == 0 & inrange(age,25,44)  [aw = famwgt], det
local limhi = `r(p95)'

forv val =  0/1 {
	sum Ftransf2val if  Ftransf2val>0 &  owner == 0 & Fowner == `val' & inrange(age,25,44) & age < first_own [aw = famwgt], det
	local mean`val'_inh = `r(mean)'
	local med`val'_inh = `r(p50)'
}

twoway ///
    (kdensity Ftransf2val if Ftransf2val<`limhi' & Ftransf2val>0 &  owner == 0 & Fowner == 0 & inrange(age,25,44) [aw = famwgt], lcolor(black) lpattern(solid)) ///
    (kdensity Ftransf2val if Ftransf2val<`limhi' & Ftransf2val>0 &  owner == 0 & Fowner == 1 & inrange(age,25,44)  [aw = famwgt] , lcolor(orange) lpattern(solid)) ///
            , legend(label(1 "Rent -> Rent") label(2 "Rent -> Own") label(3 "Own -> Own") label(4 "Own -> Rent")) ///
	    xline(`mean0_inh', lcol(black) lpattern(longdash)) ///
	    xline(`med0_inh', lcol(black) lpattern(shortdash)) ///
	    xline(`mean1_inh', lcol(orange) lpattern(longdash)) ///
	    xline(`med1_inh', lcol(orange) lpattern(shortdash)) ///
	    xtitle("Thousands of USD") ytitle("Density") /// 
    title("Distribution of Transfers, Gifts, and Inheritances") ///
    name(g1) legend(rows(2) pos(0) bplacement(neast))

    graph display, xsize(6) ysize(3.5) scale(1.8) 
graph export "tabfig/descr/PSID_transfers_renters_inheritance.pdf", replace

sum Ftotrcved_fromparent if  Ftotrcved_fromparent > 0 & inrange(age,25,44) & owner == 0  [aw = famwgt], det
local limhi = `r(p99)'
forv val =  0/1 {
	sum Ftotrcved_fromparent if Ftotrcved_fromparent > 0 &  inrange(age,25,44) & owner == 0 & Fowner == `val' & age < first_own [aw = famwgt], det
	local mean`val'_trf = `r(mean)'
	local med`val'_trf = `r(p50)'
}
twoway ///
    (kdensity Ftotrcved_fromparent if Ftotrcved_fromparent < `limhi' & Ftotrcved_fromparent > 0 & inrange(age,25,44) & owner == 0 & Fowner == 0  [aw = famwgt], lcolor(black) lpattern(solid)) ///
    (kdensity Ftotrcved_fromparent if Ftotrcved_fromparent < `limhi' & Ftotrcved_fromparent > 0 & inrange(age,25,44) & owner == 0 & Fowner == 1  [aw = famwgt], lcolor(orange) lpattern(solid)) ///
            , legend(label(1 "Rent -> Rent") label(2 "Rent -> Own")) ///
	    xline(`mean0_trf', lcol(black) lpattern(longdash)) ///
	    xline(`med0_trf', lcol(black) lpattern(shortdash)) ///
	    xline(`mean1_trf', lcol(orange) lpattern(longdash)) ///
	    xline(`med1_trf', lcol(orange) lpattern(shortdash)) ///
	    xtitle("Thousands of USD") ytitle("Density") /// 
    title("Distribution of Transfers Received from Parents") legend(off) ///
    name(g2)
   
graph display, xsize(6) ysize(3.5) scale(1.8) 
graph export "tabfig/descr/PSID_transfers_renters_transfers.pdf", replace
   
graph combine g1  g2
}

