*****************************
** AHS
*****************************
use "$basepath/data/AHS/AHS_plotdata.dta", clear

// find (approximate) standard errors using standard binomial formula (var = p(1-p))
gen var_perc_giftinherit = perc_giftinherit*(1-perc_giftinherit)
gen semean_perc_giftinherit = sqrt(var_perc_giftinherit)/sqrt(Nobs)

gen lo = (perc_giftinherit - 1.96*semean_perc_giftinherit)*100
gen hi = (perc_giftinherit + 1.96*semean_perc_giftinherit)*100
replace perc_giftinherit = perc_giftinherit*100

label var perc_giftinherit "Gift or inheritance major source of downpayment"
twoway  (line hi lo year, lcolor(gray) lpattern(dash dash)) (line perc_giftinherit year, lpattern(solid) lcolor(black)  ) (scatter perc_giftinherit year, msymbol(Dh) mcolor(black)) ///
 , xtitle("Survey Year") ytitle("%") title("Mostly Gifts/Inheritance for Downpayment")  legend(off) ylabel(2.5(1)4.5)
graph display, xsize(6) ysize(3.5) scale(1.7) 
graph export "tabfig/descr/AHS_majorsourcedown_surveyyear_paper.pdf", replace
****************************
** Main SHED plot: Time series of transfers by year of purchase
****************************
use "$basepath/data/SHED/SHED_stacked.dta", clear


drop if h_year < 0 // -1 is missing, and not clear what the -9 and -2's are in 2017 wave
keep if prevown == 0 // Only keep those reporting first-time ownership
keep if h_year>=1990 // The 2015 wave only has from year 2001
egen hyear = cut(h_year), at(1999(1)2020)
replace hgift = hgift*100 // Make it into percentages

collapse (semean) hgift_se = hgift (mean) hsave hgift hnodo   hothe  hgovt hscnd (count) Nhgift=hgift [aw = wgt], by(hyear) 

label var hyear "Year of purchase"
label var hgift "Loan or gift from family"

gen lo = hgift - 1.96*hgift_se
gen hi = hgift + 1.96*hgift_se
twoway  (line hi lo hyear, lcolor(gray) lpattern(dash dash))  (line hgift hyear, lpattern(solid) lcolor(black) xlabel(2000(4)2017))  (scatter hgift hyear , msymbol(Dh) mcolor(black) ///
	xtitle("Year of Purchase") 	ytitle("%") title("Recieved Gifts/Loan For Downpayment") legend(off)) 
graph display, xsize(6) ysize(3.5) scale(1.7) 
graph export "tabfig/descr/SHED_gift_scatter_SE_paper.pdf", replace



********************************
** Various extra plots from the SHED (discussed in text, never used)
*****************************
use "$basepath/data/SHED/SHED_stacked.dta", clear
// By regions
preserve
replace prevown = . if prevown==-1
drop if missing(prevown)
keep if prevown == 0
keep if h_year>=2001
egen hyear = cut(h_year), at(2001(1)2020)
gen region = "Rest of the US"
replace region = "Mid-North Atlantic & Pacific" if ppreg9 == 2 | ppreg9 ==1 | ppreg9 == 9

encode region, gen(region2)
drop region
rename region2 region
collapse (mean) hsave hgift hnodo   hothe  hgovt hscnd (count) Nhgift=hgift [aw = wgt], by(hyear region) 
xtset region hyear
keep hgift hyear region
label var hyear "Year of purchase"
reshape wide hgift, i(hyear) j(region)
label var hgift1 "Loan or gift from family"
twoway (qfitci hgift1 hyear, $ylabel) (scatter hgift1 hyear), name("atlanticpacific") subtitle("Pacific & Mid-North Atlantic") nodraw legend(row(1))
twoway (qfitci hgift2 hyear, $ylabel) (scatter hgift2 hyear), name("restus") subtitle("Rest of the US") nodraw
grc1leg atlanticpacific restus, ycommon 
graph display, xsize(6) ysize(3) scale(1.5)
graph export "tabfig/descr/funding_mortgage_overtime_regions.eps", replace
cap graph drop atlanticpacific restus
restore 

**********************
// Reasons for renting (discussed in text)
**********************
use "$basepath/data/SHED/SHED_stacked.dta", clear

gen rcons = (rdown == 1 | rqual ==1 & !missing(rdown) & !missing(rqual))
global ylabel "ytitle("Fraction")"
global ylabelfund "ylabel(0(0.2)0.8)"
global ylabelrent "ylabel(0(0.2)0.8)"


label var rcons "Down payment/Don't qualify"
tab year if !missing(rdown)
graph bar rcons rconv rmove  rcheap rrisk  rpref  rothe  if !missing(rdown)  [pw = wgt] ,  $ylabel   $ylabelrent ///
	legend(order (1 "`: var label rcons'" 2 "`: var label rconv'" 3 "`: var label rmove'" 4 "`: var label rcheap'" 5 "`: var label rrisk'"  6 "`: var label rpref'" 7 "`: var label rothe'"))
graph display, xsize(6) ysize(3) scale(1.5)
graph export "tabfig/descr/why_rent.eps", replace


preserve
keep if year == 2013 | year == 2014 | year == 2016
xtile incquart = totinc  [pw = wgt], nq(4)

graph bar rcons rconv rmove  rcheap rrisk  rpref  rothe  if !missing(rdown) & incquart == 1 [pw = wgt] ,  $ylabel $ylabelrent ///
	legend(order (1 "`: var label rcons'" 2 "`: var label rconv'" 3 "`: var label rmove'" 4 "`: var label rcheap'" 5 "`: var label rrisk'"  6 "`: var label rpref'" 7 "`: var label rothe'"))
graph display, xsize(6) ysize(3) scale(1.5)
graph export "tabfig/descr/why_rent_botincome.eps", replace
restore
