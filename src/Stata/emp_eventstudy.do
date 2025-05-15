
cap graph drop _all
use "$PSIDpath/panel_ind.dta" , clear
keep if id_hd!=.

xtset id_hd year
// Create necessary variables (unemployment and housing consumption growth rates)
gen unempl= ((emplstat == 3 ) & ll.emplstat == 1)  // We want the event of becoming unemployed from employment
replace unempl = . if emplstat ==.
replace unempl = . if ll.emplstat ==.



gen housing = .
replace housing = rent*12 if owner == 0 // Annualize monthly rent
replace housing = valhouse*0.05 if owner == 1 // Convert to owners equivalent rent (5% of house value)

gen lhousing = log(housing)
gen Dlhousing = lhousing - ll.lhousing

replace Dlhousing = 0 if wtrmove == 0 // Those who don't move have no change in housing consumption (same house!)


// Find percentile of parental wealth for households in age group
_pctile wealth_prnt [aw = famwgt] if inrange(age,25,44), p(75)
gen Twealthy_prnt = (wealth_prnt > `r(r1)') 
replace Twealthy_prnt = . if wealth_prnt == .


// Only keep households with one unemployment spell
bys id_hd: egen count_unemp = sum(unempl==1)
keep if count_unemp <= 1 // Only keep households unemployed at most one time

// Only keep observations +/- four years of unemployment spell
bys id_hd: egen yearofunemp = max(year*unempl)
gen yearrelunemp =  year - yearofunemp 
gen _keep = (abs(yearrelunemp) <= 4) // Only keep +/- four years of unemployment spell

// Only keep households with "simple" family changes
bysort id_hd: egen famchange_unemp= total(cond(yearrelunemp == 0, famchange, .)) // Parental wealth at unemployment
keep if famchange_unemp <= 2 // Only keep households that don't change family structure

// If wealthy at unemployment
bysort id_hd: egen wealthy_prnt= total(cond(yearrelunemp == 0, Twealthy_prnt, .)) // Parental wealth at unemployment
drop Twealthy_prnt

keep if age >= 25 & age <= 44

save "data/PSID/eventstudy.dta", replace // Used to do event study with controls

**********************************
** Without controls
**********************************
cap graph drop _all
use "data/PSID/eventstudy.dta", replace // Used to do event study with controls
keep if _keep == 1
keep if Dlhousing != .
keep if yearrelunemp != . 
keep if wealthy_prnt != .
count // observations
count if yearrelunemp == 0 //

collapse (mean) mean_Dlhousing = Dlhousing (sem) sem_Dlhousing=Dlhousing (count) id_hd [aw = famwgt] , by(yearrelunemp wealthy_prnt)
gen ub_Dlhousing= mean_Dlhousing + 1.96*sem_Dlhousing 
gen lb_Dlhousing= mean_Dlhousing - 1.96*sem_Dlhousing 
export delimited _all using "data/PSID/eventstudy.csv", replace // Creates a file with all the data that you can use in Julia

twoway  (rcap ub_Dlhousing lb_Dlhousing yearrelunemp if wealthy_prnt==0, lcolor(gray) yline(0, lpattern(solid)) )  (line mean_Dlhousing yearrelunemp if wealthy_prnt==0, lpattern(solid) lcolor(black)  ///
	xtitle("Years Relative to Unemployment") ytitle("Housing Growth Rate") ylabel(-0.15(0.1)0.25) legend(off) title("Parents in Bottom 75%")) 
graph display, xsize(6) ysize(3.5) scale(2.0) name(poor)

graph export "tabfig/descr/PSID_housinggrowthpoor_both.pdf", replace

twoway  (rcap ub_Dlhousing lb_Dlhousing yearrelunemp  if wealthy_prnt==1, lcolor(gray)  yline(0, lpattern(solid)) )  (line mean_Dlhousing yearrelunemp  if wealthy_prnt==1, lpattern(solid) lcolor(black) ///
	xtitle("Years Relative to Unemployment") ytitle("Housing Growth Rate") ylabel(-0.15(0.1)0.25) legend(off) title("Parents in Top 25%"))
graph display, xsize(6) ysize(3.5) scale(2.0) name(rich)
graph export "tabfig/descr/PSID_housinggrowthrich_both.pdf", replace


graph combine poor rich

******************************************
** with controls
******************************************
*/
use "data/PSID/eventstudy.dta", clear
cap graph drop _all
gen exp = year - yearofunemp 
replace exp = . if count_unemp == 0
recode exp (-1000/-5 = -5)
recode exp (5/1000 = 5)


gen Treat= wealthy_prnt
keep if abs(exp)<5
replace exp = exp + 5 // So that the factor variable doesn't contain negative values
xtile wealth_xtile = wealth [aw=famwgt], nq(5)
xtile income_xtile = income [aw=famwgt], nq(5)

replace income_xtile = log(income)
replace famsize = log(famsize)
winsor2 Dlhousing, replace c(1 99)

regr Dlhousing i.hs i.coll i.married i.white c.famsize i.famchange i.year i.state c.age##c.age c.age_prnt##c.age_prnt i.wealth_xtile c.income_xtile i.exp##i.Treat, vce(cluster id_hd)

mat coeff = e(b)
mat ses =  vecdiag(diag(vecdiag(e(V))))

clear
set obs  5
gen t = _n*2-1
gen Treat = .
gen Control = .

gen TreatSE = .
gen ControlSE = .
forv i = 3(2)11 {
	replace Treat = coeff[1,colnumb(coeff,"`i'.exp#1.Treat")] + coeff[1,colnumb(coeff,"1.Treat")] if t == `i'
	replace TreatSE = ses[1,colnumb(coeff,"`i'.exp#1.Treat")] if t == `i'
	replace Control = coeff[1,colnumb(coeff,"`i'.exp")] if t == `i'
	replace ControlSE = ses[1,colnumb(coeff,"`i'.exp")] if t == `i'
}

replace t = t-5
replace TreatSE = sqrt(TreatSE)
replace ControlSE = sqrt(ControlSE)


gen Control_cilo = Control - 1.96*ControlSE
gen Control_cihi = Control + 1.96*ControlSE
gen Treat_cilo = Treat - 1.96*TreatSE
gen Treat_cihi = Treat + 1.96*TreatSE

scatter Control Control_ci* t , c(l l l ) cmissing(y n n) ///
	msym(i i i) lcolor(gray gray gray) lpattern(solid dash dash) lwidth(thick medthick medthick) ///
	yline(0,  lpattern(solid)) xlabel(-4(2)4) ylabel(-0.3(0.15)0.3) xtitle("Years Relative to Unemployment") ytitle("Housing Growth Rate") legend(off) title("Parents in Bottom 75%") legend(off)  
graph display, xsize(6) ysize(3.5) scale(2.0) name(poor)
graph export "tabfig/descr/PSID_housinggrowthpoor_both_controls.pdf", replace

scatter Treat Treat_ci* t , c(l l l ) cmissing(y n n) ///
msym(i i i) lcolor(gray gray gray) lpattern(solid dash dash) lwidth(thick medthick medthick) ///
yline(0,  lpattern(solid)) xlabel(-4(2)4) ylabel(-0.3(0.15)0.3)  xtitle("Years Relative to Unemployment") ytitle("Housing Growth Rate")  legend(off) title("Parents in Top 25%") legend(off) 
graph display, xsize(6) ysize(3.5) scale(2.0) name(rich)
graph export "tabfig/descr/PSID_housinggrowthrich_both_controls.pdf", replace

graph combine poor rich
