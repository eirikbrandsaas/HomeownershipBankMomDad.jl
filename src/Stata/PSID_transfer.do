// Family file:
use "$PSIDpath/FAM13.dta", clear
rename RT13V2 fuid // Interview number in 2013 - RT13V2=ER53002
rename RT13V4 wtrwife // whether wife present
rename RT13V22 gve_prnt
rename RT13V25 rcv_prnt
rename RT13V40 gve_chld
rename RT13V43 rcv_chld
rename RT13V9 n_chld_records
rename RT13V10 n_prnt_records
rename RT13V28 wtr_hd_rcv_home
rename RT13V31 wtr_wf_rcv_home
rename RT13V48 wtr_gve_home
rename RT13V49 gve_home_amt

drop RT*
gen year = 2013

// Replace DK/NA to missing

foreach var in gve_prnt gve_chld rcv_prnt rcv_chld gve_home_amt {
	replace `var' = . if `var' == 999999999
}

sum *chld *prnt

// Make real
merge m:1 year using "data/FRED/FRED_yearly",  keepusing(year cpi_index)
sum cpi_index if year == $baseyear
replace cpi_index = cpi_index / `r(mean)'

foreach var in gve_prnt rcv_prnt gve_chld rcv_chld {
	replace `var' = `var'/cpi_index
	replace `var' = `var'/$baseunit
}
keep if _merge == 3
drop _merge cpi_index

save "$PSIDpath/transfers_family_clean.dta", replace

// Parent-child file

use "$PSIDpath/PARCHD13.dta",clear
gen id_hd = (RT13V63*1000) + RT13V64
gen id_wf = (RT13V67*1000) + RT13V68
replace id_wf = . if id_wf == 0
gen id_rcv = (RT13V74*1000) + RT13V75
replace id_rcv = . if RT13V75>=800

rename (RT13V74 RT13V75) (famnm indnm)
rename (RT13V130 RT13V131 RT13V132) (tschool thouse tother)
rename (RT13V65  RT13V69 RT13V83) (age_hd age_wf age_rcv)
rename (RT13V70 RT13V71 RT13V77 RT13V78) (rectype recparunit rel2hd rel2wf)
rename (RT13V113 RT13V114 RT13V115 RT13V116) (they_wtrhome they_incge50 they_incge75 they_incl25)
rename (RT13V124 RT13V125 RT13V126) (wtrtotgiven totgiven_amt totgiven_time)
rename (RT13V127 RT13V128 RT13V129) (wtrtotrcved totrcv_amt totrcv_time)
rename (RT13V133) (geodistance)

drop RT*
drop if id_rcv == .
foreach var in thouse tschool tother totrcv_amt totgiven_amt  {
	replace `var' = . if `var' >= 999999998
}


keep id_* thouse tschool tother *_amt rectype age_* wtr*
sort id_hd id_rcv
order id_hd id_wf id_rcv tot*

// Recoding missing/NA/DKs
replace age_hd = . if age_hd == 999
replace age_rcv = . if age_rcv == 999 | age_rcv == 998
replace age_wf = . if age_wf == 999 | age_wf == 998 | age_wf == 0
replace wtrtotgiven = 0 if wtrtotgiven == 5
replace wtrtotrcved = 0 if wtrtotrcved == 5
gen year = 2013
drop id_wf
rename id_hd id_giver
rename id_rcv id_rcver
rename (totgiven_amt totrcv_amt) (totgiven totrcved)

// Make real
merge m:1 year using "data/FRED/FRED_yearly",  keepusing(year cpi_index)
sum cpi_index if year == $baseyear
replace cpi_index = cpi_index / `r(mean)'
keep if year>=1984

foreach var in totgiven totrcved thouse tschool tother {
	replace `var' = `var'/cpi_index
	replace `var' = `var'/$baseunit
}
keep if _merge == 3
drop _merge

// Clean up
label var id_giver "pid for giver"
label var id_rcver "pid for receiver"
label var totgiven "Gave 2012 "
label var totrcved "Recieved 2012 "
label var thouse "Gift housing - total"
label var tschool "Gift school - total"
label var tother "Gift other - total"


save "$PSIDpath/transfers_clean.dta", replace

collapse (sum) t*, by(id_rcver)
rename t* sum_t*

save "$PSIDpath/transfers_clean_sum.dta", replace

