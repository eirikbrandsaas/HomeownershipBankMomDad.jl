// Step 1: Clean CSVs downloaded from AHS website (2011-2019)
forv year = 2011(2)2021 {

    import delimited "/$basepath/data/AHS/AHS TABLE13 4_15_`year'.csv", clear

    // Drop all obs before major source
    tempvar first
    gen `first' = v1 == "Major Source of Down Payment"
    replace `first' = sum(`first')
    drop if `first' == 0 


    // Drop all obs after how acquired (first after major source of down)
    tempvar last
    gen `last' = v1 == "How Acquired"
    replace `last' = sum(`last')
    drop if `last' == 1

    // basic cleaning
    drop if v1 == ""
    rename v1 how
    rename v2 value

    if _N == 0 {
        set obs 1
    }
    gen year = `year'
    drop *__* // Drop tempvar
    drop if how == "Major Source of Down Payment" // Drop header
    destring value, replace ignore(",") // Convert to numeric

    save "$basepath/data/AHS/`year'_cleaned.dta", replace
}

clear
forv year = 2011(2)2021 {
    append using "$basepath/data/AHS/`year'_cleaned.dta"
    rm "$basepath/data/AHS/`year'_cleaned.dta"
}

drop if year == 2021 // Drop 2021 as it is not does not contain source of downpayment

replace how = subinstr(how, " ", "", .)
replace how = subinstr(how, ",", "", .)
replace how = "Borrowingothr" if how == "Borrowingotherthanmortgageonthisproperty"
replace how = "Landvalue" if how == "Landwherebuildingbuiltusedforfinancing"
reshape wide value, i(year) j(how) string
rename value* *
rename Homepurchasedorbuilt Obs


gen Nobs = (Obs - Notreported - Saleofprevioushome)
gen perc_giftinherit = Inheritanceorgift/Nobs

keep year perc_giftinherit Nobs
save "$basepath/data/AHS/AHS_plotdata.dta", replace
// Step 2: Add the manually scraped data from PDFs (1991-2009)

import delimited "$basepath/data/AHS/AHS_manual_scrape.csv", clear
reshape wide value, i(year) j(how) string
rename value* *

gen Nobs = (obs - missing - prevhome)
gen perc_giftinherit = giftinherit/Nobs
keep year perc_giftinherit Nobs

append using "$basepath/data/AHS/AHS_plotdata.dta"
tsset year 
save "$basepath/data/AHS/AHS_plotdata.dta", replace
