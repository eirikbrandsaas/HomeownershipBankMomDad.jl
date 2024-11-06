************************  
** Start with yearly variables
************************

// CPI
import fred FPCPITOTLZGUSA, clear



rename FPCPITOTLZGUSA  cpi
gen year = yofd(daten)
tsset year
keep cpi year


gen cpi_index = cpi
replace cpi_index = 1 if _n==1
replace cpi_index = l.cpi_index*(1 + cpi/100) if _n>1
label var cpi_index "CPI Index (1960 = 1)"

save "$basepath/data/FRED/FRED_yearly.dta", replace
