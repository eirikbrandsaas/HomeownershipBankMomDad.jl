qui {
local datapath $datapath
local basepath $basepath 
set seed 333641 // Last 6 digits on a Highlight Pen on my desk
cap scalar drop _all
cap graph drop _all
use  "`basepath'/data/PSID/panel_winsor", clear

replace black = . if black != 1
replace black = 0 if white == 1
local lifecycle_vars = "wealth owner rent2inc hval2wealth wealthatpurchase first_own firstage mortg2inc transfsize "
local lifecycle_vars_stack "id_hd `lifecycle_vars' "

// Create necesarry variables 
gen logincome = log(income)
gen mortg = valhouse - homeequity
replace mortg = 0 if mortg < 0
replace mortg = . if owner == 0 
replace mortg = . if mortg > 1000

rename ctotrcved_fromparent  transfsize
replace transfsize = . if transfsize < $mingift

rename rcvr transfrate
replace transfrate = 0 if transfsize == . & transfrate == 1 // Replace those with gifts under threshold
gen tp2wealthp = transfsize / wealth_giver
replace tp2wealthp = . if tp2wealthp <= 0 // Can't have negative wealth and give in the model 
replace tp2wealthp = . if tp2wealthp > 1 // Avoid some outliers

gen tp2wealthk = transfsize / wealth
gen tp2inck = transfsize / (income*2)
replace tp2wealthk = . if wealth < 0 // Cant have negative wealth in the model

gen wealth_rcver = . 
gen wealth_nonrcver = . 
replace wealth_rcver = wealth if transfrate == 1
replace wealth_nonrcver = wealth if transfrate == 0

gen house2inc = valhouse/(income*2) // (Income is flow, mortg is stock so have to bake two-year annualize income)
replace house2inc = . if owner == 0
gen mortg2inc = mortg/(income*2) // (Income is flow, mortg is stock so have to bake two-year annualize income)
by id_hd: gen firsthouse = owner[1]
by id_hd: gen firstage = age[1]
by id_hd: egen first_own = min(age / (owner == 1))
replace first_own = . if firsthouse == 1
xtset id_hd year
gen wealthatpurchase = ll.wealth/owner
replace wealthatpurchase = . if age != first_own
replace wealthatpurchase = . if age >= 45
sort id_hd wealthatpurchase 

gen LTV = mortg/valhouse
gen LTVatpurchase = LTV
replace LTVatpurchase = . if age != first_own
replace LTVatpurchase = . if first_own == firstage
replace LTVatpurchase = . if LTVatpurchase > 1 & !missing(LTVatpurchase)
sort id_hd wealthatpurchase 


gen agegr = .
replace agegr = 1 if age>= 25 & age < 45
replace agegr = 2 if age>= 55 & age < 70

label define agegr 1 "young" 2 "old"
label value agegr agegr

gen wealthp_owner = wealth_prnt
replace wealthp_owner =. if owner == 0
gen wealthp_renter = wealth_prnt
replace wealthp_renter =. if owner == 1


tempfile main
rename hval2wealth h2w
rename rent2inc r2inc
replace r2inc = . if rent ==0  // Households who are renting for free are weird

gen race = ""
replace race = "All" if black == .
replace race = "Black" if black == 1
replace race = "White" if black == 0
xtset id_hd year
gen transferbuyers = .
replace transferbuyers= transfrate   if (ff.age == first_own| age == first_own |ll.age == first_own) & abs(year - 2013)<= 2 // keep households who bought around 2013 (where we see transfers)
gen transferbuyersamnt = .
replace transferbuyersamnt= transfsize  if (ff.age == first_own | age == first_own |ll.age == first_own) & abs(year - 2013)<= 2 // keep households who bought around 2013 (where we see transfers)

gen transferbuyersL = .
gen transferbuyersN = .
gen transferbuyersF = .
gen transferbuyersLamnt = .
gen transferbuyersNamnt = .
gen transferbuyersFamnt = .
replace transferbuyersL= transfrate   if (ll.age == first_own) & abs(year - 2013)<= 2 // keep households who bought around 2013 (where we see transfers)
replace transferbuyersLamnt= transfsize   if (ll.age == first_own) & abs(year - 2013)<= 2 // keep households who bought around 2013 (where we see transfers)
replace transferbuyersN= transfrate   if (age == first_own) & abs(year - 2013)<= 2 // keep households who bought around 2013 (where we see transfers)
replace transferbuyersNamnt= transfsize   if (age == first_own) & abs(year - 2013)<= 2 // keep households who bought around 2013 (where we see transfers)
replace transferbuyersF= transfrate   if (ff.age == first_own) & abs(year - 2013)<= 3 // keep households who bought around 2013 (where we see transfers)
replace transferbuyersFamnt= transfsize   if (ff.age == first_own) & abs(year - 2013)<= 2 // keep households who bought around 2013 (where we see transfers)



gen consp2consk = cons_prnt/cons
gen wealthf = wealth + wealth_prnt
// replace consp2consk = . if consp2consk > 5
// replace consp2consk = . if cons/cons_prnt > 5
local meanvars  "owner r2inc wealth h2w wealthatpurchase first_own logincome mortg2inc transfsize transfrate tp2wealthp tp2wealthk tp2inck transferbuyers transferbuyersamnt house2inc mortg LTVatpurchase wealth_rcver wealth_nonrcver transferbuyersL transferbuyersN transferbuyersF transferbuyersLamnt transferbuyersNamnt transferbuyersFamnt wealthp_owner wealthp_renter consp2consk wealthf"

replace first_own = firstage if firsthouse == 1 & firstage <27 // So that households who enter the sample as homeowners under age 27 also counts towards the age at purchase


// find overall homeownership rates (not needed for bootstrapped) (used in the end of the file, appended to the main outputs.)
preserve
keep if age >=25 & age <=74
collapse (mean) owner [aw = famwgt], by(year age)
collapse (mean) owner , by(year)
collapse (mean) owner
noi sum owner
local ownerall = `r(mean)'
restore
preserve
keep if black !=.
keep if age >=25 & age <=74
collapse (mean) owner [aw = famwgt], by(year age black)
collapse (mean) owner , by(year black)
collapse (mean) owner , by(black)
sum owner if black == 0
local ownerall_white = `r(mean)'
sum owner if black == 1
local ownerall_black = `r(mean)'
restore

keep `meanvars' agegr black famwgt age id_hd year

drop if missing(agegr)
tempfile sample
save `sample', replace // Sample the sample you bootstrap from

// bootstrapping
forv boot = 0/$Nboot {
	use `sample', clear // Reload the sample every iteration
	*************************
	** Part 2: Calculating the moments
	*************************
	
	if `boot' > 0 {
		noi bsample
		
	}
	
	
	// All households
	preserve 
	collapse (mean) `meanvars' (p50) medwealth = wealth tp2wealthk_med = tp2wealthk tp2wealthp_med = tp2wealthp transfsize_med = transfsize  medwealthp_owner = wealthp_owner medwealthp_renter = wealthp_renter medconsp2consk = consp2consk medwealthf = wealthf (sd) stdwealth = wealth [aw = famwgt] , by(agegr year) // Find the moments within years and agegroups
	collapse (mean) `meanvars' medwealth* *_med stdwealth medcons*, by(agegr) // Average moments over years
	gen medwealthpgrad = medwealthp_owner/medwealthp_renter
	gen wealthpgrad = wealthp_owner/wealthp_renter
	if `boot' ==0 { // Only do it first-time
		save `main'
	}
	else { // 
		save `main', replace
	}
	
	
	
		

	// By race (white/black)
	restore
	preserve
	drop if black == . 
	collapse (mean) `meanvars' (p50) medwealth = wealth  tp2wealthk_med = tp2wealthk tp2wealthp_med = tp2wealthp transfsize_med = transfsize medwealthp_owner = wealthp_owner medwealthp_renter = wealthp_renter  medconsp2consk = consp2consk medwealthf = wealthf (sd) stdwealth = wealth [aw = famwgt] ,  by(agegr year black)
	collapse (mean) `meanvars' medwealth* *_med stdwealth medcons*,  by(agegr black)
	gen medwealthpgrad = medwealthp_owner/medwealthp_renter
	gen wealthpgrad = wealthp_owner/wealthp_renter

	// Stack all the datasets
	if `boot' == 0 {
		append using `main'
	}
	else {
		append using `main'
	}

	// Do some cleaning, reshape etc
	gen race = ""
	replace race = "All" if black == .
	replace race = "Black" if black == 1
	replace race = "White" if black == 0
	drop black
	order agegr race

	sort race agegr
	reshape wide `meanvars' medwealth* stdwealth *_med wealthpgrad medcons* ,i(race) j(agegr) 
	rename *1 *young
	rename *2 *old
	order _all, alpha
	order race

	drop transfsizeold transfrateold tp2*old transf*old *grad*old *renterold *ownerold *consp2*old

	// convert transfers to two-year transfer rates
	if $biennialize == 1 {
		replace transfsizeyoung = transfsizeyoung * $biennializesize
		replace tp2wealthpyoung = tp2wealthpyoung * $biennializesize
		replace transfrateyoung = transfrateyoung * $biennializerate	
		foreach type in L N F {
			replace transferbuyers`type'amntyoung = transferbuyers`type'amntyoung * $biennializesize
			replace transferbuyers`type'young = transferbuyers`type'young  * $biennializerate
		}
		replace transferbuyersamntyoung = transferbuyersamntyoung * $biennializesize
		replace transferbuyersyoung = transferbuyersyoung  * $biennializerate
	}

	order race medwealthyoung medwealthold owneryoung h2wyoung r2incyoung wealthatpurchaseyoung LTVatpurchaseyoung first_ownyoung transfrateyoung transferbuyersyoung
	if `boot' == 0 {
		// Finally, add in ownerall variables (Not needed in bootstrapped)
		gen ownerall = .
		replace ownerall = `ownerall' if race == "All"
		replace ownerall = `ownerall_black' if race == "Black"
		replace ownerall = `ownerall_white' if race == "White"
		export delimited _all  using "data/moments/moments_race_age.csv", replace
		
	}
	else {
		gen boot = `boot'
		save  "data/moments/moments_race_age_boot`boot'", replace
	}
	restore
	noi disp `boot'
}
}
// Append bootstrap samples
clear
forv boot = 1/$Nboot {
	append using data/moments/moments_race_age_boot`boot'.dta
	rm data/moments/moments_race_age_boot`boot'.dta
}
export delimited _all  using "data/moments/moments_race_age_boot.csv", replace
