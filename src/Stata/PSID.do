*** Code to create a big ass dataset from the PSID, that has "all" the data you need
** Required files and stuff
* you need to install the psid tools package
* family files (famYYYY.zip)
* ind files (indYYYYer.zip)
* wealth supplements (wealthYYYY.zip)
* parent id (pidYY.zip)

*****************************

// Working the 2013 family roster and transfers
cd $basepath


// unzip the RT13.zip file
cap mkdir $PSIDpath/RT13 // Create RT13 directory, if the directory exists, nothing happens
cd $PSIDpath/RT13
unzipfile $PSIDpath/RT13.zip  // This will not replace if the contents are already there
cd $basepath


do "$PSIDpath/RT13/RT13FAM.do"
save "$PSIDpath/FAM13.dta", replace
do "$PSIDpath/RT13/RT13PARCHD.do"
save "$PSIDpath/PARCHD13.dta", replace

do "src/Stata/PSID_transfer.do" // My own file to clean and manipulate it


clear
psid install using $PSIDpath, to($PSIDpath)


psid use || valstocks		[99]ER15007 [01]ER19203 [03]ER22568 [05]ER26549 [07]ER37567 [09]ER43558 [11]ER48883 [13]ER54634 [15]ER61745 [17]ER67798 [19]ER73821 [21]ER79943  /// 
	 || valhouse 		[99]ER13041 [01]ER17044 [03]ER21043 [05]ER25029 [07]ER36029 [09]ER42030 [11]ER47330 [13]ER53030 [15]ER60031 [17]ER66031 [19]ER72031 [21]ER78032 ///
	 || houseexpend		[99]ER16515A5 [01]ER20456A5 [03]ER24138A5 [05]ER28037A5 [07]ER41027A5 [09]ER46971A5 [11]ER52395A5 [13]ER58212A5 [15]ER65414 [17]ER71491 [19]ER77520 [21]ER81847 ///
	 || rent		[99]ER13065 [01]ER17074 [03]ER21072 [05]ER25063 [07]ER36065 [09]ER42080 [11]ER47387 [13]ER53087 [15]ER60088 [17]ER66090  [19]ER72090 [21]ER78091 ///
	 || rooms		[99]ER13037 [01]ER17040 [03]ER21039 [05]ER25027 [07]ER36027 [09]ER42028 [11]ER47328 [13]ER53028 [15]ER60029 [17]ER66029 [19]ER72029 [21]ER78030 ///
	 || wealthwoequit 	[99]S416 [01]S516 [03]S616 [05]S716 [07]S816 [09]ER46968 [11]ER52392 [13]ER58209 [15]ER65406 [17]ER71483 [19]ER77509 [21]ER81836 ///
	 || wealth		[99]S417 [01]S517 [03]S617 [05]S717 [07]S817 [09]ER46970 [11]ER52394 [13]ER58211 [15]ER65408 [17]ER71485 [19]ER77511 [21]ER81838 ///
	 || homeequity 		[99]S420 [01]S520 [03]S620 [05]S720 [07]S820 [09]ER46966 [11]ER52390 [13]ER58207 [15]ER65404 [17]ER71481  [19]ER77507 [21]ER81834 ///
	 || transf		[99]ER15115 [01]ER19311 [03]ER22706 [05]ER26687 [07]ER37705 [09]ER43696 [11]ER49041 [13]ER54797 [15]ER61908 [17]ER67962 [19]ER73990 [21]ER80112 ///
	 || transfval		[99]ER15117 [01]ER19313 [03]ER22708 [05]ER26689 [07]ER37707 [09]ER43698 [11]ER49043 [13]ER54799 [15]ER61913 [17]ER67967 [19]ER73995 [21]ER80117 ///
	 || transfyear		[99]ER15116 [01]ER19312 [03]ER22707 [05]ER26688 [07]ER37706 [09]ER43697 [11]ER49042 [13]ER54798 [15]ER61910 [17]ER67964 [19]ER73992 [21]ER80114 ///
	 || income		[99]ER16462 [01]ER20456 [03]ER24099 [05]ER28037 [07]ER41027 [09]ER46935 [11]ER52343 [13]ER58152 [15]ER65349 [17]ER71426 [19]ER77448 [21]ER81775 ///
	 || state 		[99]ER13004 [01]ER17004 [03]ER21003 [05]ER25003 [07]ER36003 [09]ER42003 [11]ER47303 [13]ER53003 [15]ER60003 [17]ER66003 [19]ER72003 [21]ER78003 ///
	 || reltohead		[99]ER33503 [01]ER33603 [03]ER33703 [05]ER33803 [07]ER33903 [09]ER34003 [11]ER34103 [13]ER34203 [15]ER34303 [17]ER34503 [19]ER34703 [21]ER34903 ///
	 || sexhead 		[99]ER13011 [01]ER17014 [03]ER21018 [05]ER25018 [07]ER36018 [09]ER42018 [11]ER47318 [13]ER53018 [15]ER60018 [17]ER66018 [19]ER72018 [21]ER78018 ///
	 || educhead   		[99]ER16516 [01]ER20457 [03]ER24148 [05]ER28047 [07]ER41037 [09]ER46981 [11]ER52405 [13]ER58223 [15]ER65459 [17]ER71538 [19]ER77599 [21]ER81926 ///
	 || healthhead		[99]ER15447 [01]ER19612 [03]ER23009 [05]ER26990 [07]ER38202 [09]ER44175 [11]ER49494 [13]ER55244 [15]ER62366 [17]ER68420 [19]ER74428 [21]ER80550 ///
	 || age 		[99]ER33504 [01]ER33604 [03]ER33704 [05]ER33804 [07]ER33904 [09]ER34004 [11]ER34104 [13]ER34204 [15]ER34305 [17]ER34504 [19]ER34704 [21]ER34904 ///
	 || marstat		[99]ER16423 [01]ER20369 [03]ER24150 [05]ER28049 [07]ER41039 [09]ER46983 [11]ER52407 [13]ER58225 [15]ER65461 [17]ER71540 [19]ER77601 [21]ER81928 ///
	 || famchange		[99]ER13008A [01]ER17007 [03]ER21007 [05]ER25007 [07]ER36007 [09]ER42007 [11]ER47307 [13]ER53007 [15]ER60007 [17]ER66007 [19]ER72007 [21]ER78007 ///
	 || race		[99]ER15928 [01]ER19989 [03]ER23426 [05]ER27393 [07]ER40565 [09]ER46543 [11]ER51904 [13]ER57659 [15]ER64810 [17]ER70882  [19]ER76897 [21]ER81144 ///
	 || emplstat		[99]ER13205 [01]ER17216 [03]ER21123 [05]ER25104 [07]ER36109 [09]ER42140 [11]ER47448 [13]ER53148 [15]ER60163 [17]ER66164 [19]ER72164 [21]ER78167 ///
	 || famsize 		[99]ER13009 [01]ER17012 [03]ER21016 [05]ER25016 [07]ER36016 [09]ER42016 [11]ER47316 [13]ER53016 [15]ER60016 [17]ER66016  [19]ER72016 [21]ER78016 ///
	 || famwgt		[99]ER16518 [01]ER20394 [03]ER24179 [05]ER28078 [07]ER41069 [09]ER47012 [11]ER52436 [13]ER58257 [15]ER65492 [17]ER71570  [19]ER77631 [21]ER81958 ///
	 || expend_food		[99]ER16515A1 [01]ER20456A1 [03]ER24138A1 [05]ER28037A1 [07]ER41027A1 [09]ER46971A1 [11]ER52395A1 [13]ER58212A1 [15]ER65410 [17]ER71487 [19]ER77513 [21]ER81840 ///
	 || expend_transport	[99]ER16515B6 [01]ER20456B6 [03]ER24138B6 [05]ER28037B7 [07]ER41027B7 [09]ER46971B7 [11]ER52395B7 [13]ER58212B7 [15]ER65425 [17]ER71503 [19]ER77539 [21]ER81866 ///
	 || expend_healthcare	[99]ER16515D2 [01]ER20456D2 [03]ER24138D2 [05]ER28037D3 [07]ER41027D3 [09]ER46971D3 [11]ER52395D3 [13]ER58212D3 [15]ER65439 [17]ER71517  [19]ER77566 [21]ER81893 ///
	 || wtrmove		[99]ER13077 [01]ER17088 [03]ER21117 [05]ER25098 [07]ER36103 [09]ER42132 [11]ER47440 [13]ER53140 [15]ER60155 [17]ER66156 [19]ER72156 [21]ER78158 ///
	 || behind 									    [09]ER42052 [11]ER47359 [13]ER53059 [15]ER60060 [17]ER66062 [19]ER72062 [21]ER78063 ///
	 || pn_mom68		[]PID4 /// // These variables are here just to identify the parents of each person observed
	 || pn_mom		[]PID5 ///
	 || pn_dad68		[]PID23 ///
	 || pn_dad		[]PID24 /// 
	 || stratum		[]ER31996 ///
	 || cluster		[]ER31997 ///
	    using $PSIDpath, keepnotes design(1)

cd $basepath

psid long
svyset cluster [pweight= famwgt], strata(stratum)
rename wave year
rename x11101ll id
rename x11102 fuid
rename xsqnr sqnr
compress
save "$PSIDpath/temp17.dta", replace


*****************************************************
// Cleaning, renaming and generating new variables
use "$PSIDpath/temp17.dta", replace
gen famid = round(id/1000)*1000
gen iid = id - famid
replace famid = famid/1000

drop if sqnr == 0 // Drop those who don't have relationship to head (i.e. NA/DK and Latino sample) (See [17]ER34502)

replace pn_mom = . if pn_mom>=800 | pn_mom == 0 // (Replace with missing if inap. (=0) or parent never in any sample (>800))
replace pn_dad = . if pn_dad>=800 | pn_dad == 0

gen momid = pn_mom68*1000 + pn_mom
gen dadid = pn_dad68*1000 + pn_dad 

order id year fuid sqnr reltohead famid

gen spouse = (reltohead == 20 | reltohead == 22)
gen head = (sqnr==1 & reltohead == 10)
keep if spouse == 1 | head == 1 // Keep only heads and spouses of each household

// Now, we create a dataset of id's between heads and their spouses, and their spouses parents..
preserve 
keep if spouse == 1 | head == 1 // Keep only heads of each household
bys fuid year: drop if _N >= 3 // Drop HH with multiple spouses (Mormons?!)
keep id momid dadid fuid year head 
reshape wide id momid dadid, i(fuid year) j(head)
drop fuid
rename *0 *_wf
rename *1 *_hd
drop momid_hd dadid_hd
save "$PSIDpath/temp17_spousedata.dta", replace
restore

// Back to the main file
keep if head == 1 // Now, drop spouses
drop head spouse
rename id id_hd
merge 1:1 id_hd year using "$PSIDpath/temp17_spousedata.dta", nogen keep(3) // Only keep heads who have a single spouse

rename id_hd id


drop if famid >= 7000 & famid <= 9000 // Latino sample
// drop if famid >= 5000 & famid <= 7000 // SEO sample?

// Replace DK/NA/Refused with missing
replace educhead = . if educhead ==99
replace state =. if state == 99
replace marstat =. if marstat == 8 | marstat == 9
replace race = . if race==9 | race ==0 // Should'nt be any zeros, but there are a few so call them missing

replace valhouse = . if valhouse ==  9999998  | valhouse ==  9999999 // DK/Refused
replace valhouse = . if valhouse ==  9999997 			// topcoded
replace valhouse = . if valhouse ==  999999 & year == 1989 // topcoded in 1989
replace age = . if age==999 | age == 0
replace income = . if income == 9999999 & year == 1994 // Belongs to the latino sample
replace healthhead = . if healthhead == 8 | healthhead == 9 // replace DK/NA
replace emplstat = . if emplstat == 0 | emplstat == 22 | emplstat == 98 | emplstat == 99

replace rooms = . if rooms == 9 & year == 1984
replace rooms = . if rooms == 99 & year == 1989
replace rooms = . if rooms == 98 | rooms == 99 & year >= 1994
replace rooms = . if rooms == 98 & year >= 1994

replace rent = . if rent == 99998 | rent == 99999 

// Transfer/inheritance data is a fucking mess, lots of recoding
replace transf = . if transf == 8 | transf == 9 // replace DK/NA
replace transf = 0 if transf == 5
replace transf = 1 if transf == 1
// note: Code was changed in 2015

replace transfyear = . if transfyear == 7 | transfyear == 8 & year == 2015
replace transfyear = 2015 if transfyear == 1 & year == 2015
replace transfyear = 2014 if transfyear == 2 & year == 2015
replace transfyear = 2013 if transfyear == 3 & year == 2015 

replace transfyear = . if transfyear >= 97 & transfyear<=99 // DK/NA
replace transfyear = 1900 + transfyear if year == 1989 | year == 1984

replace transfyear = . if transfyear < 1979 // In in1984 they dont ask about last 5 years, just at any time...
replace transf = 0    if transfyear == . & year == 1984
replace transfval = 0 if transfyear == . & year == 1984
replace transfyear = . if transfyear >= 9997 & transfyear <=9999
replace transfyear = . if transfyear >= 9997 // Some years have error codes for 9997-9999, Dont know what they mean
replace transfval = . if  transfval == 9999997  | transfval == 9999998   | transfval == 9999999   | transfval == 1.0e+09


// Remake all variables into 2 year-variables to be consistent across all samples!
gen transf2 = (transfyear>= year - 2)*transf
gen transf2val = transfval*transf2
gen transf2year = transfyear*transf2
drop transf transfyear transfval

replace valstocks = . if valstocks ==  9999997 | valstocks ==  9999998 |  valstocks == 9999999
replace valstocks = . if valstocks == 1e+9
replace valstocks = . if valstocks == -1e+8

replace wtrmove = . if inlist(wtrmove,8,9) // Replace dk/na
replace wtrmove = 0 if wtrmove == 5
replace behind = . if inlist(behind,8,9) // replace dk/na
replace behind = 0 if behind == 5

// Generate new variables

sort id year
by id : gen t = _n
xtset id t

gen owner = (valhouse>0)
replace owner = . if valhouse == .  // Don't count missing as ownner
drop if owner == .
drop if wealth == .

gen white = (race == 1)
gen black = (race == 2)
drop race

gen newowner = (owner - l.owner>0    & !missing(l.owner))
replace newowner = . if l.owner == 1 |  missing(l.owner)

gen newrenter = (owner - l.owner<0 &! missing(l.owner))
replace newrenter = . if l.owner == 0 | missing(l.owner)

gen coll = (educhead == 16 | educhead == 17)
gen hs = (educhead >=12 & educhead< 16)
drop educhead

gen student = (emplstat == 7 ) 
gen retired = (emplstat == 4 ) 
gen housewife = (emplstat == 6 ) 

gen married = (marstat == 1)
drop marstat

gen alpha=. // Dont yet have these variables
gen calpha=. // Dont ye have these variables
gen alphaliq = valstocks/wealth 

gen male = (sexhead==1)
gen malehead = (sexhead == 1)
drop sex*

keep if age>=20


// Make variables real:
merge m:1 year using "$basepath/data/FRED/FRED_yearly",  keepusing(year cpi_index)

sum cpi_index if year == $baseyear

replace cpi_index = cpi_index / `r(mean)'
keep if year>=1984
drop _merge

foreach var of global realvars {
	replace `var' = `var'/cpi_index
	replace `var' = `var'/$baseunit
}


rename id id_hd
rename (momid dadid) (idmom_hd iddad_hd)
rename (momid_wf dadid_wf) (idmom_wf iddad_wf)

compress
order id_hd year id_wf id*_hd id*_wf
save "$PSIDpath/PSID17_HouseholdPanel.dta", replace
