//CREATE DATABASE WITH EMISSION SCENARIOS
//IPCC AR6 report database
cd "G:\My Drive\Google Research Projects\Stranded Assets2\Data\IPCC_AR6"
set more off
//Create one database in 'long' format for each variable
	//Raw data in 'wide' format from https://data.ece.iiasa.ac.at/ar6 
	//Remark: I manually make a variable with the SSP number in case the name of the scenario mentioned it. Reference scenario is defined as the scnario with the highest emissions in 2100 (or highest cumulative emissions?) within the scenarios of a given model (and within SSP if available). Reference scenario must at least have 2.5 °C. For groups without reference scenario, we use the mean of all other reference scenarios. 
	//In the 1.5 database, I Manually created a variable ref_family with the same number for a group of scenarios with the same reference scenario and a dummy varialbe ref_scenario in file Emissions (and drop year 2014 in file emissions)
	//
	//However, cumulative emissions as a proxy of temperature has the advantage of being available for all obs.
	foreach xxx in Emissions GDP Loss_GDP CarbonPrice Temperature {
		import excel `xxx', sheet("data") firstrow clear
		gen Model_Scenario= Model+Scenario
		drop Region Variable Unit
		local zzz 2015
		foreach yyy in F G H I J K L M N O P Q R S T U V W  {
			rename `yyy' `xxx'`zzz'
			local zzz=`zzz'+5
		}
		reshape long `xxx', i(Model_Scenario) j(year)
		save `xxx', replace
	}
//Import metadata
import excel Metadata, sheet("meta") firstrow clear
keep Model Scenario Policy_category Policy_category_name Ssp_family Climateimpacts Scenario_scope
gen  Model_Scenario= Model+Scenario
save Metadata, replace

//Merge databases containing each one variable
	use Emissions, clear
	foreach xxx in GDP Loss_GDP CarbonPrice Temperature {
		merge 1:1 Model_Scenario year using `xxx'
		if "`xxx'"!="CarbonPrice" drop if _merge==2
		drop _merge
	}
	merge m:1 Model_Scenario using Metadata 
	drop _merge
//alternative: 	//download table with all variables in the same file:	//reshape wide 2015 2020 ....2100, i(Model_Scenario) j(Variable)	//local zzz 2010	//foreach xxx in Emissions GDP ...{	//rename `zzz'`xxx' `xxx'`zzz'	//}	//reshape long Emissions GDP ... , i(Model_Scenario) j(year)	//merge 1:1 Model_Scenario year using Emissions keepusing(Ref_family)
save "G:\My Drive\Google Research Projects\Stranded Assets2\Code\Stata\IPCCAR6_long", replace

//CREATE DATABASE WITH NGFS SCENARIOS
cd "G:\My Drive\Google Research Projects\Stranded Assets2\Data\NGFS2024\"
set more off
//This database uses GDP MER, whereas for IPCC I used PPP. 
	foreach xxx in  Emissions GDP Loss_GDP CarbonPrice Temperature {  
		import excel `xxx',  firstrow clear
		gen Model_Scenario= Model+Scenario
		drop Region Variable Unit
		local zzz 2020
		foreach yyy in F G H I J K L M N O P Q R S T U V {
			rename `yyy' `xxx'`zzz'
			local zzz=`zzz'+5
		}
		destring *20*,replace
		reshape long `xxx', i(Model_Scenario) j(year)
		save `xxx', replace
	}
//Merge databases containing each one variable
	use Emissions, clear
	foreach xxx in  GDP Loss_GDP CarbonPrice {
		merge 1:1 Model_Scenario year using `xxx'
		if "`xxx'"!="CarbonPrice" drop if _merge==2
		drop _merge
	}
gen NGFS=1

save "G:\My Drive\Google Research Projects\Stranded Assets2\Code\Stata\NGFS2024_long", replace
