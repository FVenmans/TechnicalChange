//Append IPCC and NGFS DATA
cd "G:\My Drive\Google Research Projects\Stranded Assets2\Code\Stata\"
use IPCCAR6_long, replace
append using NGFS2024_long
replace NGFS=0 if NGFS!=1




//Define as panel
	egen id=group(Model_Scenario)
	gen t=(year-2015)/5
	xtset id t
	gen t2020=year-2020
//drop scenarios with zero carbon price and emissions rapidly going to zero
	drop if CarbonPrice==0 & Model =="MESSAGEix-GLOBIOM_GEI 1.0" 
//Express emissions in GtCO2, Prod in $Trillion 
	foreach xxx in Emissions GDP Loss_GDP {
	replace `xxx'=`xxx'/1000
	}
//interpolate missing data
	by id: ipolate GDP year, gen(GDPi)
	by id: ipolate Loss_GDP year, gen(Loss_GDPi)
	by id: ipolate CarbonPrice year, gen (CarbonPricei)
	by id: ipolate Emissions year, gen(Emissionsi) 	//replace Emissions=(L.Emissions +F.Emissions)/2 if Emissions==. same result 
	gen interpolated=0
	foreach xxx in GDP Loss_GDP CarbonPrice Emissions{
		replace interpolated=1 if `xxx'!=`xxx'i
		replace `xxx'=`xxx'i
		drop `xxx'i
	}
	drop if Emissions==.
//Cumulative emissions
	gen Emissions_temp=0.5*Emissions+0.5*L.Emissions if year>2020 //mean of the preceding 5 years
	bys id (t): gen E_Cum2020=sum(Emissions_temp) if year>2020 & Emissions!=. //cumulative emissions
	replace E_Cum2020=5*E_Cum2020
	replace E_Cum2020=0 if year==2020 & Emissions!=.
	by id: egen E_Cum2100=max(E_Cum2020) if Emissions[_n]!=.
	
//Use mean GDP of SSP or NGFS scenarios for missing GDP (some data has Loss_GDP but not GDP)
	bys Ssp_family year : egen GDP_meanSSP=mean(GDP) if NGFS==0
	gen GDP_calculated= GDP==. & GDP_meanSSP!=.
	replace GDP=GDP_meanSSP if GDP==.
	bys Scenario year : egen GDP_meanNGFS=mean(GDP) if NGFS==1
	replace GDP_calculated= GDP==. & GDP_meanSSP!=.
	replace GDP=GDP_meanNGFS if GDP==.
	drop GDP_mean*
//Define baseline scenario and Loss GDP	
	gen baseline= (Policy_category=="P0" |Policy_category=="P1" | Policy_category=="P1a" | Policy_category=="P0_1a" | Policy_category=="P0_1d" | Scenario=="Current Policies") & (Climateimpacts=="no"|Climateimpacts=="")  & CarbonPrice<30
	//category P0_2 is common but has abatement
	replace baseline=0 if substr(Model, 1, 5) == "POLES"  
	//The low demand scenario in NGFS scenarios has lower GDP than the Net Zero 2050 scenario despite a lower carbon price. Therefore, GDP cannot be used to find losses.
//Loss_GDP 	
	bys id: egen Total_Loss_GDP=total(Loss_GDP)
	replace Loss_GDP= -Loss_GDP if Total_Loss_GDP < 0 // both negative and positive numbers are used in the raw data. 
	replace Loss_GDP=0 if baseline==1 & Loss_GDP==. 
	//Try to construct Loss_GDP from GDP baseline observations, but this does not work
		//548 extra observations but all zero except a few values at 7% loss and zero carbonprice
		//For example there is one reference scenario made with poles, but the GDP data in 2100 are the same for all four subgroups of scenarios
		//bysort Model Ssp_family year (baseline): gen Loss_GDP_calculated= GDP[_N]!=. & GDP!=. & baseline[_N]==1 & Loss_GDP==. & NGFS==0 &  GDP_calculated==0 
		//bysort Model Ssp_family year (baseline): replace Loss_GDP= (GDP[_N]+GDP[_N-1])/2 - GDP if baseline[_N]==1 & baseline[_N-1]==1 & Loss_GDP==. & NGFS==0 &  GDP_calculated==0	
		//bysort Model Ssp_family year (baseline): replace Loss_GDP= GDP[_N] - GDP if baseline[_N]==1 & Loss_GDP==. & NGFS==0 &  GDP_calculated==0 //
	gen GDP_loss_pct=Loss_GDP/GDP 
	gen GDP_ratio=1-GDP_loss_pct
	gen lnGDP_ratio=ln(GDP_ratio)

//Truncation dummy 
	//(not used because zero for obs with only carbon price and missing GDP loss)
	gen MAC=CarbonPrice/GDP/1000
	bys year: egen GDP_ratio_5 =pctile(GDP_ratio) ,p(5)
	bys year: egen GDP_ratio_95=pctile(GDP_ratio) ,p(95)
	bys year: egen MAC_5 =pctile(MAC) ,p(5)
	bys year: egen MAC_95=pctile(MAC) ,p(95)
	gen trunc5= GDP_ratio<GDP_ratio_5 | GDP_ratio>GDP_ratio_95 | MAC<MAC_5 | MAC>MAC_95
	
	bys year: egen GDP_ratio_10=pctile(GDP_ratio) ,p(10)
	bys year: egen GDP_ratio_90=pctile(GDP_ratio) ,p(90)
	bys year: egen MAC_10=pctile(MAC) ,p(10)
	bys year: egen MAC_90=pctile(MAC) ,p(90)
	gen trunc10= GDP_ratio<GDP_ratio_10 | GDP_ratio>GDP_ratio_90 | MAC<MAC_10 | MAC>MAC_90
	
	bys year: egen GDP_ratio_25=pctile(GDP_ratio) ,p(25)
	bys year: egen GDP_ratio_75=pctile(GDP_ratio) ,p(75)
	bys year: egen MAC_25=pctile(MAC) ,p(25)
	bys year: egen MAC_75=pctile(MAC) ,p(75)
	gen trunc25= GDP_ratio<GDP_ratio_25 | GDP_ratio>GDP_ratio_75 | MAC<MAC_25 | MAC>MAC_75
//Winsorize MAC at 5% & gdp_ratio at 1%
	sum MAC, detail //5% of carbon prices which are above $1680 per tonne are set to 1680  per tonne (this is slightly more than 10% of obs after 2080)
	replace MAC=`r(p95)' if MAC>`r(p95)' & MAC!=.
	sum GDP_ratio, detail //5% of lowest GDP ratios are 95% and 1% at 90.7% of baseline GDP. 5% highest GDP ratios are at 1 (baselines)
	replace GDP_ratio=`r(p1)' if GDP_ratio<`r(p1)'
	sum lnGDP_ratio, detail
	replace lnGDP_ratio=`r(p1)' if lnGDP_ratio<`r(p1)'	
//differenced variables
	sort id t
	gen D_Emissions=d.Emissions/5
	replace D_Emissions=0 if D_Emissions>0 & D_Emissions!=.

//DEFINE PARAMETERS
	scalar delta=0.04
	scalar g=0.02
	scalar g_L=0.005
	scalar g_tau=0.015 //g+g_L-g_psi
	scalar L_0=7.403
	scalar psi_0=2
	scalar g_psi=0.01


//define stranding
	gen tau2020=1-exp(-g_tau*(year-2020)) if year>=2020
	gen r2020=(-D_Emissions+(g_tau*(1-tau2020)/(1+tau2020)-delta-g-g_L)*Emissions)/(psi_0*L_0*(1+tau2020)) 
//discount factors as weights
	gen wgt3=exp(-0.033*t) //consumption discount rate in our paper 0.011-0.005+1.35*0.02=0.033	
	gen wgt1=exp(-0.011*t) //pure time preferrence rate

//Endogenous dummy	
	gen Endo=  substr(Model, 1, 5) == "MERGE" |substr(Model, 1, 6) == "REMIND" |substr(Model, 1, 5) == "WITCH" | substr(Model, 1, 5) == "IMAGE"
//Id for Model family (fixed effects on slopes)
	gen id1 = substr(Model, 1, 3) == "AIM"
	gen id2 = substr(Model, 1, 5) == "C3IAM" | substr(Model, 1, 5) == "MERGE" | substr(Model, 1, 7) == "C-ROADS"
	gen id3 = substr(Model, 1, 4) == "GCAM"
	gen id4 = substr(Model, 1, 3) == "GEM"
	gen id5 = substr(Model, 1, 7) == "MESSAGE"
	gen id6 = substr(Model, 1, 6) == "REMIND"
	gen id7 = substr(Model, 1, 5) == "WITCH"
	gen id8 = substr(Model, 1, 5) == "POLES"
	gen id9 = substr(Model, 1, 5) == "IMAGE"
	gen id10 = substr(Model, 1, 6) == "COFFEE"
	gen id11 = substr(Model, 1, 3) == "DNE"
	
order Model_Scenario id year Emissions  D_Emissions E_Cum2020 E_Cum2100 GDP  CarbonPrice  
sort id year

save IPCC_NGFS2024, replace

//Descriptive graphs
//MAC as a function of Abatement
drop if year<2020
gen Abatement=55-Emissions
twoway 	(scatter MAC Abatement if trunc5==0  & year<2050 & MAC<.006, msize(vtiny) mcolor(red )) ///
		(lfitci  MAC Abatement if trunc5==0  & year< 2050 & MAC<.006, estopts(nocons) range(0 45) lcolor(red )) ///
		(qfit    MAC Abatement if trunc5==0  & year< 2050, estopts(nocons) range(0 45) lcolor(red) lpattern(dash)) ///
		(scatter MAC Abatement if trunc5==0 & year>=2050, msize(vtiny) mcolor(blue)) ///
		(lfitci  MAC Abatement if trunc5==0  & year>=2050, estopts(nocons) range(0 60) lcolor(blue) estopts(nocons)) ///
		(qfit    MAC Abatement if trunc5==0  & year>=2050, estopts(nocons) range(0 60) lcolor(blue) lpattern(dash)), ///
		legend(label(1 "<2050") label(2 ">=2050")) ///
		ytitle("MAC/GDP")
		graph export "G:\My Drive\Google Research Projects\Learning and negative emissions\MAC_AbatementAxis.pdf" , as(pdf) replace	
//MAC/GDP as a function of emissions
reg MAC Abatement if trunc5==0 & year <2050 , nocons 
predict lfit1
reg MAC Abatement if trunc5==0 & year >=2050 , nocons 
predict lfit2
twoway 	(scatter MAC Emissions if trunc5==0 & year< 2050  & GDP_loss_pct!=. & MAC<.006 , msize(vtiny) mcolor(red )) ///
		(qfitci  MAC Emissions if trunc5==0 & year< 2050  & GDP_loss_pct!=.,  range(10 55) lcolor(red ) fcolor(%20)) ///
		(line  lfit1 Emissions if trunc5==0 & year <2050  & GDP_loss_pct!=. & Emissions<=55 & Emissions>10,   sort lcolor(red) lpattern(dash)) ///
		(scatter MAC Emissions if trunc5==0 & year>=2050  & GDP_loss_pct!=. & MAC<.006, msize(vtiny) mcolor(blue)) ///
		(qfitci  MAC Emissions if trunc5==0 & year>=2050  & GDP_loss_pct!=. ,  range(-5 55) lcolor(blue) fcolor(%20)) ///
		(line  lfit2 Emissions if trunc5==0 & year>=2050  & GDP_loss_pct!=. & Emissions<=55 & Emissions>-5,   sort lcolor(blue) lpattern(dash)), ///
		legend(label(1 "<2050") label(2 "<2050 quadratic 95%CI") label(3 "<2050 quadratic fit") label(4 ">=2050") label(5 "<2050 linear fit") label(6 ">=2050 quadratic 95%CI") label(7 ">=2050 quadratic fit") label(8 ">=2050 linear fit")) ///
		ytitle("MAC/GDP")
graph export "G:\My Drive\Google Research Projects\Learning and negative emissions\MAC_trunc5.pdf" , as(pdf) replace	
//MAC () as a function of emissions
reg CarbonPrice Abatement if trunc5==0 & year <2050 , nocons 
predict lfit_1
reg CarbonPrice Abatement if trunc5==0 & year >=2050 , nocons 
predict lfit_2
twoway 	(scatter CarbonPrice Emissions if trunc5==0 & year< 2050  & GDP_loss_pct!=. & CarbonPrice<2000 , msize(vtiny) mcolor(red )) ///
		(qfitci  CarbonPrice Emissions if trunc5==0 & year< 2050  & GDP_loss_pct!=.,  range(10 55) lcolor(red ) fcolor(%20)) ///
		(line  lfit_1 Emissions if trunc5==0 & year <2050  & GDP_loss_pct!=. & Emissions<=55 & Emissions>10,   sort lcolor(red) lpattern(dash)) ///
		(scatter CarbonPrice Emissions if trunc5==0 & year>=2050  & GDP_loss_pct!=. & CarbonPrice<2000, msize(vtiny) mcolor(blue)) ///
		(qfitci  CarbonPrice Emissions if trunc5==0 & year>=2050  & GDP_loss_pct!=. ,  range(-5 55) lcolor(blue) fcolor(%20)) ///
		(line  lfit_2 Emissions if trunc5==0 & year>=2050  & GDP_loss_pct!=. & Emissions<=55 & Emissions>-5,   sort lcolor(blue) lpattern(dash)), ///
		legend(label(1 "<2050") label(2 "<2050 quadratic 95%CI") label(3 "<2050 quadratic fit") label(4 ">=2050") label(5 "<2050 linear fit") label(6 ">=2050 quadratic 95%CI") label(7 ">=2050 quadratic fit") label(8 ">=2050 linear fit")) ///
		ytitle("MAC ($/tCO2)")
graph export "G:\My Drive\Google Research Projects\Learning and negative emissions\MAC_Absolute_trunc5.pdf" , as(pdf) replace	



//Total Abatement Costs	as a percentage of GDP
gen Emissions_sq=Emissions^2
gen Emissions_cu=Emissions^3
reg GDP_loss_pct Emissions Emissions_sq Emissions_cu if trunc5==0  & year<2050
predict cfit1 
reg GDP_loss_pct Emissions Emissions_sq Emissions_cu if trunc5==0  & year>=2050
predict cfit2 
twoway 	(scatter GDP_loss_pct Emissions if trunc5==0 & year <2050 & GDP_loss_pct<.06, msize(vtiny) mcolor(red )) ///
		(qfitci GDP_loss_pct Emissions if trunc5==0  & year <2050, range(10 55) lcolor(red ) fcolor(%20)) ///
		(line cfit1 Emissions if trunc5==0 & year <2050 & Emissions<=55 & Emissions>10, sort lcolor(red) lpattern(dash)) ///
		(scatter GDP_loss_pct Emissions if trunc5==0 & year>=2050 & GDP_loss_pct<.06, msize(vtiny) mcolor(blue)) ///
		(qfitci GDP_loss_pct Emissions if trunc5==0  & year>=2050, range(-10 55) lcolor(blue) fcolor(%20)) ///
		(line cfit2 Emissions if trunc5==0 & year>=2050 & Emissions<=55 & Emissions>-10, sort lcolor(blue) lpattern(dash)), ///
		legend(label(1 "<2050") label(2 "<2050 quadratic 95%CI") label(3 "<2050 quadratic fit")  label(4 "<2050 cubic fit") label(5 ">=2050") label(6 ">=2050 quadratic 95%CI") label(7 ">=2050 quadratic fit") label(8 ">=2050 cubic fit")) ///
		ytitle("Total Abatement Costs/GDP")
graph export "G:\My Drive\Google Research Projects\Learning and negative emissions\TotalAbatementCostsTrunc5.pdf" , as(pdf) replace
//Total Abatement costs (absolute)
reg Loss_GDP Emissions Emissions_sq Emissions_cu if  trunc5==0  & year<2050
predict cfit_11 
reg Loss_GDP Emissions Emissions_sq Emissions_cu if trunc5==0  &  year>=2050
predict cfit_21
twoway 	(scatter Loss_GDP Emissions if trunc5==0  &  year <2050 & Loss_GDP<15, msize(.1pt) mcolor(red )) ///
		(scatter Loss_GDP Emissions if trunc5==0  &  year>=2050 & Loss_GDP<15, msize(.1pt) mcolor(blue)) ///
		(qfitci Loss_GDP Emissions if trunc5==0  &  year <2050, range(10 55) lcolor(red ) fcolor(%20)) ///
		(line cfit_1 Emissions if trunc5==0  &  year <2050 & Emissions<=55 & Emissions>10, sort lcolor(red) lpattern(dash)) ///
		(qfitci Loss_GDP Emissions if trunc5==0  &  year>=2050, range(-10 55) lcolor(blue) fcolor(%20)) ///
		(line cfit_2 Emissions if trunc5==0  &  year>=2050  & Emissions>-10 & Emissions<=55, sort lcolor(blue) lpattern(dash)), ///
		legend(label(1 "<2050") label(3 "<2050 quadratic 95%CI") label(4 "<2050 quadratic fit")  label(5 "<2050 cubic fit") label(2 ">=2050") label(6 ">=2050 quadratic 95%CI") label(7 ">=2050 quadratic fit") label(8 ">=2050 cubic fit")) ///
		ytitle("Total Abatement Costs (Trillion USD)")
graph export "G:\My Drive\Google Research Projects\Learning and negative emissions\TotalAbatementCostsAbsoluteTrunc5.pdf", as(pdf) replace