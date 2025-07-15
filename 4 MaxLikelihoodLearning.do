cd "G:\My Drive\Google Research Projects\Stranded Assets2\Code\Stata\"
set more off
use IPCC_NGFS2024, replace
//DEFINE PARAMETERS
	scalar E_Cum2100_max=1.5/0.0006 //only to use scenarios that stay below 2.7°C in 2100 (assuming 1.2°C warming in 2020)
	scalar E_BAU=55 // alternative is to take BAU emissions from the database. Fixed BAU is better because it is an assumption of our model. We fit the relationship between GDPloss and emissions, irrespective of BAU emissions. 
	scalar A_Cum0=100
	scalar weightMAC=1 //Weight of MAC equation relative to the total abatement cost (which has weight 1)
	gen A_Cum2020=t2020*E_BAU-E_Cum2020 //does not include A_Cum0
	replace A_Cum2020=0 if A_Cum2020<0
	gen Abatement=E_BAU-Emissions
	replace Abatement=0 if Abatement<0
	replace D_Emissions=0 if D_Emissions>0
//NO LEARNING
capture program drop NoLearn
program define NoLearn
	args lnf varphi //lnf is output of the function, varphi is the input
	tempvar error1 error2  	
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*`varphi' * Abatement^2)
	sum `error1'
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	
	gen `error2'= MAC      -  `varphi'*Abatement
	sum `error2'
	scalar stdev2=`r(sd)'
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)
	}
end
ml model lf NoLearn /varphi  if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. , cluster(id)  maximize init(0.00002 , copy) difficult //technique(nr 20 bhhh 20 bfgs 20 dfp 20)
ml display
scalar varphiNL=_b[varphi]
estimates store NoLearn1

//quadratic MAC
capture program drop NoLearn_sq
program define NoLearn_sq
	args lnf varphi varphi_sq theta2
	tempvar error1 error2  	
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*`varphi' * Abatement^2-0.33*`varphi_sq'*Abatement^3-`theta2'*D_Emissions^2 )
	sum `error1'
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  `varphi'*Abatement - `varphi_sq'*Abatement^2
	sum `error2'
	scalar stdev2=`r(sd)'
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)
	}
end
ml model lf NoLearn_sq /varphi /varphi_sq /theta2 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. , cluster(id)  maximize init(0.000032 0 0.0017, copy) difficult //technique(nr 20 bhhh 20 bfgs 20 dfp 20)
ml display
scalar varphi_lin=_b[varphi]
scalar varphi_sq=_b[varphi_sq]
estimates store NoLearn_sq


//add speed penalty estimated by the model
capture program drop NoLearn2
program define NoLearn2
	args lnf varphi theta2 
	tempvar error1 error2  	
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*`varphi' * Abatement^2-`theta2'*D_Emissions^2 )
	sum `error1'
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  `varphi'*Abatement
	sum `error2'
	scalar stdev2=`r(sd)'
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)
	}
end
ml model lf NoLearn2 /varphi /theta2 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=.  , cluster(id)  maximize init(0.000032 0.001, copy) difficult //technique(nr 20 bhhh 20 bfgs 20 dfp 20)
ml display
scalar varphiNL2=_b[varphi]
estimates store NoLearn2

/*//add model-specific slope and a quadratic term: 9 out of 11 slopes become negative and quadratic term is positive and significant (also if thetat2=0)
capture program drop NoLearn3
program define NoLearn3
	args lnf varphi1 varphi2 varphi3 varphi4 varphi5 varphi6 varphi7 varphi8 varphi9 varphi10 varphi11 varphi_sq  theta2
	tempvar error1 error2 
 	local varphi "id1*`varphi1' + id2*`varphi2' + id3*`varphi3' + id4*`varphi4' + id5*`varphi5' + id6*`varphi6' + id7*`varphi7' + id8*`varphi8' + id9*`varphi9' + id10*`varphi10' + id11*`varphi11' " 	
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*`varphi' * Abatement^2 -0.33*`varphi_sq' * Abatement^3 -`theta2'*D_Emissions^2 )
	sum `error1'
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  `varphi'*Abatement -`varphi_sq'*Abatement^2
	sum `error2'
	scalar stdev2=`r(sd)'
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)
	}
end
ml model lf NoLearn3 /varphi1 /varphi2 /varphi3 /varphi4 /varphi5 /varphi6 /varphi7 /varphi8 /varphi9 /varphi10 /varphi11  /varphi_sq /theta2 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.  , cluster(id)  maximize init(0.000069 0.00013 0.00011 0.000043 0.000044 0.000053 0.000072 0.00018 0.00010 0.00002 0.00019  0.00000 0, copy) difficult 
ml display
estimates store NoLearn3 */


//EXOGENOUS LEARNING////////////////////////////////////////////////////////////////////////////////////////////
//zero steady state for varphi
capture program drop ExoLearn
program define ExoLearn
	args lnf varphi g_varphi  
	tempvar error1 error2  	
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*`varphi'*exp(-`g_varphi'*t2020) * Abatement^2)   
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  `varphi'*exp(-`g_varphi'*t2020)*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)  
	}
end
ml model lf ExoLearn /varphi /g_varphi if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=.  , cluster(id)  maximize init(0.00004 0.01 , copy) difficult 
ml display
scalar varphiExo=_b[varphi]
scalar g_varphi=_b[g_varphi]
estimates store ExoLearn

//estimate steady state varphi_inf
capture program drop ExoLearn2
program define ExoLearn2
	args lnf varphi g_varphi varphi_inf 
	tempvar error1 error2  	
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*(`varphi_inf'+(`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020)) * Abatement^2)
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  (`varphi_inf'+(`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020))*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)  
	}
end
ml model lf ExoLearn2 /varphi /g_varphi /varphi_inf if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=.  , cluster(id)  maximize init(0.00005 0.02 0.00001, copy) difficult //round(year/10)==year/10
ml display
scalar varphi2=_b[varphi]
scalar varphi_inf2=_b[varphi_inf]
scalar g_varphi2=_b[g_varphi]
estimates store ExoLearn2

//speed penalty (and estimate steady state varphi_inf)
capture program drop ExoLearn3
program define ExoLearn3
	args lnf varphi g_varphi varphi_inf theta2
	tempvar error1 error2  	
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*(`varphi_inf'+(`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020)) * Abatement^2 -`theta2'*D_Emissions^2 )
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  (`varphi_inf'+(`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020))*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)  
	}
end
ml model lf ExoLearn3 /varphi /g_varphi /varphi_inf /theta2 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.  , cluster(id)  maximize init(0.00005 0.02 0.00001 0.00001, copy) difficult //round(year/10)==year/10
ml display
scalar varphi3=_b[varphi]
scalar varphi_inf3=_b[varphi_inf]
scalar g_varphi3=_b[g_varphi]
scalar theta2=_b[theta2]
estimates store ExoLearn3


//only exo models (speed penalty & steady state varphi_inf)
capture program drop ExoLearn4
program define ExoLearn4
	args lnf varphi g_varphi varphi_inf theta2
	tempvar error1 error2  	
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*(`varphi_inf'+(`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020)) * Abatement^2 -`theta2'*D_Emissions^2 )
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  (`varphi_inf'+(`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020))*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)  
	}
end
ml model lf ExoLearn4 /varphi /g_varphi /varphi_inf /theta2 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.   & Endo==0, cluster(id)  maximize init(0.00005 0.01 0.00003 0.001, copy) difficult //round(year/10)==year/10
ml display
scalar varphi4=_b[varphi]
scalar varphi_inf4=_b[varphi_inf]
scalar g_varphi4=_b[g_varphi]
estimates store ExoLearn4

//a different slope per model group (and speed penalty and varphi_inf)
capture program drop ExoLearn5
program define ExoLearn5
	args lnf varphi1 varphi2 varphi3 varphi4 varphi5 varphi6 varphi7 varphi8 varphi9 varphi10 varphi11     g_varphi varphi_inf theta2
	tempvar error1 error2 
 	local varphi "id1*`varphi1' + id2*`varphi2' + id3*`varphi3' + id4*`varphi4' + id5*`varphi5' + id6*`varphi6' + id7*`varphi7' + id8*`varphi8' + id9*`varphi9' + id10*`varphi10' + id11*`varphi11' " 
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*(`varphi_inf'+(`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020)) * Abatement^2 -`theta2'*D_Emissions^2 )
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  (`varphi_inf'+(`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020))*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)  
	}
end
ml model lf ExoLearn5 /varphi1 /varphi2 /varphi3 /varphi4 /varphi5 /varphi6 /varphi7 /varphi8 /varphi9 /varphi10 /varphi11   /g_varphi /varphi_inf /theta2 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.  , cluster(id)  maximize init(0.00007 0.00014 0.00010 0.00005 0.00005 0.00005 0.00007 0.00020 0.00010 0.00001 0.00020 0.02 0.000001 0.0001, copy) difficult 
ml display
estimates store ExoLearn5

//model-specific slope, with quadratic term (using varphi_inf and theta2 of preceding model)
estimates restore ExoLearn5
capture program drop ExoLearn6
program define ExoLearn6
	local varphi_inf= _b[varphi_inf]
	local theta2=_b[theta2]
	args lnf varphi1 varphi2 varphi3 varphi4 varphi5 varphi6 varphi7 varphi8 varphi9 varphi10 varphi11  g_varphi varphi_sq
	tempvar error1 error2 
 	local varphi "id1*`varphi1' + id2*`varphi2' + id3*`varphi3' + id4*`varphi4' + id5*`varphi5' + id6*`varphi6' + id7*`varphi7' + id8*`varphi8' + id9*`varphi9' + id10*`varphi10' + id11*`varphi11' " 
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*(`varphi_inf'+(`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020)) * Abatement^2 - 0.333*`varphi_sq'*exp(-`g_varphi'*t2020)*Abatement^3 -`theta2'*D_Emissions^2  )
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  (`varphi_inf'+(`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020))*Abatement + `varphi_sq'*exp(-`g_varphi'*t2020)*Abatement^2
	sum `error2' 
	scalar stdev2=`r(sd)'
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)  - (`varphi_sq'<0)*`varphi_sq'^2*10000    
	}
end
ml model lf ExoLearn6 /varphi1 /varphi2 /varphi3 /varphi4 /varphi5 /varphi6 /varphi7 /varphi8 /varphi9 /varphi10 /varphi11   /g_varphi  /varphi_sq if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.  , cluster(id)  maximize init(0.000069 0.00013 0.00011 0.000043 0.000044 0.000053 0.000072 0.00018 0.00010 0.00002 0.00019 0.013 0.00000, copy) difficult 
ml display
estimates store ExoLearn6

//a different slope and g_varphi per model group 
capture program drop ExoLearn7
program define ExoLearn7
	args lnf varphi1 varphi2 varphi3 varphi4 varphi5 varphi6 varphi7 varphi8 varphi9 varphi10 varphi11  g_varphi1 g_varphi2 g_varphi3 g_varphi4 g_varphi5 g_varphi6 g_varphi7 g_varphi8 g_varphi9 g_varphi10 g_varphi11  varphi_inf theta2
	tempvar error1 error2 
 	local varphi "id1*`varphi1' + id2*`varphi2' + id3*`varphi3' + id4*`varphi4' + id5*`varphi5' + id6*`varphi6' + id7*`varphi7' + id8*`varphi8' + id9*`varphi9' + id10*`varphi10' + id11*`varphi11' " 
	local g_varphi "id1*`g_varphi1' + id2*`g_varphi2' + id3*`g_varphi3' + id4*`g_varphi4' + id5*`g_varphi5' + id6*`g_varphi6' + id7*`g_varphi7' + id8*`g_varphi8' + id9*`g_varphi9' + id10*`g_varphi10' + id11*`g_varphi11' " 
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*(`varphi_inf'+(`varphi'-`varphi_inf')*exp(-(`g_varphi')*t2020)) * Abatement^2 -`theta2'*D_Emissions^2 )
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  (`varphi_inf'+(`varphi'-`varphi_inf')*exp(-(`g_varphi')*t2020))*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)  
	}
end
ml model lf ExoLearn7 /varphi1 /varphi2 /varphi3 /varphi4 /varphi5 /varphi6 /varphi7 /varphi8 /varphi9 /varphi10 /varphi11 /g_varphi1 /g_varphi2 /g_varphi3 /g_varphi4 /g_varphi5 /g_varphi6 /g_varphi7 /g_varphi8 /g_varphi9 /g_varphi10 /g_varphi11 /varphi_inf /theta2 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.  , cluster(id)  maximize init(0.00007 0.00014 0.00010 0.00005 0.00005 0.00005 0.00007 0.00020 0.00010 0.00001 0.00020 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.000001 0.0001, copy) difficult 
ml display
scalar varphi7=(_b[varphi1] + _b[varphi2] + _b[varphi3] + _b[varphi4] + _b[varphi5] + _b[varphi6] + _b[varphi7] + _b[varphi8] + _b[varphi9] + _b[varphi10] + _b[varphi11])/11
scalar g_varphi7=(_b[g_varphi1] + _b[g_varphi2] + _b[g_varphi3] + _b[g_varphi4] + _b[g_varphi5] + _b[g_varphi6] + _b[g_varphi7] + _b[g_varphi8] + _b[g_varphi9] + _b[g_varphi10] + _b[g_varphi11])/11
display g_varphi7 varphi7
estimates store ExoLearn7


/*
//Quadratic term and both model-specific varphi and varphi_inf. (Speed penalty from ExoLearn5 and varphi_inf=0) Gives the same g_varphi of 1.2%
estimates restore ExoLearn5
capture program drop ExoLearn7
program define ExoLearn7
	local varphi_inf= _b[varphi_inf]
	local theta2=_b[theta2]
	args lnf varphi ff2 ff3 ff4 ff5 ff6 ff7 ff8 ff9 ff10 ff11  varphi_sq   g_varphi 
	tempvar error1 error2 
 	local f_effect "(id1 + id2*`ff2' + id3*`ff3' + id4*`ff4' + id5*`ff5' + id6*`ff6' + id7*`ff7' + id8*`ff8' + id9*`ff9' + id10*`ff10' + id11*`ff11')" 
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*((`varphi_inf' + (`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020)) * Abatement^2 )* `f_effect' - 0.333*`varphi_sq'*exp(-`g_varphi'*t2020)*Abatement^3 -`theta2'*D_Emissions^2 ) 
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC   -  (`varphi_inf' + (`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020))*Abatement *`f_effect' + `varphi_sq'*exp(-`g_varphi'*t2020)*Abatement^2
	sum `error2' 
	scalar stdev2=`r(sd)'
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC) - (`varphi_sq'<0)*`varphi_sq'^2*1000    
	}
end
ml model lf ExoLearn7 /varphi /ff2 /ff3 /ff4 /ff5 /ff6 /ff7 /ff8 /ff9 /ff10 /ff11 /varphi_sq /g_varphi if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.  , cluster(id)  maximize init(0.000068 1.83 1.53 0.67 0.66 0.82 1.06 2.13 0.97 0.38 2.46  0  0.011 , copy) difficult 
ml display
estimates store ExoLearn7 
*/
//Model-specific slope and quadratic term (gives negative slope). (Speed penalty from ExoLearn5 and varphi_inf=0) 
/*estimates restore ExoLearn5
capture program drop ExoLearn6b
program define ExoLearn6b
	local varphi_inf= 0
	local theta2=_b[theta2]
	args lnf varphi ff1 ff2 ff3 ff4 ff5 ff6 ff7 ff8 ff9 ff10 ff11  varphi_sq   g_varphi 
	tempvar error1 error2 
 	local f_effect "(id1*`ff1' + id2*`ff2' + id3*`ff3' + id4*`ff4' + id5*`ff5' + id6*`ff6' + id7*`ff7' + id8*`ff8' + id9*`ff9' + id10*`ff10' + id11*`ff11')" 
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*(`varphi_inf' + (`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020)) * Abatement^2 - 0.333*`varphi_sq'*exp(-`g_varphi'*t2020)*Abatement^3 -`theta2'*D_Emissions^2 ) * `f_effect'
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC   -  ((`varphi_inf' + (`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020))*Abatement + `varphi_sq'*exp(-`g_varphi'*t2020)*Abatement^2) *`f_effect'
	sum `error2' 
	scalar stdev2=`r(sd)'
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC) - (`varphi_sq'<0)*`varphi_sq'^2 - (`varphi'<0)*`varphi'^2
	}
end
ml model lf ExoLearn6b /varphi /ff1 /ff2 /ff3 /ff4 /ff5 /ff6 /ff7 /ff8 /ff9 /ff10 /ff11 /varphi_sq /g_varphi if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.  , cluster(id)  maximize init(0.0001 1 1 1 1 1 1 1 1 1 1 1 0.000001  0.01 , copy) difficult 
ml display
estimates store ExoLearn6b
*/

//ENDOGENOUS LEARNING///////////////////////////////////////////////////////////////////////////////////
//impose A_Cum0=100
capture program drop EndoLearn
program define EndoLearn
	args lnf varphi chi  
	tempvar error1  error2	
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*`varphi'*((A_Cum0+A_Cum2020)/A_Cum0)^(-`chi') * Abatement^2)   
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  `varphi'*((A_Cum0+A_Cum2020)/A_Cum0)^(-`chi')*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'	
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)
	}
end
ml model lf EndoLearn /varphi /chi if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=.  , cluster(id)  maximize init(0.00005 0.3 , copy) difficult 
ml display
scalar varphiEndo=_b[varphi]
scalar chi=_b[chi]
estimates store EndoLearn 

/*//estimate A_Cum0 by the data
capture program drop EndoLearn2
program define EndoLearn2
	args lnf varphi chi A_Cum0 
	tempvar error1  error2	
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*`varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-`chi') * Abatement^2)   
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  `varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-`chi')*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'	
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)
	}
end
ml model lf EndoLearn2 /varphi /chi /A_Cum0 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=.  , cluster(id)  maximize init(0.00005 0.3 50, copy) difficult 
ml display
scalar varphiEndo2=_b[varphi]
scalar chi2=_b[chi]
scalar A_Cum02=_b[A_Cum0]
estimates store EndoLearn2 */
  
//add speed penality and estimate A_Cum0 by the data
capture program drop EndoLearn3
program define EndoLearn3
	args lnf varphi chi A_Cum0 theta2
	tempvar error1  error2	
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*`varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-`chi') * Abatement^2-`theta2'*D_Emissions^2)   
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  `varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-`chi')*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'	
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)
	}
end
ml model lf EndoLearn3 /varphi /chi /A_Cum0 /theta2 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.  , cluster(id)  maximize init(0.00005 0.3 50 0.0001, copy) difficult 
ml display
scalar varphiEndo3=_b[varphi]
scalar chi3=_b[chi]
scalar A_Cum03=_b[A_Cum0]
scalar theta2endo=_b[theta2]
estimates store EndoLearn3

//Only endo models, speed penality (estimate A_Cum0 by the data)
capture program drop EndoLearn4
program define EndoLearn4
	args lnf varphi chi A_Cum0 theta2
	tempvar error1  error2	
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*`varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-`chi') * Abatement^2-`theta2'*D_Emissions^2)   
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  `varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-`chi')*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'	
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)
	}
end
ml model lf EndoLearn4 /varphi /chi /A_Cum0 /theta2 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.   & Endo==1, cluster(id)  maximize init(0.00005 0.3 50 0.0001, copy) difficult 
ml display
estimates store EndoLearn4

//specific slope for each model family (and speed penality and estimate A_Cum0 )
capture program drop EndoLearn5
program define EndoLearn5
	args lnf varphi1 varphi2 varphi3 varphi4 varphi5 varphi6 varphi7 varphi8 varphi9 varphi10 varphi11    chi A_Cum0 theta2
	tempvar error1  error2	
	local varphi " (id1*`varphi1' + id2*`varphi2' + id3*`varphi3' + id4*`varphi4' + id5*`varphi5' + id6*`varphi6' + id7*`varphi7' + id8*`varphi8' + id9*`varphi9' + id10*`varphi10' + id11*`varphi11') " 
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*`varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-`chi') * Abatement^2 -`theta2'*D_Emissions^2)   
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  `varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-`chi')*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'	
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)
	}
end
ml model lf EndoLearn5 /varphi1 /varphi2 /varphi3 /varphi4 /varphi5 /varphi6 /varphi7 /varphi8 /varphi9 /varphi10 /varphi11 /chi /A_Cum0 /theta2 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.  , cluster(id)  maximize  init(0.00005 0.00005 0.00005 0.00005 0.00005 0.00005 0.00005 0.00005 0.00005 0.00005 0.00005 0.3 100 0.00001, copy) difficult //technique(bhhh)
ml display
estimates store EndoLearn5

//Quadratic term (learning parameters fixed)  
estimates restore EndoLearn5  //in case EndoLearn5 was not the last model that was estimated
	local A_Cum0=_b[A_Cum0]
	local theta2=_b[theta2]
capture program drop EndoLearn6
program define EndoLearn6
	args lnf varphi1 varphi2 varphi3 varphi4 varphi5 varphi6 varphi7 varphi8 varphi9 varphi10 varphi11 varphi_sq chi 
	tempvar error1  error2	
	local varphi " (id1*`varphi1' + id2*`varphi2' + id3*`varphi3' + id4*`varphi4' + id5*`varphi5' + id6*`varphi6' + id7*`varphi7' + id8*`varphi8' + id9*`varphi9' + id10*`varphi10' + id11*`varphi11') " 
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*`varphi'*((A_Cum0+A_Cum2020)/A_Cum0)^(-`chi') * Abatement^2 - 0.333*`varphi_sq'*Abatement^3 - theta2*D_Emissions^2)   
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  (`varphi'*Abatement - `varphi_sq'*Abatement^2)*((A_Cum0+A_Cum2020)/A_Cum0)^(-`chi')
	sum `error2' 
	scalar stdev2=`r(sd)'	
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)
	}
end
ml model lf EndoLearn6 /varphi1 /varphi2 /varphi3 /varphi4 /varphi5 /varphi6 /varphi7 /varphi8 /varphi9 /varphi10 /varphi11 /varphi_sq /chi if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.  , cluster(id)  maximize  init( 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5 1e-5 1e-6 0.3, copy) difficult //technique(bhhh)
ml display
estimates store EndoLearn6  

//specific slope and learning rate for each model family (and speed penality and estimate A_Cum0 )
capture program drop EndoLearn7
program define EndoLearn7
	args lnf varphi1 varphi2 varphi3 varphi4 varphi5 varphi6 varphi7 varphi8 varphi9 varphi10 varphi11  chi1 chi2 chi3 chi4 chi5 chi6 chi7 chi8 chi9 chi10 chi11 A_Cum0 theta2
	tempvar error1  error2	
	local varphi " (id1*`varphi1' + id2*`varphi2' + id3*`varphi3' + id4*`varphi4' + id5*`varphi5' + id6*`varphi6' + id7*`varphi7' + id8*`varphi8' + id9*`varphi9' + id10*`varphi10' + id11*`varphi11') " 
	local chi " (id1*`chi1' + id2*`chi2' + id3*`chi3' + id4*`chi4' + id5*`chi5' + id6*`chi6' + id7*`chi7' + id8*`chi8' + id9*`chi9' + id10*`chi10' + id11*`chi11') " 
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*`varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-(`chi')) * Abatement^2 -`theta2'*D_Emissions^2)   
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      -  `varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-(`chi'))*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'	
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC)
	}
end
ml model lf EndoLearn7 /varphi1 /varphi2 /varphi3 /varphi4 /varphi5 /varphi6 /varphi7 /varphi8 /varphi9 /varphi10 /varphi11 /chi1 /chi2 /chi3 /chi4 /chi5 /chi6 /chi7 /chi8 /chi9 /chi10 /chi11 /A_Cum0 /theta2 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.  , cluster(id)  maximize  init(0.00005 0.00005 0.00005 0.00005 0.00005 0.00005 0.00005 0.00005 0.00005 0.00005 0.00005 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 100 0.00001, copy) difficult //technique(bhhh)
ml display
scalar endovarphi7=(_b[varphi1] + _b[varphi2] + _b[varphi3] + _b[varphi4] + _b[varphi5] + _b[varphi6] + _b[varphi7] + _b[varphi8] + _b[varphi9] + _b[varphi10] + _b[varphi11])/11
scalar chi7=(_b[chi1] + _b[chi2] + _b[chi3] + _b[chi4] + _b[chi5] + _b[chi6] + _b[chi7] + _b[chi8] + _b[chi9] + _b[chi10] + _b[chi11])/11
display chi7 endovarphi7
estimates store EndoLearn7


/*// flexible power for GDP (and model-specific slope,  speed penality and estimate A_Cum0 ): No convergence (coefficient tends to be larger than 1)
//I use the log version of the TAC equation. Errors are now multiplicative, so that the algorithm puts lower weight on errors of large TAC
capture program drop EndoLearn7
program define EndoLearn7
	args lnf varphi	GDP_power   chi A_Cum0 theta2
	tempvar error1  error2	
	//local varphi " (id1*`varphi1' + id2*`varphi2' + id3*`varphi3' + id4*`varphi4' + id5*`varphi5' + id6*`varphi6' + id7*`varphi7' + id8*`varphi8' + id9*`varphi9' + id10*`varphi10' + id11*`varphi11') " 
	quietly {
	gen `error1'=ln(1-Loss_GDP/GDP^`GDP_power') - (-0.5*`varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-`chi') * Abatement^2 -`theta2'*D_Emissions^2)   
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= CarbonPrice/GDP^`GDP_power' - `varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-`chi')*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'	
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*weightMAC)/(1+weightMAC) - (`GDP_power'<0.3)*`GDP_power'^2*10 - (`chi'<0)*`chi'^2*1000
	}
end
ml model lf EndoLearn7 /varphi /GDP_power /chi /A_Cum0 /theta2 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.  , cluster(id)  maximize  init(0.00008 1 0.3 100 0.00001, copy) difficult //technique(bhhh)
ml display
estimates store EndoLearn7*/

/*//Combine endo and exo learning (does not work well, hard to get a good estimate for chi)
capture program drop EndoExoLearn
scalar omega=0.5
program define EndoExoLearn
	args lnf varphi varphi_inf g_varphi chi A_Cum0 theta2
	
	tempvar error1  error2	
	quietly {
	gen `error1'=GDP_ratio -  exp(-0.5*(omega *(`varphi_inf'+(`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020)) +(1-omega)*`varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-`chi'))* Abatement^2-`theta2'*D_Emissions^2)   
	sum `error1' 
	scalar stdev1=`r(sd)'
	replace `error1'=0 if `error1'==.
	gen `error2'= MAC      - (omega *(`varphi_inf'+(`varphi'-`varphi_inf')*exp(-`g_varphi'*t2020)) +(1-omega)* `varphi'*((`A_Cum0'+A_Cum2020)/`A_Cum0')^(-`chi'))*Abatement 
	sum `error2' 
	scalar stdev2=`r(sd)'	
	replace `lnf' = (ln(normalden(`error1',stdev1)) + ln(normalden(`error2',stdev2))*1)/(1+1) - (`chi'<0)*`chi'^2*10000 - (`chi'>0.5)*(`chi'-0.5)^2*10000 - (`A_Cum0'<0)*`A_Cum0'^2*1000  //I added non-negativity penalties to avoid a local maximum
	}
end
ml model lf EndoExoLearn /varphi /g_varphi /varphi_inf /chi /A_Cum0 /theta2 if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=. & D_Emissions!=.  , cluster(id)  maximize init(0.00005 0.01 0.00001 0.3 20 0.000001, copy) difficult 
ml display
estimates store EndoExoLearn */ 

//OUTPUT
foreach xxx in NoLearn2 ExoLearn3 ExoLearn4 ExoLearn5 ExoLearn6 ExoLearn7 EndoLearn3 EndoLearn4 EndoLearn5 EndoLearn6 EndoLearn7{
	estimates save `xxx' , replace
}
foreach xxx in NoLearn2 ExoLearn3 ExoLearn4 ExoLearn5 ExoLearn6 ExoLearn7 EndoLearn3 EndoLearn4 EndoLearn5 EndoLearn6 EndoLearn7{
	estimates use `xxx' 
}

estimates table ExoLearn7 EndoLearn7 //this is the model we use for the graphs
estimates table NoLearn2 ExoLearn3 ExoLearn4 ExoLearn5 ExoLearn6 EndoLearn3 EndoLearn4 EndoLearn5 EndoLearn6, stats(N ll BIC AIC)  b(%5.0g)  se  
estimates table NoLearn2 ExoLearn3 ExoLearn4 ExoLearn5 ExoLearn6 ExoLearn7 EndoLearn3 EndoLearn4 EndoLearn5 EndoLearn6 EndoLearn7, stats(N ll BIC AIC)  b(%5.0g) star(0.1 0.05 0.01) // se  
 


//GRAPH of MAC and Learning factor
bys t:egen A_Cum2020_Mean=mean(A_Cum2020) 	if E_Cum2100<E_Cum2100_max & MAC!=.   & (year==2020|year==2030|year==2040|year==2050|year==2060|year==2070|year==2080|year==2090|year==2100)
bys t:egen E_Mean=mean(Emissions) 			if E_Cum2100<E_Cum2100_max & MAC!=.   & (year==2020|year==2030|year==2040|year==2050|year==2060|year==2070|year==2080|year==2090|year==2100)
gen MACnoLearn=varphiNL						*(E_BAU-E_Mean) 
gen MACexo=varphiExo*exp(-g_varphi*t2020)	*(E_BAU-E_Mean)
gen MACexo2=(varphi_inf2+(varphi2-varphi_inf2)*exp(-g_varphi2*t2020))*(E_BAU-E_Mean)
gen MACexo3=(varphi_inf2+(varphi2-varphi_inf2)*exp(-g_varphi2*t2020))*(E_BAU-E_Mean)
gen MACendo=varphiEndo*((A_Cum2020_Mean+A_Cum0)/A_Cum0)^(-chi)	*(E_BAU-E_Mean)
gen MACendo2=varphiEndo2*((A_Cum2020_Mean+A_Cum02)/A_Cum02)^(-chi2)	*(E_BAU-E_Mean)
gen MACendo3=varphiEndo3*((A_Cum2020_Mean+A_Cum03)/A_Cum03)^(-chi3)	*(E_BAU-E_Mean)
gen varphiNL=varphiNL
gen varphiNL_sq=varphi_lin+2*varphi_sq*(E_BAU-E_Mean)
gen varphiNL2=varphiNL2
gen varphiExo=varphiExo*exp(-g_varphi*t2020)
gen varphiExo2=(varphi_inf2+(varphi2-varphi_inf2)*exp(-g_varphi2*t2020))
gen varphiExo3=(varphi_inf3+(varphi3-varphi_inf3)*exp(-g_varphi3*t2020))
gen varphiEndo=varphiEndo*((A_Cum2020_Mean+A_Cum0)/A_Cum0)^(-chi)
gen varphiEndo2=varphiEndo2*((A_Cum2020_Mean+A_Cum02)/A_Cum02)^(-chi2)
gen varphiEndo3=varphiEndo3*((A_Cum2020_Mean+A_Cum02)/A_Cum03)^(-chi3)

line MACnoLearn MACexo MACexo2 MACendo MACendo2 E_Mean if id==2&year>=2020 //endogenous costs depend on scenario (id)
line varphiExo varphiExo2 varphiEndo varphiEndo2 varphiNL year if id==2&year>=2020, graphregion(color(white)) ytitle("varphi(t,A)") ///
legend(label(1 "Exogenous, steady state=0") label(2 "Exogenous, fitted steady state") label(3 "Endogenous, A_Cum0=300") label(4 "Endogenous, fitted A_Cum0") label(5 "No Learning") cols(1))
graph export "C:\Users\frank\Google Drive\Google Research Projects\Learning and negative emissions\varphi.pdf" , as(pdf) replace
line varphiExo2 varphiExo3 varphiEndo2 varphiEndo3 varphiNL varphiNL2 year if id==2&year>=2020, graphregion(color(white)) ytitle("varphi(t,A)") ///
legend( label(1 "Exogenous, fitted steady state") label(2 "Exogenous, speed penalty")  label(3 "Endogenous, fitted A_Cum0") label(4 "Endogenous, speed penalty") label(5 "No Learning") label(6 "No Learning, quadratic term") cols(1))
graph export "C:\Users\frank\Google Drive\Google Research Projects\Learning and negative emissions\varphiSpeedPenalty.pdf" , as(pdf) replace


/*
//ONLY TAC or MAC fit
//EXOGENOUS
//Use only Total abatement costs
	capture program drop ExoLearnTAC //1 equation, 2 parameters
	program define ExoLearnTAC
		args lnf varphi g_varphi  
		tempvar error1 
		quietly {
		gen `error1'=GDP_ratio -  exp(-0.5*`varphi'*exp(-`g_varphi'*t2020) * Abatement^2)   
		sum `error1' 
		scalar stdev1=`r(sd)'
		replace `lnf' = ln(normalden(`error1',stdev1))  
		}
	end
	ml model lf ExoLearnTAC /varphi /g_varphi if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=.  , cluster(id) maximize init(0.00004 0.01 , copy) 
	ml display
	estimates store ExoLearnTAC
//Use only Marginal abatement costs
	capture program drop ExoLearnMAC //1 equation, 2 parameters
	program define ExoLearnMAC
		args lnf varphi g_varphi  
		tempvar error2 
		quietly {
		gen `error2'= MAC      -  `varphi'*exp(-`g_varphi'*t2020)*Abatement 
		sum `error2' 
		scalar stdev2=`r(sd)'
		replace `lnf' = ln(normalden(`error2',stdev2))  
		}
	end
	ml model lf ExoLearnMAC /varphi /g_varphi if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=.  , cluster(id) maximize init(0.00004 0.01 , copy) 
	ml display
	estimates store ExoLearnMAC 
//ENDOGENOUS
//Only total abatement cost
	capture program drop EndoLearnTAC
	program define EndoLearnTAC
		args lnf varphi chi  
		tempvar error1  	
		quietly {
		gen `error1'=GDP_ratio -  exp(-0.5*`varphi'*((A_Cum0+A_Cum2020)/A_Cum0)^(-`chi') * Abatement^2)   
		sum `error1' 
		scalar stdev1=`r(sd)'	
		replace `lnf' = ln(normalden(`error1',stdev1)) 
		}
	end
	ml model lf EndoLearnTAC /varphi /chi if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=.  , cluster(id)  maximize init(0.0005 0.05 , copy) difficult 
	ml display
	estimates store EndoLearnTAC
//Only marginal abatement cost
	capture program drop EndoLearnMAC
	program define EndoLearnMAC
		args lnf varphi chi  
		tempvar error2  	
		quietly {
		gen `error2'= MAC      -  `varphi'*((A_Cum0+A_Cum2020)/A_Cum0)^(-`chi')*Abatement 
		sum `error2' 
		scalar stdev2=`r(sd)'
		replace `lnf' = ln(normalden(`error2',stdev2))
		}
	end
	ml model lf EndoLearnMAC /varphi /chi if E_Cum2100<E_Cum2100_max & year>=2020 & MAC!=.  , cluster(id)  maximize init(0.0005 0.05 , copy) difficult 
	ml display
	estimates store EndoLearnMAC
