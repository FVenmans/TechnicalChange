This depository accompanies the paper "Exogenous and Endogenous Technical Change: making sense of the differences" published in JEEM 2025. It contains the data and code to reproduce the parameters reported in Table 2. 

Stata file 1 and 2 make the data from the exel files that were downloaded from © AR6 Scenario Explorer and Database hosted by IIASA https://data.ece.iiasa.ac.at/ar6 and © NGFS Phase 4 Scenario Explorer https://data.ene.iiasa.ac.at/ngfs. The result is file IPCC_NGFS2024.dta, which is also provided. 

The main analysis in in the code "4 MaxLikelihoodLearning", which starts with uploading the dataset IPCC_NGFS2024.dta. As the maximum likelihood estimator takes in 2 equations, this boils down to a GMM estimator with zero weight on the correlation between the errors of the marginal abatement cost equation and the total abatement cost equation.
