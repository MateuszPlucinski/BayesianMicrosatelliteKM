### This sample R script imports outcome data from antimalarial clinical efficacy trials to calculate uncorrected and corrected
### efficacy estimates using the Kaplan-Meier estimate of the survival function and a posterior probability generated from Bayesian analysis
### of genotyping data
### The input file should follow the format as specified in the "Template_DataforKaplanMeier.csv" file:
###		 SampleID: sample id
###		 Site: name of the site
###		 Drug: name of the drug
###		 DurationFollowUp: number of days of the total follow up (typically, 28 for AL/ASAQ, 42 for DP/Pyramax)
###		 DayLastVisit: last day of each participant's visit
###		 Outcome: final classification for each participant (must be one of the following):
###				ACPR: adequate clinical and parasitological response
###				ETF: early treatment failure
###				LTF: late treatment failure
###				MISSING: lost to follow up (data included up until last day of follow up, i.e. censored)
###				EXCLUDED: excluded (no data will be included for these participants)
###		 PosteriorProbRecrudescence: Posterior probability of recrudescence, as calculated from Bayesian analysis of genotyping data. 
###							        This should be blank for any participants that are not late treatment failures 

### INPUT: the user should specify the following 4 inputs
### 		 inputfile: file name with input data
### 		 nruns: number of bootstrap samples (1000 should be enough) 		
### 		 missingvalueLostToFollowup: how to analyze participants lost to follow up. 
###								A value of 0 means that these participants are not considered to have any chance of
###						 		being treatment failures. A non-zero value means that they have some imputed chance of being treatment failures.
###						 		being treatment failures.
###								Standard analysis would have this parameter set to 0. However, there are some groups that advocate differently. 
###								For a good discussion, see Dahal et al. 2019:
###									 https://malariajournal.biomedcentral.com/articles/10.1186/s12936-019-2837-4 
### 		 missingvalueLTF: similarly, this determines how patients with late treatment failure but with missing genotyping data are analyzed
### 		 						Standard analysis would have this parameter set to 0.
### 								However, some groups set this value to the average of the posterior probability of recrudescence for samples with available genotyping data
	
### OUTPUT: the script creates a csv file with the following variables for each arm: 
### 		 survival_uncorrected: estimate and 95% credible interval for the uncorrected Kaplan-Meier derivation of the efficacy
### 		 survival_corrected: estimate and 95% credible interval for the corrected Kaplan-Meier derivation of the efficacy
 		
### ~Feel free to email mplucinski@cdc.gov with any questions~

inputfile = "Template_DataforKaplanMeier.csv"
nruns= 1000
missingvalueLostToFollowup=0
missingvalueLTF=0

library(survival)

TESdata = read.csv(inputfile)

### 

TESdata$arm = paste(TESdata$Site,TESdata$Drug)
TESdata$Outcome = toupper(trimws(TESdata$Outcome ))

### exclude any participants classified as excluded
TESdata = TESdata[TESdata$Outcome != "EXCLUDED",]

armnames = unique(TESdata$arm)
arms = lapply(armnames, function (x) TESdata$arm == x)
duration = sapply(arms, function (x) TESdata$DurationFollowUp[which(x)[1]])

fitKM = function(lastvisit,outcome,trial_duration) {
	survivaldata = Surv(lastvisit,outcome)
	km_model = survfit(survivaldata~1,conf.lower = "peto")
	km_estimate = summary(km_model,times=trial_duration)$surv*100
	km_loCI = summary(km_model,times=trial_duration)$lower*100
	km_upCI = summary(km_model,times=trial_duration)$upper*100
	c(km_estimate,km_loCI,km_upCI)
}
bootstrap_KM_uncorrected = function(data,missingvalueLostToFollowup,trial_duration) {
	posteriorprob = as.numeric(as.character(data$PosteriorProbRecrudescence))
	posteriorprob[data$Outcome %in% c("MISSING")] = missingvalueLostToFollowup
	posteriorprob[data$Outcome %in% c("LTF")] = 1 ### for uncorrected, all LTFs are coded as failures
	posteriorprob[data$Outcome %in% c("ETF")] = 1 ### early treatment failures are always coded as failures
	posteriorprob[is.na(posteriorprob)] = 0 

	if (sum(posteriorprob) > 0) {
		simulated_outcome_matrix = t(sapply(posteriorprob , function (x) sample(c(0,1), nruns, replace = TRUE, prob = c(1-x,x))))
		km_estimated = sapply(1:nruns, function (x) fitKM(data$DayLastVisit,simulated_outcome_matrix[,x],trial_duration))
		result = 	c(mean(km_estimated[1,]), mean(km_estimated[2,]), mean(km_estimated[3,]))
	} else {
		result = c(100,100,100)
	}
	paste(format(result[1],digits=3)," (",format(result[2],digits=2),"-",format(result[3],digits=2),")",sep="")
}
bootstrap_KM_corrected = function(data,missingvalueLTF,missingvalueLostToFollowup,trial_duration) {
	posteriorprob = as.numeric(as.character(data$PosteriorProbRecrudescence))
	posteriorprob[data$Outcome %in% c("MISSING")] = missingvalueLostToFollowup*missingvalueLTF
	posteriorprob[data$Outcome %in% c("LTF") & is.na(posteriorprob)] = missingvalueLTF
	posteriorprob[data$Outcome %in% c("ETF")] = 1 ### early treatment failures are always coded as failures
	posteriorprob[is.na(posteriorprob)] = 0 

	if (sum(posteriorprob) > 0) {
		simulated_outcome_matrix = t(sapply(posteriorprob , function (x) sample(c(0,1), nruns, replace = TRUE, prob = c(1-x,x))))
		km_estimated = sapply(1:nruns, function (x) fitKM(data$DayLastVisit,simulated_outcome_matrix[,x],trial_duration))
		result = 	c(mean(km_estimated[1,]), mean(km_estimated[2,]), mean(km_estimated[3,]))
	} else {
		result = c(100,100,100)
	}
	paste(format(result[1],digits=3)," (",format(result[2],digits=2),"-",format(result[3],digits=2),")",sep="")
}

survival_uncorrected = sapply(1:length(arms),function (x) bootstrap_KM_uncorrected (TESdata[arms[[x]],],missingvalueLostToFollowup,duration[x]))
survival_corrected = sapply(1:length(arms),function (x) bootstrap_KM_corrected(TESdata[arms[[x]],],missingvalueLTF,missingvalueLostToFollowup,duration[x]))

write.csv(cbind(armnames,survival_uncorrected,survival_corrected),"KM_estimates.csv")
