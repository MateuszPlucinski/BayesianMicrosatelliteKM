setwd("\\\\cdc.gov\\private\\M326\\wif7\\Microsatellite\\Bayesian KM")

library(survival)

TESdata = read.csv("DRCDataforKaplanMeier.csv")
TESdata$arm = paste(TESdata$Site,TESdata$Arm)

armnames = unique(TESdata$arm)
arms = lapply(armnames, function (x) TESdata$arm == x)
duration = rep(28,length(armnames))
duration[grep("DP",armnames)] = 42
outcome_uncorrected=rep(0,length(TESdata$Classification))
outcome_uncorrected[TESdata$Classification %in% c("ECT","EPT")] = 1

lastvisit_uncorrected = as.numeric(sapply(1:dim(TESdata)[1], function (x) strsplit(as.character(TESdata$Jour.de.la.classification[x]),"our")[[1]][2]))

TESdata$outcome = outcome_uncorrected
TESdata$lastvisit = lastvisit_uncorrected



fitKM = function(lastvisit,outcome,trial_duration) {
	survivaldata = Surv(lastvisit,outcome)
	km_model = survfit(survivaldata~1,conf.lower = "peto")
	km_estimate = summary(km_model,times=trial_duration)$surv*100
	km_loCI = summary(km_model,times=trial_duration)$lower*100
	km_upCI = summary(km_model,times=trial_duration)$upper*100
	c(km_estimate,km_loCI,km_upCI)
}
bootstrap_KM = function(data,missingvalue,trial_duration) {
	posteriorprob = trimws(data$Prob_Recr)
	posteriorprob[posteriorprob == "" & data$Classification %in% c("PDV")] = 0 
	posteriorprob = as.numeric(as.character(posteriorprob ))
	posteriorprob[is.na(posteriorprob)] = missingvalue 
	if (sum(posteriorprob) > 0) {
		simulated_outcome_matrix = t(sapply(posteriorprob , function (x) sample(c(0,1), nruns, replace = TRUE, prob = c(1-x,x))))
		km_estimated = sapply(1:nruns, function (x) fitKM(data$lastvisit,simulated_outcome_matrix[,x],trial_duration))
		result = 	c(mean(km_estimated[1,]), mean(km_estimated[2,]), mean(km_estimated[3,]))
	} else {
		result = c(100,100,100)
	}
	paste(format(result[1],digits=3)," (",format(result[2],digits=2),"-",format(result[3],digits=2),")",sep="")
}

survivaldata_uncorrected = sapply(arms, function (x) Surv(lastvisit_uncorrected[x],outcome_uncorrected[x]))
km_uncorrected = lapply(survivaldata_uncorrected,function (x) survfit(x~1,conf.lower = "peto"))
survival_uncorrected = sapply(1:length(arms),function (x) paste(format(summary(km_uncorrected[[x]],times=duration[x])$surv*100, digits=3)," (",format(summary(km_uncorrected[[x]],times=duration[x])$lower*100, digits=2),"–",format(summary(km_uncorrected[[x]],times=duration[x])$upper*100, digits=2),")",sep=""))

nruns= 1000
mean_posteriorprob = mean(as.numeric(as.character(TESdata$Prob_Recr[TESdata$Classification %in% c("ECT","EPT")])),na.rm=TRUE)
survival_corrected50 = sapply(1:length(arms),function (x) bootstrap_KM(TESdata[arms[[x]],],mean_posteriorprob,duration[x]))
survival_corrected0 = sapply(1:length(arms),function (x) bootstrap_KM(TESdata[arms[[x]],],0,duration[x]))

write.csv(cbind(armnames,survival_uncorrected,survival_corrected0,survival_corrected50),"KM_estimates.csv")
