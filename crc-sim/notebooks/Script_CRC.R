	## Load Library 'dga' ##
# install.packages('dga')
# install.packages('chron')
# install.packages('drpop')
library(drpop)
library(dga)
library(chron)

	## Load Library 'drpop' ##
data = simuldata(500, l = 3, K = 3)$data 
str(data)

	# TMLE using logistic regression #
qhat = popsize(data = data, funcname = c("logit"), nfolds = 2, margin = 0.005) # core function #1
psin_estimate <- popsize(data = data, getnuis = qhat$nuis, idfold = qhat$idfold) # core function #2

# Sample size estimates using Doubly Robust estimation
estimated_N_logit_DR <- psin_estimate$result$n[1] # Estimated population size
estimated_N_logit_DR_LB <- psin_estimate$result$cin.l[1] # 95% lower limit of the estimated population size 
estimated_N_logit_DR_UB <- psin_estimate$result$cin.u[1] # 95% upper limit of the estimated population size

# Sample size estimates using Plug-in estimation
estimated_N_logit_PI <- psin_estimate$result$n[2] 
estimated_N_logit_PI_LB <- psin_estimate$result$cin.l[2]
estimated_N_logit_PI_UB <- psin_estimate$result$cin.u[2]

# Sample size estimates using TML estimation
estimated_N_logit_TMLE <- psin_estimate$result$n[3] 
estimated_N_logit_TMLE_LB <- psin_estimate$result$cin.l[3]
estimated_N_logit_TMLE_UB <- psin_estimate$result$cin.u[3]

	# TMLE using generalized additive model #
qhat = popsize(data = data, funcname = c("gam"), nfolds = 2, margin = 0.005) 
psin_estimate <- popsize(data = data, getnuis = qhat$nuis, idfold = qhat$idfold) 
estimated_N_gam_DR <- psin_estimate$result$n[1] 
estimated_N_gam_DR_LB <- psin_estimate$result$cin.l[1]
estimated_N_gam_DR_UB <- psin_estimate$result$cin.u[1]

estimated_N_gam_PI <- psin_estimate$result$n[2] 
estimated_N_gam_PI_LB <- psin_estimate$result$cin.l[2]
estimated_N_gam_PI_UB <- psin_estimate$result$cin.u[2]

estimated_N_gam_TMLE <- psin_estimate$result$n[3] 
estimated_N_gam_TMLE_LB <- psin_estimate$result$cin.l[3]
estimated_N_gam_TMLE_UB <- psin_estimate$result$cin.u[3]


	# TMLE using "rangerlogit" #
qhat = popsize(data = data, funcname = c("rangerlogit"), nfolds = 2, margin = 0.005) 
psin_estimate <- popsize(data = data, getnuis = qhat$nuis, idfold = qhat$idfold) 
estimated_N_rangerlogit_DR <- psin_estimate$result$n[1] 
estimated_N_rangerlogit_DR_LB <- psin_estimate$result$cin.l[1]
estimated_N_rangerlogit_DR_UB <- psin_estimate$result$cin.u[1]

estimated_N_rangerlogit_PI <- psin_estimate$result$n[2] 
estimated_N_rangerlogit_PI_LB <- psin_estimate$result$cin.l[2]
estimated_N_rangerlogit_PI_UB <- psin_estimate$result$cin.u[2]

estimated_N_rangerlogit_TMLE <- psin_estimate$result$n[3] 
estimated_N_rangerlogit_TMLE_LB <- psin_estimate$result$cin.l[3]
estimated_N_rangerlogit_TMLE_UB <- psin_estimate$result$cin.u[3]

	# TMLE using "ranger" #
qhat = popsize(data = data, funcname = c("ranger"), nfolds = 2, margin = 0.005) 
psin_estimate <- popsize(data = data, getnuis = qhat$nuis, idfold = qhat$idfold) 
estimated_N_ranger_DR <- psin_estimate$result$n[1] 
estimated_N_ranger_DR_LB <- psin_estimate$result$cin.l[1]
estimated_N_ranger_DR_UB <- psin_estimate$result$cin.u[1]

estimated_N_ranger_PI <- psin_estimate$result$n[2] 
estimated_N_ranger_PI_LB <- psin_estimate$result$cin.l[2]
estimated_N_ranger_PI_UB <- psin_estimate$result$cin.u[2]

estimated_N_ranger_TMLE <- psin_estimate$result$n[3] 
estimated_N_ranger_TMLE_LB <- psin_estimate$result$cin.l[3]
estimated_N_ranger_TMLE_UB <- psin_estimate$result$cin.u[3]
	save.image("Script_CRC.RData")




