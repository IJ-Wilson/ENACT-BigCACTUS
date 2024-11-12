################################################################################
## Annotated code to calculate key parameters from the CACTUS pilot trial     ##
## Laura Flight                                                               ##
## 31Aug2022                                                                  ##
################################################################################

#####################
## Setup workspace ##
#####################

# Load required packages -------------------------------------------------------

rm(list=ls())
require(devtools)

install_version("heemod", version = "0.13.0", repos = "http://cran.us.r-project.org")

pkgs<-c("dplyr","ggplot2","reshape2","foreach","doParallel","xtable", "heemod")
inst<-lapply(pkgs,library,character.only=TRUE)

# Set working directory --------------------------------------------------------

setwd("xxxxxx/Value-Based Sequential 2-Arm Design")

# Load functions (should be saved in the working directory) --------------------

source("Programs/Functions/Unadjusted_pilot_data_analysis.R")
source("Programs/Functions/CACTUS_HE_model.R")

################################################################################

#####################
## Load pilot data ##
#####################

# Data has been created from the original health economic model analysis -------

data.cactus<-read.csv("Data/Raw/pilot_data_set.csv")

# Set trt to be 0=Control and 1=Intervention -----------------------------------

data.cactus<-data.cactus%>%
  mutate(trt=ifelse(trt=="Intervention",1,0))

################################################################################

########################
## Analyse pilot data ##
########################

pilot_param<-EVSI_analysis_pilot(terminal_analysis=1,
                                 
                                 terminal_design='FIX',
                                 
                                 terminal_analysis_number=1,
                                 
                                 terminal_n_intervention=data.cactus %>%
                                   filter(trt==1) %>%
                                   nrow(),
                                 
                                 terminal_n_control=data.cactus %>%
                                   filter(trt==0) %>%
                                   nrow(),
                                 
                                 terminal_dataset=data.cactus, 
                                 
                                 alpha=0.05,
                                 
                                 beta=0.1,
                                 
                                 rho_time=0.5)

################################################################################

#############################################
## Fit health economic model to pilot data ##
#############################################

pilot_model_fit<-HE_model_func(in_age=67.91667,
                               in_prob_death_disease=0.0111503621,
                               in_prob_relapse_month_int=unlist(pilot_param$unadjusted_list$p_relapse_month_unadjusted),
                               in_prob_good_response_int=unlist(pilot_param$unadjusted_list$p_good_6_intervention_unadjusted),
                               in_cost_intervention= 769.25,
                               in_cost_cont=unlist(pilot_param$unadjusted_list$resource_cost_control_unadjusted),
                               in_cost_int_resp=unlist(pilot_param$unadjusted_list$resource_cost_intervention_unadjusted), #cost not including intervention cost
                               in_cost_int_non_resp=unlist(pilot_param$unadjusted_list$resource_cost_intervention_unadjusted),
                               in_dis_no_resp_cont=1-unlist(pilot_param$unadjusted_list$baseline_utility_control_unadjusted),
                               in_utility_change_int=unlist(pilot_param$unadjusted_list$utility_increment_unadjusted),
                               threshold=20000)

################################################################################

##################################################################
## Bootstrapping to generate probabilistic sensitivity analysis ##
##################################################################

B <- 1000
SEED <- 0
PSA_mat <- matrix(NA,ncol=44,nrow=B)
cl <- makeCluster(3)
registerDoParallel(cl)

PSA_mat <- foreach(b=1:B,.combine='rbind',.packages=c("dplyr","tidyr")) %dopar% {
  RNGversion("3.4.3")
  set.seed(12+b+SEED)

  boot_dat<-rbind(data.frame(sample_n(data.cactus%>%
                                        filter(trt==0),
                                      data.cactus%>%
                                        filter(trt==0)%>%
                                        nrow(),
                                      replace=TRUE)),
                  data.frame(sample_n(data.cactus%>%
                                        filter(trt==1),
                                      data.cactus%>%
                                        filter(trt==1)%>%
                                        nrow(),
                                      replace=TRUE)))
  
  boot_param<-EVSI_analysis_pilot(terminal_analysis=1,
                                  terminal_design='FIX',
                                  terminal_analysis_number=1,
                                  terminal_n_intervention=boot_dat%>%
                                    filter(trt==1)%>%
                                    nrow(),
                                  terminal_n_control=boot_dat%>%
                                    filter(trt==0)%>%
                                    nrow(),
                                  terminal_dataset=boot_dat, 
                                  alpha=0.05,
                                  beta=0.1 ,
                                  rho_time=0.5)
  
  # model parameters
  boot_model_param<-boot_param$data_param
  
  # data parameters
  boot_data_param<-boot_param$unadjusted_list
  
  boot_row<-matrix(unlist(boot_param),nrow=1)
}
stopCluster(cl)
registerDoSEQ()

PSA_mat <- data.frame(PSA_mat)
colnames(PSA_mat) <- names(unlist(pilot_param))

# Fit the model to the bootstrapped values -------------------------------------

model_boot_run <- mapply(HE_model_func,in_age=67.91667,
                       in_prob_death_disease=0.0111503621,
                       in_prob_relapse_month_int=PSA_mat$unadjusted_list.p_relapse_month_unadjusted,
                       in_prob_good_response_int=PSA_mat$unadjusted_list.p_good_6_intervention_unadjusted,
                       in_cost_intervention= 769.25,
                       in_cost_cont=PSA_mat$unadjusted_list.resource_cost_control_unadjusted.mean,
                       in_cost_int_resp=PSA_mat$unadjusted_list.resource_cost_intervention_unadjusted.mean, #cost not including intervention cost
                       in_cost_int_non_resp=PSA_mat$unadjusted_list.resource_cost_intervention_unadjusted.mean,
                       in_dis_no_resp_cont=1-PSA_mat$unadjusted_list.baseline_utility_control_unadjusted.mean,
                       in_utility_change_int=PSA_mat$unadjusted_list.utility_increment_unadjusted.mean,
                       threshold=20000)

model_boot_run <- (t(model_boot_run))

# Mean INB based on probabilistic health economic model ------------------------

prob_mean_INB <- mean(unlist(model_boot_run[,3]))

# SD INB based on probabilistic health economic model --------------------------

prob_sd_INB <- sd(unlist(model_boot_run[,3]))
SX <- prob_sd_INB*sqrt(7)

################################################################################

##################
## VOI analysis ##
##################

# Calculating EVSI and EVPI ----------------------------------------------------

EVSI_norm<-function(m,m0,n0,sx){
  s0=sqrt(sx^2/n0)
  s0.=sqrt((sx^2/n0)*(m/(n0+m)))
  u0.=abs(m0)/s0.
  EVSI<-s0.*(dnorm(u0.) - u0. * pnorm(u0., lower.tail=FALSE))
  u0=abs(m0)/s0
  EVPI<-s0*(dnorm(u0) - u0 * pnorm(u0, lower.tail=FALSE))
  out=list(EVSI=EVSI,EVPI=EVPI)
  out
}

EVSI<-data.frame(unlist(mapply(EVSI_norm,m=seq(1,1000,1),
                               m0=prob_mean_INB,
                               n0=7,
                               sx=SX)[1,]),unlist(mapply(EVSI_norm,
                                                         m=seq(1,1000,1),
                                                         m0=prob_mean_INB,
                                                         n0=7,
                                                         sx=SX)[2,]))

EVSI_dat<-data.frame(cbind(EVSI,seq(1,1000,1)))
colnames(EVSI_dat)<-c("EVSI","EVPI","m")

fixed_cost<-0
variable_cost<-4706
population_benefit<-215378

ENBS_dat<-EVSI_dat%>%
  mutate(cost_sampling=fixed_cost+m*variable_cost)%>%
  mutate(popEVSI=EVSI*population_benefit)%>%
  mutate(ENBS=popEVSI-cost_sampling)

ENBS_dat%>%
  filter(m==435)
