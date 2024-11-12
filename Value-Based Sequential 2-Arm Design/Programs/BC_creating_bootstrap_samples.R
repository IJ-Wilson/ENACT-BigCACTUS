################################################################################
## Annotated code to bootstrap the Big CACTUS data and calculate prior/       ##
## posterior parameters                                                       ##
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

pkgs<-c("dplyr","ggplot2","reshape2","tidyr","tibble","xtable","heemod","foreach","doParallel","scales")
inst<-lapply(pkgs,library,character.only=TRUE)
  
# Set working directory --------------------------------------------------------

setwd("xxxxxx/Value-Based Sequential 2-Arm Design")

################################################################################

# Create function --------------------------------------------------------------

boot_func_block<-function(B, # number of boots
                          a, # additional seed number
                          Tmax, # sample size of trial
                          block){
  
  plot.dat.mat<-NA
  cl<-makeCluster(7)
  registerDoParallel(cl)
  
  plot.dat.mat<-foreach(b=1:B,.packages=c("dplyr","heemod")) %dopar% {
    set.seed(12+b+a)
    
    # Load data ----------------------------------------------------------------         
    DATA<-read.csv("Data/Raw/data_set.csv")
    
    # Load functions (should be saved in the working directory) ----------------
    source("Programs/Functions/Big_CACTUS_R_model.R")
    source("Programs/Functions/BC_model_param_function.R")

    # Set parameters -----------------------------------------------------------
    n0<-7 # prior sample size
    tau<-ceiling(95*(Tmax/95)/((638.751*(Tmax/95))/365.25) ) # delay
    
    # Bootstrap the data with replacement in each of the three trial arms ------
    set.seed(1+b+a)
    boot.dat<-rbind(data.frame(sample_n(DATA%>%
                                          filter(random_group==1),
                                        ceiling((Tmax/95)*95),
                                        replace=TRUE)),
                    data.frame(sample_n(DATA%>%
                                          filter(random_group==2),
                                        ceiling((Tmax/95)*95),
                                        replace=TRUE)),
                    data.frame(sample_n(DATA%>%
                                          filter(random_group==3),
                                        ceiling((Tmax/95)*95),
                                        replace=TRUE)))
    
    set.seed(1+b+a)
    data1<-boot.dat%>%
      filter(random_group==1)%>%
      mutate(ordering=seq(1:ceiling((Tmax/95)*95)))
    
    data2<-boot.dat%>%
      filter(random_group==2)%>%
      mutate(ordering=seq(1:ceiling((Tmax/95)*95)))
    
    data3<-boot.dat%>%
      filter(random_group==3)%>%
      mutate(ordering=seq(1:ceiling((Tmax/95)*95)))
    
    # Data set is now sequential so patients allocated to 1,2,3,1,2,3,1,2,3... -
    paired.data <- rbind(data1,data2,data3)
    paired.data <- paired.data[order(paired.data$ordering),]
    
    # FIXED DESIGNS ------------------------------------------------------------
    
    # Just want to fit the model to the full dataset 
    fixed_param<-param_func(paired.data)
    fixed_model<-BC_HE_model_func(in.age=65,
                               in.prob.death.disease=0.024690, # 3-month probability
                               in.p_relapse_6_9_CSLT=unlist(fixed_param$relapse_9_CSLT),
                               in.p_relapse_6_9_UC=unlist(fixed_param$relapse_9_UC),
                               in.p_relapse_6_9_AC=unlist(fixed_param$relapse_9_AC),
                               # Probability of relapse (9–12 months)
                               in.p_relapse_9_12_CSLT=unlist(fixed_param$relapse_12_CSLT),
                               in.p_relapse_9_12_UC=unlist(fixed_param$relapse_12_UC),
                               in.p_relapse_9_12_AC=unlist(fixed_param$relapse_12_AC),
                               # Probability of relapse (12 months onwards)
                               in.p_relapse_12_CSLT=unlist(fixed_param$relapse_12_CSLT),
                               in.p_relapse_12_UC=unlist(fixed_param$relapse_12_UC),
                               in.p_relapse_12_AC= unlist(fixed_param$relapse_12_AC),
                               # Probability of good response (0–6 months)
                               in.p_good_0_6_CSLT=1-exp(-(-(log(1-unlist(fixed_param$response_6_CSLT))/6))*3),
                               in.p_good_0_6_UC=1-exp(-(-(log(1-unlist(fixed_param$response_6_UC))/6))*3),
                               in.p_good_0_6_AC=1-exp(-(-(log(1-unlist(fixed_param$response_6_AC))/6))*3),
                               # Probability of new good response (6–9 months)
                               in.p_good_6_9_CSLT=unlist(fixed_param$response_9_CSLT),
                               in.p_good_6_9_UC=unlist(fixed_param$response_9_UC),
                               in.p_good_6_9_AC=unlist(fixed_param$response_9_AC),
                               # Probability of new good responseor renewed response in people who responded at 6 months and relapsed at 9 months (9–12 months)
                               in.p_good_12_CSLT=unlist(fixed_param$response_12_CSLT),
                               in.p_good_12_UC=unlist(fixed_param$response_12_UC),
                               in.p_good_12_AC=unlist(fixed_param$response_12_AC),
                               # Costs
                               in.cost_CSLT=unlist(fixed_param$cost_CSLT),
                               in.cost_UC=0,
                               in.cost_AC=unlist(fixed_param$cost_AC),
                               in.u_aphasia=fixed_param$utility0,
                               in.u_response_6=fixed_param$utility0+fixed_param$u_response_6,
                               in.u_response_9=fixed_param$utility0+fixed_param$u_response_9,
                               in.u_response_12=fixed_param$utility0+fixed_param$u_response_12,
                               in.cycle.length=3,
                               in.threshold=20000 ,
                               in.DR=0.035)
    
    fixed_mu<-fixed_model$INB_CSLT_UC
    fixed_mu_n<-rep(NA,2) # +1 because there is the prior and then the trial allocations
    fixed_mu_n[1]<-3190.78
    fixed_m<-c(7,7+Tmax)
 
    fixed_mu_n[2]<-((fixed_m[1]*fixed_mu_n[1])+((Tmax)*fixed_mu[1]))/(fixed_m[1]+(Tmax))
   
    # Number to include in the analysis (block sizes plus the incomplete block if not a multiple)
    n_seq<-unique(c(seq(block*3,ceiling(95*(Tmax/95))*3,block*3),
             ceiling(95*(Tmax/95))*3) )
  
    # SEQUENTIAL DESIGNS -------------------------------------------------------
    
    param<-list()
    for(i in n_seq){
      param[[i]]<- param_func(paired.data[1:i,])
    }
    
    param2<-param[-which(sapply(param, is.null))]
    INMB_CSLT_UC<-rep(NA)
    INMB_CSLT_AC<-rep(NA)
    n_CSLT<-rep(NA)
    n_UC<-rep(NA)
    n_AC<-rep(NA)
    QALY.diff_CSLT_UC<-rep(NA)
    cost.diff_CSLT_UC<-rep(NA)
    param<-matrix(NA,nrow=length(n_seq),ncol=30)
  # ptm <- proc.time()
    
    for(i in 1:length(n_seq)){
      X<-BC_HE_model_func(in.age=65,
                       in.prob.death.disease=0.024690, # 3-month probability
                       in.p_relapse_6_9_CSLT=unlist(param2[[i]]$relapse_9_CSLT),
                       in.p_relapse_6_9_UC=unlist(param2[[i]]$relapse_9_UC),
                       in.p_relapse_6_9_AC=unlist(param2[[i]]$relapse_9_AC),
                       # Probability of relapse (9–12 months)
                       in.p_relapse_9_12_CSLT=unlist(param2[[i]]$relapse_12_CSLT),
                       in.p_relapse_9_12_UC=unlist(param2[[i]]$relapse_12_UC),
                       in.p_relapse_9_12_AC=unlist(param2[[i]]$relapse_12_AC),
                       # Probability of relapse (12 months onwards)
                       in.p_relapse_12_CSLT=unlist(param2[[i]]$relapse_12_CSLT),
                       in.p_relapse_12_UC=unlist(param2[[i]]$relapse_12_UC),
                       in.p_relapse_12_AC= unlist(param2[[i]]$relapse_12_AC),
                       # Probability of good response (0–6 months)
                       in.p_good_0_6_CSLT=1-exp(-(-(log(1-unlist(param2[[i]]$response_6_CSLT))/6))*3),
                       in.p_good_0_6_UC=1-exp(-(-(log(1-unlist(param2[[i]]$response_6_UC))/6))*3),
                       in.p_good_0_6_AC=1-exp(-(-(log(1-unlist(param2[[i]]$response_6_AC))/6))*3),
                       # Probability of new good response (6–9 months)
                       in.p_good_6_9_CSLT=unlist(param2[[i]]$response_9_CSLT),
                       in.p_good_6_9_UC=unlist(param2[[i]]$response_9_UC),
                       in.p_good_6_9_AC=unlist(param2[[i]]$response_9_AC),
                       # Probability of new good response renewed response in people who responded at 6 months and relapsed at 9 months (9–12 months)
                       in.p_good_12_CSLT=unlist(param2[[i]]$response_12_CSLT),
                       in.p_good_12_UC=unlist(param2[[i]]$response_12_UC),
                       in.p_good_12_AC=unlist(param2[[i]]$response_12_AC),
                       # Costs
                       in.cost_CSLT=unlist(param2[[i]]$cost_CSLT),
                       in.cost_UC=0,
                       in.cost_AC=unlist(param2[[i]]$cost_AC),
                       in.u_aphasia=param2[[i]]$utility0,
                       in.u_response_6=param2[[i]]$utility0+param2[[i]]$u_response_6,
                       in.u_response_9=param2[[i]]$utility0+param2[[i]]$u_response_9,
                       in.u_response_12=param2[[i]]$utility0+param2[[i]]$u_response_12,
                       in.cycle.length=3,
                       in.threshold=20000,
                       in.DR=0.035)
      
      param[i,]<-(X$param)
      INMB_CSLT_UC[i]<-X$INB_CSLT_UC
      # INMB_CSLT_AC[i]<-X$INB_CSLT_AC
      QALY.diff_CSLT_UC[i]<-X$QALY.diff_CSLT_UC
      cost.diff_CSLT_UC[i]<-X$cost.diff_CSLT_UC
      
      n_CSLT[i]<-(param2[[i]]$n_CSLT)
      n_UC[i]<-param2[[i]]$n_UC
      # n_AC[i]<-param2[[i]]$n_AC
    }
    
    # Storing the data to create a pliot of the INB observed during the trial
    plot.dat<-data.frame(cbind(na.omit(unlist(n_CSLT)),
                               # na.omit(unlist(n_AC)),
                               na.omit(unlist(n_UC)),
                               na.omit(param),
                               na.omit(INMB_CSLT_UC),
                               # na.omit(INMB_CSLT_AC),
                               na.omit(QALY.diff_CSLT_UC),
                               na.omit(cost.diff_CSLT_UC),
                               unique(c(seq(block*2,ceiling(95*(Tmax/95))*2,block*2),
                                        ceiling(95*(Tmax/95))*2) )))
    colnames(plot.dat)<-c("n_CSLT",
                        # "n_AC",
                          "n_UC",
                          "in.age",
                          "in.prob.death.disease", # 3-month probability
                          "in.p_relapse_6_9_CSLT",
                          "in.p_relapse_6_9_UC",
                          "in.p_relapse_6_9_AC",
                          # Probability of relapse (9–12 months)
                          "in.p_relapse_9_12_CSLT",
                          "in.p_relapse_9_12_UC",
                          "in.p_relapse_9_12_AC" ,
                          # Probability of relapse (12 months onwards)
                          "in.p_relapse_12_CSLT",
                          "in.p_relapse_12_UC",
                          "in.p_relapse_12_AC",
                          # Probability of good response (0–6 months)
                          "in.p_good_0_6_CSLT",
                          "in.p_good_0_6_UC",
                          "in.p_good_0_6_AC",
                          # Probability of new good response (6–9 months)
                          "in.p_good_6_9_CSLT",
                          "in.p_good_6_9_UC",
                          "in.p_good_6_9_AC",
                          # Probability of new good response, or renewed response in people who responded at 6 months and relapsed at 9 months (9–12 months)
                          "in.p_good_12_CSLT",
                          "in.p_good_12_UC",
                          "in.p_good_12_AC",
                          # Costs
                          "in.cost_CSLT",
                          "in.cost_UC",
                          "in.cost_AC",
                          "in.u_aphasia",
                          "in.u_response_6",
                          "in.u_response_9",
                          "in.u_response_12",
                          "in.cycle.length",
                          "in.threshold",
                          "in.DR",
                          "INMB_CSLT_UC",
                          # "INMB_CSLT_AC",
                          "QALY_CSLT_UC",
                          "Cost_CSLT_UC",
                          "patients")
    
    # Calculating the prior/posterior INB
    mu<-c(plot.dat$INMB_CSLT_UC)
    mu_n<-rep(NA,length(  unique(c(seq(block*2,ceiling(95*(Tmax/95))*2,block*2),
                                   ceiling(95*(Tmax/95))*2) ))+1) # +1 because there is the prior and then the trial allocations
    
    mu_n[1]<-3190.78
    
    m<-c(7, 7+unique(c(seq(block,ceiling(95*(Tmax/95)),block),
                      ceiling(95*(Tmax/95))) ))
    t<-c(seq(0,Tmax,block),Tmax)
    
    for(i in 2:(length(unique(c(seq(block*2,ceiling(95*(Tmax/95))*2,block*2),
                                         ceiling(95*(Tmax/95))*2) ))+1)){
      mu_n[i]<-((m[1]*mu_n[1])+((t[i])*mu[(i-1)]))/(m[1]+(t[i]))
    }
    
    list(seq=data.frame(c(0,mu),mu_n,(m),m+round(tau),b),
         fixed=data.frame(fixed_mu_n,(fixed_m),fixed_m+round(tau),b))
  }
  stopCluster(cl)
  registerDoSEQ()
  
  colname<-c("cum_INMB","post_INMB","m","m_delay","label","fixed_INMB","fixed_m","fixed_m_delay")
  plot.dat.mat.list<-lapply(lapply(plot.dat.mat, bind_rows), setNames, colname)
  
  plot.dat.mat.frame<-bind_rows(plot.dat.mat.list, .id = 'label2')
  
  write.csv(plot.dat.mat.frame,paste("Data/Derived/Bootstraps/block5_CSLTvsUC_Tmax_",Tmax,"_attempt2","_B_",B,"a_",a,".csv",sep="" ))
}

# Uncomment the scenario you want to run ---------------------------------------

## Scenario 1 ------------------------------------------------------------------
ptm <- proc.time()
boot_func_block(B=2500,
                 a=0,
                 Tmax=95,
                 block=5)
proc.time()-ptm

## Scenario 2 ------------------------------------------------------------------
ptm <- proc.time()
boot_func_block(B=2500,
                a=0,
                Tmax=435,
                block=5)
proc.time()-ptm

## Scenario 3 ------------------------------------------------------------------
ptm <- proc.time()
boot_func_block(B=2500,
                a=2500,
                Tmax=95,
                block=5)
proc.time()-ptm

## Scenario 4 ------------------------------------------------------------------
ptm <- proc.time()
boot_func_block(B=2500,
                a=2500,
                Tmax=435,
                block=5)
proc.time()-ptm
